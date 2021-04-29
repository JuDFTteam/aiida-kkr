#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the workflow that does a decimation calculation for half-infinite crystals
"""

from aiida_kkr.tools import neworder_potential_wf
from aiida_kkr.calculations import VoronoiCalculation
from aiida_kkr.calculations import KkrCalculation
from aiida.engine import WorkChain, ToContext, calcfunction
from aiida.orm import Code, Dict, Int, Float, RemoteData, KpointsData, XyData, StructureData, FolderData
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode
from aiida_kkr.tools import kkrparams
import numpy as np
from masci_tools.io.common_functions import get_Ry2eV

__copyright__ = (u'Copyright (c), 2020, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = u'Philipp Rüßmann'

_eV2Ry = 1.0 / get_Ry2eV()


class kkr_decimation_wc(WorkChain):
    """
    Workchain a decimation calculation with KKR starting from a thick slab (center potential should be bulk-like).

    The workchain
     - creates the structure nodes of the substrate and decimation region from thick slab structure
     - creates starting potentials of the slab structure
     - runs auxiliary voronoi steps to get starting setup correctly for the KKR calculations
     - runs the deci-out step in serial
     - runs the decimation step

     The workflow starts either from a converged thick film of from a previous decimation calculation (skips the structure and starting potential setup and the voronoi steps).

     The workflow parameters input can be:
          {'nkz' : 30,   # number of k-points in z-direction for substrate
           'nprinc': 4,  # number of layer in principle layer
           'nplayer': 4, # number of principle layers (naez deci: nprinc*nplayer)
           'dosmode': False, # run DOS calculation
           'dos_params': {'emin_EF': -5.0, # EMIN-EF in eV
                          'emax_EF':  3.0, # EMAX-EF in eV
                          'nepts': 96,     # number of points in contour
                          'tempr': 100,    # smearing temperature
                          'kmesh': [50, 50, 50]}, # k-mesh used in dos calculation
           }

    :param wf_parameters: (Dict); Workchain specifications
    :param options: (Dict); specifications for the computer (used in decimation step only)
    :param remote_data: (RemoteData), mandatory; either parent slab or previous decimation calculation
    :param kkr: (Code), mandatory; KKR code for running deci-out and decimation steps
    :param voronoi: (Code), mandatory if starting from slab calculation; voronoi code for auxiliary calculations

    :return results: (Dict), Information of workflow results
        like Success, last result node, list with convergence behavior
    :return deci_calc: (RemoteData), Remote data of decimation calculation, used to reuse decimation setup if calculation is continued (e.g. for DOS after scf)
    """

    _workflowversion = __version__
    _wf_label = 'kkr_decimation_wc'
    _wf_description = 'Workflow for a decimation calculation starting from a thick slab or a previous decimation calculation.'
    _wf_default = {
        'nkz': 30,  # number of k-points in z-direction for substrate
        'nprinc': 4,  # number of layer in principle layer
        'nplayer': 4,  # number of principle layers (naez deci: nprinc*nplayer)
        'dosmode': False,  # run DOS calculation
        'dos_params': {
            'emin_EF': -5.0,  # EMIN-EF in eV
            'emax_EF': 3.0,  # EMAX-EF in eV
            'nepts': 96,  # number of points in contour
            'tempr': 100,  # smearing temperature
            'kmesh': [50, 50, 50]
        },  # k-mesh used in dos calculation
    }
    _options_default = {
        'resources': {
            'tot_num_mpiprocs': 1,
            'num_machines': 1
        },  # resources to allowcate for the job
        'max_wallclock_seconds': 60 * 60,  # walltime after which the job gets killed (gets parsed to KKR)
    }

    _keys2d = [
        'INTERFACE',
        '<NLBASIS>',
        '<RBLEFT>',
        'ZPERIODL',
        '<NRBASIS>',
        '<RBRIGHT>',
        'ZPERIODR',  # standard names
        'NLBASIS',
        'RBLEFT',
        'NRBASIS',
        'RBRIGHT'  # version of keywords without brackets
    ]

    # intended to guide user interactively in setting up a valid wf_params node
    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """
        if not silent:
            print('Version of workflow: {}'.format(self._workflowversion))
        return self._wf_default

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_decimation_wc, cls).define(spec)
        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._wf_default),
            help='parameters for decimation setup (used only if not started from previous decimation calculation).'
        )
        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._wf_default),
            help=
            'Computer options used in the deicmation step (voronoi and deci-out steps run serially but use the walltime given here).'
        )
        spec.input(
            'remote_data',
            valid_type=RemoteData,
            required=True,
            help=
            'Parent calculation (either previous decimation calculation, then voronoi steps are skipped or slab calculation).'
        )
        spec.input('kkr', valid_type=Code, required=True, help='KKRhost code.')
        spec.input(
            'voronoi',
            valid_type=Code,
            required=False,
            help='Voronoi code. Only needed if remote_data is slab claculation and not a previous decimation run.'
        )
        spec.input(
            'kpoints',
            valid_type=KpointsData,
            required=False,
            help='If given this triggers a bandstructure (i.e. qdos) calculation.'
        )
        spec.input(
            'calc_parameters',
            valid_type=Dict,
            required=False,
            help='If given overwrites KKR parameters starting from slab params (can be used to run DOS for instance).'
        )

        # define outputs
        spec.output(
            'structure_decimate', valid_type=StructureData, required=True, help='Structure of decimation region.'
        )
        spec.output(
            'structure_substrate',
            valid_type=StructureData,
            required=True,
            help='Structure of substrate lattice continuation.'
        )
        spec.output('out_params_calc_deci_out', valid_type=Dict, help='Output parameter node of deci-out calculation.')
        spec.output(
            'out_params_calc_decimate', valid_type=Dict, help='Output parameter node of decimation calculation.'
        )
        spec.output('out_remote_calc_decimate', valid_type=RemoteData, help='Remote folder of decimation calculation.')
        spec.output(
            'out_retrieved_calc_decimate', valid_type=FolderData, help='Retrieved folder of decimation calculation.'
        )
        # optional output nodes for DOS mode
        spec.output(
            'dos_data',
            valid_type=XyData,
            required=False,
            help='DOS data with finite imaginary part in the energy contour.'
        )
        spec.output(
            'dos_data_interpol', valid_type=XyData, required=False, help='interpolated DOS data onto the real axis.'
        )

        # Here the structure of the workflow is defined
        spec.outline(
            # initialize workflow
            cls.start,
            # validate input
            cls.validate_input,
            # prepare or extract structure and startpot
            cls.prepare_deci_from_slab,
            # run voroaus steps if needed
            cls.run_voroaux,
            # run deci-out step
            cls.run_deciout,
            # run decimatino step
            cls.run_decimation,
            #  collect results and get output nodes
            cls.return_results
        )

        # definition of exit code in case something goes wrong in this workflow
        spec.exit_code(
            300, 'ERROR_INVALID_INPUT_REMOTE_DATA',
            'Given remote_data is not correct (needs to be a slab or decimation KKR calculation)'
        )
        spec.exit_code(
            301, 'ERROR_KKRCODE_NOT_CORRECT', 'The code you provided for kkr does not use the plugin kkr.kkr'
        )
        spec.exit_code(
            302, 'ERROR_VORONOICODE_NOT_CORRECT', 'The code you provided for voronoi does not use the plugin kkr.voro'
        )

    ###################################################################################################################
    # functions from outline

    def start(self):
        """
        init context and some parameters
        """
        self.report('INFO: started KKR decimation workflow version {}' ''.format(self._workflowversion))

        ####### init    #######

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()
        options_dict = self.inputs.options.get_dict()

        # add missing values from defaults
        for k, v in self._wf_default.items():
            if k not in wf_dict:
                wf_dict[k] = v
                self.report('INFO: using default parameter {}:{}'.format(k, v))
        for k, v in self._options_default.items():
            if k not in options_dict:
                options_dict[k] = v
                self.report('INFO: using default option {}:{}'.format(k, v))

        # get basic input parameters for setup of calculation
        self.ctx.nprinc = wf_dict.get('nprinc', self._wf_default['nprinc'])
        self.ctx.nplayer = wf_dict.get('nplayer', self._wf_default['nplayer'])
        self.ctx.nkz = wf_dict.get('nkz', self._wf_default['nkz'])
        self.ctx.dosmode = wf_dict.get('dosmode', self._wf_default.get('dosmode'))
        self.ctx.dos_params = wf_dict.get('dos_params', self._wf_default.get('dos_params'))

        # set values, or defaults
        self.ctx.options = options_dict

        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []

    def validate_input(self):
        """
        # validate input and find out which path (starting from slab of previous decimation calc) to take
        """
        inputs = self.inputs

        # extract correct remote folder of last calculation if input remote_folder node is not from KkrCalculation but kkr_scf_wc workflow
        input_remote = self.inputs.remote_data
        # check if input_remote has single KkrCalculation parent
        parents = input_remote.get_incoming(node_class=KkrCalculation)
        nparents = len(parents.all_link_labels())
        if nparents != 1:
            return self.exit_codes.ERROR_INVALID_INPUT_REMOTE_DATA
        parent_calc = parents.first().node
        # check if parent_calc is decimation calculation
        self.ctx.parent_params = parent_calc.inputs.parameters.get_dict()
        runopts = self.ctx.parent_params.get('RUNOPT')
        if runopts is None:
            runopts = []
        if 'DECIMATE' in runopts:
            self.ctx.parent_is_deci = True
        else:
            self.ctx.parent_is_deci = False
            self.ctx.slab_calc = parent_calc

        # check codes
        try:
            test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT
        if 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_VORONOICODE_NOT_CORRECT

    def prepare_deci_from_slab(self):
        """
        prepare deci-out and decimation steps (creates structure and startpot)
        """
        if not self.ctx.parent_is_deci:
            # create aiida structures for substrate and decimation region
            self._create_structures()
            # copy correct starting potentials from slab calculation
            self._create_startpots()
        else:
            # extract startpots and structures from parent calculation
            deci_remote = self.inputs.remote_data
            deci_calc = deci_remote.get_incoming(node_class=KkrCalculation).first().node

            struc_substrate, voro_substrate = VoronoiCalculation.find_parent_structure(deci_calc.inputs.deciout_parent)
            self.ctx.struc_substrate, self.ctx.voroaux_substrate = struc_substrate, voro_substrate

            struc_deci, voro_deci = VoronoiCalculation.find_parent_structure(deci_remote)
            self.ctx.struc_decimation, self.ctx.voroaux_decimation = struc_deci, voro_deci

            self.ctx.startpot_substrate = voro_substrate.inputs.potential_overwrite
            self.ctx.startpot_decimate = voro_deci.inputs.potential_overwrite

            startpot_calcfun = struc_deci.get_incoming(link_label_filter='result').first().node
            slab_calc_remote = startpot_calcfun.inputs.slab_calc
            self.ctx.slab_calc = slab_calc_remote.get_incoming(node_class=KkrCalculation).first().node

        # get input params for deci-out and decimation steps
        self.ctx.params_overwrite = None
        if 'calc_parameters' in self.inputs:
            self.ctx.params_overwrite = self.inputs.calc_parameters

        # take care of qdos or dos modes
        qdos_mode = 0
        if 'kpoints' in self.inputs:
            qdos_mode = 1
        if qdos_mode == 1 or self.ctx.dosmode:
            if self.ctx.params_overwrite is None:
                params_new = {}
            else:
                params_new = self.ctx.params_overwrite.get_dict()
            # set dos contour from dos_params
            params_new['NPOL'] = 0
            params_new['NPT1'] = 0
            params_new['NPT2'] = self.ctx.dos_params['nepts']
            params_new['IEMXD'] = self.ctx.dos_params['nepts']
            params_new['NPT3'] = 0
            params_new['TEMPR'] = self.ctx.dos_params['tempr']
            ef_ry = self.ctx.slab_calc.outputs.output_parameters['fermi_energy']
            params_new['EMIN'] = ef_ry + self.ctx.dos_params['emin_EF'] * _eV2Ry
            params_new['EMAX'] = ef_ry + self.ctx.dos_params['emax_EF'] * _eV2Ry
            params_new['BZDIVIDE'] = self.ctx.dos_params['kmesh']
            if qdos_mode == 1:
                if params_new['TEMPR'] > 100.:
                    params_new['TEMPR'] = 50.
                params_new['BZDIVIDE'] = [1, 1, 1]
            self.ctx.params_overwrite = Dict(dict=params_new)

        alat_slab = self.ctx.slab_calc.outputs.output_parameters['alat_internal']

        out = make_decimation_param_nodes(
            self.ctx.slab_calc.inputs.parameters, Float(alat_slab), self.ctx.struc_decimation, Int(self.ctx.nkz),
            self.ctx.params_overwrite
        )

        self.ctx.dsubstrate = out['dsubstrate']
        self.ctx.ddecimation = out['ddeci']

        # extract angles from slab calculation
        if 'initial_noco_angles' in self.ctx.slab_calc.inputs:
            init_noco_slab = self.ctx.slab_calc.inputs.initial_noco_angles
            struc_decimation = self.ctx.struc_decimation
            struc_substrate = self.ctx.struc_substrate
            noco_angles = get_noco_angles_deci(init_noco_slab, struc_decimation, struc_substrate)
            self.ctx.noco_angles_substrate = noco_angles['initial_noco_angles_substrate']
            self.ctx.noco_angles_decimation = noco_angles['initial_noco_angles_decimation']

    def run_voroaux(self):
        """
        run auxiliary voronoi steps if needed
        """
        # only needed if parent is not already a decimation calculation
        if not self.ctx.parent_is_deci:
            # set up voronoi calculation for substrate
            builder = VoronoiCalculation.get_builder()
            builder.code = self.inputs.voronoi
            builder.parameters = self.ctx.dsubstrate
            builder.structure = self.ctx.struc_substrate
            builder.metadata.label = 'auxiliary_voronoi_substrate'
            builder.metadata.options = self.ctx.options
            builder.metadata.options['resources'] = {'tot_num_mpiprocs': 1, 'num_machines': 1}
            builder.potential_overwrite = self.ctx.startpot_substrate

            # submit voroaux for substrate calculation
            future_substrate = self.submit(builder)
            self.report('INFO: running voroaux for substrate (pk: {})'.format(future_substrate.pk))

            # set up voronoi calculation for decimation
            builder = VoronoiCalculation.get_builder()
            builder.code = self.inputs.voronoi
            builder.parameters = self.ctx.ddecimation
            builder.structure = self.ctx.struc_decimation
            builder.metadata.label = 'auxiliary_voronoi_decimation'
            builder.metadata.options = self.ctx.options
            builder.metadata.options['resources'] = {'tot_num_mpiprocs': 1, 'num_machines': 1}
            builder.potential_overwrite = self.ctx.startpot_decimation

            # submit voroaux for substrate calculation
            future_decimation = self.submit(builder)
            self.report('INFO: running voroaux for decimation region (pk: {})'.format(future_decimation.pk))

            return ToContext(voroaux_substrate=future_substrate, voroaux_decimation=future_decimation)

        else:
            self.report('Skip voroaux steps due to previous decimation input')

    def run_deciout(self):
        """
        run KKR calculation for deci-out step
        """
        builder = KkrCalculation.get_builder()
        builder.code = self.inputs.kkr
        builder.parameters = self.ctx.dsubstrate
        builder.metadata.options = self.ctx.options
        builder.metadata.options['resources'] = {'tot_num_mpiprocs': 1, 'num_machines': 1}  # force serial run
        builder.metadata.label = 'deci-out'
        builder.parent_folder = self.ctx.voroaux_substrate.outputs.remote_folder
        # create and set initial nonco_angles if needed
        if 'initial_noco_angles' in self.ctx.slab_calc.inputs:
            builder.initial_noco_angles = self.ctx.noco_angles_substrate

        future = self.submit(builder)
        self.report('INFO: running deci-out step (pk: {})'.format(future.pk))

        return ToContext(deciout_calc=future)

    def run_decimation(self):
        """
        run KKR calculation for decimation step
        """
        builder = KkrCalculation.get_builder()
        builder.code = self.inputs.kkr
        builder.parameters = self.ctx.ddecimation
        builder.metadata.options = self.ctx.options
        builder.metadata.label = 'decimation'
        builder.parent_folder = self.ctx.voroaux_decimation.outputs.remote_folder
        builder.deciout_parent = self.ctx.deciout_calc.outputs.remote_folder
        if 'kpoints' in self.inputs:
            builder.kpoints = self.inputs.kpoints
            self.report('INFO: detected kpoints input: run qdos calculation')
        # create and set initial nonco_angles if needed
        if 'initial_noco_angles' in self.ctx.slab_calc.inputs:
            builder.initial_noco_angles = self.ctx.noco_angles_decimation

        future = self.submit(builder)
        self.report('INFO: running deci-out step (pk: {})'.format(future.pk))

        return ToContext(decimation_calc=future)

    def return_results(self):
        """
        check if the calculation was successful and return the result nodes
        """
        # add output nodes
        self.out('structure_decimate', self.ctx.struc_decimation)
        self.out('structure_substrate', self.ctx.struc_substrate)
        self.out('out_params_calc_deci_out', self.ctx.deciout_calc.outputs.output_parameters)
        self.out('out_params_calc_decimate', self.ctx.decimation_calc.outputs.output_parameters)
        self.out('out_remote_calc_decimate', self.ctx.decimation_calc.outputs.remote_folder)
        self.out('out_retrieved_calc_decimate', self.ctx.decimation_calc.outputs.retrieved)

    ###################################################################################################################
    # helper functions

    def _create_structures(self):
        """
        create substrate and slab structures for deci-out and decimation steps
        """
        struc_deci = get_deci_structure(
            Int(self.ctx.nprinc), Int(self.ctx.nplayer), self.ctx.slab_calc.outputs.remote_folder
        )
        struc_substrate = get_substrate_structure(
            Int(self.ctx.nprinc), Int(self.ctx.nplayer), self.ctx.slab_calc.outputs.remote_folder
        )

        self.ctx.struc_decimation = struc_deci
        self.ctx.struc_substrate = struc_substrate

    def _create_startpots(self):
        """
        create substrate and slab potentials from slab calculation
        """
        scf_slab_remote = self.inputs.remote_data

        # create decimation region startpot
        Nlayer_deci = len(self.ctx.struc_decimation.sites)
        new_pot_indices = list(range(Nlayer_deci))

        settings = Dict(
            dict={
                'pot1': 'out_potential',
                'out_pot': 'startpot_decimate',
                'neworder': new_pot_indices,
                'label': 'potential_decimation'
            }
        )
        startpot_deci = neworder_potential_wf(settings_node=settings, parent_calc_folder=scf_slab_remote)

        # create substrate startpot
        nrbasis = self.ctx.struc_decimation.extras['kkr_settings']['NRBASIS']
        Nlayer_deci = len(self.ctx.struc_decimation.sites)
        new_pot_indices = list(range(Nlayer_deci, Nlayer_deci + nrbasis))

        settings = Dict(
            dict={
                'pot1': 'out_potential',
                'out_pot': 'startpot_substrate',
                'neworder': new_pot_indices,
                'label': 'potential_substrate'
            }
        )
        startpot_substrate = neworder_potential_wf(settings_node=settings, parent_calc_folder=scf_slab_remote)

        self.ctx.startpot_substrate = startpot_substrate
        self.ctx.startpot_decimation = startpot_deci


###################################################################################################################


@calcfunction
def get_deci_structure(nprinc, nplayer, slab_calc):
    """
    calcfunction that creates the decimation structure
    """
    nprinc, nplayer = nprinc.value, nplayer.value
    struc, voro_calc = VoronoiCalculation.find_parent_structure(slab_calc)
    voro_params = voro_calc.inputs.parameters.get_dict()
    nrbasis = voro_params.get('<NRBASIS>', voro_params.get('NRBASIS'))
    zperiodr = voro_params.get('ZPERIODR')

    # create decimation structure
    cell = struc.cell
    # adapt third bravias vector
    cell[2] = list(np.array(cell[2]) * (nplayer * nprinc) / len(struc.sites))
    struc_deci = StructureData(cell=cell)
    # add layers that are included in decimation region
    for i in range(nplayer * nprinc):
        struc_deci.append_atom(position=struc.sites[i].position, symbols=struc.sites[i].kind_name)
    # 2D periodic boundary conditions
    struc_deci.pbc = (True, True, False)

    # create settings for KKR lattice continuation (ATTENTION: RBRIGHT etc needs to match substrate calculation)
    kkr_settings = {}
    for k in kkr_decimation_wc._keys2d:
        if k in voro_params:
            kkr_settings[k] = voro_params[k]
    kkr_settings['NRBASIS'] = nrbasis
    kkr_settings['<RBRIGHT>'] = list([list(struc.sites[nplayer * nprinc + i].position) for i in range(nrbasis)])
    if len(kkr_settings['<RBRIGHT>']) == 1:
        kkr_settings['<RBRIGHT>'] = kkr_settings['<RBRIGHT>'][0]
    kkr_settings['NPRINCD'] = nprinc
    # save as extras to decimation structure
    struc_deci.extras['kkr_settings'] = kkr_settings

    return struc_deci


@calcfunction
def get_substrate_structure(nprinc, nplayer, slab_calc):
    """
    calcfunction that creates the decimation structure
    """
    nprinc, nplayer = nprinc.value, nplayer.value
    struc, voro_calc = VoronoiCalculation.find_parent_structure(slab_calc)
    voro_params = voro_calc.inputs.parameters.get_dict()
    nrbasis = voro_params.get('<NRBASIS>', voro_params.get('NRBASIS'))
    zperiodr = voro_params.get('ZPERIODR')

    # make continuation structure (i.e. the substrate)
    cell = struc.cell
    cell[2] = zperiodr
    struc_substrate = StructureData(cell=cell)
    for i in range(nrbasis):
        struc_substrate.append_atom(
            position=struc.sites[nplayer * nprinc + i].position, symbols=struc.sites[nplayer * nprinc + i].kind_name
        )
    struc_substrate.pbc = (True, True, True)

    return struc_substrate


@calcfunction
def make_decimation_param_nodes(slab_calc_params, slab_alat, struc_deci, nkz, params_overwrite=None):
    """
    Create parameter nodes for deci-out and
    """
    # prepare decimation params
    d = {k: v for k, v in slab_calc_params.get_dict().items() if v is not None}

    # make kkr params for substrate calculation (i.e. deci-out mode, needs to run in serial!)
    dsubstrate = d.copy()
    for k in kkr_decimation_wc._keys2d:
        if k in dsubstrate:
            dsubstrate.pop(k)
    dsubstrate = kkrparams(**dsubstrate)
    dsubstrate.set_value('NSTEPS', 1)
    # add deci-out option
    runopts = dsubstrate.get_value('RUNOPT')
    if runopts is None:
        runopts = []
    runopts = [i for i in runopts if i != 'DECIMATE']  # remove decimate flag
    runopts += ['deci-out']
    dsubstrate.set_value('RUNOPT', runopts)
    # increase BZDIVIDE[2] if needed
    bzdiv = dsubstrate.get_value('BZDIVIDE')
    if bzdiv is not None:
        bzdiv[2] = nkz.value
        dsubstrate.set_value('BZDIVIDE', bzdiv)
    # overwrite params from input node
    if params_overwrite is not None:
        for k, v in params_overwrite.get_dict().items():
            dsubstrate.set_value(k, v)
    # use alat from slab calculation
    dsubstrate = {
        k: v for k, v in dsubstrate.get_dict().items() if v is not None
    }  # clean up removing None values from dict
    dsubstrate['ALATBASIS'] = slab_alat.value
    dsubstrate['use_input_alat'] = True
    # now create Dict node with params
    dsubstrate = Dict(dict=dsubstrate)

    # set kkr params for substrate writeout
    ddeci = kkrparams(**d)
    ddeci.set_value('NSTEPS', 1)
    ddeci.set_value('DECIFILES', ['vacuum', 'decifile'])
    ddeci.set_multiple_values(**struc_deci.extras['kkr_settings'])
    # add decimate runopt
    runopts = ddeci.get_value('RUNOPT')
    if runopts is None:
        runopts = []
    runopts = [i for i in runopts if i != 'deci-out']  # remove deci-out flag
    runopts += ['DECIMATE']
    ddeci.set_value('RUNOPT', runopts)
    if params_overwrite is not None:
        for k, v in params_overwrite.get_dict().items():
            ddeci.set_value(k, v)
    # use alat from slab calculation
    ddeci = {k: v for k, v in ddeci.get_dict().items() if v is not None}  # clean up removing None values from dict
    ddeci['ALATBASIS'] = slab_alat.value
    ddeci['use_input_alat'] = True
    # now create Dict node with params
    ddeci = Dict(dict=ddeci)

    return {'dsubstrate': dsubstrate, 'ddeci': ddeci}


@calcfunction
def get_noco_angles_deci(init_noco_slab, struc_decimation, struc_substrate):
    """
    create noco angles for substrate and decimation regions from initial angles of slab calc
    """

    theta_slab = init_noco_slab['theta']
    phi_slab = init_noco_slab['phi']
    fix_dir_slab = init_noco_slab['fix_dir']
    # get number of layers in decimation and substrate regions
    Nlayer_deci = len(struc_decimation.sites)
    nrbasis = len(struc_substrate.sites)

    # set correct noco angles for substrate or decimation region

    # deci-out
    theta = theta_slab[Nlayer_deci:Nlayer_deci + nrbasis]
    phi = phi_slab[Nlayer_deci:Nlayer_deci + nrbasis]
    fix_dir = fix_dir_slab[Nlayer_deci:Nlayer_deci + nrbasis]
    initial_noco_angles_substrate = Dict(dict={'theta': theta, 'phi': phi, 'fix_dir': fix_dir})

    # decimation
    theta = theta_slab[:Nlayer_deci]
    phi = phi_slab[:Nlayer_deci]
    fix_dir = fix_dir_slab[:Nlayer_deci]
    initial_noco_angles_decimation = Dict(dict={'theta': theta, 'phi': phi, 'fix_dir': fix_dir})

    return {
        'initial_noco_angles_substrate': initial_noco_angles_substrate,
        'initial_noco_angles_decimation': initial_noco_angles_decimation
    }
