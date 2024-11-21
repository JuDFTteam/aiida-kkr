#!/usr/bin/env python
# coding: utf-8
"""
This module contains the band structure workflow for KKR which is done by calculating the k-resolved spectral density
also known as Bloch spectral function.
"""

from aiida.orm import Code, Dict, RemoteData, StructureData, Float, Str, WorkChainNode, load_node, CalcJobNode, ArrayData, KpointsData
from aiida.engine import WorkChain, ToContext, calcfunction
from aiida.tools.data.array.kpoints import get_explicit_kpoints_path
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode, update_params_wf, get_inputs_kkr
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.save_output_nodes import create_out_dict_node
from aiida_kkr.tools.extract_kkrhost_noco_angles import extract_noco_angles
from masci_tools.io.kkr_params import kkrparams
from masci_tools.io.common_functions import get_Ry2eV
import numpy as np

__copyright__ = (u'Copyright (c), 2020, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.8'
__contributors__ = (u'Rubel Mozumder', u'Philipp Rüßmann')


class kkr_bs_wc(WorkChain):
    """
    Workchain for BandStructure calculation, starting from RemoteFolderData of the previous converged KKR calculation remote folder data

    inputs:
    :param wf_parameters: (Dict), (optional); Workchain Specifications, contains nepts, tempr, emin (in eV relative to EF), emax (in eV),
                          and RCLUSTZ (can be used to increase the screening cluster radius) keys.

    :param options: (Dict), (optional); Computer Specifications, scheduler command, parallel or serial
    :param kpoints: (KpointsData),(optional); Kpoints data type from the structure,
                                   but not mendatory as it can be extracted from structure internaly from the remote data
    :param remote_data: (RemoteData)(mendaory); From the previous kkr-converged calculation.
    :param kkr: (Code)(mendaory); KKR code specifiaction
    :param label: (Str) (optional) ; label for WC but will be found in the 'result_wf' output
                                     Dict as 'BS_wf_label' key
    :param description: (Str) (optional) : description for WC but will be found in the 'result_wf' output
                                     Dict as 'BS_wf_description' key


    returns:
    :out BS_Data : (ArrayData) ; Consist of BlochSpectralFunction, k_points (list), energy_points (list), special_kpoints(dict)
    :out result_wf: (Dict); work_chain_specifications node, BS_data node, remote_folder node
    """

    _wf_version = __version__
    _wf_label = 'kkr_BandStructure_wc'
    _wf_description = """Workflow for a bandstructure calculation starting eithe from a structure with automatic voronoi'
                        calculation or a valid RemoteData of a previous calculation."""
    _wf_default = {
        'emin': -10.0,  # start of the energy range in eV, relative to the Fermi energy
        'emax': 5.0,  # end of the energy range in eV, relative to the Fermi energy
        'nepts': 96,  # number of energy points
        'RCLUSTZ': None,  # can be used to increase the cluster radius if a value is set here
        'tempr': 50.,  # smearing temperature in K
        'kmesh': None,  # k-point integration mesh, only useful for CPA calculation
    }

    _options_default = {
        'max_wallclock_seconds': 36000,
        'resources': {
            'num_machines': 1
        },
        'withmpi': True,
        'queue_name': '',
        'prepend_text': '',
        'append_text': '',
        'additional_retrieve_list': None
    }

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Return the default values of the workflow parameters (wf_parameters input node)
        """
        if not silent:
            print(f'Version of the kkr_bs_wc workflow: {self._wf_version}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(kkr_bs_wc, cls).define(spec)
        # here inputs are defined
        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._wf_default),
            help='Parameters of the bandstructure workflow (see output of kkr_bs_wc.get_wf_default() for more details).'
        )
        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help='Computer options (walltime etc.) passed onto KkrCalculation'
        )
        spec.input(
            'remote_data',
            valid_type=RemoteData,
            required=True,
            help='Parent folder of previously converged KkrCalculation'
        )
        spec.input('kkr', valid_type=Code, required=True, help='KKRhost code, needed to run the qdos KkrCalculation')
        spec.input(
            'kpoints',
            valid_type=KpointsData,
            required=False,
            help=
            'K-points data for the calculation. If not given the seekpath library is used to find the irreducable k-points of a structure.'
        )
        spec.input('label', valid_type=Str, required=False, help='label for the workflow')
        spec.input('description', valid_type=Str, required=False, help='description for the workflow')
        spec.input(
            'initial_noco_angles',
            valid_type=Dict,
            required=False,
            help="""Initial non-collinear angles for the magnetic moments. See KkrCalculation for details.
            If this is found in the input potentially extracted nonco angles from the parent calulation are overwritten!"""
        )
        # maybe overwrite some settings from the KKRhost convergence run
        spec.input(
            'params_kkr_overwrite',
            valid_type=Dict,
            required=False,
            help='Overwrite some input parameters of the parent KKR calculation.'
        )
        # expose LDAU input node
        spec.input(
            'settings_LDAU',
            valid_type=Dict,
            required=False,
            help="""
Settings for running a LDA+U calculation. The Dict node should be of the form
    settings_LDAU = Dict(dict={'iatom=0':{
        'L': 3,         # l-block which gets U correction (1: p, 2: d, 3: f-electrons)
        'U': 7.,        # U value in eV
        'J': 0.75,      # J value in eV
        'Eref_EF': 0.,  # reference energy in eV relative to the Fermi energy. This is the energy where the projector wavefunctions are calculated (should be close in energy where the states that are shifted lie (e.g. for Eu use the Fermi energy))
    }})
    Note: you can add multiple entries like the one for iatom==0 in this example. The atom index refers to the corresponding atom in the impurity cluster.
"""
        )

        # Here outputs are defined
        spec.output('results_wf', valid_type=Dict, required=True)
        spec.output('BS_Data', valid_type=ArrayData, required=True)

        # Here outlines are being specified
        spec.outline(
            # For initialiging workflow
            cls.start,
            cls.validate_input,
            cls.set_params_BS,
            cls.get_BS,
            cls.return_results
        )
        # definition of exit code in case something goes wrong in this workflow
        spec.exit_code(161, 'ERROR_NO_INPUT_REMOTE_DATA', 'No remote_data was provided as Input')
        spec.exit_code(
            162, 'ERROR_KKRCODE_NOT_CORRECT', 'The code you provided for kkr does not use the plugin kkr.kkr'
        )
        spec.exit_code(
            163, 'ERROR_CALC_PARAMETERS_INVALID',
            'calc_parameters given are not consistent! Hint: did you give an unknown keyword?'
        )
        spec.exit_code(164, 'ERROR_CALC_PARAMETERS_INCOMPLETE', 'calc_parameters not complete')
        spec.exit_code(165, 'ERROR_BS_CALC_FAILED', 'KKR Band Structure calculation failed')
        spec.exit_code(166, 'ERROR_NO_KPOINTS_EXTRACTED', 'No K-POINTS can be extracted from the structure data')
        spec.exit_code(
            167, 'ERROR_INCORRECT_KPOINTS_EXTRACTED',
            'No K-POINTS can be extracted from the primtive structure data rather conventional structure data'
        )
        spec.exit_code(
            168, 'ERROR_INVALID_REMOTE_DATA_TPYE',
            'Input remote_data node neither output of a KKR/voronoi calculation nor of kkr_scf_wc workflow'
        )

    def start(self):
        """
        set up context of the workflow
        """
        self.report(f'INFO: started KKR Band Structure workflow version {self._wf_version}')
        wf_dict = self.inputs.wf_parameters.get_dict()
        # Count energy points only once
        if 'NPT2' in wf_dict.keys():
            npt2 = wf_dict.pop('NPT2', None)
            wf_dict['nepts'] = npt2
        # add missing default values
        for key, val in self._wf_default.items():
            if ((key not in wf_dict.keys()) and (key.swapcase() not in wf_dict.keys()) and (val is not None)):

                self.report(f'INFO: Using default wf parameter {key}: {val}')
                wf_dict[key] = val

        options_dict = self.inputs.options.get_dict()
        if options_dict == {}:
            self.report('INFO: Using default wf Options')
            options_dict = self._options_default
        self.ctx.append_text = options_dict.get('append_text', self._options_default['append_text'])
        self.ctx.prepend_text = options_dict.get('prepend_text', self._options_default['prepend_text'])
        self.ctx.additional_retrieve_list = options_dict.get(
            'additional_retrieve_list', self._options_default['additional_retrieve_list']
        )
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds', self._options_default['max_wallclock_seconds']
        )
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', '')
        self.ctx.BS_params_dict = wf_dict
        self.ctx.BS_kkrparams = None  # is set in set_params_BS
        self.ctx.BS_kpoints = None
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)
        self.report(
            'INFO: use the following parameter:\n'
            'withmpi: {}\n'
            'Resources: {}\n'
            'Walltime (s): {}\n'
            'queue name: {}\n'
            'scheduler command: {}\n'
            'description_wf: {}\n'
            'label_wf: {}\n'
            'BS_params: {}\n'.format(
                self.ctx.withmpi, self.ctx.resources, self.ctx.max_wallclock_seconds, self.ctx.queue,
                self.ctx.custom_scheduler_commands, self.ctx.description_wf, self.ctx.label_wf, self.ctx.BS_params_dict
            )
        )

        self.ctx.successful = True
        self.ctx.errors = []

    def validate_input(self):
        """
        validate input and find out which path ( converged kkr calc or wf ) to take
        return True means run voronoi if false run kkr directly
        """
        inputs = self.inputs
        if 'remote_data' in inputs:
            input_ok = True
        else:
            input_ok = False
            return self.exit_codes.ERROR_NO_INPUT_REMOTE_DATA  # pylint: disable=no-member
        input_remote = self.inputs.remote_data
        parents = input_remote.get_incoming(node_class=CalcJobNode)
        nparents = len(parents.all_link_labels())

        if nparents != 1:
            # extract parent workflow and get uuid of last calc from output node
            parent_workflow = input_remote.inputs.last_RemoteData
            if not isinstance(parent_workflow, WorkChainNode):
                return self.exit_codes.ERROR_INVALID_REMOTE_DATA_TPYE  # pylint: disable=no-member

            parent_workflow_out = parent_workflow.outputs.output_kkr_scf_wc_ParameterResults
            uuid_last_calc = parent_workflow_out.get_dict().get('last_calc_nodeinfo').get('uuid')
            last_calc = load_node(uuid_last_calc)

            if not isinstance(last_calc, KkrCalculation) and not isinstance(last_calc, VoronoiCalculation):
                return self.exit_code.ERROR_INVALID_REMOTE_DATA_TPYE

            # overwrite remote_data node with extracted remote folder
            output_remote = last_calc.outputs.remote_folder

            self.inputs.remote_data = output_remote

        # extract structure
        struc_kkr, _ = VoronoiCalculation.find_parent_structure(self.inputs.remote_data)
        # save if structure is an alloy
        self.ctx.struc_is_alloy = struc_kkr.is_alloy

        # To validate for kpoints
        if 'kpoints' in inputs:
            self.ctx.BS_kpoints = inputs.kpoints
            input_ok = True
            self.ctx.structure_data = 'None (kpoints taken from input)'
        else:
            #create an auxiliary structure with unique kind_names, this leads to using the input structure in the seekpath method instead of finding the primitive one
            cell = np.array(struc_kkr.cell)
            if not struc_kkr.pbc[2]:
                # 2D structure, make sure the third bravais vector points along z
                cell[2] = np.cross(cell[0], cell[1])
            saux = StructureData(cell=cell)
            for isite, site in enumerate(struc_kkr.sites):
                kind = struc_kkr.get_kind(site.kind_name)
                symbols = kind.symbols
                weights = kind.weights
                # replace old way of empty site with new 'X' kind
                if site.kind_name == 'HX' and kind.weights[0] < 1e-8:
                    symbols = 'X'
                    weights = 1.0
                saux.append_atom(
                    name='atom' + str(isite) + ':' + site.kind_name,
                    symbols=symbols,
                    weights=weights,
                    position=site.position
                )
            # use auxiliary structure inside k-point generator
            output = get_explicit_kpoints_path(saux)
            primitive_struc = output['primitive_structure']
            conventional_struc = output['conv_structure']
            kpoints_ok = True

            #check if primitive_structure and input structure are identical:
            maxdiff_cell = sum(abs(np.array(primitive_struc.cell) - np.array(saux.cell))).max()

            if maxdiff_cell > 3 * 10**-9:
                self.report(f'Error in cell : {maxdiff_cell}')
                self.report(
                    'WARNING : The structure data from the voronoi calc is not the primitive structure type and in come cases it is medatory'
                )
                self.report(f'prim: {primitive_struc.cell} {primitive_struc.sites}')
                self.report(f'conv: {conventional_struc.cell} {conventional_struc.sites}')
                self.ctx.structure_data = 'conventional_unit_cell '
            else:
                self.ctx.structure_data = 'primitive_unit_cell'

            if not kpoints_ok:
                return self.exit_codes.ERROR_INCORRECT_KPOINTS_EXTRACTED  # pylint: disable=no-member
            else:
                kpts = output['explicit_kpoints']

            self.ctx.BS_kpoints = kpts
            if isinstance(KpointsData(), type(kpts)):
                input_ok = True
            else:
                input_ok = False
                return self.exit_codes.ERROR_NO_KPOINTS_EXTRACTED  # pylint: disable=no-member

        # To validate for kkr
        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                input_ok = False
                return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT  # pylint: disable=no-member

        # set self.ctx.input_params_KKR
        self.ctx.input_params_KKR = get_parent_paranode(self.inputs.remote_data)
        self.report(f'The validation input_ok {input_ok}')

    def set_params_BS(self):
        """
        set kkr parameters for the bandstructure (i.e. qdos) calculation
        """
        params = self.ctx.input_params_KKR

        # maybe overwrite some inputs
        if 'params_kkr_overwrite' in self.inputs:
            self.report(f'found params_kkr_overwrite: {self.inputs.params_kkr_overwrite.get_dict()}')
            updatenode = self.inputs.params_kkr_overwrite
            updatenode.label = 'params overwrite'
            params = update_params_wf(params, updatenode)

        input_dict = params.get_dict()
        para_check = kkrparams()
        try:
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)
        except:
            return self.exit_codes.ERROR_CALC_PARAMETERS_INVALID  # pylint: disable=no-member
        label = ''
        descr = f'(pk - {self.inputs.remote_data.pk}, and uuid - {self.inputs.remote_data.uuid})'

        missing_list = para_check.get_missing_keys(use_aiida=True)

        if missing_list != []:
            kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
            kkrdefaults_updated = []
            for key_default, val_default in list(kkrdefaults.items()):
                if key_default in missing_list:
                    para_check.set_value(key_default, val_default)
                    kkrdefaults_updated.append(key_default)
                    missing_list.remove(key_default)
            if len(missing_list) > 0:
                self.report(f'ERROR: calc_parameters misses keys: {missing_list}')
                return self.exit_codes.ERROR_CALC_PARAMETERS_INCOMPLETE  # pylint: disable=no-member
            else:
                self.report(f'updated KKR parameter node with default values: {kkrdefaults_updated}')
                label = 'add_defaults_'
                descr = 'added missing default keys, '
        ##+++ Starts to add the NTP2, EMAX and EMIN from the
        econt_new = self.ctx.BS_params_dict
        if self.ctx.struc_is_alloy:
            if econt_new.get('kmesh', None) is None:
                econt_new['kmesh'] = [1, 1, 1]  # overwrite kmesh since the kpoints are used from the input
        kkr_calc = self.inputs.remote_data.get_incoming(node_class=KkrCalculation).first().node
        ef = kkr_calc.outputs.output_parameters.get_dict()['fermi_energy']  # unit in Ry
        self.ctx.fermi_energy = ef  ## in Ry unit

        # Set BS params
        try:
            para_check = set_energy_params(econt_new, ef, para_check)
        except:
            return self.exit_codes.ERROR_CALC_PARAMETERS_INVALID  # pylint: disable=no-member

        para_check.set_multiple_values(
            NPT1=0,
            NPT3=0,
            NPOL=0,
        )

        updatenode = Dict(para_check.get_dict())
        updatenode.label = label + 'KKRparam_BS'
        updatenode.description = 'KKR parameter node extracted from remote_folder' + descr + ' as well as wf_parameter input node.'

        paranode_BS = update_params_wf(self.ctx.input_params_KKR, updatenode)
        self.ctx.BS_kkrparams = paranode_BS

    def get_BS(self):
        """
        submit the KkrCalcultion with the qdos settings for a bandstructure calculation
        """
        label = 'KKR BS calc.'
        BS_dict = self.ctx.BS_params_dict
        key_list = list(BS_dict)
        description = 'User defined BandStructure parameters '
        for key in key_list:
            description += f'{key}= {BS_dict[key]} ,'

        code = self.inputs.kkr
        remote = self.inputs.remote_data
        params = self.ctx.BS_kkrparams
        kpoints = self.ctx.BS_kpoints
        options = {
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'resources': self.ctx.resources,
            'queue_name': self.ctx.queue,
        }
        if self.ctx.custom_scheduler_commands:
            options['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        if self.ctx.append_text:
            options['append_text'] = self.ctx.append_text
        if self.ctx.prepend_text:
            options['prepend_text'] = self.ctx.prepend_text
        if self.ctx.additional_retrieve_list:
            options['additional_retrieve_list'] = self.ctx.additional_retrieve_list

        # get inputs for band structure calculation
        inputs = get_inputs_kkr(
            code, remote, options, label, description, parameters=params, serial=(not self.ctx.withmpi)
        )
        inputs.kpoints = kpoints

        # add nonco angles if found in the parent calculation or in the input
        if 'initial_noco_angles' in self.inputs:
            # overwrite nonco_angles from the input if given
            inputs['initial_noco_angles'] = self.inputs.initial_noco_angles
            self.report('used nonco angles from input to workflow')
        else:
            # extract from the parent calculation
            parent_calc = remote.get_incoming(node_class=KkrCalculation).first().node
            if 'initial_noco_angles' in parent_calc.inputs:
                noco_angles = extract_noco_angles(
                    fix_dir_threshold=Float(1e-6),  # make small enough
                    old_noco_angles=parent_calc.inputs.initial_noco_angles,
                    last_retrieved=parent_calc.outputs.retrieved
                )
                # set nonco angles (either from input or from output if it was updated)
                if noco_angles == {}:
                    noco_angles = parent_calc.inputs.initial_noco_angles
                self.report(f'extract nonco angles and use from parent ({noco_angles})')

        # LDA+U settings
        if 'settings_LDAU' in self.inputs:
            self.report('Add settings_LDAU input node')
            inputs.settings_LDAU = self.inputs.settings_LDAU

        BS_run = self.submit(KkrCalculation, **inputs)
        self.ctx.last_calc = BS_run

        return ToContext(BS_run=BS_run)

    def return_results(self):
        """
        Collect results, parse BS_calc output and link output nodes to workflow node
        """

        caching_info = f'INFO: cache_source of BS calc node: {self.ctx.BS_run.get_cache_source}'
        self.report(caching_info)

        if not self.ctx.BS_run.is_finished_ok:
            self.ctx.successful = False
            error = f'ERROR BS calculation failed somehow it is in state {self.ctx.BS_run.process_state}'
            self.report(error)
            self.ctx.errors.append(error)
            return self.exit_codes.ERROR_BS_CALC_FAILED  # pylint: disable=no-member

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._wf_version
        outputnode_dict['withmpi'] = self.ctx.withmpi
        outputnode_dict['resources'] = self.ctx.resources
        outputnode_dict['max_wallclock_seconds'] = self.ctx.max_wallclock_seconds
        outputnode_dict['queue_name'] = self.ctx.queue
        outputnode_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        outputnode_dict['BS_params'] = self.ctx.BS_params_dict
        if 'kpoints' not in self.inputs:
            outputnode_dict['structure_type'] = self.ctx.structure_data
        outputnode_dict['BS_wf_description'] = self.ctx.description_wf
        outputnode_dict['BS_wf_label'] = self.ctx.label_wf
        try:
            outputnode_dict['nspin'] = self.ctx.BS_run.res.nspin
        except:
            error = 'ERROR: nspin not extracted'
            self.report(error)
            self.ctx.successful = False
            self.ctx.errors.append(error)
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['list_of_errors'] = self.ctx.errors

        # create output node with data-provenance
        outputnode = Dict(outputnode_dict)
        outputnode.label = 'kkr_BS_wc_results'
        outputnode.description = 'Contains the info of the WC'

        self.report('INFO: create Banstructure results nodes')
        try:
            self.report(
                f'INFO: create Bandstructure results nodes. BS calc retrieved node={self.ctx.BS_run.outputs.retrieved}'
            )
            has_BS_run = True
        except AttributeError as e:
            self.report('ERROR: No Bandstructure calc retrieved node found')
            self.report(f'Caught AttributeError {e}')
            return self.exit_codes.ERROR_BS_CALC_FAILED  # pylint: disable=no-member

        if has_BS_run:
            BS_retrieved = self.ctx.BS_run.outputs.retrieved

        ef = self.ctx.fermi_energy  # in Ry unit
        kpoints = self.ctx.BS_kpoints

        # Here outdict dictionary has been created to set the Dict result_wf, BS_data
        # to the output(spec.output) of the wf
        outdict = {}
        if has_BS_run:
            ArraData = parse_BS_data(BS_retrieved, Float(ef), kpoints)
            outdict['BS_Data'] = ArraData['BS_Data']

        # link to the BS output nodes
        link_nodes = outdict.copy()

        outdict['results_wf'] = create_out_dict_node(outputnode, **link_nodes)

        # create links to output nodes
        for link_name, node in outdict.items():
            self.out(link_name, node)

        self.report('INFO: done with BS_workflow!\n')


def set_energy_params(econt_new, ef, para_check):
    """
    set energy contour values to para_check
    internally convert from relative eV units to absolute Ry units
    """
    evscal = get_Ry2eV()

    for key, val in econt_new.items():
        if key in ['kmesh', 'BZDIVIDE', 'KMESH', 'bzdivide']:
            key = 'BZDIVIDE'
        elif key in ['nepts', 'NPT2']:
            key = 'NPT2'
            # also add IEMXD which has to be big enough
            para_check.set_value('IEMXD', val, silent=True)
        elif key in ['emin', 'EMIN']:
            key = 'EMIN'
            val = (ef + val / evscal)  # converting the Energy value to Ry while the fermi_energy in Ry
        elif key in ['emax', 'EMAX']:
            key = 'EMAX'
            val = (ef + val / evscal)  # Converting to the Ry (unit of the energy)
        elif key in ['tempr', 'TEMPR']:
            key = 'TEMPR'
        elif key in ['RCLUSTZ', 'rclustz']:
            key = 'RCLUSTZ'
        para_check.set_value(key, val, silent=True)

    # set the rest of the DOS contour
    para_check.set_multiple_values(
        NPT1=0,
        NPT3=0,
        NPOL=0,
        use_semi_circle_contour=False,  # this is needed to get a DOS contour
    )

    # set KPOIBZ to match BZDIVIDE setting
    # this is only done if 'KPOIBZ' is not given already in the input
    bzdiv = para_check.get_value('BZDIVIDE')
    if bzdiv is not None and 'KPOIBZ' not in econt_new:
        para_check.set_value('KPOIBZ', np.prod(bzdiv), silent=True)
    # we need to make sure to deactivate the semi-circle contour, otherwise DOS contour is not used
    para_check.set_value('<USE_SEMI_CIRCLE_CONTOUR>', False, silent=True)

    return para_check


@calcfunction
def parse_BS_data(retrieved_folder, fermi_level, kpoints):
    """
    parse the qdos files from the retreived folderand save as ArrayData
    """
    # conversion factor from Ry to eV
    eVscale = get_Ry2eV()

    retrieved_list = retrieved_folder.list_object_names()
    qdos_file_list = [i for i in retrieved_list if 'qdos.' in i]
    q_vec_file = 'qvec.dat'

    with retrieved_folder.open(q_vec_file) as file_opened:
        q_vec = np.loadtxt(file_opened, skiprows=1)

    for icount, fname in enumerate(qdos_file_list):
        with retrieved_folder.open(fname) as _f:
            loaded_file = np.loadtxt(_f)
            if icount == 0:
                total_qdos = loaded_file
            else:
                total_qdos[:, 5:] += loaded_file[:, 5:]

    ef = fermi_level.value  # in Ry unit
    total_qdos[:, 0] = (total_qdos[:, 0] - ef) * eVscale
    eng_points = set(total_qdos[:, 0])
    eng_points = np.sort(list(eng_points))
    no_eng_points = len(eng_points)

    qdos_intensity = np.ndarray(shape=(no_eng_points, len(q_vec)))
    for ne in range(np.shape(qdos_intensity)[0]):
        nk = np.shape(qdos_intensity)[1]
        # sum up all l-channels (5 is only the s-channel!)
        qdos_intensity[ne, :] = np.sum(total_qdos[ne * nk:(ne + 1) * nk, 5:], axis=1) / eVscale

    qdos_intensity = qdos_intensity.T  # setting eng-kpts corresponds to x-y asix
    q_vec = np.asarray(q_vec)  # converting q_vec into array
    eng_points = (np.asarray(eng_points))  # converting eng_popints into array in Ry unit

    # To save into the ArrayData
    array = ArrayData()
    array.set_array('BlochSpectralFunction', qdos_intensity)
    array.set_array('Kpts', q_vec)
    array.set_array('energy_points', eng_points)
    if kpoints.labels is not None:
        klbl_dict = dict(kpoints.labels)  # Special k-points
        array.extras['k-labels'] = klbl_dict

    return {'BS_Data': array}
