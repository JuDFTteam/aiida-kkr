#!/usr/bin/env python
# coding: utf-8
"""
This module contains the workflow that can be used to calculate the exchange coupling constants
"""

import numpy as np
from aiida.engine import WorkChain, ToContext, calcfunction
from aiida.orm import Dict, RemoteData, StructureData, ArrayData, CalcJobNode, Code
from aiida_kkr.tools.find_parent import get_calc_from_remote, find_parent_structure
from aiida_kkr.tools.jij_tools import parse_jij_calc, get_sites
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode, update_params_wf, get_inputs_kkr
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.tools.save_output_nodes import create_out_dict_node
from masci_tools.io.common_functions import get_alat_from_bravais
from masci_tools.io.kkr_params import kkrparams
from masci_tools.util.constants import BOHR_A

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.2'
__contributors__ = (u'Philipp Rüßmann')


class kkr_jij_wc(WorkChain):
    """
    Workchain for calculation of exchange coupling constants Jij and Dij if parent calculation used the SOC solver.

    inputs::

        :param wf_parameters: optional Dict node of workchain specifications, contains settings like Jij radius cutoff,
                              selection of sites for i and j and numerical cutoffs. None values in the accuracy sub-dict
                              means that values from parent calculation are coptied.
        :param remote_data: mandatory RemoteData node of parent (i.e. converged) KkrCalculation
        :param kkr: optional Code for KKRhost executable (if not given the same as in the parent calculation is used)
        :param options: optional Dict computer options like scheduler command or parallelization

    returns::

        :return jij_data: ArrayData with the arrays 'Jij_expanded' (Table of all Jij and Dij pairs) and 'positions_expanded' (positions of all ij pairs)
        :return structure_jij_sites: StructureData
    """

    _wf_version = __version__
    _wf_default = {
        'jijrad_ang': 5.0, # default cutoff radius
        'jijsite_i': None, # use all sites by default
        'jijsite_j': None, # use all sites by default
        'accuracy': { # accuracy settings, typically set to larger values than in scf run
            'NATOMIMPD': 500,
            'NSHELD': 2000,
            'TEMPR': None,
            'RCLUSTZ': None,
            'kmesh': None,
        },
    }
    _options_default = {
        'max_wallclock_seconds': 36000,
        'resources': {
            'num_machines': 1
        },
        'withmpi': True,
        'queue_name': ''
    }

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Return the default values of the workflow parameters (wf_parameters input node)
        """
        if not silent:
            print(f'Version of the kkr_jij_wc workflow: {self._wf_version}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(kkr_jij_wc, cls).define(spec)

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
            help=
            'Computer options (walltime etc.) passed onto KkrCalculation, fall back to settings from parent calculation if not given'
        )

        spec.input(
            'remote_data',
            valid_type=RemoteData,
            required=True,
            help='Parent folder of previously converged KkrCalculation'
        )

        spec.input('kkr', valid_type=Code, required=True, help='KKRhost code, needed to run the Jij KkrCalculation')

        spec.input(
            'params_kkr_overwrite',
            valid_type=Dict,
            required=False,
            help='Overwrite some input parameters of the parent KKR calculation.'
        )

        # Here outputs are defined
        spec.output('results_wf', valid_type=Dict, required=True)
        spec.output('jij_data', valid_type=ArrayData, required=True)
        spec.output('structure_jij_sites', valid_type=StructureData, required=True)

        # Here outlines are being specified
        spec.outline(
            # For initialiging workflow
            cls.start,
            cls.validate_input,
            cls.set_jij_params,
            cls.submit_Jij_calcs,
            cls.return_results
        )
        # definition of exit code in case something goes wrong in this workflow
        spec.exit_code(
            160, 'ERROR_KKRCODE_NOT_CORRECT', 'The code you provided for kkr does not use the plugin kkr.kkr'
        )
        spec.exit_code(161, 'ERROR_INVALID_PARENT', 'Parent calculation is not valid')
        spec.exit_code(162, 'ERROR_CALC_FAILED', 'KKR Band Structure calculation failed')
        spec.exit_code(163, 'ERROR_PARSING_FAILED', 'Parsing of Jij calculations failed')

    def start(self):
        """
        set up context of the workflow
        """
        self.report(f'INFO: started KKR Jij workflow version {self._wf_version}')
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
        else:
            wf_dict = {}
        # add missing default valuesi
        for key, val in self._wf_default.items():
            if ((key not in wf_dict.keys()) and (key.swapcase() not in wf_dict.keys()) and (val is not None)):
                self.report(f'INFO: Using default wf parameter {key}: {val}')
                wf_dict[key] = val

        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
        else:
            options_dict = self._options_default
        self.ctx.options = options_dict

        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds', self._options_default['max_wallclock_seconds']
        )
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', '')
        self.ctx.wf_dict = wf_dict
        self.report(
            'INFO: use the following parameter:\n'
            'withmpi: {}\n'
            'Resources: {}\n'
            'Walltime (s): {}\n'
            'queue name: {}\n'
            'scheduler command: {}\n'
            'Workflow parameters: {}\n'.format(
                self.ctx.withmpi, self.ctx.resources, self.ctx.max_wallclock_seconds, self.ctx.queue,
                self.ctx.custom_scheduler_commands, wf_dict
            )
        )

    def validate_input(self):
        """
        validate inputs
        """

        # save parent calculation
        input_remote = self.inputs.remote_data
        parents = input_remote.get_incoming(node_class=CalcJobNode).all()
        if len(parents) != 1:
            # check if parent is unique
            return self.exit_codes.ERROR_INVALID_PARENT  # pylint: disable=no-member
        self.ctx.parent_calc = get_calc_from_remote(input_remote)

        # validate for kkrcode
        try:
            test_and_get_codenode(self.inputs.kkr, 'kkr.kkr', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT  # pylint: disable=no-member

        # get parent input parameter
        self.ctx.input_params_KKR = get_parent_paranode(self.inputs.remote_data)

    def set_jij_params(self):
        """
        set kkr parameters for the Jij calculation
        """

        # input parameters from parent
        params = self.ctx.input_params_KKR

        # maybe overwrite some inputs
        if 'params_kkr_overwrite' in self.inputs:
            self.report(f'found params_kkr_overwrite: {self.inputs.params_kkr_overwrite.get_dict()}')
            updatenode = self.inputs.params_kkr_overwrite
            updatenode.label = 'params overwrite'
            params = update_params_wf(params, updatenode)

        # set Jij parameters
        para_jij, runopts = self._get_para_jij(params)
        updatenode = Dict(para_jij.get_dict())
        updatenode.label = 'Jij params'
        paranode_jij = update_params_wf(params, updatenode)
        self.ctx.jij_params = paranode_jij

        # find out if we have a calculation with or without SOC (then no DMI is calculated)
        self.ctx.noSOC = True
        if 'NEWSOSOL' in runopts or para_jij.get_value('<USE_CHEBYCHEV_SOLVER>'):
            self.ctx.noSOC = False

    def submit_Jij_calcs(self):
        """
        submit the KkrCalcultion with the Jij settings
        """

        # get inputs for band structure calculation
        inputs = get_inputs_kkr(
            self.inputs.kkr,
            self.inputs.remote_data,
            self.ctx.options,
            label='Jij_calc',
            description='',
            parameters=self.ctx.jij_params,
            serial=(not self.ctx.withmpi)
        )

        # noSOC, only m||z
        if self.ctx.noSOC:
            inputs.metadata.label = 'jij_calc_z'  # pylint: disable=no-member
            jij_calc_z = self.submit(KkrCalculation, **inputs)
            self.ctx.jij_calc_z = jij_calc_z
            self.ctx.jij_calc_x = None
            self.ctx.jij_calc_y = None
        else:
            # create nonco angles for the three calculations
            nonco_angles = _make_nonco_angles(parent_remote=self.inputs.remote_data)
            init_angles_x, init_angles_y, init_angles_z = nonco_angles['init_angles_x'], nonco_angles[
                'init_angles_y'], nonco_angles['init_angles_z']

            # submit m||z calculation
            inputs.initial_noco_angles = init_angles_z
            inputs.metadata.label = 'jij_calc_z'  # pylint: disable=no-member
            jij_calc_z = self.submit(KkrCalculation, **inputs)
            self.ctx.jij_calc_z = jij_calc_z

            # submit m||x calculation
            inputs.initial_noco_angles = init_angles_x
            inputs.metadata.label = 'jij_calc_x'  # pylint: disable=no-member
            jij_calc_x = self.submit(KkrCalculation, **inputs)
            self.ctx.jij_calc_x = jij_calc_x

            # submit m||y calculation
            inputs.initial_noco_angles = init_angles_y
            inputs.metadata.label = 'jij_calc_y'  # pylint: disable=no-member
            jij_calc_y = self.submit(KkrCalculation, **inputs)
            self.ctx.jij_calc_y = jij_calc_y

        # add to context (needed to tell aiida to wait for processes to finish)
        futures = {'jij_calc_z': self.ctx.jij_calc_z}
        if self.ctx.jij_calc_x is not None:
            futures['jij_calc_x'] = self.ctx.jij_calc_x
        if self.ctx.jij_calc_y is not None:
            futures['jij_calc_y'] = self.ctx.jij_calc_y
        return ToContext(**futures)

    def return_results(self):
        """
        Collect results, parse Jij output and link output nodes to workflow node
        """

        # check if calculations finished ok
        success = True
        if self.ctx.jij_calc_z is not None and not self.ctx.jij_calc_z.is_finished_ok:
            success = False
        if self.ctx.jij_calc_x is not None and not self.ctx.jij_calc_x.is_finished_ok:
            success = False
        if self.ctx.jij_calc_y is not None and not self.ctx.jij_calc_y.is_finished_ok:
            success = False
        if not success:
            self.ctx.successful = False
            error = f'ERROR Jij calculation failed somehow it is in state {self.ctx.jij_calc_z.process_state}'
            if self.ctx.jij_calc_x is not None:
                error += f'; {self.ctx.jij_calc_x.process_state} (x)'
            if self.ctx.jij_calc_y is not None:
                error += f'; {self.ctx.jij_calc_y.process_state} (y)'
            self.report(error)
            return self.exit_codes.ERROR_CALC_FAILED  # pylint: disable=no-member

        # now parse calculation output
        try:
            jij_data, structure_jij_sites = parse_jij_calc(
                self.ctx.jij_calc_z, jij_calc_x=self.ctx.jij_calc_x, jij_calc_y=self.ctx.jij_calc_y, verbose=False
            )
        except Exception as err:
            self.report(f'Caught error when trying to parse Jij output:{err}')
            return self.exit_codes.ERROR_PARSING_FAILED  # pylint: disable=no-member

        # collect output nodes
        outdict = {'jij_data': jij_data, 'structure_jij_sites': structure_jij_sites}

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._wf_version
        outputnode_dict['successful'] = success

        # create output node with data-provenance
        outputnode = Dict(outputnode_dict)

        # link to the output nodes
        link_nodes = outdict.copy()

        outdict['results_wf'] = create_out_dict_node(outputnode, **link_nodes)

        # create links to output nodes
        for link_name, node in outdict.items():
            self.out(link_name, node)

        self.report('INFO: done with Jij workflow!')

    # Helper functions

    def _get_para_jij(self, params):
        """
        Set the Jij parameters from the input.
        Returns a kkrparams instance with the set values
        """

        # get input parameters
        input_dict = params.get_dict()
        para_jij = kkrparams(**input_dict)

        # set Jij parameters
        # add 'XCPL' runopt to list of runopts (activates Jij calculation)
        runopts = input_dict.get('RUNOPT')
        if runopts is None:
            runopts = []
        runopts.append('XCPL    ')
        para_jij.set_value('RUNOPT', runopts)
        para_jij.set_value('NSTEPS', 1)  # one-shot run

        # accuracy settings
        tempr = self.ctx.wf_dict.get('accuracy', {}).get('TEMPR')
        if tempr is not None:
            para_jij.set_value('TEMPR', tempr)  # slightly reduce temperature

        rclustz = self.ctx.wf_dict.get('accuracy', {}).get('RCLUSTZ')
        if rclustz is not None:
            para_jij.set_value('RCLUSTZ', rclustz)  # increase cluster radius

        kmesh = self.ctx.wf_dict.get('accuracy', {}).get('kmesh')
        if kmesh is not None:
            para_jij.set_value('BZDIVIDE', kmesh)  # increase k-points
            para_jij.set_value('KPOIBZ', np.product(kmesh))  # array dimension

        # array dimensions
        NATOMIMPD = self.ctx.wf_dict.get('accuracy', {}).get('NATOMIMPD')
        if NATOMIMPD is not None:
            para_jij.set_value('NATOMIMPD', NATOMIMPD)  # array dimension

        NSHELD = self.ctx.wf_dict.get('accuracy', {}).get('NSHELD')
        if NSHELD is not None:
            para_jij.set_value('NSHELD', NSHELD)  # array dimension

        # Jij settings
        jijrad = self._get_jijrad()
        if jijrad is not None:
            para_jij.set_value('JIJRAD', jijrad)  # radius in lattice constants up to which the Jijs are calculated

        # set optional Jij parameters
        # i and j index for Jij calculation in internal units
        # uses site index (i.e. needs to be <=10)
        JIJSITEI = self.ctx.wf_dict.get('jijsite_i')
        JIJSITEJ = self.ctx.wf_dict.get('jijsite_j')
        if JIJSITEI is not None:
            if JIJSITEJ is None:
                JIJSITEJ = JIJSITEI
            para_jij.set_value('JIJSITEI', [len(JIJSITEI)] + [i + 1 for i in JIJSITEI])
        if JIJSITEJ is not None:
            para_jij.set_value('JIJSITEJ', [len(JIJSITEJ)] + [i + 1 for i in JIJSITEJ])

        return para_jij, runopts

    def _get_jijrad(self):
        """
        get Jij radius convert from Ang to internal alat units
        """

        # Jij radius in Ang.
        jijrad_ang = self.ctx.wf_dict.get('jijrad_ang')

        # find structure from calculation
        struc, _ = find_parent_structure(self.ctx.parent_calc)
        self.ctx.structure = struc

        # get alat from structure
        alat_ang = get_alat_from_bravais(np.array(struc.cell), struc.pbc[2])

        # maybe use value provided in input instead
        para = {k.lower(): v for k, v in self.ctx.parent_calc.inputs.parameters.get_dict().items() if v is not None}
        if para.get('use_alat_input', False) or para.get('use_input_alat', False):
            alat_ang = para.get('alatbasis') * BOHR_A

        # now have Jij radius in alat units
        jijrad = jijrad_ang / alat_ang

        return jijrad


@calcfunction
def _make_nonco_angles(parent_remote):
    """
    Create nonco angles for the 3 directions (x y, z)
    """
    # find structure to count number of sites
    structure = find_parent_structure(parent_remote)[0]
    Nsites = len(get_sites(structure))

    # create nonco angles for m||z
    init_angles_z = Dict({
        'fix_dir': [True for i in range(Nsites)],
        'theta': [0.0 for i in range(Nsites)],
        'phi': [0.0 for i in range(Nsites)],
    })

    # create nonco angles for m||x
    init_angles_x = Dict({
        'fix_dir': [True for i in range(Nsites)],
        'theta': [90.0 for i in range(Nsites)],
        'phi': [0.0 for i in range(Nsites)],
    })

    # create nonco angles for m||y
    init_angles_y = Dict({
        'fix_dir': [True for i in range(Nsites)],
        'theta': [90.0 for i in range(Nsites)],
        'phi': [90.0 for i in range(Nsites)],
    })

    return {'init_angles_x': init_angles_x, 'init_angles_y': init_angles_y, 'init_angles_z': init_angles_z}
