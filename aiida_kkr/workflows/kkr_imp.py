# -*- coding: utf-8 -*-
"""
In this module you find the total workflow for a kkr impurity calculation
and some helper methods to do so with AiiDA
"""

from aiida.orm import Code, load_node, RemoteData, StructureData, Dict, SinglefileData, FolderData
from aiida.orm import CalcJobNode
from aiida.engine import WorkChain, ToContext, if_
from aiida.engine import calcfunction
from aiida_kkr.calculations.voro import VoronoiCalculation
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.tools import test_and_get_codenode, neworder_potential_wf, update_params_wf
from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
import numpy as np
from aiida_kkr.tools.save_output_nodes import create_out_dict_node

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.9.3'
__contributors__ = (u'Fabian Bertoldo', u'Philipp Rüßmann')
#TODO: generalize workflow to multiple impurities
#TODO: add additional checks for the input
#TODO: maybe work on a clearer outputnode structure


class kkr_imp_wc(WorkChain):
    """
    Workchain of a kkrimp calculation starting either from scratch (with a structure
    and impurity_info node), or with a converged host potential and impurity
    startpotentials, ... to calculate the converged host-impurity potential of the system.

    :param options: (Dict), Workchain specifications
    :param wf_parameters: (Dict), specifications for the kkr impurity workflow
    :param voro_aux_parameters: (Dict), specification for the auxiliary voronoi calculation for the impurity
    :param kkrimp: (Code), mandatory: KKRimp code converging the host-imp-potential
    :param kkr: (Code), mandatory: KKR code for calculation the host potential
    :param voronoi: (Code), mandatory: Voronoi code to generate the impurity startpot
    :param remote_data_gf: (RemoteData): remote folder of a previous kkrflex
                                         calculation containing the flexfiles ...
    :param remote_data_host: (RemoteData): remote folder of a converged KKR
                                           host calculation

    :return workflow_info: (Dict), Information of workflow results
    :return last_calc_output_parameters: (Dict), output parameters of
                                         the last called calculation
    :return last_calc_info: (Dict), information of the last called calculation
    """

    _workflowversion = __version__
    _wf_label = 'kkr_imp_wc'
    _wf_description = 'Workflow for a KKRimp calculation'

    _options_default = {
        'queue_name': '',  # Queue name to submit jobs too
        'resources': {
            'num_machines': 1  # resources to allowcate for the job
        },
        'max_wallclock_seconds': 60 * 60,  # walltime after which the job gets killed (gets parsed to KKR)}
        'custom_scheduler_commands': '',  # some additional scheduler commands
        'withmpi': True  # execute KKR with mpi or without
    }

    # settings for sub workflow (impurity convergence)
    _wf_default = kkr_imp_sub_wc.get_wf_defaults(silent=True)
    # add control to retrieve kkrflex files to repository or leave on remote computer only
    _wf_default['retrieve_kkrflex'] = True

    # settings for vorostart workflow, used to generate starting potential
    _voro_aux_default = kkr_startpot_wc.get_wf_defaults(silent=True)

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create
        set of wf_parameters.

        returns _wf_defaults
        """
        if not silent:
            print(f'Version of workflow: {self._workflowversion}')
        return self._options_default.copy(), self._wf_default.copy(), self._voro_aux_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """

        super(kkr_imp_wc, cls).define(spec)

        # expose inputs of sub workflow
        # TODO also expose the other inputs in next release, but put deprecation warnings in first
        spec.expose_inputs(
            kkr_imp_sub_wc,
            namespace='scf',
            include=(
                # 'kkrimp',
                'options',
                # 'wf_parameters',
                'params_overwrite',
                'initial_noco_angles',
                'rimpshift',
                'settings_LDAU',
            )
        )

        # define the inputs of the workflow
        spec.input('kkr', valid_type=Code, required=False, help='KKRhost code used to run GF writeout step.')
        spec.input(
            'voronoi',
            valid_type=Code,
            required=True,
            help='Voronoi code used to create the impurity starting potential.'
        )
        spec.input(
            'kkrimp', valid_type=Code, required=True, help='KKRimp code used to converge the impurity calculation'
        )
        spec.input(
            'impurity_info',
            valid_type=Dict,
            required=True,
            help='Information of the impurity like position in the unit cell, screening cluster, atom type.'
        )
        spec.input(
            'remote_data_host',
            valid_type=RemoteData,
            required=False,
            help='RemoteData node of the converged host calculation. Used to write out the host Green function.'
        )
        spec.input(
            'remote_data_gf',
            valid_type=RemoteData,
            required=False,
            help='RemoteData node of precomputed host Green function.'
        )
        spec.input(
            'remote_data_gf_Efshift',
            valid_type=RemoteData,
            required=False,
            help=
            'RemoteData node of precomputed host Green function with Fermi level shift (overwrite kkrflex_green and tmat files from first remote_data_gf node.'
        )
        spec.input('options', valid_type=Dict, required=False, help='Options for running the codes (walltime etc.).')
        spec.input(
            'options_voronoi',
            valid_type=Dict,
            required=False,
            help='Options for running the Voronoi code (if differing from general `options` node)'
        )
        spec.input(
            'voro_aux_parameters',
            valid_type=Dict,
            required=False,
            help='Parameters for the auxiliary voronoi starting potential workflow.'
        )
        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            help='Parameters for the KKRimp selfconsistency workflow.'
        )
        spec.input(
            'voro_params_overwrite',
            valid_type=Dict,
            required=False,
            help=
            'If given, overwrite the some parameters used as input for auxiliary voronoi calculation of starting potential.'
        )
        spec.input(
            'params_kkr_overwrite',
            valid_type=Dict,
            required=False,
            help='Set some input parameters of the KKR calculation for the GF writeout step.'
        )
        spec.input(
            'startpot',
            valid_type=SinglefileData,
            required=False,
            help='Set starting potential (e.g. from preconverged calculation'
        )
        spec.expose_inputs(kkr_flex_wc, namespace='gf_writeout', include=('params_kkr_overwrite', 'options', 'kkr'))

        # structure of the workflow
        spec.outline(
            cls.start,                                                          # initialize workflow
            if_(cls.validate_input)(                                            # validate the input (if true, run_gf_writeout, else skip)
                cls.run_gf_writeout),                                           # write out the host GF
            if_(cls.has_starting_potential_input)(                              # check if strarting potential exists in input already (otherwise create it)
                cls.run_voroaux,                                                  # calculate the auxiliary impurity potentials
                cls.construct_startpot),                                          # construct the host-impurity startpotential
            cls.run_kkrimp_scf,                                                 # run the kkrimp_sub workflow to converge the host-imp startpot
            cls.return_results,                                                 # check if the calculation was successful and return the result nodes
            cls.error_handler)

        # define the possible exit codes
        spec.exit_code(
            141,
            'ERROR_INVALID_INPUT_CODE',
            message='ERROR: one or more of the codes you provided do not '
            'use the necessary plugins: kkr.voro, kkr.kkr, kkr.kkrimp'
        )
        spec.exit_code(
            142,
            'ERROR_MISSING_KKRCODE',
            message='ERROR: since GF writeout step has to be conducted, '
            "'kkrcode' is needed as an input"
        )
        spec.exit_code(
            143,
            'ERROR_MISSING_REMOTE',
            message='ERROR: neither converged host remote nor GF writeout '
            'remote is given as an input. One of them is needed to '
            'proceed with this workflow!'
        )
        spec.exit_code(
            144, 'ERROR_KKRIMP_SUB_WORKFLOW_FAILURE', message='ERROR: sub-workflow for KKRimp convergence failed'
        )
        spec.exit_code(
            145,
            'ERROR_KKRSTARTPOT_WORKFLOW_FAILURE',
            message='ERROR: sub-workflow Kkr_startpot failed (look for failure of voronoi calculation).'
        )

        # define the outputs of the workflow
        spec.output('workflow_info', valid_type=Dict)
        spec.output('last_calc_output_parameters', valid_type=Dict)
        spec.output('last_calc_info', valid_type=Dict)
        spec.output('converged_potential', valid_type=SinglefileData, required=False)
        spec.output('remote_data_gf', valid_type=RemoteData)

    def start(self):
        """
        Init context and some parameters
        """

        self.report(f'INFO: started KKR impurity workflow version {self._workflowversion}')

        # get input parameters
        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
        else:
            options_dict = self._options_default
            self.report('INFO: using default options')
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
        else:
            wf_dict = self._wf_default
            self.report('INFO: using default workflow parameters for KKRimp scf cycle')
        if 'voro_aux_parameters' in self.inputs:
            voro_aux_dict = self.inputs.voro_aux_parameters.get_dict()
        else:
            voro_aux_dict = self._voro_aux_default
            self.report('INFO: using default workflow parameters for auxiliary voronoi calculation')

        # get parameters that will be added/overwritten for voronoi run
        if 'voro_params_overwrite' in self.inputs:
            self.ctx.change_voro_params = self.inputs.voro_params_overwrite.get_dict()
        else:
            self.ctx.change_voro_params = {}

        # set option parameters from input, or defaults
        self.ctx.exit_code = None  # collect errors here which are passed on to the error handler in the end
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds', self._options_default['max_wallclock_seconds']
        )
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get(
            'custom_scheduler_commands', self._options_default['custom_scheduler_commands']
        )
        self.ctx.options_params_dict = Dict({
            'withmpi': self.ctx.withmpi,
            'resources': self.ctx.resources,
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'queue_name': self.ctx.queue,
            'custom_scheduler_commands': self.ctx.custom_scheduler_commands
        })
        if 'options_voronoi' in self.inputs:
            self.ctx.options_params_dict_voronoi = self.inputs.options_voronoi.get_dict()
            self.report(f'INFO: Use different options for voronoi code ({self.ctx.options_params_dict_voronoi})')
        else:
            self.ctx.options_params_dict_voronoi = self.ctx.options_params_dict.get_dict()

        # set label and description of the workflow
        self.ctx.description_wf = self.inputs.get(
            'description', 'Workflow for a KKR impurity calculation starting from a host-impurity potential'
        )
        self.ctx.label_wf = self.inputs.get('label', 'kkr_imp_sub_wc')

        # set parameters for the auxiliary voronoi calculation
        self.ctx.voro_dos_params = voro_aux_dict.get('dos_params', self._voro_aux_default['dos_params'])
        self.ctx.voro_num_rerun = voro_aux_dict.get('num_rerun', self._voro_aux_default['num_rerun'])
        self.ctx.voro_fac_cls_increase = voro_aux_dict.get(
            'fac_cls_increase', self._voro_aux_default['fac_cls_increase']
        )
        self.ctx.voro_natom_in_cls_min = voro_aux_dict.get(
            'natom_in_cls_min', self._voro_aux_default['natom_in_cls_min']
        )
        self.ctx.voro_delta_e_min = voro_aux_dict.get('delta_e_min', self._voro_aux_default['delta_e_min'])
        self.ctx.voro_threshold_dos_zero = voro_aux_dict.get(
            'threshold_dos_zero', self._voro_aux_default['threshold_dos_zero']
        )
        self.ctx.voro_check_dos = voro_aux_dict.get('check_dos', self._voro_aux_default['check_dos'])
        #self.ctx.voro_delta_e_min_core_states = voro_aux_dict.get(
        #    'delta_e_min_core_states', self._voro_aux_default['delta_e_min_core_states']
        #)
        # set up new parameter dict to pass to voronoi subworkflow later
        self.ctx.voro_params_dict = Dict({
            'queue_name': self.ctx.queue,
            'resources': self.ctx.resources,
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'withmpi': self.ctx.withmpi,
            'custom_scheduler_commands': self.ctx.custom_scheduler_commands,
            'dos_params': self.ctx.voro_dos_params,
            'num_rerun': self.ctx.voro_num_rerun,
            'fac_cls_increase': self.ctx.voro_fac_cls_increase,
            'natom_in_cls_min': self.ctx.voro_natom_in_cls_min,
            'delta_e_min': self.ctx.voro_delta_e_min,
            'threshold_dos_zero': self.ctx.voro_threshold_dos_zero,
            'check_dos': self.ctx.voro_check_dos,
        })

        # set workflow parameters for the KKR impurity calculation
        self.ctx.kkr_runmax = wf_dict.get('kkr_runmax', self._wf_default['kkr_runmax'])
        self.ctx.convergence_criterion = wf_dict.get('convergence_criterion', self._wf_default['convergence_criterion'])
        self.ctx.mixreduce = wf_dict.get('mixreduce', self._wf_default['mixreduce'])
        self.ctx.threshold_aggressive_mixing = wf_dict.get(
            'threshold_aggressive_mixing', self._wf_default['threshold_aggressive_mixing']
        )
        self.ctx.strmix = wf_dict.get('strmix', self._wf_default['strmix'])
        self.ctx.nsteps = wf_dict.get('nsteps', self._wf_default['nsteps'])
        self.ctx.aggressive_mix = wf_dict.get('aggressive_mix', self._wf_default['aggressive_mix'])
        self.ctx.aggrmix = wf_dict.get('aggrmix', self._wf_default['aggrmix'])
        self.ctx.broyden_number = wf_dict.get('broyden-number', self._wf_default['broyden-number'])
        self.ctx.mag_init = wf_dict.get('mag_init', self._wf_default['mag_init'])
        self.ctx.hfield = wf_dict.get('hfield', self._wf_default['hfield'])
        self.ctx.init_pos = wf_dict.get('init_pos', self._wf_default['init_pos'])
        self.ctx.accuracy_params = wf_dict.get('accuracy_params', self._wf_default['accuracy_params'])
        # set up new parameter dict to pass to kkrimp subworkflow later
        self.ctx.kkrimp_params_dict = Dict({
            'nsteps': self.ctx.nsteps,
            'kkr_runmax': self.ctx.kkr_runmax,
            'threshold_aggressive_mixing': self.ctx.threshold_aggressive_mixing,
            'convergence_criterion': self.ctx.convergence_criterion,
            'mixreduce': self.ctx.mixreduce,
            'strmix': self.ctx.strmix,
            'aggressive_mix': self.ctx.aggressive_mix,
            'aggrmix': self.ctx.aggrmix,
            'broyden-number': self.ctx.broyden_number,
            'mag_init': self.ctx.mag_init,
            'hfield': self.ctx.hfield,
            'init_pos': self.ctx.init_pos,
            'accuracy_params': self.ctx.accuracy_params,
        })

        # retrieve option for kkrlfex files
        self.ctx.retrieve_kkrflex = wf_dict.get('retrieve_kkrflex', self._wf_default['retrieve_kkrflex'])

        # report the chosen parameters to the user
        self.report(
            'INFO: use the following parameter:\n'
            '\nGeneral settings\n'
            'use mpi: {}\n'
            'resources: {}\n'
            'walltime (s): {}\n'
            'queue name: {}\n'
            'scheduler command: {}\n'
            'description: {}\n'
            'label: {}\n'
            'parameters for the voroaux calculation: {}\n'
            'parameters for the kkrimp scf: {}\n'
            ''.format(
                self.ctx.withmpi, self.ctx.resources, self.ctx.max_wallclock_seconds, self.ctx.queue,
                self.ctx.custom_scheduler_commands, self.ctx.description_wf, self.ctx.label_wf,
                self.ctx.voro_params_dict.get_dict(), self.ctx.kkrimp_params_dict.get_dict()
            )
        )

    def validate_input(self):
        """
        Validate the input and catch possible errors from the input
        """

        inputs = self.inputs
        inputs_ok = True

        if 'kkrimp' and 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.kkrimp, 'kkr.kkrimp', use_exceptions=True)
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                inputs_ok = False
                self.report(self.exit_codes.ERROR_INVALID_INPUT_CODE)  # pylint: disable=no-member
                self.ctx.exit_code = self.exit_codes.ERROR_INVALID_INPUT_CODE  # pylint: disable=no-member
                return False
        elif 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                inputs_ok = False
                self.report(self.exit_codes.ERROR_INVALID_INPUT_CODE)  # pylint: disable=no-member
                self.ctx.exit_code = self.exit_codes.ERROR_INVALID_INPUT_CODE  # pylint: disable=no-member
                return False

        if 'impurity_info' in inputs:
            self.report(f'INFO: found the following impurity info node in input: {inputs.impurity_info.get_dict()}')

        if 'remote_data_gf' in inputs and 'remote_data_host' in inputs:
            self.report(
                'INFO: both converged host remote (pid: {}) and GF writeout remote (pid: {}) found in input. '
                'Converged host remote will not be used. Skip GF writeout step and '
                'start workflow with auxiliary voronoi calculations.'.format(
                    inputs.remote_data_host.pk, inputs.remote_data_gf.pk
                )
            )
            do_gf_calc = False
        elif 'remote_data_host' in inputs:
            self.report(
                'INFO: found converged host remote (pid: {}) in input. '
                'Start workflow by calculating the host GF.'.format(inputs.remote_data_host.pk)
            )
            if 'kkr' in inputs:
                do_gf_calc = True
            else:
                inputs_ok = False
                self.report(self.exit_codes.ERROR_MISSING_KKRCODE)  # pylint: disable=no-member
                self.ctx.exit_code = self.exit_codes.ERROR_MISSING_KKRCODE  # pylint: disable=no-member
                return False
        elif 'remote_data_gf' in inputs:
            remote_data_gf_node = load_node(inputs.remote_data_gf.pk)
            pk_kkrflex_writeoutcalc = remote_data_gf_node.get_incoming(link_label_filter=u'remote_folder'
                                                                       ).first().node.pk
            self.report(
                'INFO: found remote_data node (pid: {}) from previous KKRFLEX calculation (pid: {}) in input. '
                'Skip GF writeout step and start workflow by auxiliary voronoi calculations.'.format(
                    inputs.remote_data_gf.pk, pk_kkrflex_writeoutcalc
                )
            )
            do_gf_calc = False
            # check if second remote_data_gf node is given (used to overwrite Fermi level)
            if 'remote_data_gf_Efshift' in inputs:
                self.report(
                    'INFO: found remote_data_gf_Efshift (pid: {}) used to overwrite Fermi level.'.format(
                        inputs.remote_data_gf_Efshift.pk
                    )
                )
        else:
            inputs_ok = False
            self.report(self.exit_codes.ERROR_MISSING_REMOTE)  # pylint: disable=no-member
            self.ctx.exit_code = self.exit_codes.ERROR_MISSING_REMOTE  # pylint: disable=no-member
            return False

        self.ctx.do_gf_calc = do_gf_calc
        self.report(f'INFO: validated input successfully: {inputs_ok}. Do GF writeout calc: {self.ctx.do_gf_calc}.')

        return do_gf_calc

    def run_gf_writeout(self):
        """
        Run the gf_writeout workflow to calculate the host Green's function and the
        KKR flexfiles using the converged host remote folder and the impurity info node
        """

        # collect inputs
        kkrcode = self.inputs.kkr
        imp_info = self.inputs.impurity_info
        converged_host_remote = self.inputs.remote_data_host
        options = self.ctx.options_params_dict

        # set label and description of the calc
        sub_label = f'GF writeout (conv. host pid: {converged_host_remote.pk}, imp_info pid: {imp_info.pk})'
        sub_description = 'GF writeout sub workflow for kkrimp_wc using converged host remote data (pid: {}) and impurity_info node (pid: {})'.format(
            converged_host_remote.pk, imp_info.pk
        )

        builder = kkr_flex_wc.get_builder()

        builder.metadata.label = sub_label  # pylint: disable=no-member
        builder.metadata.description = sub_description  # pylint: disable=no-member
        builder.kkr = kkrcode
        builder.options = options
        builder.remote_data = converged_host_remote
        builder.impurity_info = imp_info

        if 'params_kkr_overwrite' in self.inputs:
            builder.params_kkr_overwrite = self.inputs.params_kkr_overwrite

        if 'gf_writeout' in self.inputs:
            if 'options' in self.inputs.gf_writeout:
                builder.options = self.inputs.gf_writeout.options
            if 'params_kkr_overwrite' in self.inputs.gf_writeout:
                builder.params_kkr_overwrite = self.inputs.gf_writeout.params_kkr_overwrite

        # maybe set kkrflex_retrieve
        wf_params_gf = {}
        if not self.ctx.retrieve_kkrflex:
            wf_params_gf['retrieve_kkrflex'] = self.ctx.retrieve_kkrflex
        wf_params_gf = Dict(wf_params_gf)
        builder.wf_parameters = wf_params_gf

        future = self.submit(builder)

        self.report(f'INFO: running GF writeout (pk: {future.pk})')

        return ToContext(gf_writeout=future, last_calc_gf=future)

    def has_starting_potential_input(self):
        """
        check whether or not a starting potential needs to be created
        """
        # initialize
        self.ctx.create_startpot = True

        # check if startpot exists in input
        # TODO maybe implement some consistency checks
        if 'startpot' in self.inputs:
            self.ctx.startpot_kkrimp = self.inputs.startpot
            self.ctx.create_startpot = False

        if self.ctx.exit_code is not None:
            # skip creation of starting potential by overwriting with True return value
            return True

        return self.ctx.create_startpot

    def run_voroaux(self):
        """
        Perform a voronoi calculation for every impurity charge using the structure
        from the converged KKR host calculation
        """
        # TODO: generalize to multiple impurities

        # collect inputs
        vorocode = self.inputs.voronoi
        kkrcode = self.inputs.kkr
        imp_info = self.inputs.impurity_info
        voro_params = self.ctx.voro_params_dict

        if self.ctx.do_gf_calc:
            self.report('INFO: get converged host remote from inputs to extract structure for Voronoi calculation')
            converged_host_remote = self.inputs.remote_data_host
        else:
            self.report(
                'INFO: get converged host remote from GF_host_calc and graph to extract structure for Voronoi calculation'
            )
            remote_data_gf_node = load_node(self.inputs.remote_data_gf.pk)
            GF_host_calc = remote_data_gf_node.get_incoming(link_label_filter=u'remote_folder').first().node
            converged_host_remote = GF_host_calc.inputs.parent_folder

        # find host structure
        structure_host, voro_calc = VoronoiCalculation.find_parent_structure(converged_host_remote)

        # get previous kkr parameters following remote_folder->calc->parameters links
        prev_kkrparams = converged_host_remote.get_incoming(link_label_filter='remote_folder'
                                                            ).first().node.get_incoming(link_label_filter='parameters'
                                                                                        ).first().node
        calc_params = prev_kkrparams

        # set Fermi level for auxiliary impurity potential correctly (extract from EF that is used in impurity calc)
        set_efermi = self.get_ef_from_parent()
        self.report(f'INFO: set Fermi level in jellium starting potential to {set_efermi}')
        # change voronoi parameters
        updatenode = Dict({'ef_set': set_efermi, 'add_direct': True})
        updatenode.label = 'Added Fermi energy'
        voro_params = update_params_wf(voro_params, updatenode)

        # add or overwrite some parameters (e.g. things that are only used by voronoi)
        calc_params_dict = calc_params.get_dict()
        # add some voronoi specific parameters automatically if found (RMTREF should also set RMTCORE to the same value)
        if calc_params_dict.get('<RMTREF>', None) is not None:
            self.report('INFO: add rmtcore to voro params')
            self.ctx.change_voro_params['<RMTCORE>'] = calc_params_dict['<RMTREF>']
            self.report(self.ctx.change_voro_params)

        # add some voronoi-specific settings starting from host's voronoi run (might have been overwritten in between)
        # this is necessary to make sure that voronoi creates the same radial mesh for the impurity potential, otherwise KKRimp will fail
        if voro_calc.inputs.parameters.get_dict().get('RUNOPT', None) is not None:
            self.report("INFO: copy runopt from host's voronoi run")
            runopt = self.ctx.change_voro_params.get('RUNOPT', None)
            if runopt is None:
                runopt = []
            runopt += voro_calc.inputs.parameters['RUNOPT']
            self.ctx.change_voro_params['RUNOPT'] = runopt
        if voro_calc.inputs.parameters.get_dict().get('RCLUSTZ', None) is not None:
            self.report("INFO: copy RCLUSTZ from host's voronoi run")
            self.ctx.change_voro_params['RCLUSTZ'] = voro_calc.inputs.parameters['RCLUSTZ']

        changed_params = False
        for key, val in self.ctx.change_voro_params.items():
            if key in ['RUNOPT', 'TESTOPT']:
                opt_old = calc_params_dict.get(key, [])
                if type(val) != list:
                    val = [val]
                val = opt_old + val
            calc_params_dict[key] = val
            changed_params = True
        if changed_params:
            updatenode = Dict(calc_params_dict)
            updatenode.label = f'Changed params for voroaux: {list(self.ctx.change_voro_params.keys())}'
            updatenode.description = 'Overwritten voronoi input parameter from kkr_imp_wc input.'
            calc_params = update_params_wf(calc_params, updatenode)

        # for every impurity, generate a structure and launch the voronoi workflow
        # to get the auxiliary impurity startpotentials
        self.ctx.voro_calcs = {}
        inter_struc = change_struc_imp_aux_wf(structure_host, imp_info)
        sub_label = f"voroaux calc for Zimp: {imp_info.get_dict().get('Zimp')} in host-struc"
        sub_description = 'Auxiliary voronoi calculation for an impurity with charge '
        sub_description += '{} in the host structure from pid: {}'.format(
            imp_info.get_dict().get('Zimp'), converged_host_remote.pk
        )

        builder = kkr_startpot_wc.get_builder()
        builder.metadata.label = sub_label  # pylint: disable=no-member
        builder.metadata.description = sub_description  # pylint: disable=no-member
        builder.structure = inter_struc
        builder.voronoi = vorocode
        builder.kkr = kkrcode
        builder.wf_parameters = voro_params
        builder.calc_parameters = calc_params
        builder.options = Dict(self.ctx.options_params_dict_voronoi)
        future = self.submit(builder)

        tmp_calcname = f'voro_aux_{1}'
        self.ctx.voro_calcs[tmp_calcname] = future
        self.report(f"INFO: running voro aux (Zimp= {imp_info.get_dict().get('Zimp')}, pid: {future.pk})")

        return ToContext(last_voro_calc=future)

    def get_ef_from_parent(self):
        """
        Extract Fermi level in Ry to which starting potential is set
        """
        # first choose calculation to start from (3 possibilities)
        if self.ctx.do_gf_calc:
            parent_remote = self.inputs.remote_data_host
        elif 'remote_data_gf_Efshift' in self.inputs:
            parent_remote = load_node(self.inputs.remote_data_gf_Efshift.pk)
        else:
            parent_remote = load_node(self.inputs.remote_data_gf.pk)

        # now extract output parameters
        parent_calc = parent_remote.get_incoming(link_label_filter='remote_folder').first().node
        output_params = parent_calc.outputs.output_parameters.get_dict()

        # get fermi energy in Ry from output of KkrCalculation and return result
        set_efermi = output_params.get('fermi_energy')

        return set_efermi

    def construct_startpot(self):
        """
        Take the output of GF writeout and the converged host potential as well as the
        auxiliary startpotentials for the impurity to construct the startpotential for the
        KKR impurity sub workflow
        """

        if not self.ctx.last_voro_calc.is_finished_ok:
            self.report(self.exit_codes.ERROR_KKRSTARTPOT_WORKFLOW_FAILURE)  # pylint: disable=no-member
            self.ctx.exit_code = self.exit_codes.ERROR_KKRSTARTPOT_WORKFLOW_FAILURE  # pylint: disable=no-member

        # collect all nodes necessary to construct the startpotential
        if self.ctx.do_gf_calc:
            GF_host_calc_pk = self.ctx.gf_writeout.outputs.workflow_info.get_dict().get('pk_flexcalc')
            self.report(f'GF_host_calc_pk: {GF_host_calc_pk}')
            GF_host_calc = load_node(GF_host_calc_pk)
            converged_host_remote = self.inputs.remote_data_host
        else:
            remote_data_gf_node = load_node(self.inputs.remote_data_gf.pk)
            GF_host_calc = remote_data_gf_node.get_incoming(link_label_filter=u'remote_folder').first().node
            self.report(f'GF_host_calc_pk: {GF_host_calc.pk}')
            # follow parent_folder link up to get remote folder
            converged_host_remote = GF_host_calc.get_incoming(link_label_filter='parent_folder').first().node

        # get remote folder of last voronoi calculation (i.e. the one from where we take the starting potential)
        print(self.ctx.last_voro_calc)
        all_nodes = self.ctx.last_voro_calc.get_outgoing(node_class=CalcJobNode).all()
        print(all_nodes)
        pk_last_voronoi = max([i.node.pk for i in all_nodes])
        print(pk_last_voronoi)
        voro_calc_remote = load_node(pk_last_voronoi).outputs.remote_folder
        print(voro_calc_remote)

        print(load_node(pk_last_voronoi).outputs.retrieved.list_object_names())
        print(GF_host_calc)
        print(GF_host_calc.outputs.retrieved.list_object_names())

        # check wether or not calculation was taked from cached node
        caching_info = f'INFO: cache_source of GF_host_calc node: {GF_host_calc.get_cache_source()}'
        print(caching_info)
        self.report(caching_info)
        caching_info = f'INFO: cache_source of voronoi node: {load_node(pk_last_voronoi).get_cache_source()}'
        print(caching_info)
        self.report(caching_info)

        imp_info = self.inputs.impurity_info
        nspin = GF_host_calc.outputs.output_parameters.get_dict().get('nspin')

        ilayer_cent = imp_info.get_dict().get('ilayer_center', 0)  # defaults to layer 0

        # prepare settings dict
        potname_converged = 'out_potential'
        potname_impvorostart = 'output.pot'
        potname_imp = 'potential_imp'

        if nspin < 2:
            replacelist_pot2 = [[0, ilayer_cent]]
        else:
            replacelist_pot2 = [[0, 2 * ilayer_cent], [1, 2 * ilayer_cent + 1]]
        try:
            with GF_host_calc.outputs.retrieved.open('scoef') as _f:
                neworder_pot1 = [int(i) for i in np.loadtxt(_f, skiprows=1)[:, 3] - 1]
        except:
            with GF_host_calc.outputs.retrieved.open('scoef') as _f:
                neworder_pot1 = [int(np.loadtxt(_f, skiprows=1)[3] - 1)]

        settings_label = f'startpot_KKRimp for imp_info node {imp_info.pk}'
        settings_description = f'starting potential for impurity info: {imp_info}'

        settings = Dict({
            'out_pot': potname_imp,
            'neworder': neworder_pot1,
            'replace_newpos': replacelist_pot2,
            'label': settings_label,
            'description': settings_description
        })
        print('startpot_kkrimp construction:', settings, converged_host_remote, voro_calc_remote)
        startpot_kkrimp = neworder_potential_wf(
            settings_node=settings, parent_calc_folder=converged_host_remote, parent_calc_folder2=voro_calc_remote
        )

        # add starting potential for kkrimp calculation to context
        self.ctx.startpot_kkrimp = startpot_kkrimp

        self.report(
            'INFO: created startpotential (pid: {}) for the impurity calculation '
            'by using information of the GF host calculation (pid: {}), the potential of the '
            'converged host system (remote pid: {}) and the potential of the auxiliary voronoi '
            'calculation (remote pid: {})'.format(
                startpot_kkrimp.pk, GF_host_calc.pk, converged_host_remote.pk, self.ctx.last_voro_calc.pk
            )
        )

    def run_kkrimp_scf(self):
        """
        Uses both the previously generated host-impurity startpotential and the output from
        the GF writeout workflow as inputs to run the kkrimp_sub workflow in order to
        converge the host-impurity potential
        """

        # collect all necessary input nodes
        kkrimpcode = self.inputs.kkrimp
        startpot = self.ctx.startpot_kkrimp
        kkrimp_params = self.ctx.kkrimp_params_dict
        options = self.ctx.options_params_dict
        imp_info = self.inputs.impurity_info
        if self.ctx.do_gf_calc:
            self.report(f'INFO: get GF remote from gf_writeout sub wf (pid: {self.ctx.gf_writeout.pk})')
            gf_remote = self.ctx.gf_writeout.outputs.GF_host_remote
        else:
            self.report(f'INFO: get GF remote from input node (pid: {self.inputs.remote_data_gf.pk})')
            gf_remote = self.inputs.remote_data_gf
        # save in context to return as output node
        self.ctx.gf_remote = gf_remote

        # set label and description
        sub_label = f'kkrimp_sub scf wf (GF host remote: {gf_remote.pk}, imp_info: {self.inputs.impurity_info.pk})'
        sub_description = 'convergence of the host-impurity potential (pk: {}) using GF remote (pk: {})'.format(
            startpot.pk, gf_remote.pk
        )

        builder = kkr_imp_sub_wc.get_builder()
        builder.metadata.label = sub_label  # pylint: disable=no-member
        builder.metadata.description = sub_description  # pylint: disable=no-member
        builder.kkrimp = kkrimpcode
        builder.options = options
        builder.impurity_info = imp_info
        builder.host_imp_startpot = startpot
        builder.remote_data = gf_remote
        if 'remote_data_gf_Efshift' in self.inputs:
            builder.remote_data_Efshift = self.inputs.remote_data_gf_Efshift
        if 'scf' in self.inputs:
            if 'params_overwrite' in self.inputs.scf:
                builder.params_overwrite = self.inputs.scf.params_overwrite
            if 'options' in self.inputs.scf:
                builder.options = self.inputs.scf.options
            if 'initial_noco_angles' in self.inputs.scf:
                builder.initial_noco_angles = self.inputs.scf.initial_noco_angles
            if 'rimpshift' in self.inputs.scf:
                builder.rimpshift = self.inputs.scf.rimpshift
            if 'settings_LDAU' in self.inputs.scf:
                builder.settings_LDAU = self.inputs.scf.settings_LDAU
        builder.wf_parameters = kkrimp_params
        future = self.submit(builder)

        self.report(
            f'INFO: running kkrimp_sub_wf (startpot: {startpot.pk}, GF_remote: {gf_remote.pk}, wf pid: {future.pk})'
        )

        return ToContext(kkrimp_scf_sub=future)

    def return_results(self):
        """
        Return the results and create all of the output nodes
        """

        self.report('INFO: creating output nodes for the KKR impurity workflow ...')

        if self.ctx.kkrimp_scf_sub.is_finished_ok:
            link_nodes = {'results_scf_workflow': self.ctx.kkrimp_scf_sub.outputs.workflow_info}

            last_calc_pk = self.ctx.kkrimp_scf_sub.outputs.workflow_info.get_dict().get('last_calc_nodeinfo')['pk']
            last_calc_output_params = load_node(last_calc_pk).outputs.output_parameters
            last_calc_info = self.ctx.kkrimp_scf_sub.outputs.workflow_info
            outputnode_dict = {}
            outputnode_dict['workflow_name'] = self.__class__.__name__
            outputnode_dict['workflow_version'] = self._workflowversion
            if self.ctx.do_gf_calc:
                outputnode_dict['used_subworkflows'] = {
                    'gf_writeout': self.ctx.gf_writeout.pk,
                    'kkr_imp_sub': self.ctx.kkrimp_scf_sub.pk
                }
                outputnode_dict['gf_wc_success'] = self.ctx.gf_writeout.outputs.workflow_info.get_dict(
                ).get('successful')
                link_nodes['gf_writeout'] = self.ctx.gf_writeout.outputs.workflow_info
            else:
                outputnode_dict['used_subworkflows'] = {'kkr_imp_sub': self.ctx.kkrimp_scf_sub.pk}
            if self.ctx.create_startpot:
                outputnode_dict['used_subworkflows']['auxiliary_voronoi'] = self.ctx.last_voro_calc.pk
                res_voro_info = self.ctx.last_voro_calc.outputs.results_vorostart_wc
                outputnode_dict['voro_wc_success'] = res_voro_info.get_dict().get('successful')
                link_nodes['results_startpot_workflow'] = self.ctx.last_voro_calc.outputs.results_vorostart_wc
            outputnode_dict['converged'] = last_calc_info.get_dict().get('convergence_reached')
            outputnode_dict['number_of_rms_steps'] = len(last_calc_info.get_dict().get('convergence_values_all_steps'))
            outputnode_dict['convergence_values_all_steps'] = last_calc_info.get_dict(
            ).get('convergence_values_all_steps')
            outputnode_dict['impurity_info'] = self.inputs.impurity_info.get_dict()
            outputnode_dict['kkrimp_wc_success'] = last_calc_info.get_dict().get('successful')
            outputnode_dict['last_calculation_uuid'] = load_node(last_calc_pk).uuid

            # create results node and link all sub-workflow output nodes
            outputnode_t = create_out_dict_node(Dict(dict=outputnode_dict), **link_nodes)
            outputnode_t.label = 'kkrimp_wc_inform'
            outputnode_t.description = 'Contains information for workflow'
            self.report(f'INFO: workflow_info node: {outputnode_t.uuid}')

            self.out('workflow_info', outputnode_t)
            self.out('last_calc_output_parameters', last_calc_output_params)
            self.out('last_calc_info', last_calc_info)
            self.out('converged_potential', self.ctx.kkrimp_scf_sub.outputs.host_imp_pot)
            self.out('remote_data_gf', self.ctx.gf_remote)

            # print final message before exiting
            self.report('INFO: created 3 output nodes for the KKR impurity workflow.')
            self.report(
                '\n'
                '|------------------------------------------------------------------------------------------------------------------|\n'
                '|-------------------------------------| Done with the KKR impurity workflow! |-------------------------------------|\n'
                '|------------------------------------------------------------------------------------------------------------------|'
            )
        else:
            self.report(self.exit_codes.ERROR_KKRIMP_SUB_WORKFLOW_FAILURE)  # pylint: disable=no-member
            return self.exit_codes.ERROR_KKRIMP_SUB_WORKFLOW_FAILURE  # pylint: disable=no-member

    def error_handler(self):
        """Capture errors raised in validate_input"""
        if self.ctx.exit_code is not None:
            return self.ctx.exit_code


@calcfunction
def change_struc_imp_aux_wf(struc, imp_info):  # Note: works for single imp at center only!
    from aiida.common.constants import elements as PeriodicTableElements
    _atomic_numbers = {data['symbol']: num for num, data in PeriodicTableElements.items()}

    new_struc = StructureData(cell=struc.cell)
    new_struc.pbc = struc.pbc  # take also pbc values from parent struc
    isite = 0
    for site in struc.sites:
        sname = site.kind_name
        kind = struc.get_kind(sname)
        pos = site.position
        # intermediate fix to avoid crash for old structures with vacuum:'{H0.00X1.00}'
        # use atom kind='X' in the future for new structures
        if kind.get_symbols_string() == '{H0.00X1.00}':
            zatom = 0
        else:
            zatom = _atomic_numbers[kind.get_symbols_string()]
        if isite == imp_info.get_dict().get('ilayer_center', 0):
            zatom = imp_info.get_dict().get('Zimp')
            if type(zatom) == list:
                zatom = zatom[0]  # here this works for single impurity only!
        symbol = PeriodicTableElements.get(zatom).get('symbol')
        new_struc.append_atom(position=pos, symbols=symbol)
        isite += 1

    return new_struc
