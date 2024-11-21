# -*- coding: utf-8 -*-
"""
In this module you find the sub workflow for the kkrimp self consistency cycle
and some helper methods to do so with AiiDA
"""

from aiida import __version__ as aiida_core_version
from aiida.orm import Float, Code, CalcJobNode, RemoteData, StructureData, Dict, SinglefileData, FolderData, Bool
from aiida.engine import WorkChain, ToContext, while_, if_, calcfunction
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.tools import test_and_get_codenode, get_inputs_kkrimp, kick_out_corestates_wf, get_ldaumatrices, get_LDAU_initmatrices_dict
from aiida_kkr.calculations.kkrimp import KkrimpCalculation
from numpy import array
import tarfile, os
from aiida_kkr.tools.save_output_nodes import create_out_dict_node

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.10.4'
__contributors__ = (u'Fabian Bertoldo', u'Philipp Ruessmann', u'David Antognini Silva')

#TODO: work on return results function
#TODO: edit inspect_kkrimp function
#TODO: get rid of create_scf_result node and create output nodes differently
#TOTO: check if calculation parameters from previous calculation have to be
#      loaded (in validate input, compare to kkr workflow)
#TODO: maybe add decrease mixing factor option as in kkr_scf wc
#TODO: add option to check if the convergence is on track


class kkr_imp_sub_wc(WorkChain):
    """
    Workchain of a kkrimp self consistency calculation starting from the
    host-impurity potential of the system. (Not the entire kkr_imp workflow!)

    :param options: (Dict), Workchain specifications
    :param wf_parameters: (Dict), specifications for the calculation
    :param host_imp_startpot: (RemoteData), mandatory; input host-impurity potential
    :param kkrimp: (Code), mandatory; KKRimp code converging the host-imp-potential
    :param remote_data: (RemoteData), mandatory; remote folder of a previous
                           kkrflex calculation containing the flexfiles ...
    :param kkrimp_remote: (RemoteData), remote folder of a previous kkrimp calculation
    :param impurity_info: (Dict), Parameter node with information
                          about the impurity cluster

    :return workflow_info: (Dict), Information of workflow results
                                   like success, last result node, list with
                                   convergence behavior
    :return host_imp_pot: (SinglefileData), output potential of the sytem
    """

    _workflowversion = __version__
    _wf_label = 'kkr_imp_sub_wc'
    _wf_description = 'Workflow for a KKRimp self consistency calculation to converge a given host-impurity potential'

    _options_default = {
        'queue_name': '',  # Queue name to submit jobs too
        'resources': {
            'num_machines': 1
        },  # resources to allowcate for the job
        'max_wallclock_seconds': 60 * 60,  # walltime after which the job gets killed (gets parsed to KKR)}
        'custom_scheduler_commands': '',  # some additional scheduler commands
        'withmpi': True
    }  # execute KKR with mpi or without

    _wf_default = {
        'kkr_runmax': 5,  # Maximum number of kkr jobs/starts (defauld iterations per start)
        'convergence_criterion': 1e-7,  # Stop if charge denisty is converged below this value
        'mixreduce': 0.5,  # reduce mixing factor by this factor if calculaito fails due to too large mixing
        'threshold_aggressive_mixing': 1e-2,  # threshold after which agressive mixing is used
        'strmix': 0.01,  # mixing factor of simple mixing
        'nsteps': 50,  # number of iterations done per KKR calculation
        'aggressive_mix': 5,  # type of aggressive mixing (3: broyden's 1st, 4: broyden's 2nd, 5: generalized anderson)
        'aggrmix': 0.01,  # mixing factor of aggressive mixing
        'broyden-number': 20,  # number of potentials to 'remember' for Broyden's mixing
        'nsimplemixfirst': 0,  # number of simple mixing step at the beginning of Broyden mixing
        'mag_init': False,  # initialize and converge magnetic calculation
        'hfield': [0.02, 5],  # Ry                     # external magnetic field used in initialization step
        'init_pos': None,  # position in unit cell where magnetic field is applied [default (None) means apply to all]
        'dos_run':
        False,  # specify if DOS should be calculated (!KKRFLEXFILES with energy contour necessary as GF_remote_data!)
        'lmdos': True,  # specify if DOS calculation should calculate l-resolved or l and m resolved output
        'jij_run': False,  # specify if Jijs should be calculated (!changes behavior of the code!!!)
        #                   # Some parameter for direct solver (if None, use the same as in host code, otherwise overwrite)
        'accuracy_params': {
            'RADIUS_LOGPANELS': None,  # where to set change of logarithmic to linear radial mesh
            'NPAN_LOG': None,  # number of panels in log mesh
            'NPAN_EQ': None,  # number of panels in linear mesh
            'NCHEB': None
        }  # number of chebychev polynomials in each panel (total number of points in radial mesh NCHEB*(NPAN_LOG+NPAN_EQ))
    }

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create
        set of wf_parameters.

        returns _wf_defaults
        """
        if not silent:
            print(f'Version of workflow: {self._workflowversion}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """

        super(kkr_imp_sub_wc, cls).define(spec)

        # Define the inputs of the workflow
        spec.input('kkrimp', valid_type=Code, required=True)
        spec.input('host_imp_startpot', valid_type=SinglefileData, required=False)
        spec.input('remote_data', valid_type=RemoteData, required=False)
        spec.input('remote_data_Efshift', valid_type=RemoteData, required=False)
        spec.input('kkrimp_remote', valid_type=RemoteData, required=False)
        spec.input('impurity_info', valid_type=Dict, required=False)
        spec.input('options', valid_type=Dict, required=False)
        spec.input('wf_parameters', valid_type=Dict, required=False, default=lambda: Dict(dict=cls._wf_default))
        spec.input(
            'settings_LDAU', valid_type=Dict, required=False, help='LDA+U settings. See KKRimpCalculation for details.'
        )
        spec.input(
            'params_overwrite',
            valid_type=Dict,
            required=False,
            help='Dict of parameters that are given to the KKRimpCalculation. Overwrites automatically set values!'
        )

        spec.expose_inputs(KkrimpCalculation, include=('initial_noco_angles', 'rimpshift'))

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            if_(cls.validate_input)(
                while_(cls.condition)(
                    cls.update_kkrimp_params,
                    # TODO: encapsulate this in restarting mechanism (should be a base class of workflows that start calculations)
                    # i.e. use base_restart_calc workchain as parent
                    cls.run_kkrimp,
                    cls.inspect_kkrimp
                ),
                cls.return_results
            ),
            cls.error_handler
        )

        # exit codes
        spec.exit_code(
            121,
            'ERROR_HOST_IMP_POT_GF',
            message='ERROR: Not both host-impurity potential and GF remote '
            'found in the inputs. Provide either both of them or a '
            'RemoteData from a previous kkrimp calculation.'
        )
        spec.exit_code(
            122,
            'ERROR_INVALID_INPUT_KKRIMP',
            message='ERROR: The code you provided for KKRimp does not '
            'use the plugin kkr.kkrimp'
        )
        spec.exit_code(
            123,
            'ERROR_INVALID_HOST_IMP_POT',
            message='ERROR: Unable to extract parent paremeter node of '
            'input remote folder'
        )
        # probably not necessary
        spec.exit_code(124, 'ERROR_NO_CALC_PARAMS', message='ERROR: No calculation parameters provided')

        spec.exit_code(125, 'ERROR_SUB_FAILURE', message='ERROR: Last KKRcalc in SUBMISSIONFAILED state!\nstopping now')
        spec.exit_code(
            126, 'ERROR_MAX_STEPS_REACHED', message='ERROR: Maximal number of KKR restarts reached. Exiting now!'
        )
        spec.exit_code(
            127,
            'ERROR_SETTING_LAST_REMOTE',
            message='ERROR: Last_remote could not be set to a previous succesful calculation'
        )
        spec.exit_code(128, 'ERROR_MISSING_PARAMS', message='ERROR: There are still missing calculation parameters')
        spec.exit_code(129, 'ERROR_PARAMETER_UPDATE', message='ERROR: Parameters could not be updated')
        spec.exit_code(
            130, 'ERROR_LAST_CALC_NOT_FINISHED_OK', message='ERROR: Last calculation is not in finished state'
        )
        spec.exit_code(
            131,
            'ERROR_NO_CALC_FOUND_FOR_REMOTE_DATA',
            message='The input `remote_data` node has no valid calculation parent.'
        )
        spec.exit_code(
            132,
            'ERROR_REMOTE_DATA_CALC_UNSUCCESFUL',
            message='The parent calculation of the input `remote_data` node was not succesful.'
        )
        spec.exit_code(
            133,
            'ERROR_NO_OUTPUT_POT_FROM_LAST_CALC',
            message='ERROR: Last calculation does not have an output potential.'
        )

        # Define the outputs of the workflow
        spec.output('workflow_info', valid_type=Dict)
        spec.output('host_imp_pot', valid_type=SinglefileData, required=False)

    def start(self):
        """
        init context and some parameters
        """
        message = f'INFO: started KKR impurity convergence workflow version {self._workflowversion}'
        self.report(message)

        ####### init #######

        # internal para /control para
        self.ctx.loop_count = 0
        self.ctx.last_mixing_scheme = 0
        self.ctx.calcs = []
        self.ctx.exit_code = None
        # flags used internally to check whether the individual steps were successful
        self.ctx.kkr_converged = False
        self.ctx.kkr_step_success = False
        self.ctx.kkr_higher_accuracy = False
        # links to previous calculations
        self.ctx.last_calc = None
        self.ctx.last_params = None
        self.ctx.last_remote = None
        # link to previous host impurity potential
        self.ctx.last_pot = None
        # convergence info about rms etc. (used to determine convergence behavior)
        self.ctx.last_rms_all = []
        self.ctx.rms_all_steps = []
        self.ctx.last_neutr_all = []
        self.ctx.neutr_all_steps = []
        # LDA+U settings, either from input port or extracted from kkrimp_remote
        self.ctx.settings_LDAU = None

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()
        options_dict = self.inputs.options.get_dict()

        if options_dict == {}:
            options_dict = self._options_default
            message = 'INFO: using default options'
            self.report(message)

        if wf_dict == {}:
            wf_dict = self._wf_default
            message = 'INFO: using default wf parameter'
            self.report(message)

        # set option parameters from input, or defaults
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds', self._options_default['max_wallclock_seconds']
        )
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get(
            'custom_scheduler_commands', self._options_default['custom_scheduler_commands']
        )

        # set workflow parameters from input, or defaults
        self.ctx.max_number_runs = wf_dict.get('kkr_runmax', self._wf_default['kkr_runmax'])
        self.ctx.description_wf = self.inputs.get(
            'description', 'Workflow for '
            'a KKR impurity calculation'
            'starting from a host-impurity'
            'potential'
        )
        self.ctx.label_wf = self.inputs.get('label', 'kkr_imp_sub_wc')
        self.ctx.strmix = wf_dict.get('strmix', self._wf_default['strmix'])
        self.ctx.convergence_criterion = wf_dict.get('convergence_criterion', self._wf_default['convergence_criterion'])
        self.ctx.mixreduce = wf_dict.get('mixreduce', self._wf_default['mixreduce'])
        self.ctx.threshold_aggressive_mixing = wf_dict.get(
            'threshold_aggressive_mixing', self._wf_default['threshold_aggressive_mixing']
        )
        self.ctx.type_aggressive_mixing = wf_dict.get('aggressive_mix', self._wf_default['aggressive_mix'])
        self.ctx.aggrmix = wf_dict.get('aggrmix', self._wf_default['aggrmix'])
        self.ctx.nsteps = wf_dict.get('nsteps', self._wf_default['nsteps'])
        self.ctx.broyden_num = wf_dict.get('broyden-number', self._wf_default['broyden-number'])
        self.ctx.nsimplemixfirst = wf_dict.get('nsimplemixfirst', self._wf_default['nsimplemixfirst'])
        self.ctx.mesh_params = wf_dict.get('accuracy_params', {})

        # initial magnetization
        self.ctx.mag_init = wf_dict.get('mag_init', self._wf_default['mag_init'])
        self.ctx.hfield = wf_dict.get('hfield', self._wf_default['hfield'])
        self.ctx.xinit = wf_dict.get('init_pos', self._wf_default['init_pos'])
        self.ctx.mag_init_step_success = False

        # DOS
        self.ctx.dos_run = wf_dict.get('dos_run', self._wf_default['dos_run'])
        self.ctx.lmdos = wf_dict.get('lmdos', self._wf_default['lmdos'])
        # Jij
        self.ctx.jij_run = wf_dict.get('jij_run', self._wf_default['jij_run'])

        self.report(
            'INFO: use the following parameter:\n'
            '\nGeneral settings\n'
            'use mpi: {}\n'
            'max number of KKR runs: {}\n'
            'Resources: {}\n'
            'Walltime (s): {}\n'
            'queue name: {}\n'
            'scheduler command: {}\n'
            'description: {}\n'
            'label: {}\n'
            '\nMixing parameter\n'
            'Straight mixing factor: {}\n'
            'Nsteps scf cycle: {}\n'
            'threshold_aggressive_mixing: {}\n'
            'Aggressive mixing technique: {}\n'
            'Aggressive mixing factor: {}\n'
            'Mixing decrease factor if convergence fails: {}\n'
            'Convergence criterion: {}\n'
            '\nAdditional parameter\n'
            'init magnetism in first step: {}\n'
            'init magnetism, hfield: {}\n'
            'init magnetism, init_pos: {}\n'
            ''.format(
                self.ctx.withmpi, self.ctx.max_number_runs, self.ctx.resources, self.ctx.max_wallclock_seconds,
                self.ctx.queue, self.ctx.custom_scheduler_commands, self.ctx.description_wf, self.ctx.label_wf,
                self.ctx.strmix, self.ctx.nsteps, self.ctx.threshold_aggressive_mixing, self.ctx.type_aggressive_mixing,
                self.ctx.aggrmix, self.ctx.mixreduce, self.ctx.convergence_criterion, self.ctx.mag_init,
                self.ctx.hfield, self.ctx.xinit
            )
        )

        # return para/vars
        self.ctx.successful = False
        self.ctx.rms = []
        self.ctx.neutr = []
        self.ctx.warnings = []
        self.ctx.formula = ''

        # for results table each list gets one entry per iteration that has been performed
        self.ctx.KKR_steps_stats = {}
        # later contains these keys:
        # 'success', 'isteps', 'imix', 'mixfac', 'qbound', 'high_sett', 'first_rms', 'last_rms'
        # 'first_neutr', 'last_neutr', 'pk', 'uuid'

    def validate_input(self):
        """
        validate input and catch possible errors from the input
        """

        inputs = self.inputs
        inputs_ok = True

        if not 'kkrimp_remote' in inputs:
            if not ('host_imp_startpot' in inputs and 'remote_data' in inputs):
                inputs_ok = False
                self.ctx.exit_code = self.exit_codes.ERROR_HOST_IMP_POT_GF  # pylint: disable=maybe-no-member

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkrimp', use_exceptions=True)
            except ValueError:
                inputs_ok = False
                self.ctx.exit_code = self.exit_codes.ERROR_INVALID_INPUT_KKRIMP  # pylint: disable=maybe-no-member

        # check if LDA+U settings should be set from input port
        if 'settings_LDAU' in inputs:
            self.ctx.settings_LDAU = inputs.settings_LDAU

        if 'kkrimp_remote' in inputs:
            self.ctx.start_from_imp_remote = True
            kkrimp_remote = inputs.kkrimp_remote
            self.ctx.last_remote = kkrimp_remote
            # check if LDA+U settings should be set kkrimp parent calculation
            if 'settings_LDAU' not in inputs:
                # check if kkrimp parent calculation has LDA+U input
                parent_kkrimp_calc = kkrimp_remote.get_incoming(node_class=CalcJobNode).first().node
                if 'settings_LDAU' in parent_kkrimp_calc.inputs:
                    self.ctx.settings_LDAU = parent_kkrimp_calc.inputs.settings_LDAU

        # check if input remote_data node is fine
        if 'remote_data' in inputs:
            if len(inputs.remote_data.get_incoming(link_label_filter='remote_folder').all()) < 1:
                self.ctx.exit_code = self.exit_codes.ERROR_NO_CALC_FOUND_FOR_REMOTE_DATA  # pylint: disable=maybe-no-member
            else:
                if not inputs.remote_data.get_incoming(link_label_filter='remote_folder').first().node.is_finished_ok:
                    self.ctx.exit_code = self.exit_codes.ERROR_REMOTE_DATA_CALC_UNSUCCESFUL  # pylint: disable=maybe-no-member

        # set starting potential
        if 'host_imp_startpot' in inputs:
            self.ctx.last_pot = inputs.host_imp_startpot

        # TBD!!!
        if 'wf_parameters' in inputs:
            self.ctx.last_params = inputs.wf_parameters
        else:
            inputs_ok = False
            self.ctx.exit_code = self.exit_codes.ERROR_NO_CALC_PARAMS  # pylint: disable=maybe-no-member

        message = f'INFO: validated input successfully: {inputs_ok}'
        self.report(message)
        if not inputs_ok:
            message = f'Exit code: {self.exit_codes.ERROR_NO_CALC_PARAMS}'  # pylint: disable=maybe-no-member
            self.report(message)

        return inputs_ok

    def condition(self):
        """
        check convergence condition
        """

        do_kkr_step = True
        stopreason = ''

        #increment KKR runs loop counter
        self.ctx.loop_count += 1

        # check if previous calculation reached convergence criterion
        if self.ctx.kkr_converged:
            if not self.ctx.kkr_higher_accuracy:
                do_kkr_step = do_kkr_step & True
            else:
                stopreason = 'KKR converged'
                self.ctx.successful = True
                do_kkr_step = False
        else:
            do_kkr_step = do_kkr_step & True

        # check if previous calculation was successful
        if self.ctx.loop_count > 1 and not self.ctx.last_calc.is_finished_ok:
            message = 'ERROR: last calc not finished_ok'
            self.report(message)
            return self.exit_codes.ERROR_SUB_FAILURE  # pylint: disable=no-member

        # next check only needed if another iteration should be done after validating convergence etc. (previous checks)
        if do_kkr_step:
            # check if maximal number of iterations has been reached
            if self.ctx.loop_count <= self.ctx.max_number_runs:
                do_kkr_step = do_kkr_step & True
            else:
                do_kkr_step = False

        message = f'INFO: done checking condition for kkr step (result={do_kkr_step})'
        self.report(message)

        if not do_kkr_step:
            message = f'INFO: Stopreason={stopreason}'
            self.report(message)

        return do_kkr_step

    def update_kkrimp_params(self):
        """
        update set of KKR parameters (check for reduced mixing, change of
        mixing strategy, change of accuracy setting)
        """

        decrease_mixing_fac = False
        switch_agressive_mixing = False
        switch_higher_accuracy = False
        initial_settings = False

        # only do something other than simple mixing after first kkr run
        if self.ctx.loop_count != 1:
            # first determine if previous step was successful (otherwise try to find some rms value and decrease mixing to try again)
            if not self.ctx.kkr_step_success:
                decrease_mixing_fac = True
                message = 'INFO: last KKR calculation failed. Trying decreasing mixfac'
                self.report(message)

            convergence_on_track = self.convergence_on_track()

            # check if calculation was on its way to converge
            if not convergence_on_track:
                decrease_mixing_fac = True
                message = 'INFO: Last KKR did not converge. Trying decreasing mixfac'
                self.report(message)
                # reset last_remote to last successful calculation
                last_calcs_list = list(range(len(self.ctx.calcs)))  # needs to be list to support slicing
                if len(last_calcs_list) > 1:
                    last_calcs_list = array(last_calcs_list)[::-1]  # make sure to go from latest calculation backwards
                for icalc in last_calcs_list:
                    message = f"INFO: last calc success? {icalc} {self.ctx.KKR_steps_stats['success'][icalc]}"
                    self.report(message)
                    if self.ctx.KKR_steps_stats['success'][icalc]:
                        if self.ctx.KKR_steps_stats['last_rms'][icalc] < self.ctx.KKR_steps_stats['first_rms'][icalc]:
                            self.ctx.last_remote = self.ctx.calcs[icalc].outputs.remote_folder
                            break  # exit loop if last_remote was found successfully
                        else:
                            self.ctx.last_remote = None
                    else:
                        self.ctx.last_remote = None
                # now cover case when last_remote needs to be set to initial remote folder (from input)
                if self.ctx.last_remote is None:
                    if 'kkrimp_remote' in self.inputs:
                        messager = 'INFO: no successful and converging calculation to take RemoteData from. Reuse RemoteData from input instead.'
                        self.report(message)
                        self.ctx.last_remote = self.inputs.kkrimp_remote
                    elif 'impurity_info' in self.inputs or 'remote_data' in self.inputs:
                        self.ctx.last_remote = None
                # check if last_remote has finally been set and abort if this is not the case
                if self.ctx.last_remote is None:
                    messager = 'ERROR: last remote not found'
                    self.report(message)
                    return self.exit_codes.ERROR_SETTING_LAST_REMOTE  # pylint: disable=no-member

            # check if mixing strategy should be changed
            last_mixing_scheme = self.ctx.last_params.get_dict()['IMIX']
            if last_mixing_scheme is None:
                last_mixing_scheme = 0

            if convergence_on_track:
                last_rms = self.ctx.last_rms_all[-1]
                if last_rms < self.ctx.threshold_aggressive_mixing and last_mixing_scheme == 0:
                    switch_agressive_mixing = True
                    message = 'INFO: rms low enough, switch to agressive mixing'
                    self.report(message)

                # check if switch to higher accuracy should be done
                if not self.ctx.kkr_higher_accuracy:
                    if self.ctx.kkr_converged:  # or last_rms < self.ctx.threshold_switch_high_accuracy:
                        switch_higher_accuracy = True


#                        self.report("INFO: rms low enough, switch to higher accuracy settings")
        else:
            initial_settings = True
            self.ctx.kkr_step_success = True

        if self.ctx.loop_count > 1:
            last_rms = self.ctx.last_rms_all[-1]

        # extract values from host calculation
        host_GF_calc = self.inputs.remote_data.get_incoming(node_class=CalcJobNode).first().node
        host_GF_outparams = host_GF_calc.outputs.output_parameters.get_dict()
        host_GF_inparams = host_GF_calc.inputs.parameters.get_dict()
        nspin = host_GF_outparams.get('nspin')
        non_spherical = host_GF_inparams.get('INS')
        if non_spherical is None:
            non_spherical = kkrparams.get_KKRcalc_parameter_defaults()[0].get('INS')
        self.ctx.spinorbit = host_GF_outparams.get('use_newsosol')

        # if needed update parameters
        if decrease_mixing_fac or switch_agressive_mixing or switch_higher_accuracy or initial_settings or self.ctx.mag_init:
            if initial_settings:
                label = 'initial KKR scf parameters'
                description = 'initial parameter set for scf calculation'
            else:
                label = ''
                description = ''

            # step 1: extract info from last input parameters and check consistency
            para_check = kkrparams(params_type='kkrimp')
            para_check.get_all_mandatory()
            message = 'INFO: get kkrimp keywords'
            self.report(message)

            # init new_params dict where updated params are collected
            new_params = {}

            # step 1.2: check if all mandatory keys are there and add defaults if missing
            missing_list = para_check.get_missing_keys(use_aiida=True)
            if missing_list != []:
                kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
                kkrdefaults_updated = []
                for key_default, val_default in list(kkrdefaults.items()):
                    if key_default in missing_list:
                        new_params[key_default] = kkrdefaults.get(key_default)
                        kkrdefaults_updated.append(key_default)
                if len(kkrdefaults_updated) > 0:
                    self.report('ERROR: no default param found')
                    return self.exit_codes.ERROR_MISSING_PARAMS  # pylint: disable=no-member
                else:
                    message = f'updated KKR parameter node with default values: {kkrdefaults_updated}'
                    self.report(message)

            # step 2: change parameter (contained in new_params dictionary)
            last_mixing_scheme = para_check.get_value('IMIX')
            if last_mixing_scheme is None:
                last_mixing_scheme = 0

            strmixfac = self.ctx.strmix
            aggrmixfac = self.ctx.aggrmix
            nsteps = self.ctx.nsteps

            # TODO: maybe add decrease mixing factor option as in kkr_scf wc
            # step 2.1 fill new_params dict with values to be updated
            if decrease_mixing_fac:
                if last_mixing_scheme == 0:
                    self.report(f'(strmixfax, mixreduce)= ({strmixfac}, {self.ctx.mixreduce})')
                    self.report(f'type(strmixfax, mixreduce)= {type(strmixfac)} {type(self.ctx.mixreduce)}')
                    strmixfac = strmixfac * self.ctx.mixreduce
                    self.ctx.strmix = strmixfac
                    label += f'decreased_mix_fac_str (step {self.ctx.loop_count})'
                    description += f'decreased STRMIX factor by {self.ctx.mixreduce}'
                else:
                    self.report(f'(aggrmixfax, mixreduce)= ({aggrmixfac}, {self.ctx.mixreduce})')
                    self.report(f'type(aggrmixfax, mixreduce)= {type(aggrmixfac)} {type(self.ctx.mixreduce)}')
                    aggrmixfac = aggrmixfac * self.ctx.mixreduce
                    self.ctx.aggrmix = aggrmixfac
                    label += 'decreased_mix_fac_bry'
                    description += f'decreased AGGRMIX factor by {self.ctx.mixreduce}'

            if switch_agressive_mixing:
                last_mixing_scheme = self.ctx.type_aggressive_mixing
                label += ' switched_to_agressive_mixing'
                description += f' switched to agressive mixing scheme (IMIX={last_mixing_scheme})'

            # add number of scf steps, spin
            new_params['SCFSTEPS'] = nsteps
            new_params['NSPIN'] = nspin
            new_params['INS'] = non_spherical

            # add ldos runoption if dos_run = True
            if self.ctx.dos_run:
                if self.ctx.lmdos:
                    runflags = new_params.get('RUNFLAG', []) + ['lmdos']
                else:
                    runflags = new_params.get('RUNFLAG', []) + ['ldos']
                new_params['RUNFLAG'] = runflags
                new_params['SCFSTEPS'] = 1

            # turn on Jij calculation if jij_run == True
            if self.ctx.jij_run:
                new_params['CALCJIJMAT'] = 1

            # add newsosol
            if self.ctx.spinorbit:
                testflags = new_params.get('TESTFLAG', []) + ['tmatnew']
                new_params['TESTFLAG'] = testflags
                new_params['SPINORBIT'] = 1
                new_params['NCOLL'] = 1
                # TODO add deprecation warning and remove these lines (can be set with params_overwrite instead)
                if self.ctx.mesh_params.get('RADIUS_LOGPANELS', None) is not None:
                    new_params['RADIUS_LOGPANELS'] = self.ctx.mesh_params['RADIUS_LOGPANELS']
                if self.ctx.mesh_params.get('NCHEB', None) is not None:
                    new_params['NCHEB'] = self.ctx.mesh_params['NCHEB']
                if self.ctx.mesh_params.get('NPAN_LOG', None) is not None:
                    new_params['NPAN_LOG'] = self.ctx.mesh_params['NPAN_LOG']
                if self.ctx.mesh_params.get('NPAN_EQ', None) is not None:
                    new_params['NPAN_EQ'] = self.ctx.mesh_params['NPAN_EQ']
                new_params['CALCORBITALMOMENT'] = 1
            else:
                new_params['SPINORBIT'] = 0
                new_params['NCOLL'] = 0
                new_params['CALCORBITALMOMENT'] = 0
                new_params['TESTFLAG'] = []

            # set mixing schemes and factors
            if last_mixing_scheme > 2:
                new_params['ITDBRY'] = self.ctx.broyden_num
                new_params['IMIX'] = last_mixing_scheme
                new_params['MIXFAC'] = aggrmixfac
                new_params['NSIMPLEMIXFIRST'] = self.ctx.nsimplemixfirst
            elif last_mixing_scheme == 0:
                new_params['IMIX'] = last_mixing_scheme
                new_params['MIXFAC'] = strmixfac

            # add mixing scheme to context
            self.ctx.last_mixing_scheme = last_mixing_scheme

            if switch_higher_accuracy:
                self.ctx.kkr_higher_accuracy = True

            # add convergence settings
            if self.ctx.loop_count == 1 or self.ctx.last_mixing_scheme == 0:
                new_params['QBOUND'] = self.ctx.threshold_aggressive_mixing
            else:
                new_params['QBOUND'] = self.ctx.convergence_criterion

            # initial magnetization
            if initial_settings and self.ctx.mag_init:
                if self.ctx.hfield[0] <= 0.0 or self.ctx.hfield[1] == 0:
                    self.report(
                        '\nWARNING: magnetization initialization chosen but hfield is zero. Automatically change back to default value (hfield={})\n'
                        .format(self._wf_default['hfield'])
                    )
                    self.ctx.hfield = self._wf_default['hfield']
                new_params['HFIELD'] = self.ctx.hfield
            elif self.ctx.mag_init and self.ctx.mag_init_step_success:  # turn off initialization after first (successful) iteration
                new_params['HFIELD'] = [0.0, 0]
            elif not self.ctx.mag_init:
                self.report("INFO: mag_init is False. Overwrite 'HFIELD' to '0.0' and 'LINIPOL' to 'False'.")
                # reset mag init to avoid resinitializing
                new_params['HFIELD'] = [0.0, 0]

            # set nspin to 2 if mag_init is used
            if self.ctx.mag_init:
                nspin_in = nspin
                if nspin_in is None:
                    nspin_in = 1
                if nspin_in < 2:
                    self.report('WARNING: found NSPIN=1 but for maginit needs NPIN=2. Overwrite this automatically')
                    new_params['NSPIN'] = 2
            message = f'new_params: {new_params}'
            self.report(message)

            # overwrite values from additional input node
            if 'params_overwrite' in self.inputs:
                print('use params_overwrite', self.inputs.params_overwrite.get_dict())
                self._overwrite_parameters_from_input(new_params)

            # step 2.2 update values
            try:
                for key, val in new_params.items():
                    para_check.set_value(key, val, silent=True)
            except:
                message = 'ERROR: failed to set some parameters'
                self.report(message)
                return self.exit_codes.ERROR_PARAMETER_UPDATE  # pylint: disable=no-member

            # step 3:
            message = f'INFO: update parameters to: {para_check.get_set_values()}'
            self.report(message)
            updatenode = Dict(para_check.get_dict())
            updatenode.label = label
            updatenode.description = description
            paranode_new = updatenode  #update_params_wf(self.ctx.last_params, updatenode)
            self.ctx.last_params = paranode_new
        else:
            message = 'INFO: reuse old settings'
            self.report(message)

        message = 'INFO: done updating kkr param step'
        self.report(message)

    def run_kkrimp(self):
        """
        submit a KKR impurity calculation
        """
        message = f'INFO: setting up kkrimp calculation step {self.ctx.loop_count}'
        self.report(message)

        label = f'KKRimp calculation step {self.ctx.loop_count} (IMIX={self.ctx.last_mixing_scheme})'
        description = 'KKRimp calculation of step {}, using mixing scheme {}'.format(
            self.ctx.loop_count, self.ctx.last_mixing_scheme
        )
        code = self.inputs.kkrimp
        params = self.ctx.last_params
        host_GF = self.inputs.remote_data
        imp_pot = self.ctx.last_pot
        last_remote = self.ctx.last_remote
        if 'remote_data_Efshift' in self.inputs:
            host_GF_Efshift = self.inputs.remote_data_Efshift
        else:
            host_GF_Efshift = None

        options = {
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'resources': self.ctx.resources,
            'queue_name': self.ctx.queue
        }
        if self.ctx.custom_scheduler_commands:
            options['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands

        if last_remote is None:
            # make sure no core states are in energy contour
            # extract emin from output of GF host calculation
            GF_out_params = host_GF.get_incoming(link_label_filter='remote_folder'
                                                 ).first().node.outputs.output_parameters
            emin = GF_out_params.get_dict().get('energy_contour_group').get('emin')
            # then use this value to get rid of all core states that are lower than emin (return the same input potential if no states have been removed
            imp_pot = kick_out_corestates_wf(imp_pot, Float(emin))
            if 'impurity_info' in self.inputs:
                message = 'INFO: using impurity_info node as input for kkrimp calculation'
                self.report(message)
                imp_info = self.inputs.impurity_info
                label = 'KKRimp calculation step {} (IMIX={}, Zimp: {})'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme,
                    imp_info.get_dict().get('Zimp')
                )
                description = 'KKRimp calculation of step {}, using mixing scheme {}'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme
                )
                inputs = get_inputs_kkrimp(
                    code,
                    options,
                    label,
                    description,
                    params,
                    not self.ctx.withmpi,
                    imp_info=imp_info,
                    host_GF=host_GF,
                    imp_pot=imp_pot,
                    host_GF_Efshift=host_GF_Efshift
                )
            else:
                message = 'INFO: getting impurity_info node from previous GF calculation'
                self.report(message)
                label = 'KKRimp calculation step {} (IMIX={}, GF_remote: {})'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme, host_GF.pk
                )
                description = 'KKRimp calculation of step {}, using mixing scheme {}'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme
                )
                inputs = get_inputs_kkrimp(
                    code,
                    options,
                    label,
                    description,
                    params,
                    not self.ctx.withmpi,
                    host_GF=host_GF,
                    imp_pot=imp_pot,
                    host_GF_Efshift=host_GF_Efshift
                )
        elif last_remote is not None:
            # fix to get Zimp properly
            if 'impurity_info' in self.inputs:
                message = 'INFO: using RemoteData from previous kkrimp calculation and impurity_info node as input'
                self.report(message)
                imp_info = self.inputs.impurity_info
                label = 'KKRimp calculation step {} (IMIX={}, Zimp: {})'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme,
                    imp_info.get_dict().get('Zimp')
                )
                description = 'KKRimp calculation of step {}, using mixing scheme {}'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme
                )
                inputs = get_inputs_kkrimp(
                    code,
                    options,
                    label,
                    description,
                    params,
                    not self.ctx.withmpi,
                    imp_info=imp_info,
                    host_GF=host_GF,
                    kkrimp_remote=last_remote,
                    host_GF_Efshift=host_GF_Efshift
                )
            else:
                message = 'INFO: using RemoteData from previous kkrimp calculation as input'
                self.report(message)
                label = 'KKRimp calculation step {} (IMIX={}, Zimp: {})'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme, None
                )
                description = 'KKRimp calculation of step {}, using mixing scheme {}'.format(
                    self.ctx.loop_count, self.ctx.last_mixing_scheme
                )
                inputs = get_inputs_kkrimp(
                    code,
                    options,
                    label,
                    description,
                    params,
                    not self.ctx.withmpi,
                    host_GF=host_GF,
                    kkrimp_remote=last_remote,
                    host_GF_Efshift=host_GF_Efshift
                )

        # add LDA+U input node if it was set in parent calculation of last kkrimp_remote or from input port
        if self.ctx.settings_LDAU is not None:
            inputs['settings_LDAU'] = self.ctx.settings_LDAU  # pylint: disable=possibly-used-before-assignment

        # set nonco angles if given
        if 'initial_noco_angles' in self.inputs:
            inputs['initial_noco_angles'] = self.inputs.initial_noco_angles
        if 'rimpshift' in self.inputs:
            inputs['rimpshift'] = self.inputs.rimpshift

        # run the KKR calculation
        message = 'INFO: doing calculation'
        self.report(message)
        kkrimp_run = self.submit(KkrimpCalculation, **inputs)
        print('caching_info KKRimpCalc:', kkrimp_run.get_cache_source())
        print('hash: ', kkrimp_run.get_hash())
        print('_get_objects_to_hash: ', kkrimp_run._get_objects_to_hash())

        return ToContext(kkr=kkrimp_run, last_calc=kkrimp_run)

    def inspect_kkrimp(self):
        """
        check for convergence and store some of the results of the last calculation to context
        """

        self.ctx.calcs.append(self.ctx.last_calc)
        self.ctx.kkrimp_step_success = True

        # check calculation state
        if not self.ctx.last_calc.is_finished_ok:
            self.ctx.kkrimp_step_success = False
            message = 'ERROR: last calc not finished_ok'
            self.report(message)
            return self.exit_codes.ERROR_LAST_CALC_NOT_FINISHED_OK  # pylint: disable=no-member

        message = f'INFO: kkrimp_step_success: {self.ctx.kkrimp_step_success}'
        self.report(message)

        # get potential from last calculation
        try:
            retrieved_folder = self.ctx.kkr.outputs.retrieved
            imp_pot_sfd = extract_imp_pot_sfd(retrieved_folder)
            self.ctx.last_pot = imp_pot_sfd
            print('use potfile sfd:', self.ctx.last_pot)
        except:
            message = 'ERROR: no output potential found'
            self.report(message)
            return self.exit_codes.ERROR_NO_OUTPUT_POT_FROM_LAST_CALC  # pylint: disable=no-member

        # extract convergence info about rms etc. (used to determine convergence behavior)
        try:
            message = f'INFO: trying to find output of last_calc: {self.ctx.last_calc}'
            self.report(message)
            last_calc_output = self.ctx.last_calc.outputs.output_parameters.get_dict()
            found_last_calc_output = True
        except:
            found_last_calc_output = False
        message = f'INFO: found_last_calc_output: {found_last_calc_output}'
        self.report(message)

        # try to extract remote folder
        try:
            self.ctx.last_remote = self.ctx.kkr.outputs.remote_folder
        except:
            self.ctx.last_remote = None
            self.ctx.kkrimp_step_success = False

        message = f'INFO: last_remote: {self.ctx.last_remote}'
        self.report(message)

        if self.ctx.kkrimp_step_success and found_last_calc_output:
            # check convergence
            self.ctx.kkr_converged = last_calc_output['convergence_group']['calculation_converged']
            # check rms
            self.ctx.rms.append(last_calc_output['convergence_group']['rms'])
            rms_all_iter_last_calc = list(last_calc_output['convergence_group']['rms_all_iterations'])
            # check rms of LDAU pot (if LDAU is set)
            try:
                rms_LDAU = last_calc_output['convergence_group']['rms_LDAU']
                self.ctx.rms_LDAU = rms_LDAU
                if rms_LDAU != 0.0:
                    rms_LDAU_all_iter_last_calc = list(last_calc_output['convergence_group']['rms_LDAU_all_iterations'])
                    self.ctx.last_rms_LDAU_all = rms_LDAU_all_iter_last_calc
            except:
                pass

            # add lists of last iterations
            self.ctx.last_rms_all = rms_all_iter_last_calc
            if self.ctx.kkrimp_step_success and self.convergence_on_track():
                self.ctx.rms_all_steps += rms_all_iter_last_calc
        else:
            self.ctx.kkr_converged = False

        message = f'INFO: kkr_converged: {self.ctx.kkr_converged}'
        self.report(message)
        message = f'INFO: rms: {self.ctx.rms}'
        self.report(message)
        message = f'INFO: last_rms_all: {self.ctx.last_rms_all}'
        self.report(message)

        # turn off initial magnetization once one step was successful (update_kkr_params) used in
        if self.ctx.mag_init and self.convergence_on_track():  # and self.ctx.kkrimp_step_success:
            self.ctx.mag_init_step_success = True
        else:
            self.ctx.mag_init_step_success = False

        # store some statistics used to print table in the end of the report
        tmplist = self.ctx.KKR_steps_stats.get('success', [])
        message = f'INFO: append kkr_step_success {tmplist}, {self.ctx.kkr_step_success}'
        self.report(message)
        tmplist.append(self.ctx.kkr_step_success)
        self.ctx.KKR_steps_stats['success'] = tmplist
        try:
            isteps = self.ctx.last_calc.outputs.output_parameters.get_dict(
            )['convergence_group']['number_of_iterations']
        except:
            self.ctx.warnings.append('cound not set isteps in KKR_steps_stats dict')
            isteps = -1

        try:
            first_rms = self.ctx.last_rms_all[0]
            last_rms = self.ctx.last_rms_all[-1]
        except:
            self.ctx.warnings.append('cound not set first_rms, last_rms in KKR_steps_stats dict')
            first_rms = -1
            last_rms = -1

        if self.ctx.last_mixing_scheme == 0:
            mixfac = self.ctx.strmix
        elif self.ctx.last_mixing_scheme > 2:
            mixfac = self.ctx.aggrmix

        if self.ctx.kkr_higher_accuracy:
            qbound = self.ctx.convergence_criterion
        else:
            qbound = self.ctx.threshold_aggressive_mixing

        # store some values in self.ctx.KKR_steps_stats
        for name, val in {
            'isteps': isteps,
            'imix': self.ctx.last_mixing_scheme,
            'mixfac': mixfac,  # pylint: disable=possibly-used-before-assignment
            'qbound': qbound,
            'high_sett': self.ctx.kkr_higher_accuracy,
            'first_rms': first_rms,
            'last_rms': last_rms,
            'pk': self.ctx.last_calc.pk,
            'uuid': self.ctx.last_calc.uuid
        }.items():
            tmplist = self.ctx.KKR_steps_stats.get(name, [])
            tmplist.append(val)
            self.ctx.KKR_steps_stats[name] = tmplist

        message = 'INFO: done inspecting kkrimp results step'
        self.report(message)

    def convergence_on_track(self):
        """
        Check if convergence behavior of the last calculation is on track (i.e. going down)
        """

        on_track = True
        threshold = 5.  # used to check condition if at least one of charnge_neutrality, rms-error goes down fast enough

        # first check if previous calculation was stopped due to reaching the QBOUND limit
        try:
            calc_reached_qbound = self.ctx.last_calc.outputs.output_parameters.get_dict(
            )['convergence_group']['calculation_converged']
        except AttributeError:  # captures error when last_calc dies not have an output node
            calc_reached_qbound = False
        except KeyError:  # captures
            calc_reached_qbound = False

        if self.ctx.kkrimp_step_success and not calc_reached_qbound:
            first_rms = self.ctx.last_rms_all[0]
            # skip first if this is the initial LDA+U iteration because there we see the original non-LDAU convergence value
            if 'settings_LDAU' in self.inputs and self.ctx.loop_count < 2 and len(self.ctx.last_rms_all) > 1:
                first_rms = self.ctx.last_rms_all[1]
            last_rms = self.ctx.last_rms_all[-1]
            # use this trick to avoid division by zero
            if last_rms == 0:
                last_rms = 10**-16
            r = last_rms / first_rms
            message = f'INFO: convergence check: first/last rms {first_rms}, {last_rms}'
            self.report(message)
            if r < 1:
                message = 'INFO: convergence check: rms goes down'
                self.report(message)
                on_track = True
            elif r > threshold:
                message = 'INFO: convergence check: rms goes up too fast, convergence is not expected'
                self.report(message)
                on_track = False
                if self.ctx.rms_LDAU != 0:
                    if self.ctx.last_rms_LDAU_all[-1] < self.ctx.last_rms_LDAU_all[0]:
                        message = 'INFO: convergence check: rms LDAU potential goes down, convergence could still be expected in the next KkrimpCalulation run'
                        self.report(message)
                        on_track = True
            elif len(self.ctx.last_rms_all) == 1:
                message = 'INFO: convergence check: already converged after single iteration'
                self.report(message)
                on_track = True
            else:
                message = 'INFO: convergence check: rms does not shrink fast enough, convergence is not expected'
                self.report(message)
                on_track = False
                if self.ctx.rms_LDAU != 0:
                    if self.ctx.last_rms_LDAU_all[-1] < self.ctx.last_rms_LDAU_all[0]:
                        message = 'INFO: convergence check: rms LDAU potential goes down, convergence could still be expected in the next KkrimpCalulation run'
                        self.report(message)
                        on_track = True
        elif calc_reached_qbound:
            message = 'INFO: convergence check: calculation reached QBOUND'
            self.report(message)
            on_track = True
        else:
            message = 'INFO: convergence check: calculation unsuccessful'
            self.report(message)
            on_track = False

        message = f'INFO: convergence check result: {on_track}'
        self.report(message)

        return on_track

    def return_results(self):
        """
        Return the results of the calculations
        This should run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """

        message = 'INFO: entering return_results'
        self.report(message)

        # try/except to capture as mnuch as possible (everything that is there even when workflow exits unsuccessfully)
        # capture pk and uuids of last calc, params and remote
        try:
            last_calc_uuid = self.ctx.last_calc.uuid
            last_calc_pk = self.ctx.last_calc.pk
            last_params_uuid = self.ctx.last_params.uuid
            last_params_pk = self.ctx.last_params.pk
            last_remote_uuid = self.ctx.last_remote.uuid
            last_remote_pk = self.ctx.last_remote.pk
        except:
            last_calc_uuid = None
            last_calc_pk = None
            last_params_uuid = None
            last_params_pk = None
            last_remote_uuid = None
            last_remote_pk = None

        all_pks = []
        for calc in self.ctx.calcs:
            try:
                all_pks.append(calc.pk)
            except:
                self.ctx.warnings.append(f'cound not get pk of calc {calc}')

        # capture links to last parameter, calcualtion and output
        try:
            last_calc_out = self.ctx.kkr.out['output_parameters']
            last_calc_out_dict = last_calc_out.get_dict()
            last_RemoteData = self.ctx.last_remote
            last_InputParameters = self.ctx.last_params
        except:
            last_InputParameters = None
            last_RemoteData = None
            last_calc_out = None
            last_calc_out_dict = {}

        # capture convergence info
        try:
            last_rms = self.ctx.rms[-1]
        except:
            last_rms = None

        # now collect results saved in results node of workflow
        message = 'INFO: collect outputnode_dict'
        self.report(message)
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['material'] = self.ctx.formula
        outputnode_dict['loop_count'] = self.ctx.loop_count
        outputnode_dict['warnings'] = self.ctx.warnings
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['last_params_nodeinfo'] = {'uuid': last_params_uuid, 'pk': last_params_pk}
        outputnode_dict['last_remote_nodeinfo'] = {'uuid': last_remote_uuid, 'pk': last_remote_pk}
        outputnode_dict['last_calc_nodeinfo'] = {'uuid': last_calc_uuid, 'pk': last_calc_pk}
        outputnode_dict['pks_all_calcs'] = all_pks
        outputnode_dict['convergence_value'] = last_rms
        outputnode_dict['convergence_values_all_steps'] = array(self.ctx.rms_all_steps)
        outputnode_dict['convergence_values_last_step'] = array(self.ctx.last_rms_all)
        outputnode_dict['convergence_reached'] = self.ctx.kkr_converged
        outputnode_dict['kkr_step_success'] = self.ctx.kkr_step_success
        outputnode_dict['used_higher_accuracy'] = self.ctx.kkr_higher_accuracy

        # report the status
        if self.ctx.successful:
            self.report(
                'STATUS: Done, the convergence criteria are reached.\n'
                'INFO: The charge density of the KKR calculation pk= {} '
                'converged after {} KKR runs and {} iterations to {} \n'
                ''.format(
                    last_calc_pk, self.ctx.loop_count - 1, sum(self.ctx.KKR_steps_stats.get('isteps', [])),
                    self.ctx.last_rms_all[-1]
                )
            )
        else:  # Termination ok, but not converged yet...
            self.report(
                'STATUS/WARNING: Done, the maximum number of runs '
                'was reached or something failed.\n INFO: The '
                'charge density of the KKR calculation pk= '
                'after {} KKR runs and {} iterations is {} "me/bohr^3"\n'
                ''.format(
                    self.ctx.loop_count - 1, sum(self.ctx.KKR_steps_stats.get('isteps', [])), self.ctx.last_rms_all[-1]
                )
            )

        # create results node and link all calculations
        message = 'INFO: create results nodes'
        self.report(message)
        link_nodes = {}
        icalc = 0
        for calc in self.ctx.calcs:
            link_nodes[f'KkrimpCalc{icalc}'] = calc.outputs.remote_folder
            icalc += 1
        if not self.ctx.dos_run:
            link_nodes['final_imp_potential'] = self.ctx.last_pot
        outputnode_t = create_out_dict_node(Dict(dict=outputnode_dict), **link_nodes)
        outputnode_t.label = 'kkr_scf_wc_results'
        outputnode_t.description = 'Contains results of workflow (e.g. workflow version number, info about success of wf, lis tof warnings that occured during execution, ...)'

        self.out('workflow_info', outputnode_t)
        # store out_potential as SingleFileData only if this was no DOS run
        if not self.ctx.dos_run:
            self.out('host_imp_pot', self.ctx.last_pot)

        # print results table for overview
        # table layout:
        message = 'INFO: overview of the result:\n\n'
        message += '|------|---------|--------|------|--------|---------|-----------------|---------------------------------------------|\n'
        message += '| irun | success | isteps | imix | mixfac | qbound  |       rms       |                pk and uuid                  |\n'
        message += '|      |         |        |      |        |         | first  |  last  |                                             |\n'
        message += '|------|---------|--------|------|--------|---------|--------|--------|---------------------------------------------|\n'
        KKR_steps_stats = self.ctx.KKR_steps_stats
        for irun in range(len(KKR_steps_stats.get('success', []))):
            message += '|%6i|%9s|%8i|%6i|%.2e|%.3e|%.2e|%.2e|' % (
                irun + 1, KKR_steps_stats.get('success')[irun], KKR_steps_stats.get('isteps')[irun],
                KKR_steps_stats.get('imix')[irun], KKR_steps_stats.get('mixfac')[irun],
                KKR_steps_stats.get('qbound')[irun], KKR_steps_stats.get('first_rms')[irun],
                KKR_steps_stats.get('last_rms')[irun]
            )
            message += f" {KKR_steps_stats.get('pk')[irun]} | {KKR_steps_stats.get('uuid')[irun]}|\n"
            message += '|------|---------|--------|------|--------|---------|-----------------|---------------------------------------------|\n'
            """
            message += "#|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".format(irun+1,
                          KKR_steps_stats.get('success')[irun], KKR_steps_stats.get('isteps')[irun],
                          KKR_steps_stats.get('imix')[irun], KKR_steps_stats.get('mixfac')[irun],
                          KKR_steps_stats.get('qbound')[irun],
                          KKR_steps_stats.get('first_rms')[irun], KKR_steps_stats.get('last_rms')[irun])
            """
        self.report(message)

        self.report('INFO: done with kkr_scf workflow!\n')

    def error_handler(self):
        """Capture errors raised in validate_input"""
        if self.ctx.exit_code is not None:
            return self.ctx.exit_code

    def _overwrite_parameters_from_input(self, new_params):
        """Overwrite input parameters for KKRimpCalculation if found in input"""
        params_overwrite = self.inputs.params_overwrite.get_dict()
        for key, val in params_overwrite.items():
            if key in new_params:
                self.report('ATTENTION: overwriting parameter from "params_overwrite" input')
                self.report(f'key: {key}')
                self.report(f'old value: {new_params[key]}')
                self.report(f'overwritten value: {val}')
            new_params[key] = val


@calcfunction
def extract_imp_pot_sfd(retrieved_folder):
    """
    Extract potential file from retrieved folder and save as SingleFileData
    """
    # take output potential file either from tarfile or directy from output folder

    if KkrimpCalculation._FILENAME_TAR in retrieved_folder.list_object_names():
        print('take potfile from tar file of retrieved')
        # take potfile after extracting tar file
        # get full filename
        with retrieved_folder.open(KkrimpCalculation._FILENAME_TAR) as tar_file:
            tarfilename = tar_file.name
        print('tarfile name:', tarfilename)
        # open tarfile and extract potfile
        with tarfile.open(tarfilename) as tar_file:
            print('extract potfile:', KkrimpCalculation._OUT_POTENTIAL)
            tar_file.extract(KkrimpCalculation._OUT_POTENTIAL, os.path.dirname(tarfilename))
            with retrieved_folder.open(KkrimpCalculation._OUT_POTENTIAL, 'rb') as pot_file:
                print('get potfile sfd:', pot_file)
                imp_pot_sfd = SinglefileData(file=pot_file)

        # delete extracted potfile again
        print('delete potfile from outfile:', KkrimpCalculation._OUT_POTENTIAL)
        retrieved_folder.delete_object(KkrimpCalculation._OUT_POTENTIAL, force=True)
    else:
        # take potfile directly from output
        with retrieved_folder.open(KkrimpCalculation._OUT_POTENTIAL, 'rb') as pot_file:
            imp_pot_sfd = SinglefileData(file=pot_file)

    return imp_pot_sfd
