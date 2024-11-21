#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for converging a kkr calculation and
some helper methods to do so with AiiDA
"""
import numpy as np
from masci_tools.io.kkr_params import kkrparams
from masci_tools.io.common_functions import get_Ry2eV, get_ef_from_potfile
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistentAttributeError
from aiida.engine import (
    WorkChain,
    while_,
    if_,
    ToContext,
    workfunction,
    calcfunction,
)
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.extract_kkrhost_noco_angles import extract_noco_angles
from aiida_kkr.tools.common_workfunctions import (
    test_and_get_codenode,
    get_inputs_kkr,
    get_parent_paranode,
    update_params_wf,
)
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.dos import kkr_dos_wc

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.11.2'
__contributors__ = (u'Jens Broeder', u'Philipp Rüßmann', u'David Antognini Silva')

eV2Ry = 1.0 / get_Ry2eV()

# TODO: magnetism (init and converge magnetic state)
# TODO: check convergence (RMAX, GMAX etc.)
# TODO: save timing info of the steps
# TODO: switch to LLOYD
# TODO: emin-emax setting
# TODO: restart from workflow output instead of calculation output
# TODO: add warnings
# TODO: maybe define the energy point density instead of a fixed number as in input?
# TODO: overwrite defaults from parent if parent is previous kkr_scf run
# TODO: retrieve DOS within scf run


class kkr_scf_wc(WorkChain):
    """
    Workchain for converging a KKR calculation (SCF).

    It converges the charge potential.
    Two paths are possible:

    (1) Start from a structure and run a voronoi calculation first,
    optional with calc_parameters
    (2) Start from an existing Voronoi or KKR calculation, with a remoteData

    :param wf_parameters: (Dict), Workchain Specifications
    :param options: (Dict); specifications for the computer
    :param structure: (StructureData), Crystal structure
    :param calc_parameters: (Dict), Voronoi/Kkr Parameters
    :param remote_data: (RemoteData), from a KKR, or Voronoi calculation
    :param voronoi: (Code)
    :param kkr: (Code)

    :return output_kkr_scf_wc_para: (Dict), Information of workflow results
        like Success, last result node, list with convergence behavior

    minimum input example:
    1. Code1, Code2, Structure, (Parameters), (wf_parameters)
    2. Code2, remote_data, (Parameters), (wf_parameters)

    maximum input example:
    1. Code1, Code2, Structure, Parameters
        wf_parameters: {'queue_name' : String,
                        'resources' : dict({"num_machines": int, "num_mpiprocs_per_machine" : int})
                        'walltime' : int}
    2. Code2, (remote-data), wf_parameters as in 1.

    Hints:
    1. This workflow does not work with local codes!
    """

    _workflowversion = __version__

    # default workflow parameters
    _wf_default = AttributeDict()
    # Maximum number of kkr jobs/starts (defauld iterations per start)
    _wf_default.kkr_runmax = 5
    # Stop if charge density is converged below this value
    _wf_default.convergence_criterion = 10**-8
    # reduce mixing factor by this factor if calculation fails due to too large mixing
    _wf_default.mixreduce = 0.5
    # temperature increasing steps if workflow fails
    _wf_default.tempr_increase = 50.
    # threshold after which agressive mixing is used
    _wf_default.threshold_aggressive_mixing = 8 * 10**-3
    _wf_default.strmix = 0.03  # mixing factor of simple mixing
    _wf_default.brymix = 0.05  # mixing factor of aggressive mixing
    # number of iterations done per KKR calculation
    _wf_default.nsteps = 50
    # wether or not pre-convergnce on a coarse energy and k-point grid should be done
    _wf_default.coarse_preconvergence = True
    # convergence settings for coarse preconvergence
    _wf_default.convergence_setting_coarse = AttributeDict()  # setting of the coarse preconvergence
    _wf_default.convergence_setting_coarse.npol = 7
    _wf_default.convergence_setting_coarse.n1 = 3
    _wf_default.convergence_setting_coarse.n2 = 11
    _wf_default.convergence_setting_coarse.n3 = 3
    _wf_default.convergence_setting_coarse.tempr = 1000.0
    _wf_default.convergence_setting_coarse.kmesh = [10, 10, 10]
    # threshold after which final conversion settings are used
    _wf_default.threshold_switch_high_accuracy = 10**-3
    # setting of the final convergence (lower tempr, 48 epts, denser k-mesh)
    _wf_default.convergence_setting_fine = AttributeDict()
    _wf_default.convergence_setting_fine.npol = 5
    _wf_default.convergence_setting_fine.n1 = 7
    _wf_default.convergence_setting_fine.n2 = 29
    _wf_default.convergence_setting_fine.n3 = 7
    _wf_default.convergence_setting_fine.tempr = 600.0
    _wf_default.convergence_setting_fine.kmesh = [30, 30, 30]
    # initialize and converge magnetic calculation
    _wf_default.mag_init = False
    # external magnetic field used in initialization step in Ry
    _wf_default.hfield = 0.02
    # position in unit cell where magnetic field is applied [default (None) means apply to all]
    _wf_default.init_pos = None
    # fix direction of magnetic moment if the direction changes less
    # than this value in degrees (calculated as sqrt((delta theta)**2 + (delta phi)**2))
    _wf_default.fix_dir_threshold = 1.0
    # add DOS to testopts and retrieve dos.atom files in each scf run
    _wf_default.retreive_dos_data_scf_run = False
    # minimal distance of start of the energy contour to highest lying core state in Ry
    _wf_default.delta_e_min_core_states = 0.2  # Ry
    # set these keys from defaults in kkr_startpot workflow since they are only passed onto that workflow
    for key, value in kkr_startpot_wc.get_wf_defaults(silent=True).items():
        if key in [
            'dos_params',
            'fac_cls_increase',
            'natom_in_cls_min',
            'delta_e_min',
            'threshold_dos_zero',
            'check_dos',
        ]:
            _wf_default[key] = value

    # default options
    _options_default = AttributeDict()
    _options_default.queue_name = ''  # Queue name to submit jobs too
    # resources to allowcate for the job
    _options_default.resources = AttributeDict()
    _options_default.resources.num_machines = 1
    # walltime after which the job gets killed (gets parsed to KKR)
    _options_default.max_wallclock_seconds = 60 * 60
    _options_default.withmpi = True  # execute KKR with mpi or without
    _options_default.custom_scheduler_commands = ''  # some additional scheduler commands

    # intended to guide user interactively in setting up a valid wf_params node

    @classmethod
    def get_wf_defaults(cls, silent=False):
        """Print and return _wf_default dictionary.

        Can be used to easily create set of wf_parameters.
        returns _wf_default, _options_default
        """
        if not silent:
            print(f'Version of workflow: {cls._workflowversion}')
        return cls._wf_default, cls._options_default.copy()

    @classmethod
    def define(cls, spec):
        """Defines the outline of the workflow."""
        # Take input of the workflow or use defaults defined above
        super(kkr_scf_wc, cls).define(spec)
        spec.input(
            'wf_parameters',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
            help="""
            Settings for the workflow. Use `KkrCalculation.get_wf_defaults()`
            to get the default values and default options.
            """
        )
        spec.input(
            'options',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
            help="""
            Computer settings used by the calculations in the workflow
            (see also help string of wf_parameters).
            """
        )
        spec.input(
            'structure',
            valid_type=orm.StructureData,
            required=False,
            help="""
            Input structure for which a calculation is started with a
            VoronoiCalculation.
            Can be skipped if a previous KkrCalculation is given with the
            `remote_data` input node.
            """
        )
        spec.input(
            'calc_parameters',
            valid_type=orm.Dict,
            required=False,
            help="""
            KKR-specific calculation parameters (LMAX etc.),
            usually set up with the help of the `kkrparams` class.
            """
        )
        spec.input(
            'remote_data',
            valid_type=orm.RemoteData,
            required=False,
            help="""
            RemoteFolder node of a preconverged calculation.
            Can be used as a starting point to skip the Voronoi step.
            """
        )
        spec.input(
            'voronoi',
            valid_type=orm.Code,
            required=False,
            help="""
            Voronoi code node, needed only if `structure` input node is given.
            """
        )
        spec.input(
            'kkr',
            valid_type=orm.Code,
            required=True,
            help='KKRhost code node which will run the KkrCalculations',
        )
        spec.input(
            'startpot_overwrite',
            valid_type=orm.SinglefileData,
            required=False,
            help="""
            Potential SinglefileData, can be used to overwrite
            the starting potential from Voronoi
            (the shapefun will be used though and thus needs to be compatible).
            This can be used to construct a better starting potential from a
            preconverged calculation (e.g. in a smaller unit cell).
            """
        )
        spec.input(
            'initial_noco_angles',
            valid_type=orm.Dict,
            required=False,
            help="""
            Initial non-collinear angles for the magnetic moments of the
            impurities. See KkrCalculation for details.
            """
        )
        spec.input(
            'params_kkr_overwrite',
            valid_type=orm.Dict,
            required=False,
            help='Set some input parameters of the KKR calculation.'
        )

        # define output nodes
        spec.output(
            'output_kkr_scf_wc_ParameterResults',
            valid_type=orm.Dict,
            required=True,
        )
        spec.output(
            'last_calc_out',
            valid_type=orm.Dict,
            required=False,
        )
        spec.output(
            'last_RemoteData',
            valid_type=orm.RemoteData,
            required=False,
        )
        spec.output(
            'last_InputParameters',
            valid_type=orm.Dict,
            required=False,
        )
        spec.output(
            'results_vorostart',
            valid_type=orm.Dict,
            required=False,
        )
        spec.output(
            'starting_dosdata_interpol',
            valid_type=orm.XyData,
            required=False,
        )
        spec.output(
            'final_dosdata_interpol',
            valid_type=orm.XyData,
            required=False,
        )
        spec.output(
            'last_noco_angles',
            valid_type=orm.Dict,
            required=False,
        )

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            # check if voronoi run needed, otherwise skip this step
            if_(cls.validate_input)(
                # run kkr_startpot workflow (sets up voronoi input, runs voro calc, does some consistency checks)
                cls.run_voronoi,
                # check output of run_voronoi and determine if calculation has to be terminated here
                cls.check_voronoi
            ),
            # while loop for KKR run(s), first simple mixing
            # then Anderson with corase pre-convergence settings
            # finally convergence step with higher accuracy
            while_(cls.condition)(
                # update parameters for kkr step using previous output(s)
                cls.update_kkr_params,
                # run kkr step
                # TODO: encapsulate this in restarting mechanism
                # (should be a base class of workflows that start calculations)
                # i.e. use base_restart_calc workchain as parent
                cls.run_kkr,
                # check results for convergence and collect some intermediate results
                cls.inspect_kkr
            ),
            # compute final dos if check_dos is True
            cls.get_dos,
            cls.check_dos,
            # finalize calculation and create output nodes
            cls.return_results
        )

        # definition of exit codes if the workflow needs to be terminated
        spec.exit_code(
            221,
            'ERROR_NO_PARENT_PARAMS_FOUND',
            message='Unable to extract parent paremeter node of input remote folder',
        )
        spec.exit_code(
            222,
            'ERROR_INVALID_KKR_CODE',
            message='The code you provided for kkr does not use the plugin kkr.kkr',
        )
        spec.exit_code(
            223,
            'ERROR_INVALID_VORONOI_CODE',
            message='The code you provided for voronoi does not use the plugin kkr.voro',
        )
        spec.exit_code(
            224,
            'ERROR_NO_VORONOI_CODE_GIVEN',
            message='ERROR: StructureData was provided, but no voronoi code was provided',
        )
        spec.exit_code(
            225,
            'ERROR_NOT_ENOUGH_INPUTS',
            message='ERROR: No StructureData nor remote_data was provided as Input',
        )
        spec.exit_code(
            226,
            'ERROR_KKR_STARTPOT_FAILED',
            message='ERROR: kkr_startpot_wc step failed!',
        )
        spec.exit_code(
            227,
            'ERROR_DOS_RUN_UNSUCCESSFUL',
            message='DOS run unsuccessful. Check inputs.',
        )
        spec.exit_code(
            228,
            'ERROR_CALC_PARAMETERS_INCOMPLETE',
            message='ERROR: calc_parameters given are not consistent! Missing mandatory keys',
        )
        spec.exit_code(
            229,
            'ERROR_CALC_PARAMTERS_INCONSISTENT',
            message='ERROR: calc_parameters given are not consistent! Hint: did you give an unknown keyword?',
        )
        spec.exit_code(
            230,
            'ERROR_NO_CALC_PARAMETERS_GIVEN',
            message='ERROR: calc_parameters not given as input but are needed!',
        )
        spec.exit_code(
            231,
            'ERROR_PARAM_UPDATE_FAILED',
            message='ERROR: parameter update unsuccessful: some key, value pair not valid!',
        )
        spec.exit_code(
            232,
            'ERROR_CALC_PARAMTERS_INCOMPLETE',
            message='ERROR: calc_parameters misses keys',
        )
        spec.exit_code(
            233,
            'ERROR_LAST_REMOTE_NOT_FOUND',
            message='ERROR: last_remote could not be set to a previous successful calculation',
        )
        spec.exit_code(
            234,
            'ERROR_MAX_KKR_RESTARTS_REACHED',
            message='ERROR: maximal number of KKR restarts reached. Exiting now!',
        )
        spec.exit_code(
            235,
            'ERROR_CALC_SUBMISSION_FAILED',
            message='ERROR: last KKRcalc in SUBMISSIONFAILED state',
        )

    def start(self):
        """Initialize context and some parameters."""
        self.report(f'INFO: started KKR convergence workflow version {self._workflowversion}')

        ####### init #######

        # internal para /control para
        self.ctx.loop_count = 0
        self.ctx.last_mixing_scheme = 0
        self.ctx.calcs = []
        self.ctx.abort = False
        # flags used internally to check whether the individual steps were successful
        self.ctx.dos_ok = False
        self.ctx.voro_step_success = False
        self.ctx.kkr_step_success = False
        self.ctx.kkr_converged = False
        self.ctx.kkr_higher_accuracy = False
        # links to previous calculations
        self.ctx.last_calc = None
        self.ctx.last_params = None
        self.ctx.last_remote = None
        # convergence info about rms etc. (used to determine convergence behavior)
        self.ctx.last_rms_all = []
        self.ctx.rms_all_steps = []
        self.ctx.last_neutr_all = []
        self.ctx.neutr_all_steps = []

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()
        options_dict = self.inputs.options.get_dict()

        if wf_dict == {}:
            wf_dict = self._wf_default
            self.report('INFO: using default wf parameter')
        if options_dict == {}:
            options_dict = self._options_default
            self.report('INFO: using default options')

        # set values from input, or defaults
        self.ctx.withmpi = options_dict.get(
            'withmpi',
            self._options_default['withmpi'],
        )
        self.ctx.resources = options_dict.get(
            'resources',
            self._options_default['resources'],
        )
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds',
            self._options_default['max_wallclock_seconds'],
        )
        self.ctx.queue = options_dict.get(
            'queue_name',
            self._options_default['queue_name'],
        )
        self.ctx.custom_scheduler_commands = options_dict.get(
            'custom_scheduler_commands',
            self._options_default['custom_scheduler_commands'],
        )
        self.ctx.options_params_dict = orm.Dict(
            dict={
                'withmpi': self.ctx.withmpi,
                'resources': self.ctx.resources,
                'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
                'queue_name': self.ctx.queue,
                'custom_scheduler_commands': self.ctx.custom_scheduler_commands
            }
        )

        # set label and description
        self.ctx.description_wf = self.inputs.get(
            'description', 'Workflow for '
            'a KKR scf calculation starting '
            'either from a structure with '
            'automatic voronoi calculation '
            'or a valid RemoteData node of '
            'a previous calculation'
        )
        self.ctx.label_wf = self.inputs.get('label', 'kkr_scf_wc')

        # set workflow parameters
        self.ctx.max_number_runs = wf_dict.get(
            'kkr_runmax',
            self._wf_default['kkr_runmax'],
        )
        self.ctx.strmix = wf_dict.get(
            'strmix',
            self._wf_default['strmix'],
        )
        self.ctx.brymix = wf_dict.get(
            'brymix',
            self._wf_default['brymix'],
        )
        self.ctx.check_dos = wf_dict.get(
            'check_dos',
            self._wf_default['check_dos'],
        )
        self.ctx.dos_params = wf_dict.get(
            'dos_params',
            self._wf_default['dos_params'],
        )
        self.ctx.convergence_criterion = wf_dict.get(
            'convergence_criterion',
            self._wf_default['convergence_criterion'],
        )
        self.ctx.coarse_preconvergence = wf_dict.get(
            'coarse_preconvergence',
            self._wf_default['coarse_preconvergence'],
        )
        self.ctx.convergence_setting_coarse = wf_dict.get(
            'convergence_setting_coarse',
            self._wf_default['convergence_setting_coarse'],
        )
        self.ctx.convergence_setting_fine = wf_dict.get(
            'convergence_setting_fine',
            self._wf_default['convergence_setting_fine'],
        )
        self.ctx.mixreduce = wf_dict.get(
            'mixreduce',
            self._wf_default['mixreduce'],
        )
        self.ctx.tempr_increase = wf_dict.get(
            'tempr_increase',
            self._wf_default['tempr_increase'],
        )
        self.ctx.nsteps = wf_dict.get(
            'nsteps',
            self._wf_default['nsteps'],
        )
        self.ctx.threshold_aggressive_mixing = wf_dict.get(
            'threshold_aggressive_mixing',
            self._wf_default['threshold_aggressive_mixing'],
        )
        self.ctx.threshold_switch_high_accuracy = wf_dict.get(
            'threshold_switch_high_accuracy',
            self._wf_default['threshold_switch_high_accuracy'],
        )

        # initial magnetization
        self.ctx.mag_init = wf_dict.get(
            'mag_init',
            self._wf_default['mag_init'],
        )
        self.ctx.hfield = wf_dict.get(
            'hfield',
            self._wf_default['hfield'],
        )
        self.ctx.xinit = wf_dict.get(
            'init_pos',
            self._wf_default['init_pos'],
        )
        self.ctx.mag_init_step_success = False

        # difference in eV to emin (e_fermi) if emin (emax) are larger
        # (smaller) than emin (e_fermi)
        self.ctx.delta_e = wf_dict.get(
            'delta_e_min',
            self._wf_default['delta_e_min'],
        )
        # threshold for dos comparison (comparison of dos at emin)
        self.ctx.threshold_dos_zero = wf_dict.get(
            'threshold_dos_zero',
            self._wf_default['threshold_dos_zero'],
        )
        # distance of EMIN to highest core state (used to modify EMIN start)
        self.ctx.delta_e_min_core_states = wf_dict.get(
            'delta_e_min_core_states',
            self._wf_default['delta_e_min_core_states'],
        )

        # this will be used to store the Fermi level
        self.ctx.efermi = None

        # set starting noco angles, gets updated in between the KKR runs if
        # fix_dir == False for all atoms
        if 'initial_noco_angles' in self.inputs:
            self.ctx.initial_noco_angles = self.inputs.initial_noco_angles
            self.ctx.fix_dir_threshold = wf_dict.get(
                'fix_dir_threshold',
                self._wf_default['fix_dir_threshold'],
            )

        # retreive dos data in each scf run
        self.ctx.scf_dosdata = wf_dict.get(
            'retreive_dos_data_scf_run',
            self._wf_default['retreive_dos_data_scf_run'],
        )

        self.report(
            'INFO: use the following parameter:\n'
            '\nGeneral settings\n'
            f'use mpi: {self.ctx.withmpi}\n'
            f'max number of KKR runs: {self.ctx.max_number_runs}\n'
            f'Resources: {self.ctx.resources}\n'
            f'Walltime (s): {self.ctx.max_wallclock_seconds}\n'
            f'queue name: {self.ctx.queue}\n'
            f'scheduler command: {self.ctx.custom_scheduler_commands}\n'
            f'description: {self.ctx.description_wf}\n'
            f'label: {self.ctx.label_wf}\n'
            f'\nMixing parameter\n'
            f'Straight mixing factor: {self.ctx.strmix}\n'
            f'Anderson mixing factor: {self.ctx.brymix}\n'
            f'Nsteps scf cycle: {self.ctx.nsteps}\n'
            f'Convergence criterion: {self.ctx.convergence_criterion}\n'
            f'threshold_aggressive_mixing: {self.ctx.threshold_aggressive_mixing}\n'
            f'threshold_switch_high_accuracy: {self.ctx.threshold_switch_high_accuracy}\n'
            f'use coarse preconvergence: {self.ctx.coarse_preconvergence}\n'
            f'convergence_setting_coarse: {self.ctx.convergence_setting_coarse}\n'
            f'convergence_setting_fine: {self.ctx.convergence_setting_fine}\n'
            f'factor reduced mixing if failing calculation: {self.ctx.mixreduce}\n'
            f'temperature increasing value if failing calculation: {self.ctx.tempr_increase}\n'
            f'\nAdditional parameter\n'
            f'check DOS between runs: {self.ctx.check_dos}\n'
            f'DOS parameters: {self.ctx.dos_params}\n'
            f'init magnetism in first step: {self.ctx.mag_init}\n'
            f'init magnetism, hfield: {self.ctx.hfield}\n'
            f'init magnetism, init_pos: {self.ctx.xinit}\n'
        )

        self.ctx.successful = True
        self.ctx.rms = []
        self.ctx.neutr = []
        self.ctx.warnings = []
        self.ctx.errors = []
        self.ctx.formula = ''

        # for results table each list gets one entry per iteration that has been performed
        self.ctx.KKR_steps_stats = {
            'success': [],
            'isteps': [],
            'imix': [],
            'mixfac': [],
            'qbound': [],
            'high_sett': [],
            'first_rms': [],
            'last_rms': [],
            'first_neutr': [],
            'last_neutr': [],
            'pk': [],
            'uuid': []
        }

    def validate_input(self):
        """Validate input and find out which path (1, or 2) to take

        return True means run voronoi if false run kkr directly
        """
        run_voronoi = True
        inputs = self.inputs

        if 'structure' in inputs:
            self.report('INFO: Found structure in input. Start with Voronoi calculation.')
            if not 'voronoi' in inputs:
                return self.exit_codes.ERROR_NO_VORONOI_CODE_GIVEN  # pylint: disable=no-member
        elif 'remote_data' in inputs:
            self.report('INFO: Found remote_data in input. Continue calculation without running voronoi step.')
            run_voronoi = False
        else:
            return self.exit_codes.ERROR_NOT_ENOUGH_INPUTS  # pylint: disable=no-member

        if 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_INVALID_VORONOI_CODE  # pylint: disable=no-member

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_INVALID_KKR_CODE  # pylint: disable=no-member

        # set params and remote folder to input if voronoi step is skipped
        if not run_voronoi:
            self.ctx.last_remote = inputs.remote_data
            num_parents = len(self.ctx.last_remote.get_incoming(node_class=orm.CalcJobNode).all_link_labels())
            if num_parents == 0:
                pk_last_remote = self.ctx.last_remote.inputs.last_RemoteData.outputs.output_kkr_scf_wc_ParameterResults.get_dict(
                ).get('last_calc_nodeinfo').get('pk')
                last_calc = orm.load_node(pk_last_remote)
                self.ctx.last_remote = last_calc.outputs.remote_folder
            try:  # first try parent of remote data output of a previous calc.
                parent_params = get_parent_paranode(inputs.remote_data)
            except AttributeError:
                try:  # next try to extract parameter from previous kkr_scf_wc output
                    parent_params = inputs.remote_data.inputs.last_RemoteData.inputs.calc_parameters
                except AttributeError:
                    return self.exit_codes.ERROR_NO_PARENT_PARAMS_FOUND  # pylint: disable=no-member
            if 'calc_parameters' in inputs:
                self.ctx.last_params = inputs.calc_parameters
                # TODO: check last_params consistency against parent_params
                #parent_params_dict = parent_params.get_dict()
                # for key, val in self.ctx.last_params.get_dict().iteritems():
                #    if key in parent_params_dict.keys():
                #        if val != parent_params_dict[key]:
                #
                #    else:
                #
            else:
                self.ctx.last_params = parent_params
            self.ctx.voro_step_success = True

        return run_voronoi

    def run_voronoi(self):
        """Run the voronoi step calling voro_start workflow."""

        # collects inputs
        structure = self.inputs.structure
        self.ctx.formula = structure.get_formula()
        voronoicode = self.inputs.voronoi
        kkrcode = self.inputs.kkr

        # set KKR parameters if any are given, otherwise use defaults
        if 'calc_parameters' in self.inputs:
            params = self.inputs.calc_parameters
        else:
            params = None

        # check if default values are missing and set appropriately from defaults
        defaults, version = kkrparams.get_KKRcalc_parameter_defaults()
        if params is None:
            params = orm.Dict(dict=defaults)
            params.label = 'default values'
        newparams = {}
        for key, val in defaults.items():
            # this one is automatically set by kkr_startpot_wc so we skip it here
            if params.get_dict().get(key, None) is None and key != 'RCLUSTZ':
                newparams[key] = val
                self.report(f'INFO: Automatically added default values to KKR parameters: {key} {val}')
        if newparams != {}:
            for key, val in params.get_dict().items():
                if val is not None:
                    newparams[key] = val
            updatenode = orm.Dict(dict=newparams)
            updatenode.label = 'added defaults to KKR input parameter'
            updatenode.description = 'Overwritten KKR input parameter to correct missing default values automatically'
            params = update_params_wf(params, updatenode)

        # set nspin to 2 if mag_init is used
        if self.ctx.mag_init:
            input_dict = params.get_dict()
            para_check = kkrparams()
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)
            nspin_in = para_check.get_value('NSPIN')
            if nspin_in is None:
                nspin_in = 1
            if nspin_in < 2:
                self.report('WARNING: found NSPIN=1 but for maginit needs NPIN=2. Overwrite this automatically')
                para_check.set_value('NSPIN', 2, silent=True)
                self.report(f'INFO: update parameters to: {para_check.get_set_values()}')
                updatenode = orm.Dict(dict=para_check.get_dict())
                updatenode.label = 'overwritten KKR input parameters'
                updatenode.description = 'Overwritten KKR input parameter to correct NSPIN to 2'
                params = update_params_wf(params, updatenode)

        # check consistency of input parameters before running calculation
        self.check_input_params(params, is_voronoi=True)

        # set parameters of voro_start sub workflow
        sub_wf_params_dict = kkr_startpot_wc.get_wf_defaults(silent=True)
        label, description = 'voro_start_default_params', 'workflow parameters for voro_start'
        if 'wf_parameters' in self.inputs:
            wf_params_input = self.inputs.wf_parameters.get_dict()
            num_updated = 0
            for key in list(sub_wf_params_dict.keys()):
                if key in list(wf_params_input.keys()):
                    val = wf_params_input[key]
                    sub_wf_params_dict[key] = val
                    num_updated += 1
            if num_updated > 0:
                label = 'voro_start_updated_params'
        sub_wf_params = orm.Dict(dict=sub_wf_params_dict)
        sub_wf_params.label = label
        sub_wf_params.description = description

        self.report('INFO: run voronoi step')
        self.report(f'INFO: using calc_params ({params})')
        self.report(f'INFO: using wf_parameters ({sub_wf_params})')

        wf_label = 'kkr_startpot (voronoi)'
        wf_desc = 'subworkflow to set up the input of a KKR calculation'
        builder = kkr_startpot_wc.get_builder()
        builder.kkr = kkrcode
        builder.voronoi = voronoicode
        builder.calc_parameters = params
        builder.wf_parameters = sub_wf_params
        builder.structure = structure
        builder.metadata.label = wf_label  # pylint: disable=no-member
        builder.metadata.description = wf_desc  # pylint: disable=no-member
        builder.options = self.ctx.options_params_dict
        if 'startpot_overwrite' in self.inputs:
            builder.startpot_overwrite = self.inputs.startpot_overwrite
        future = self.submit(builder)

        return ToContext(voronoi=future, last_calc=future)

    def check_voronoi(self):
        """Check output of kkr_startpot_wc workflow.

        It checks the starting potential, shapefun etc.
        """

        self.report('INFO: checking voronoi output')
        voro_step_ok = False

        # check some output
        kkrstartpot_results = self.ctx.voronoi.outputs.results_vorostart_wc.get_dict()
        if kkrstartpot_results['successful']:
            voro_step_ok = True
        # initialize last_remote and last_params (gets updated in loop of KKR calculations)
        self.ctx.last_params = self.ctx.voronoi.outputs.last_params_voronoi
        self.ctx.last_remote = self.ctx.voronoi.outputs.last_voronoi_remote

        # store result in context
        self.ctx.voro_step_success = voro_step_ok

        # abort calculation if something failed in voro_start step
        if not voro_step_ok:
            return self.exit_codes.ERROR_KKR_STARTPOT_FAILED  # pylint: disable=no-member

        self.report('INFO: done checking voronoi output')
        return None

    def condition(self):
        """
        check convergence condition
        """
        self.report('INFO: checking condition for kkr step')
        do_kkr_step = True
        stopreason = ''

        # increment KKR runs loop counter
        self.ctx.loop_count += 1

        # check if initial step was successful
        if not self.ctx.voro_step_success:
            stopreason = 'voronoi step unsucessful'
            do_kkr_step = False
            return do_kkr_step

        # check if previous calculation reached convergence criterion
        if self.ctx.kkr_converged:
            if not self.ctx.kkr_higher_accuracy:
                do_kkr_step = do_kkr_step & True
            else:
                stopreason = 'KKR converged'
                do_kkr_step = False
        else:
            do_kkr_step = do_kkr_step & True

        # next check only needed if another iteration should be done after validating convergence etc. (previous checks)
        if do_kkr_step:
            # check if maximal number of iterations has been reached
            if self.ctx.loop_count <= self.ctx.max_number_runs:
                do_kkr_step = do_kkr_step & True
            else:
                do_kkr_step = False
                # TODO do this differently, return of exit code needs to be done outside of while loop!
                # return self.exit_codes.ERROR_MAX_KKR_RESTARTS_REACHED
                return do_kkr_step

        self.report(f'INFO: done checking condition for kkr step (result={do_kkr_step})')

        if not do_kkr_step:
            self.report(f'INFO: stopreason={stopreason}')

        return do_kkr_step

    def update_kkr_params(self):
        """Update set of KKR parameters.

        (check for reduced mixing, change of mixing strategy, change of accuracy setting)
        """
        self.report('INFO: updating kkr param step')
        decrease_mixing_fac = False
        switch_agressive_mixing = False
        switch_higher_accuracy = False
        initial_settings = False

        # only do something other than sample mixing after first kkr run
        if self.ctx.loop_count != 1:
            # first determine if previous step was successful
            # (otherwise try to find some rms value and decrease mixing to try again)
            if not self.ctx.kkr_step_success:
                try:
                    # check if calculation did start (maybe cluster had some hiccup)
                    calc = self.ctx.calcs[-1]
                    has_output_node = len(calc.get_outgoing(link_label_filter='output_parameters').all()) > 0
                    self.report(
                        'INFO: last KKR calculation failed. '
                        'Probably because the cluster had some issue. '
                        'Try to resubmit the same calculation'
                    )
                except:
                    # otherwise try to decrease the mixing factor
                    decrease_mixing_fac = True
                    self.report('INFO: last KKR calculation failed. Trying to decrease mixfac')

            convergence_on_track = self.convergence_on_track()

            # check if calculation was on its way to converge
            if not convergence_on_track:
                decrease_mixing_fac = True
                self.report('INFO: last KKR did not converge. trying decreasing mixfac')
                self.report(f'INFO: ctx.calcs: {self.ctx.calcs} {type(self.ctx.calcs)}')
                # reset last_remote to last successful calculation
                # needs to be list because `(x)range` does not support slicing
                calclist = list(range(len(self.ctx.calcs)))
                if len(calclist) > 1:
                    calclist = calclist[::-1]  # go backwards through list
                for icalc in calclist:
                    self.report(f'INFO: last calc success? {icalc} {self.ctx.KKR_steps_stats["success"][icalc]}')
                    if self.ctx.KKR_steps_stats['success'][icalc]:
                        self.ctx.last_remote = self.ctx.calcs[icalc].outputs.remote_folder
                        break  # exit loop if last_remote was found successfully
                    self.ctx.last_remote = None
                # if no previous calculation was successful take voronoi output
                # or remote data from input (depending on the inputs)
                self.report(f'INFO: last_remote is None? {self.ctx.last_remote is None} {"structure" in self.inputs}')
                if self.ctx.last_remote is None:
                    if 'structure' in self.inputs:
                        self.ctx.last_remote = self.ctx.voronoi.outputs.last_voronoi_remote
                    else:
                        self.ctx.last_remote = self.inputs.remote_data
                # check if last_remote has finally been set and abort if this is not the case
                self.report(f'INFO: last_remote is still None? {self.ctx.last_remote is None}')
                if self.ctx.last_remote is None:
                    return self.exit_codes.ERROR_LAST_REMOTE_NOT_FOUND  # pylint: disable=no-member

            # check if mixing strategy should be changed
            last_mixing_scheme = self.ctx.last_params.get_dict()['IMIX']
            if last_mixing_scheme is None:
                last_mixing_scheme = 0

            if convergence_on_track:
                last_rms = self.ctx.last_rms_all[-1]

                if (last_rms < self.ctx.threshold_aggressive_mixing and last_mixing_scheme == 0):
                    switch_agressive_mixing = True
                    self.report('INFO: rms low enough, switch to agressive mixing')

                # check if switch to higher accuracy should be done
                if not self.ctx.kkr_higher_accuracy and last_rms < self.ctx.threshold_switch_high_accuracy:
                    switch_higher_accuracy = True
                    self.report('INFO: rms low enough, switch to higher accuracy settings')
        else:
            initial_settings = True

        # check if coarse pre-convergence is done
        if not self.ctx.coarse_preconvergence:
            switch_higher_accuracy = True

        # if needed update parameters
        if (
            decrease_mixing_fac or switch_agressive_mixing or switch_higher_accuracy or initial_settings or
            self.ctx.mag_init
        ):

            if initial_settings:
                label = 'initial KKR scf parameters'
                description = 'initial parameter set for scf calculation'
            else:
                label = ''
                description = ''

            # step 1: extract info from last input parameters and check consistency
            params = self.ctx.last_params
            input_dict = params.get_dict()
            para_check = kkrparams()

            # step 1.1: try to fill keywords
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)

            # step 1.2: check if all mandatory keys are there and add defaults if missing
            new_params = self.initialize_params(para_check)

            # step 2: change parameter (contained in new_params dictionary)

            # adapt EMIN if core states are too close
            new_params = self.adapt_emin(new_params)

            # adapt mixing and convergence settings
            new_params, label, description = self.change_conv_para(
                new_params, para_check, initial_settings, decrease_mixing_fac, switch_agressive_mixing,
                switch_higher_accuracy, label, description
            )

            # initial magnetization
            new_params = self.initial_mag(new_params, initial_settings)

            # step 2.2 update values
            try:
                for key, val in new_params.items():
                    para_check.set_value(key, val, silent=True)
            except:
                return self.exit_codes.ERROR_PARAM_UPDATE_FAILED  # pylint: disable=no-member

            # step 3:
            self.report(f'INFO: update parameters to: {para_check.get_set_values()}')
            updatenode = orm.Dict(dict=para_check.get_dict())
            updatenode.label = label
            updatenode.description = description

            paranode_new = update_params_wf(self.ctx.last_params, updatenode)
            self.ctx.last_params = paranode_new
        else:
            self.report('INFO: reuse old settings')

        self.report('INFO: done updating kkr param step')

        return None

    def initialize_params(self, para_check):
        """Initialize new_params with missing defaults"""
        # init new_params dict where updated params are collected
        new_params = {}

        # find and add missing defaults
        missing_list = para_check.get_missing_keys(use_aiida=True)
        if missing_list != []:
            kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
            kkrdefaults_updated = []
            for key_default, val_default in list(kkrdefaults.items()):
                if key_default in missing_list:
                    new_params[key_default] = kkrdefaults.get(key_default)
                    kkrdefaults_updated.append(key_default)
            if len(kkrdefaults_updated) > 0:
                self.report(f'ERROR: calc_parameters misses keys: {missing_list}')
                return self.exit_codes.ERROR_CALC_PARAMTERS_INCOMPLETE  # pylint: disable=no-member
            self.report(f'updated KKR parameter node with default values: {kkrdefaults_updated}')

        return new_params

    def adapt_emin(self, new_params):
        """Change EMIN if core states are too close."""

        # get parent calculation (either Voronoi or KkrCalculation)
        parent_calc = self.ctx.last_remote.get_incoming(node_class=orm.CalcJobNode).first().node

        # get highest lying core state
        parent_calc_out = parent_calc.outputs.output_parameters
        ecores = parent_calc_out['core_states_group']['energy_highest_lying_core_state_per_atom']
        ecore_max = max(list(set([ec for ec in ecores if ec is not None])))

        # get emin from Voronoi or KkrCalculation
        if parent_calc.process_class == VoronoiCalculation:
            emin_old = parent_calc.outputs.output_parameters['emin']
        else:
            emin_old = parent_calc.outputs.output_parameters['energy_contour_group']['emin']

        # find new emin if distance to core states is too small
        ecore_mindist = self.ctx.delta_e_min_core_states  # in Ry
        if abs(ecore_max - emin_old) < ecore_mindist:
            # move emin lower by delta_e (converted to Ry units)
            emin_new = ecore_max - self.ctx.delta_e * eV2Ry
            self.report(
                f'INFO: Core states too close to start of contour.'
                f'\nChanging EMIN to {emin_new} (ecore_max={ecore_max}, emin_old={emin_old})'
            )
            new_params['EMIN'] = emin_new

        return new_params

    def change_conv_para(
        self, new_params, para_check, initial_settings, decrease_mixing_fac, switch_agressive_mixing,
        switch_higher_accuracy, label, description
    ):
        """Adapt the kkr parameters to change the convergence settings and the mixing"""

        if initial_settings and 'structure' in self.inputs:
            # make sure to ignore IMIX from input node
            # (start with simple mixing even if IMIX is set otherwise)
            # this is enforced whenever voronoi step is starting point
            # (otherwise you may want to continue a preconverged calculation)
            last_mixing_scheme = None
        else:
            last_mixing_scheme = para_check.get_value('IMIX')
        if last_mixing_scheme is None:
            last_mixing_scheme = 0

        strmixfac = self.ctx.strmix
        brymixfac = self.ctx.brymix
        nsteps = self.ctx.nsteps

        # add number of scf steps
        new_params['NSTEPS'] = nsteps

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
                self.report(f'(brymixfax, mixreduce)= ({brymixfac}, {self.ctx.mixreduce})')
                self.report(f'type(brymixfax, mixreduce)= {type(brymixfac)} {type(self.ctx.mixreduce)}')
                brymixfac = brymixfac * self.ctx.mixreduce
                self.ctx.brymix = brymixfac
                label += 'decreased_mix_fac_bry'
                description += f'decreased BRYMIX factor by {self.ctx.mixreduce}'

        # add mixing factor
        new_params['STRMIX'] = strmixfac
        new_params['BRYMIX'] = brymixfac

        if switch_agressive_mixing:
            last_mixing_scheme = 5
            label += ' switched_to_agressive_mixing'
            description += f' switched to agressive mixing scheme (IMIX={last_mixing_scheme})'

        # add mixing scheme
        new_params['IMIX'] = last_mixing_scheme
        self.ctx.last_mixing_scheme = last_mixing_scheme

        if switch_higher_accuracy:
            self.ctx.kkr_higher_accuracy = True
            convergence_settings = self.ctx.convergence_setting_fine
            label += ' use_higher_accuracy'
            description += ' using higher accuracy settings goven in convergence_setting_fine'
        else:
            convergence_settings = self.ctx.convergence_setting_coarse

        # slightly increase temperature if previous calculation was
        # unsuccessful for the second time
        if decrease_mixing_fac and not self.convergence_on_track():
            self.report(
                'INFO: last calculation did not converge and convergence not on track. Try to increase temperature by 50K.'
            )
            convergence_settings['tempr'] += self.ctx.tempr_increase
            label += ' TEMPR_increased'
            description += f' with increased temperature of {self.ctx.tempr_increase}K'

        # add convegence settings
        if self.ctx.loop_count == 1 or self.ctx.last_mixing_scheme == 0:
            new_params['QBOUND'] = self.ctx.threshold_aggressive_mixing
        else:
            new_params['QBOUND'] = self.ctx.convergence_criterion
        new_params['NPOL'] = convergence_settings['npol']
        new_params['NPT1'] = convergence_settings['n1']
        new_params['NPT2'] = convergence_settings['n2']
        new_params['NPT3'] = convergence_settings['n3']
        new_params['TEMPR'] = convergence_settings['tempr']
        new_params['BZDIVIDE'] = convergence_settings['kmesh']

        return new_params, label, description

    def initial_mag(self, new_params, initial_settings):
        """Add settings for initial magnetization"""
        if initial_settings and self.ctx.mag_init:
            if self.ctx.hfield <= 0:
                self.report(
                    f'\nWARNING: magnetization initialization chosen but hfield is zero. '
                    f'Automatically change back to default value (hfield={self._wf_default["hfield"]})\n'
                )
                self.ctx.hfield = self._wf_default['hfield']
            xinipol = self.ctx.xinit
            if xinipol is None:
                # find structure to determine needed length on xinipol
                if 'structure' in self.inputs:
                    struc = self.inputs.structure
                else:
                    struc, voro_parent = VoronoiCalculation.find_parent_structure(self.ctx.last_remote)
                natom = len(get_site_symbols(struc))
                xinipol = np.ones(natom)
            new_params['LINIPOL'] = True
            new_params['HFIELD'] = self.ctx.hfield
            new_params['XINIPOL'] = xinipol
        # turn off initialization after first (successful) iteration
        elif self.ctx.mag_init and self.ctx.mag_init_step_success:
            new_params['LINIPOL'] = False
            new_params['HFIELD'] = 0.0
        elif not self.ctx.mag_init:
            self.report("INFO: mag_init is False. Overwrite 'HFIELD' to '0.0' and 'LINIPOL' to 'False'.")
            # reset mag init to avoid reinitializing
            new_params['HFIELD'] = 0.0
            new_params['LINIPOL'] = False

        # set nspin to 2 if mag_init is used
        if self.ctx.mag_init:
            nspin_in = new_params.get('NSPIN')
            if nspin_in is None:
                nspin_in = 1
            if nspin_in < 2:
                self.report('WARNING: found NSPIN=1 but for maginit needs NPIN=2. Overwrite this automatically')
                new_params['NSPIN'] = 2

        return new_params

    def run_kkr(self):
        """
        submit a KKR calculation
        """
        self.report(f'INFO: setting up kkr calculation step {self.ctx.loop_count}')

        label = f'KKR calculation step {self.ctx.loop_count} (IMIX={self.ctx.last_mixing_scheme})'
        description = f'KKR calculation of step {self.ctx.loop_count}, using mixing scheme {self.ctx.last_mixing_scheme}'
        code = self.inputs.kkr
        remote = self.ctx.last_remote
        params = self.ctx.last_params.get_dict()

        # overwrite some parameters of the KKR calculation by hand before setting mandatory keys
        if 'params_kkr_overwrite' in self.inputs:
            for key, val in self.inputs.params_kkr_overwrite.get_dict().items():
                params[key] = val
                self.report(f'INFO: overwriting KKR parameter: {key} with {val} from params_kkr_overwrite input node')

        options = {
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'resources': self.ctx.resources,
            'queue_name': self.ctx.queue
        }
        if self.ctx.custom_scheduler_commands:
            options['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands

        inputs = get_inputs_kkr(
            code,
            remote,
            options,
            label,
            description,
            parameters=orm.Dict(params),
            serial=(not self.ctx.withmpi),
        )

        # pass nonco angles setting to KkrCalculation
        if 'initial_noco_angles' in self.inputs:
            inputs['initial_noco_angles'] = self.ctx.initial_noco_angles

        # run the KKR calculation
        self.report('INFO: doing calculation')
        kkr_run = self.submit(KkrCalculation, **inputs)

        return ToContext(kkr=kkr_run, last_calc=kkr_run)

    def inspect_kkr(self):
        """Check for convergence and store the results of the last calculation."""
        self.report('INFO: inspecting kkr results step')

        self.report(f'Caching info: {self.ctx.last_calc.get_cache_source()}')

        self.ctx.calcs.append(self.ctx.last_calc)
        self.ctx.kkr_step_success = True

        # check calculation state
        if not self.ctx.last_calc.is_finished_ok:
            self.ctx.kkr_step_success = False
            self.report('ERROR: last calculation not finished correctly')

        self.report(f'INFO: kkr_step_success: {self.ctx.kkr_step_success}')

        # extract convergence info about rms etc. (used to determine convergence behavior)
        try:
            self.report(f'INFO: trying to find output of last_calc: {self.ctx.last_calc}')
            last_calc_output = self.ctx.last_calc.outputs.output_parameters.get_dict()
            found_last_calc_output = True
        except:
            found_last_calc_output = False
        self.report(f'INFO: found_last_calc_output: {found_last_calc_output}')

        # try to extract remote folder
        try:
            self.ctx.last_remote = self.ctx.kkr.outputs.remote_folder
        except:
            self.ctx.last_remote = None
            self.ctx.kkr_step_success = False

        self.report(f'INFO: last_remote: {self.ctx.last_remote}')

        if self.ctx.kkr_step_success and found_last_calc_output:
            # check rms, compare spin and charge values and take bigger one
            rms_charge = last_calc_output['convergence_group']['rms']
            # returning 0 if not found allows to reuse older verisons (e.g. in caching)
            rms_spin = last_calc_output['convergence_group'].get('rms_spin', 0)
            if rms_spin is None:
                rms_spin = 0  # this happens for NSPIN==1
            if rms_charge >= rms_spin:
                rms_max = rms_charge
                use_rms_charge = True
            else:
                rms_max = rms_spin
                use_rms_charge = False
            self.ctx.rms.append(rms_max)
            if use_rms_charge:
                rms_all_iter_last_calc = list(last_calc_output['convergence_group']['rms_all_iterations'])
            else:
                rms_all_iter_last_calc = list(last_calc_output['convergence_group']['rms_spin_all_iterations'])
            # check charge neutrality
            self.ctx.neutr.append(last_calc_output['convergence_group']['charge_neutrality'])
            neutr_all_iter_last_calc = list(last_calc_output['convergence_group']['charge_neutrality_all_iterations'])

            # add lists of last iterations
            self.ctx.last_rms_all = rms_all_iter_last_calc
            self.ctx.last_neutr_all = neutr_all_iter_last_calc
            if self.ctx.kkr_step_success and self.convergence_on_track():
                self.ctx.rms_all_steps += rms_all_iter_last_calc
                self.ctx.neutr_all_steps += neutr_all_iter_last_calc

            # check if calculation converged
            self.ctx.kkr_converged = True
            if rms_max > self.ctx.convergence_criterion:
                self.ctx.kkr_converged = False

        else:
            # if last step did not succeed we know the calculation did not converge
            self.ctx.kkr_converged = False

        self.report(f'INFO: kkr_converged: {self.ctx.kkr_converged}')
        self.report(f'INFO: rms: {self.ctx.rms}')
        self.report(f'INFO: last_rms_all: {self.ctx.last_rms_all}')
        #self.report("INFO: rms_all_steps: {}".format(self.ctx.rms_all_steps))
        self.report(f'INFO: charge_neutrality: {self.ctx.neutr}')
        self.report(f'INFO: last_neutr_all: {self.ctx.last_neutr_all}')
        #self.report("INFO: neutr_all_steps: {}".format(self.ctx.neutr_all_steps))

        # turn off initial magnetization once one step was successful (update_kkr_params) used in
        if self.ctx.mag_init and self.ctx.kkr_step_success:
            self.ctx.mag_init_step_success = True

        # TODO: extract something else (maybe total energy, charge neutrality, magnetisation)?

        # store some statistics used to print table in the end of the report
        tmplist = self.ctx.KKR_steps_stats.get('success', [])
        self.report(f'INFO: append kkr_step_success {tmplist}, {self.ctx.kkr_step_success}')
        tmplist.append(self.ctx.kkr_step_success)
        self.ctx.KKR_steps_stats['success'] = tmplist
        try:
            isteps = self.ctx.last_calc.outputs.output_parameters.get_dict(
            )['convergence_group']['number_of_iterations']
        except (AttributeError, NotExistentAttributeError):
            self.ctx.warnings.append('cound not set isteps in KKR_steps_stats dict')
            isteps = -1

        try:
            first_rms = self.ctx.last_rms_all[0]
            last_rms = self.ctx.last_rms_all[-1]
        except:
            self.ctx.warnings.append('cound not set first_rms, last_rms in KKR_steps_stats dict')
            first_rms = -1
            last_rms = -1

        try:
            first_neutr = self.ctx.last_neutr_all[0]
            last_neutr = self.ctx.last_neutr_all[-1]
        except:
            self.ctx.warnings.append('cound not set first_neutr, last_neutr in KKR_steps_stats dict')
            first_neutr = -999
            last_neutr = -999

        if self.ctx.last_mixing_scheme == 0:
            mixfac = self.ctx.strmix
        else:
            mixfac = self.ctx.brymix

        if self.ctx.kkr_higher_accuracy:
            qbound = self.ctx.convergence_criterion
        else:
            qbound = self.ctx.threshold_switch_high_accuracy

        # store results in KKR_steps_stats dict
        for key, val in {
            'isteps': isteps,
            'imix': self.ctx.last_mixing_scheme,
            'mixfac': mixfac,
            'qbound': qbound,
            'high_sett': self.ctx.kkr_higher_accuracy,
            'first_rms': first_rms,
            'last_rms': last_rms,
            'first_neutr': first_neutr,
            'last_neutr': last_neutr,
            'pk': self.ctx.last_calc.pk,
            'uuid': self.ctx.last_calc.uuid
        }.items():
            tmplist = self.ctx.KKR_steps_stats.get(key, [])
            tmplist.append(val)
            self.ctx.KKR_steps_stats[key] = tmplist

        # update noco angles
        if 'initial_noco_angles' in self.inputs and (self.ctx.kkr_step_success and found_last_calc_output):
            self._get_new_noco_angles()

        self.report('INFO: done inspecting kkr results step')

    def convergence_on_track(self):
        """Check if convergence behavior of the last calculation is on track.

        (i.e. going down)
        """
        on_track = True
        threshold = 5.  # used to check condition if at least one of charnge_neutrality, rms-error goes down fast enough
        eps_neutr = 1e-4
        eps_rms = 1e-6

        # first check if previous calculation was stopped due to reaching the QBOUND limit
        try:
            calc_reached_qbound = self.ctx.last_calc.outputs.output_parameters.get_dict(
            )['convergence_group']['calculation_converged']
        except:  # captures error when last_calc dies not have an output node
            calc_reached_qbound = False

        if self.ctx.kkr_step_success and not calc_reached_qbound:
            first_rms = self.ctx.last_rms_all[0]
            first_neutr = abs(self.ctx.last_neutr_all[0])
            last_rms = self.ctx.last_rms_all[-1]
            last_neutr = abs(self.ctx.last_neutr_all[-1])
            # use this trick to avoid division by zero
            if first_neutr == 0:
                first_neutr = 10**-16
            if first_rms == 0:
                first_rms = 10**-16
            r, n = last_rms / first_rms, last_neutr / first_neutr

            # deal with small differences
            if abs(first_neutr - last_neutr) < eps_neutr:
                n = 0
            if abs(first_rms - last_rms) > eps_rms:
                r = 0

            self.report(
                f'INFO convergence check: first/last rms {first_rms},'
                f' {last_rms}; first/last neutrality {first_neutr}, {last_neutr}'
            )
            if r < 1 and n < 1:
                self.report('INFO convergence check: both rms and neutrality go down')
                on_track = True
            elif n > threshold or r > threshold:
                self.report('INFO convergence check: rms or neutrality goes up too fast, convergence is not expected')
                on_track = False
            elif n * r < 1:
                self.report('INFO convergence check: either rms goes up and neutrality goes down or vice versa')
                self.report('INFO convergence check: but product goes down fast enough')
                on_track = True
            elif len(self.ctx.last_rms_all) == 1:
                self.report('INFO convergence check: already converged after single iteration')
                on_track = True
            else:
                self.report(
                    'INFO convergence check: rms or neutrality do not shrink fast enough, convergence is not expected'
                )
                on_track = False

        elif calc_reached_qbound:
            self.report('INFO convergence check: calculation reached QBOUND')
            on_track = True
        else:
            self.report('INFO convergence check: calculation unsuccessful')
            on_track = False

        self.report(f'INFO convergence check result: {on_track}')

        return on_track

    def return_results(self):
        """Return the results of the calculations.

        This should run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """

        self.report('INFO: entering return_results')

        # try/except to capture as mnuch as possible
        # (everything that is there even when workflow exits unsuccessfully)
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
            last_calc_out = self.ctx.kkr.outputs.output_parameters
            last_calc_out_dict = last_calc_out.get_dict()
            self.report('Found last_calc_out')
            last_RemoteData = self.ctx.last_remote
            self.report('Found last_remote')
            last_InputParameters = self.ctx.last_params
            self.report('Found last_params')
        except:
            self.report('Error in finding last_calc_out etc.')
            last_InputParameters = None
            last_RemoteData = None
            last_calc_out = None
            last_calc_out_dict = {}

        # capture convergence info
        try:
            last_rms = self.ctx.rms[-1]
            last_neutr = self.ctx.neutr[-1]
        except:
            last_rms = None
            last_neutr = None

        # capture result of vorostart sub-workflow
        try:
            results_vorostart = self.ctx.voronoi.outputs.results_vorostart_wc
        except:
            results_vorostart = None
        try:
            starting_dosdata_interpol = self.ctx.voronoi.outputs.last_doscal_dosdata_interpol
        except:
            starting_dosdata_interpol = None

        try:
            final_dosdata_interpol = self.ctx.doscal.outputs.dos_data_interpol
        except:
            final_dosdata_interpol = None

        # now collect results saved in results node of workflow
        self.report('INFO: collect outputnode_dict')
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
        outputnode_dict['errors'] = self.ctx.errors
        outputnode_dict['convergence_value'] = last_rms
        outputnode_dict['convergence_values_all_steps'] = np.array(self.ctx.rms_all_steps)
        outputnode_dict['convergence_values_last_step'] = np.array(self.ctx.last_rms_all)
        outputnode_dict['charge_neutrality'] = last_neutr
        outputnode_dict['charge_neutrality_all_steps'] = np.array(self.ctx.neutr_all_steps)
        outputnode_dict['charge_neutrality_last_step'] = np.array(self.ctx.last_neutr_all)
        outputnode_dict['dos_check_ok'] = self.ctx.dos_ok
        outputnode_dict['convergence_reached'] = self.ctx.kkr_converged
        outputnode_dict['voronoi_step_success'] = self.ctx.voro_step_success
        outputnode_dict['kkr_step_success'] = self.ctx.kkr_step_success
        outputnode_dict['used_higher_accuracy'] = self.ctx.kkr_higher_accuracy

        # report the status
        if self.ctx.successful:
            self.report(
                'STATUS: Done, the convergence criteria are reached.\n'
                f'INFO: The charge density of the KKR calculation pk= {last_calc_pk} '
                f'converged after {self.ctx.loop_count} KKR runs and '
                f'{self.ctx.loop_count} iterations to {last_rms} \n'
            )
        else:  # Termination ok, but not converged yet...
            if self.ctx.abort:  # some error occured, donot use the output.
                self.report('STATUS/ERROR: I abort, see logs and '
                            'erros/warning/hints in output_kkr_scf_wc_para')
            else:
                self.report(
                    'STATUS/WARNING: Done, the maximum number of runs '
                    'was reached or something failed.\n INFO: The '
                    'charge density of the KKR calculation pk= '
                    f'after {self.ctx.loop_count} KKR runs and '
                    f"{sum(self.ctx.KKR_steps_stats.get('isteps')),} "
                    f"iterations is {last_rms} 'me/bohr^3'\n"
                )

        # create results  node
        # : {}".format(outputnode_dict))
        self.report('INFO: create results node')
        outputnode_t = orm.Dict(dict=outputnode_dict)
        outputnode_t.label = 'kkr_scf_wc_results'
        outputnode_t.description = (
            'Contains results of workflow'
            ' (e.g. workflow version number, info about success of wf,'
            ' lis tof warnings that occured during execution, ...)'
        )

        # collect nodes in outputs dictionary
        out_nodes = {'outpara': outputnode_t}
        if last_calc_out is not None:
            out_nodes['last_calc_out'] = last_calc_out
        if last_RemoteData is not None:
            out_nodes['last_RemoteData'] = last_RemoteData
        if last_InputParameters is not None:
            out_nodes['last_InputParameters'] = last_InputParameters
        if final_dosdata_interpol is not None:
            out_nodes['final_dosdata_interpol'] = final_dosdata_interpol
        if starting_dosdata_interpol is not None:
            out_nodes['starting_dosdata_interpol'] = starting_dosdata_interpol
        if results_vorostart is not None:
            out_nodes['results_vorostart'] = results_vorostart
        if 'initial_noco_angles' in self.inputs and not all(self.ctx.initial_noco_angles['fix_dir']):
            # was updated in inspect_kkr to last noco angles output
            out_nodes['last_noco_angles'] = self.ctx.initial_noco_angles

        # call helper function to create output nodes in correct AiiDA graph structure
        outdict = create_scf_result_node(**out_nodes)

        for link_name, node in outdict.items():
            #self.report("INFO: storing node {} {} with linkname {}".format(type(node), node, link_name))
            self.out(link_name, node)

        # print results table for overview
        # table layout:
        message = 'INFO: overview of the result:\n\n'
        message += '|------|---------|--------|------|--------|-------------------|-----------------|-----------------|-------------\n'
        message += '| irun | success | isteps | imix | mixfac | accuracy settings |       rms       | abs(neutrality) | pk and uuid \n'
        message += '|      |         |        |      |        | qbound  | higher? | first  |  last  | first  |  last  |             \n'
        message += '|------|---------|--------|------|--------|---------|---------|--------|--------|--------|--------|-------------\n'
        # | %6i  | %9s     | %8i    | %6i  | %.2e   | %.3e    | %9s     | %.2e   |  %.2e  |  %.2e  |  %.2e  |
        KKR_steps_stats = self.ctx.KKR_steps_stats
        for irun in range(len(KKR_steps_stats.get('success'))):
            KKR_steps_stats.get('first_neutr')[irun] = abs(KKR_steps_stats.get('first_neutr')[irun])
            KKR_steps_stats.get('last_neutr')[irun] = abs(KKR_steps_stats.get('last_neutr')[irun])
            message += (
                '|{:6d}|{:9s}|{:8d}|{:6d}|{:.2e}|{:.3e}|{:9s}|{:.2e}|{:.2e}|{:.2e}|{:.2e}|'.format(
                    irun + 1, str(KKR_steps_stats.get('success')[irun]),
                    KKR_steps_stats.get('isteps')[irun],
                    KKR_steps_stats.get('imix')[irun],
                    KKR_steps_stats.get('mixfac')[irun],
                    KKR_steps_stats.get('qbound')[irun], str(KKR_steps_stats.get('high_sett')[irun]),
                    KKR_steps_stats.get('first_rms')[irun],
                    KKR_steps_stats.get('last_rms')[irun],
                    KKR_steps_stats.get('first_neutr')[irun],
                    KKR_steps_stats.get('last_neutr')[irun]
                )
            )
            message += f" {KKR_steps_stats.get('pk')[irun]} | {KKR_steps_stats.get('uuid')[irun]}\n"
            """
            message += "#|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".format(irun+1,
                          KKR_steps_stats.get('success')[irun], KKR_steps_stats.get('isteps')[irun],
                          KKR_steps_stats.get('imix')[irun], KKR_steps_stats.get('mixfac')[irun],
                          KKR_steps_stats.get('qbound')[irun], KKR_steps_stats.get('high_sett')[irun],
                          KKR_steps_stats.get('first_rms')[irun], KKR_steps_stats.get('last_rms')[irun],
                          KKR_steps_stats.get('first_neutr')[irun], KKR_steps_stats.get('last_neutr')[irun])
            """
        self.report(message)

        self.report('\nINFO: done with kkr_scf workflow!\n')

    def check_input_params(self, params, is_voronoi=False):
        """Check the input parameter consistency and aborts wf if check fails."""
        if params is None:
            return self.exit_codes.ERROR_NO_CALC_PARAMETERS_GIVEN  # pylint: disable=no-member
        input_dict = params.get_dict()
        if is_voronoi:
            para_check = kkrparams(params_type='voronoi')
        else:
            para_check = kkrparams()

        # step 1 try to fill keywords
        try:
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)
        except:
            return self.exit_codes.ERROR_CALC_PARAMTERS_INCONSISTENT  # pylint: disable=no-member

        # step 2: check if all mandatory keys are there
        missing_list = para_check.get_missing_keys(use_aiida=True)
        if missing_list != []:
            all_defaults = True
            for key in missing_list:
                if key not in kkrparams.get_KKRcalc_parameter_defaults()[0]:
                    all_defaults = False
            if not all_defaults:
                self.report(
                    f'ERROR: calc_parameters given are not consistent! '
                    f'Missing mandatory keys: {missing_list}'
                )
                return self.exit_codes.ERROR_CALC_PARAMETERS_INCOMPLETE  # pylint: disable=no-member
        return None

    def get_dos(self):
        """Call the dos sub workflow pass the input and submit the calculation."""
        if self.ctx.check_dos:
            self.report('INFO: Doing final DOS calculation')

            # fix emin/emax to include emin, ef of scf contour
            # remember: efermi, emin and emax are in internal units (Ry) but delta_e is in eV!
            emin = self.ctx.dos_params['emin']  # from dos params
            # from contour (sets limit of dos emin!)
            emin_cont = self.ctx.last_calc.outputs.output_parameters.get_dict().get('energy_contour_group').get('emin')
            if emin_cont - self.ctx.delta_e * eV2Ry < emin:
                self.ctx.dos_params['emin'] = (emin_cont - self.ctx.delta_e * eV2Ry)
                self.report(
                    f'INFO: emin ({emin_cont} Ry) - '
                    f'delta_e ({self.ctx.delta_e * eV2Ry} Ry) smaller than '
                    f'emin ({emin} Ry) of voronoi output. '
                    f'Setting automatically to {emin_cont - self.ctx.delta_e * eV2Ry}Ry'
                )
            with self.ctx.last_calc.outputs.retrieved.open('out_potential') as _file:
                self.ctx.efermi = get_ef_from_potfile(_file)
            emax = self.ctx.dos_params['emax']
            if emax < self.ctx.efermi + self.ctx.delta_e * eV2Ry:
                self.ctx.dos_params['emax'] = self.ctx.efermi + self.ctx.delta_e * eV2Ry
                self.report(
                    f'INFO: self.ctx.efermi ({self.ctx.efermi} Ry) '
                    f'+ delta_e ({self.ctx.delta_e * eV2Ry} Ry) larger than '
                    f'emax ({emax} Ry). Setting automatically to '
                    f'{self.ctx.efermi + self.ctx.delta_e * eV2Ry}Ry'
                )

            wfdospara_node = orm.Dict(dict=self.ctx.dos_params)
            wfdospara_node.label = 'DOS params'
            wfdospara_node.description = 'DOS parameter set for final DOS calculation of kkr_scf_wc'

            code = self.inputs.kkr
            remote = self.ctx.last_calc.outputs.remote_folder
            wf_label = ' final DOS calculation'
            wf_desc = ' subworkflow of a DOS calculation'

            builder = kkr_dos_wc.get_builder()
            builder.metadata.description = wf_desc  # pylint: disable=no-member
            builder.metadata.label = wf_label  # pylint: disable=no-member
            builder.kkr = code
            builder.wf_parameters = wfdospara_node
            builder.options = self.ctx.options_params_dict
            builder.remote_data = remote

            future = self.submit(builder)

            return ToContext(doscal=future)

    def check_dos(self):
        """Checks if dos of final potential is ok."""
        # initialize dos_ok variable
        self.ctx.dos_ok = True

        # first check if dos should be checked or if the test is skipped
        if not self.ctx.check_dos:
            self.report('INFO: skipping DOS check')
            return

        self.report('INFO: checking DOS for consistency (EMIN position, negative DOS, etc.)')

        # check parser output
        doscal = self.ctx.doscal
        try:
            dos_outdict = doscal.outputs.results_wf.get_dict()
            if not dos_outdict['successful']:
                self.report('ERROR: DOS workflow unsuccessful')
                self.ctx.dos_ok = False
                return self.exit_codes.ERROR_DOS_RUN_UNSUCCESSFUL  # pylint: disable=no-member

            if dos_outdict['list_of_errors'] != []:
                self.report(f'ERROR: DOS wf output contains errors: {dos_outdict["list_of_errors"]}')
                self.ctx.dos_ok = False
                return self.exit_codes.ERROR_DOS_RUN_UNSUCCESSFUL  # pylint: disable=no-member
        except AttributeError:
            self.ctx.dos_ok = False
            return self.exit_codes.ERROR_DOS_RUN_UNSUCCESSFUL  # pylint: disable=no-member

        # check for negative DOS
        try:
            dosdata = doscal.outputs.dos_data
            natom = self.ctx.last_calc.outputs.output_parameters.get_dict()['number_of_atoms_in_unit_cell']
            nspin = dos_outdict['nspin']

            ener = dosdata.get_x()[1]
            totdos = dosdata.get_y()[0][1]

            if len(ener) != nspin * natom:
                self.report(
                    f'ERROR: DOS output shape does not fit nspin, natom '
                    f'information: len(energies)={len(ener)}, natom={natom}, nspin={nspin}'
                )
                self.ctx.doscheck_ok = False
                return self.exit_codes.ERROR_DOS_RUN_UNSUCCESSFUL  # pylint: disable=no-member

            # deal with snpin==1 or 2 cases and check negtive DOS
            for iatom in range(natom // nspin):
                for ispin in range(nspin):
                    x, y = ener[iatom * nspin + ispin], totdos[iatom * nspin + ispin]
                    if nspin == 2 and ispin == 0:
                        y = -y
                    if y.min() < 0:
                        self.report(
                            f'INFO: negative DOS value found in '
                            f'(atom, spin)=({iatom},{ispin}) at iteration '
                            f'{self.ctx.loop_count}: {y.min()}'
                        )
                        self.ctx.dos_ok = False

            # check starting EMIN
            dosdata_interpol = doscal.outputs.dos_data_interpol

            ener = dosdata_interpol.get_x()[1]  # shape= natom*nspin, nept
            totdos = dosdata_interpol.get_y()[0][1]  # shape= natom*nspin, nept
            Ry2eV = get_Ry2eV()

            for iatom in range(natom // nspin):
                for ispin in range(nspin):
                    x, y = ener[iatom * nspin + ispin], totdos[iatom * nspin + ispin]
                    xrel = abs(x - (self.ctx.dos_params_dict['emin'] - self.ctx.efermi) * Ry2eV)
                    mask_emin = np.where(xrel == xrel.min())
                    ymin = abs(y[mask_emin])
                    if ymin > self.ctx.threshold_dos_zero:
                        self.report(f'INFO: DOS at emin not zero! {ymin}>{self.ctx.threshold_dos_zero}')
                        self.ctx.dos_ok = False
        except AttributeError:
            self.ctx.dos_ok = False

    def _get_new_noco_angles(self):
        """Extract nonco angles from output of calculation.

        If fix_dir is True we skip this and leave the initial angles unchanged
        Here we update self.ctx.initial_noco_angles with the new values
        """
        # first check if we need to update the angles
        if all(self.ctx.initial_noco_angles['fix_dir']):
            return

        # extract output angles from last calculation and update initial_noco_angles in context
        self.ctx.initial_noco_angles = extract_noco_angles(
            fix_dir_threshold=orm.Float(self.ctx.fix_dir_threshold),
            old_noco_angles=self.ctx.initial_noco_angles,
            last_retrieved=self.ctx.last_calc.outputs.retrieved
        )
        if self.ctx.initial_noco_angles == {}:
            self.ctx.initial_noco_angles = self.ctx.initial_noco_angles


@workfunction
def create_scf_result_node(**kwargs):
    """This is a pseudo wf, to create the right graph structure of AiiDA.

    This workfunction will create the output node in the database.
    It also connects the output_node to all nodes the information commes from.
    So far it is just also parsed in as argument, because so far we are to lazy
    to put most of the code overworked from return_results in here.
    """

    has_last_outpara = False
    has_last_calc_out_dict = False
    has_last_RemoteData = False
    has_vorostart_output = False
    has_starting_dos = False
    has_final_dos = False
    has_last_InputParameters = False

    for key, val in kwargs.items():
        if key == 'outpara':  # should always be there
            outpara = val
            has_last_outpara = True
        elif key == 'last_calc_out':
            has_last_calc_out_dict = True
            last_calc_out_dict = val
        elif key == 'last_RemoteData':
            last_RemoteData_dict = val
            has_last_RemoteData = True
        elif key == 'last_InputParameters':
            last_InputParameters_dict = val
            has_last_InputParameters = True
        elif key == 'results_vorostart':
            has_vorostart_output = True
            vorostart_output_dict = val
        elif key == 'starting_dosdata_interpol':
            has_starting_dos = True
            start_dosdata_interpol_dict = val
        elif key == 'final_dosdata_interpol':
            has_final_dos = True
            final_dosdata_interpol_dict = val

    outdict = {}

    if has_last_outpara:
        outputnode = outpara  # pylint: disable=possibly-used-before-assignment
        outputnode.label = 'workflow_Results'
        outputnode.description = ('Contains self-consistency results and '
                                  'information of an kkr_scf_wc run.')
        outdict['output_kkr_scf_wc_ParameterResults'] = outputnode

    if has_last_calc_out_dict:
        outputnode = last_calc_out_dict  # pylint: disable=used-before-assignment
        outputnode.label = 'last_calc_out'
        outputnode.description = (
            'Contains the Results Parameter node from the output '
            'of the last calculation done in the workflow.'
        )
        outdict['last_calc_out'] = outputnode

    if has_last_RemoteData:
        outputnode = last_RemoteData_dict  # pylint: disable=used-before-assignment
        outputnode.label = 'last_RemoteData'
        outputnode.description = (
            'Contains a link to the latest remote data node '
            'where the output of the calculation can be accessed.'
        )
        outdict['last_RemoteData'] = outputnode

    if has_last_InputParameters:
        outputnode = last_InputParameters_dict  # pylint: disable=used-before-assignment
        outputnode.label = 'last_InputParameters'
        outputnode.description = (
            'Contains the latest parameter data node '
            'where the input of the last calculation can be found.'
        )
        outdict['last_InputParameters'] = outputnode

    if has_vorostart_output:
        outputnode = vorostart_output_dict  # pylint: disable=used-before-assignment
        outputnode.label = 'results_vorostart'
        outputnode.description = (
            'Contains the results parameter data node '
            'of the vorostart sub-workflow (sets up starting portentials).'
        )
        outdict['results_vorostart'] = outputnode

    if has_starting_dos:
        outputnode = start_dosdata_interpol_dict  # pylint: disable=used-before-assignment
        outputnode.label = 'starting_dosdata_interpol'
        outputnode.description = ('Contains the interpolated DOS data note, computed '
                                  'from the starting portential.')
        outdict['starting_dosdata_interpol'] = outputnode

    if has_final_dos:
        outputnode = final_dosdata_interpol_dict  # pylint: disable=used-before-assignment
        outputnode.label = 'final_dosdata_interpol'
        outputnode.description = ('Contains the interpolated DOS data note, computed '
                                  'from the converged potential.')
        outdict['final_dosdata_interpol'] = outputnode

    return outdict


def get_site_symbols(structure):
    """
    extract the site number taking into account a possible CPA structure
    """
    sites = structure.sites
    sitelist = []
    for isite, site in enumerate(sites):
        sitekind = structure.get_kind(site.kind_name)
        for ikind in range(len(sitekind.symbols)):
            site_symbol = sitekind.symbols[ikind]
            sitelist.append([isite, site_symbol])

    return sitelist
