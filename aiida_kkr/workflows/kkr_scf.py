#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for converging a kkr calculation and
some helper methods to do so with AiiDA
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from aiida.orm import Code, load_node
from aiida.plugins import DataFactory
from aiida.engine import WorkChain, while_, if_, ToContext
from aiida.engine import workfunction as wf
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import (test_and_get_codenode, get_inputs_kkr,
                                                  get_parent_paranode, update_params_wf)
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.dos import kkr_dos_wc
from masci_tools.io.common_functions import get_Ry2eV, get_ef_from_potfile
from numpy import array, where, ones
from six.moves import range

__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.9.2"
__contributors__ = (u"Jens Broeder", u"Philipp Rüßmann")

#TODO: magnetism (init and converge magnetic state)
#TODO: check convergence (RMAX, GMAX etc.)
#TODO: save timing info of the steps
#TODO: switch to LLOYD
#TODO: emin-emax setting
#TODO: restart from workflow output instead of calculation output
#TODO: add warnings
#TODO: maybe define the energy point density instead of a fixed number as in input?
#TODO: overwrite defaults from parent if parent is previous kkr_scf run
#TODO: retrieve DOS within scf run

RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
Dict = DataFactory('dict')

class kkr_scf_wc(WorkChain):
    """
    Workchain for converging a KKR calculation (SCF).

    It converges the charge potential.
    Two paths are possible:

    (1) Start from a structure and run a voronoi calculation first,
    optional with calc_parameters
    (2) Start from an existing Voronoi or KKR calculation, with a remoteData

    :param wf_parameters: (Dict), Workchain Spezifications
    :param options: (Dict); specifications for the computer
    :param structure: (StructureData), Crystal structure
    :param calc_parameters: (Dict), Vornoi/Kkr Parameters
    :param remote_data: (RemoteData), from a KKR, or Vornoi calculation
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
    _wf_default = {'kkr_runmax': 5,                           # Maximum number of kkr jobs/starts (defauld iterations per start)
                   'convergence_criterion' : 10**-8,          # Stop if charge denisty is converged below this value
                   'mixreduce': 0.5,                          # reduce mixing factor by this factor if calculaito fails due to too large mixing
                   'threshold_aggressive_mixing': 8*10**-3,   # threshold after which agressive mixing is used
                   'strmix': 0.03,                            # mixing factor of simple mixing
                   'brymix': 0.05,                            # mixing factor of aggressive mixing
                   'nsteps': 50,                              # number of iterations done per KKR calculation
                   'convergence_setting_coarse': {            # setting of the coarse preconvergence
                        'npol': 7,
                        'n1': 3,
                        'n2': 11,
                        'n3': 3,
                        'tempr': 1000.0,
                        'kmesh': [10, 10, 10]},
                   'threshold_switch_high_accuracy': 10**-3,  # threshold after which final conversion settings are used
                   'convergence_setting_fine': {              # setting of the final convergence (lower tempr, 48 epts, denser k-mesh)
                        'npol': 5,
                        'n1': 7,
                        'n2': 29,
                        'n3': 7,
                        'tempr': 600.0,
                        'kmesh': [30, 30, 30]},
                   'mag_init' : False,                        # initialize and converge magnetic calculation
                   'hfield' : 0.02, # Ry                      # external magnetic field used in initialization step
                   'init_pos' : None,                         # position in unit cell where magnetic field is applied [default (None) means apply to all]
                   'retreive_dos_data_scf_run' : False,       # add DOS to testopts and retrieve dos.atom files in each scf run
                   }
    _options_default = {'queue_name' : '',                         # Queue name to submit jobs too
                        'resources': {"num_machines": 1},          # resources to allowcate for the job
                        'max_wallclock_seconds' : 60*60,           # walltime after which the job gets killed (gets parsed to KKR)
                        'use_mpi' : True,                          # execute KKR with mpi or without
                        'custom_scheduler_commands' : ''           # some additional scheduler commands
                        }
    # set these keys from defaults in kkr_startpot workflow since they are only passed onto that workflow
    for key, value in kkr_startpot_wc.get_wf_defaults(silent=True).items():
        if key in ['dos_params', 'fac_cls_increase', 'r_cls', 'natom_in_cls_min', 'delta_e_min', 'threshold_dos_zero', 'check_dos']:
            _wf_default[key] = value


    # intended to guide user interactively in setting up a valid wf_params node
    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """
        if not silent: print('Version of workflow: {}'.format(self._workflowversion))
        return self._wf_default


    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_scf_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=Dict, required=False,
                   default=Dict(dict=cls._wf_default))
        spec.input("options", valid_type=Dict, required=False,
                   default=Dict(dict=cls._wf_default))
        spec.input("structure", valid_type=StructureData, required=False)
        spec.input("calc_parameters", valid_type=Dict, required=False)
        spec.input("remote_data", valid_type=RemoteData, required=False)
        spec.input("voronoi", valid_type=Code, required=False)
        spec.input("kkr", valid_type=Code, required=True)

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            # check if voronoi run needed, otherwise skip this step
            if_(cls.validate_input)(
                    # run kkr_startpot workflow (sets up voronoi input, runs voro calc, does some consistency checks)
                    cls.run_voronoi,
                    # check output of run_voronoi and determine if calculation has to be terminated here
                    cls.check_voronoi),
            # while loop for KKR run(s), first simple mixing
            # then Anderson with corase pre-convergence settings
            # finally convergence step with higher accuracy
            while_(cls.condition)(
                # update parameters for kkr step using previous output(s)
                cls.update_kkr_params,
                # run kkr step
                cls.run_kkr,
                # check results for convergence and collect some intermediate results
                cls.inspect_kkr),
            # compute final dos if check_dos is True
            cls.get_dos,
            cls.check_dos,
            # finalize calculation and create output nodes
            cls.return_results
        )

        # definition of exit codes if the workflow needs to be terminated
        spec.exit_code(221, 'ERROR_NO_PARENT_PARAMS_FOUND',
          message='Unable to extract parent paremeter node of input remote folder')
        spec.exit_code(222, 'ERROR_INVALID_KKR_CODE',
          message='The code you provided for kkr does not use the plugin kkr.kkr')
        spec.exit_code(223, 'ERROR_INVALID_VORONOI_CODE',
          message='The code you provided for voronoi does not use the plugin kkr.voro')
        spec.exit_code(224, 'ERROR_NO_VORONOI_CODE_GIVEN',
          message='ERROR: StructureData was provided, but no voronoi code was provided')
        spec.exit_code(225, 'ERROR_NOT_ENOUGH_INPUTS',
          message='ERROR: No StructureData nor remote_data was provided as Input')
        spec.exit_code(226, 'ERROR_KKR_STARTPOT_FAILED',
          message='ERROR: kkr_startpot_wc step failed!')
        spec.exit_code(227, 'ERROR_DOS_RUN_UNSUCCESFUL',
          message='DOS run unsuccessful. Check inputs.')
        spec.exit_code(228, 'ERROR_CALC_PARAMETERS_INCOMPLETE',
          message='ERROR: calc_parameters given are not consistent! Missing mandatory keys')
        spec.exit_code(229, 'ERROR_CALC_PARAMTERS_INCONSISTENT',
          message='ERROR: calc_parameters given are not consistent! Hint: did you give an unknown keyword?')
        spec.exit_code(230, 'ERROR_NO_CALC_PARAMETERS_GIVEN',
          message='ERROR: calc_parameters not given as input but are needed!')
        spec.exit_code(231, 'ERROR_PARAM_UPDATE_FAILED',
          message='ERROR: parameter update unsuccessful: some key, value pair not valid!')
        spec.exit_code(232, 'ERROR_CALC_PARAMTERS_INCOMPLETE',
          message='ERROR: calc_parameters misses keys')
        spec.exit_code(233, 'ERROR_LAST_REMOTE_NOT_FOUND',
          message='ERROR: last_remote could not be set to a previous succesful calculation')
        spec.exit_code(234, 'ERROR_MAX_KKR_RESTARTS_REACHED',
          message='ERROR: maximal number of KKR restarts reached. Exiting now!')
        spec.exit_code(235, 'ERROR_CALC_SUBMISSION_FAILED',
          message='ERROR: last KKRcalc in SUBMISSIONFAILED state')


    def start(self):
        """
        init context and some parameters
        """
        self.report('INFO: started KKR convergence workflow version {}'
                    ''.format(self._workflowversion))

        ####### init #######

        # internal para /control para
        self.ctx.loop_count = 0
        self.ctx.last_mixing_scheme = 0
        self.ctx.calcs = []
        self.ctx.abort = False
        # flags used internally to check whether the individual steps were successful
        self.ctx.kkr_converged = False
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
        self.ctx.use_mpi = options_dict.get('use_mpi', self._options_default['use_mpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.walltime_sec = options_dict.get('max_wallclock_seconds', self._options_default['max_wallclock_seconds'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])
        self.ctx.options_params_dict = Dict(dict={'use_mpi': self.ctx.use_mpi, 'resources': self.ctx.resources,
                                                           'max_wallclock_seconds': self.ctx.walltime_sec, 'queue_name': self.ctx.queue,
                                                           'custom_scheduler_commands': self.ctx.custom_scheduler_commands})

        # set label and description
        self.ctx.description_wf = self.inputs.get('description', 'Workflow for '
                                                  'a KKR scf calculation starting '
                                                  'either from a structure with '
                                                  'automatic voronoi calculation '
                                                  'or a valid RemoteData node of '
                                                  'a previous calculation')
        self.ctx.label_wf = self.inputs.get('label', 'kkr_scf_wc')

        # set workflow parameters
        self.ctx.max_number_runs = wf_dict.get('kkr_runmax', self._wf_default['kkr_runmax'])
        self.ctx.strmix = wf_dict.get('strmix', self._wf_default['strmix'])
        self.ctx.brymix = wf_dict.get('brymix', self._wf_default['brymix'])
        self.ctx.check_dos = wf_dict.get('check_dos', self._wf_default['check_dos'])
        self.ctx.dos_params = wf_dict.get('dos_params', self._wf_default['dos_params'])
        self.ctx.convergence_criterion = wf_dict.get('convergence_criterion', self._wf_default['convergence_criterion'])
        self.ctx.convergence_setting_coarse = wf_dict.get('convergence_setting_coarse', self._wf_default['convergence_setting_coarse'])
        self.ctx.convergence_setting_fine = wf_dict.get('convergence_setting_fine', self._wf_default['convergence_setting_fine'])
        self.ctx.mixreduce = wf_dict.get('mixreduce', self._wf_default['mixreduce'])
        self.ctx.nsteps = wf_dict.get('nsteps', self._wf_default['nsteps'])
        self.ctx.threshold_aggressive_mixing = wf_dict.get('threshold_aggressive_mixing', self._wf_default['threshold_aggressive_mixing'])
        self.ctx.threshold_switch_high_accuracy = wf_dict.get('threshold_switch_high_accuracy', self._wf_default['threshold_switch_high_accuracy'])

        # initial magnetization
        self.ctx.mag_init = wf_dict.get('mag_init', self._wf_default['mag_init'])
        self.ctx.hfield = wf_dict.get('hfield', self._wf_default['hfield'])
        self.ctx.xinit = wf_dict.get('init_pos', self._wf_default['init_pos'])
        self.ctx.mag_init_step_success = False

        # difference in eV to emin (e_fermi) if emin (emax) are larger (smaller) than emin (e_fermi)
        self.ctx.delta_e = wf_dict.get('delta_e_min', self._wf_default['delta_e_min'])
        # threshold for dos comparison (comparison of dos at emin)
        self.ctx.threshold_dos_zero = wf_dict.get('threshold_dos_zero', self._wf_default['threshold_dos_zero'])
        self.ctx.efermi = None

        # retreive dos data in each scf run
        self.ctx.scf_dosdata = wf_dict.get('retreive_dos_data_scf_run', self._wf_default['retreive_dos_data_scf_run'])

        self.report('INFO: use the following parameter:\n'
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
                    'Anderson mixing factor: {}\n'
                    'Nsteps scf cycle: {}\n'
                    'Convergence criterion: {}\n'
                    'threshold_aggressive_mixing: {}\n'
                    'threshold_switch_high_accuracy: {}\n'
                    'convergence_setting_coarse: {}\n'
                    'convergence_setting_fine: {}\n'
                    'factor reduced mixing if failing calculation: {}\n'
                    '\nAdditional parameter\n'
                    'check DOS between runs: {}\n'
                    'DOS parameters: {}\n'
                    'init magnetism in first step: {}\n'
                    'init magnetism, hfield: {}\n'
                    'init magnetism, init_pos: {}\n'
                    ''.format(self.ctx.use_mpi, self.ctx.max_number_runs,
                                self.ctx.resources, self.ctx.walltime_sec,
                                self.ctx.queue, self.ctx.custom_scheduler_commands,
                                self.ctx.description_wf, self.ctx.label_wf,
                                self.ctx.strmix, self.ctx.brymix, self.ctx.nsteps,
                                self.ctx.convergence_criterion,
                                self.ctx.threshold_aggressive_mixing,
                                self.ctx.threshold_switch_high_accuracy,
                                self.ctx.convergence_setting_coarse,
                                self.ctx.convergence_setting_fine,
                                self.ctx.mixreduce, self.ctx.check_dos,
                                self.ctx.dos_params, self.ctx.mag_init,
                                self.ctx.hfield, self.ctx.xinit)
                    )

        # return para/vars
        self.ctx.successful = True
        self.ctx.rms = []
        self.ctx.neutr = []
        self.ctx.warnings = []
        self.ctx.errors = []
        self.ctx.formula = ''

        # for results table each list gets one entry per iteration that has been performed
        self.ctx.KKR_steps_stats = {'success':[],
                                    'isteps':[],
                                    'imix':[],
                                    'mixfac':[],
                                    'qbound':[],
                                    'high_sett':[],
                                    'first_rms':[],
                                    'last_rms':[],
                                    'first_neutr':[],
                                    'last_neutr':[],
                                    'pk':[],
                                    'uuid':[]}

    def validate_input(self):
        """
        # validate input and find out which path (1, or 2) to take
        # return True means run voronoi if false run kkr directly
        """
        run_voronoi = True
        inputs = self.inputs

        if 'structure' in inputs:
            self.report('INFO: Found structure in input. Start with Voronoi calculation.')
            if not 'voronoi' in inputs:
                return self.exit_codes.ERROR_NO_VORONOI_CODE_GIVEN
        elif 'remote_data' in inputs:
            self.report('INFO: Found remote_data in input. Continue calculation without running voronoi step.')
            run_voronoi = False
        else:
            return self.exit_codes.ERROR_NOT_ENOUGH_INPUTS

        if 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_INVALID_VORONOI_CODE

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_INVALID_KKR_CODE

        # set params and remote folder to input if voronoi step is skipped
        if not run_voronoi:
            self.ctx.last_remote = inputs.remote_data
            num_parents = len(self.ctx.last_remote.get_incoming(node_class=CalcJobNode).all_link_labels())
            if num_parents == 0:
                pk_last_remote = self.ctx.last_remote.inputs.last_RemoteData.outputs.output_kkr_scf_wc_ParameterResults.get_dict().get('last_calc_nodeinfo').get('pk')
                last_calc = load_node(pk_last_remote)
                self.ctx.last_remote = last_calc.outputs.remote_folder
            try: # first try parent of remote data output of a previous calc.
                parent_params = get_parent_paranode(inputs.remote_data)
            except AttributeError:
                try: # next try to extract parameter from previous kkr_scf_wc output
                    parent_params = inputs.remote_data.inputs.last_RemoteData.inputs.calc_parameters
                except AttributeError:
                    return self.exit_codes.ERROR_NO_PARENT_PARAMS_FOUND
            if  'calc_parameters' in inputs:
                self.ctx.last_params = inputs.calc_parameters
                # TODO: check last_params consistency against parent_params
                #parent_params_dict = parent_params.get_dict()
                #for key, val in self.ctx.last_params.get_dict().iteritems():
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
        """
        run the voronoi step calling voro_start workflow
        """

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
            params = Dict(dict=defaults)
            params.label = 'default values'
        newparams = {}
        for key, val in defaults.items():
            if params.get_dict().get(key, None) is None:
                if key != 'RCLUSTZ': # this one is automatically set by kkr_startpot_wc so we skip it here
                    newparams[key] = val
                    self.report("INFO: Automatically added default values to KKR parameters: {} {}".format(key, val))
        if newparams!={}:
            for key, val in params.get_dict().items():
                if val is not None:
                    newparams[key] = val
            updatenode = Dict(dict=newparams)
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
                self.report("INFO: update parameters to: {}".format(para_check.get_set_values()))
                updatenode = Dict(dict=para_check.get_dict())
                updatenode.label = 'overwritten KKR input parameters'
                updatenode.description = 'Overwritten KKR input parameter to correct NSPIN to 2'
                params = update_params_wf(params, updatenode)

        # check consistency of input parameters before running calculation
        self.check_input_params(params, is_voronoi=True)

        # set parameters of voro_start sub workflow
        sub_wf_params_dict = kkr_startpot_wc.get_wf_defaults(silent=True)
        label, description = "voro_start_default_params", "workflow parameters for voro_start"
        if 'wf_parameters' in self.inputs:
            wf_params_input = self.inputs.wf_parameters.get_dict()
            num_updated = 0
            for key in list(sub_wf_params_dict.keys()):
                if key in list(wf_params_input.keys()):
                    val = wf_params_input[key]
                    sub_wf_params_dict[key] = val
                    num_updated += 1
            if num_updated > 0:
                label = "voro_start_updated_params"
        sub_wf_params = Dict(dict=sub_wf_params_dict)
        sub_wf_params.label = label
        sub_wf_params.description = description

        self.report('INFO: run voronoi step')
        self.report('INFO: using calc_params ({}): {}'.format(params, params.get_dict()))
        self.report('INFO: using wf_parameters ({}): {}'.format(sub_wf_params, sub_wf_params.get_dict()))

        wf_label= 'kkr_startpot (voronoi)'
        wf_desc = 'subworkflow to set up the input of a KKR calculation'
        future = self.submit(kkr_startpot_wc, kkr=kkrcode, voronoi=voronoicode,
                             calc_parameters=params, wf_parameters=sub_wf_params,
                             structure=structure, label=wf_label, description=wf_desc,
                             options=self.ctx.options_params_dict)

        return ToContext(voronoi=future, last_calc=future)


    def check_voronoi(self):
        """
        check output of kkr_startpot_wc workflow that creates starting potential, shapefun etc.
        """

        self.report("INFO: checking voronoi output")
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
            return self.exit_codes.ERROR_KKR_STARTPOT_FAILED

        self.report("INFO: done checking voronoi output")


    def condition(self):
        """
        check convergence condition
        """
        self.report("INFO: checking condition for kkr step")
        do_kkr_step = True
        stopreason = ''

        #increment KKR runs loop counter
        self.ctx.loop_count += 1

        # check if initial step was succesful
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
                #return self.exit_codes.ERROR_MAX_KKR_RESTARTS_REACHED
                return do_kkr_step

        self.report("INFO: done checking condition for kkr step (result={})".format(do_kkr_step))

        if not do_kkr_step:
            self.report("INFO: stopreason={}".format(stopreason))

        return do_kkr_step


    def update_kkr_params(self):
        """
        update set of KKR parameters (check for reduced mixing, change of mixing strategy, change of accuracy setting)
        """
        self.report("INFO: updating kkr param step")
        decrease_mixing_fac = False
        switch_agressive_mixing = False
        switch_higher_accuracy= False
        initial_settings = False

        # only do something other than somple mixing after first kkr run
        if self.ctx.loop_count != 1:
            # first determine if previous step was successful (otherwise try to find some rms value and decrease mixing to try again)
            if not self.ctx.kkr_step_success:
                decrease_mixing_fac = True
                self.report("INFO: last KKR calculation failed. trying decreasing mixfac")

            convergence_on_track = self.convergence_on_track()

            # check if calculation was on its way to converge
            if not convergence_on_track:
                decrease_mixing_fac = True
                self.report("INFO: last KKR did not converge. trying decreasing mixfac")
                # reset last_remote to last successful calculation
                for icalc in range(len(self.ctx.calcs))[::-1]:
                    self.report("INFO: last calc success? {} {}".format(icalc, self.ctx.KKR_steps_stats['success'][icalc]))
                    if self.ctx.KKR_steps_stats['success'][icalc]:
                        self.ctx.last_remote = self.ctx.calcs[icalc].outputs.remote_folder
                        break # exit loop if last_remote was found successfully
                    else:
                        self.ctx.last_remote = None
                # if no previous calculation was succesful take voronoi output or remote data from input (depending on the inputs)
                self.report("INFO: last_remote is None? {} {}".format(self.ctx.last_remote is None, 'structure' in self.inputs))
                if self.ctx.last_remote is None:
                    if 'structure' in self.inputs:
                        self.ctx.voronoi.outputs.last_voronoi_remote
                    else:
                        self.ctx.last_remote = self.inputs.remote_data
                # check if last_remote has finally been set and abort if this is not the case
                self.report("INFO: last_remote is still None? {}".format(self.ctx.last_remote is None))
                if self.ctx.last_remote is None:
                    return self.exit_codes.ERROR_LAST_REMOTE_NOT_FOUND

            # check if mixing strategy should be changed
            last_mixing_scheme = self.ctx.last_params.get_dict()['IMIX']
            if last_mixing_scheme is None:
                last_mixing_scheme = 0

            if convergence_on_track:
                last_rms = self.ctx.last_rms_all[-1]

                if last_rms < self.ctx.threshold_aggressive_mixing and last_mixing_scheme == 0:
                    switch_agressive_mixing = True
                    self.report("INFO: rms low enough, switch to agressive mixing")

                # check if switch to higher accuracy should be done
                if not self.ctx.kkr_higher_accuracy:
                    if self.ctx.kkr_converged or last_rms < self.ctx.threshold_switch_high_accuracy:
                        switch_higher_accuracy = True
                        self.report("INFO: rms low enough, switch to higher accuracy settings")
        else:
            initial_settings = True

        # if needed update parameters
        if decrease_mixing_fac or switch_agressive_mixing or switch_higher_accuracy or initial_settings or self.ctx.mag_init:

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
                if len(kkrdefaults_updated)>0:
                    self.report('ERROR: calc_parameters misses keys: {}'.format(missing_list))
                    return self.exit_codes.ERROR_CALC_PARAMTERS_INCOMPLETE
                else:
                    self.report('updated KKR parameter node with default values: {}'.format(kkrdefaults_updated))

            # step 2: change parameter (contained in new_params dictionary)
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
                    self.report('(strmixfax, mixreduce)= ({}, {})'.format(strmixfac, self.ctx.mixreduce))
                    self.report('type(strmixfax, mixreduce)= {} {}'.format(type(strmixfac), type(self.ctx.mixreduce)))
                    strmixfac = strmixfac * self.ctx.mixreduce
                    self.ctx.strmix = strmixfac
                    label += 'decreased_mix_fac_str (step {})'.format(self.ctx.loop_count)
                    description += 'decreased STRMIX factor by {}'.format(self.ctx.mixreduce)
                else:
                    self.report('(brymixfax, mixreduce)= ({}, {})'.format(brymixfac, self.ctx.mixreduce))
                    self.report('type(brymixfax, mixreduce)= {} {}'.format(type(brymixfac), type(self.ctx.mixreduce)))
                    brymixfac = brymixfac * self.ctx.mixreduce
                    self.ctx.brymix = brymixfac
                    label += 'decreased_mix_fac_bry'
                    description += 'decreased BRYMIX factor by {}'.format(self.ctx.mixreduce)

            #add mixing factor
            new_params['STRMIX'] = strmixfac
            new_params['BRYMIX'] = brymixfac

            if switch_agressive_mixing:
                last_mixing_scheme = 5
                label += ' switched_to_agressive_mixing'
                description += ' switched to agressive mixing scheme (IMIX={})'.format(last_mixing_scheme)

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

            # slightly increase temperature if previous calculation was unsuccessful for the second time
            if decrease_mixing_fac and not self.convergence_on_track():
                self.report('INFO: last calculation did not converge and convergence not on track. Try to increase temperature by 50K.')
                convergence_settings['tempr'] += 50.
                label += ' TEMPR+50K'
                description += ' with increased temperature of 50K'

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

            # initial magnetization
            if initial_settings and self.ctx.mag_init:
                if self.ctx.hfield <= 0:
                    self.report('\nWARNING: magnetization initialization chosen but hfield is zero. Automatically change back to default value (hfield={})\n'.format(self._wf_default['hfield']))
                    self.ctx.hfield = self._wf_default['hfield']
                xinipol = self.ctx.xinit
                if xinipol is None:
                    # find structure to determine needed length on xinipol
                    if 'structure' in self.inputs:
                        struc = self.inputs.structure
                    else:
                        struc, voro_parent = VoronoiCalculation.find_parent_structure(self.ctx.last_remote)
                    natom = len(struc.sites)
                    xinipol = ones(natom)
                new_params['LINIPOL'] = True
                new_params['HFIELD'] = self.ctx.hfield
                new_params['XINIPOL'] = xinipol
            elif self.ctx.mag_init and self.ctx.mag_init_step_success: # turn off initialization after first (successful) iteration
                new_params['LINIPOL'] = False
                new_params['HFIELD'] = 0.0
            elif not self.ctx.mag_init:
                self.report("INFO: mag_init is False. Overwrite 'HFIELD' to '0.0' and 'LINIPOL' to 'False'.")
                # reset mag init to avoid reinitializing
                new_params['HFIELD'] = 0.0
                new_params['LINIPOL'] = False

            # set nspin to 2 if mag_init is used
            if self.ctx.mag_init:
                nspin_in = para_check.get_value('NSPIN')
                if nspin_in is None:
                    nspin_in = 1
                if nspin_in < 2:
                    self.report('WARNING: found NSPIN=1 but for maginit needs NPIN=2. Overwrite this automatically')
                    new_params['NSPIN'] = 2

            # step 2.2 update values
            try:
                for key, val in new_params.items():
                    para_check.set_value(key, val, silent=True)
            except:
                return self.exit_codes.ERROR_PARAM_UPDATE_FAILED

            # step 3:
            self.report("INFO: update parameters to: {}".format(para_check.get_set_values()))
            updatenode = Dict(dict=para_check.get_dict())
            updatenode.label = label
            updatenode.description = description

            paranode_new = update_params_wf(self.ctx.last_params, updatenode)
            self.ctx.last_params = paranode_new
        else:
            self.report("INFO: reuse old settings")

        self.report("INFO: done updating kkr param step")


    def run_kkr(self):
        """
        submit a KKR calculation
        """
        self.report("INFO: setting up kkr calculation step {}".format(self.ctx.loop_count))


        label = 'KKR calculation step {} (IMIX={})'.format(self.ctx.loop_count, self.ctx.last_mixing_scheme)
        description = 'KKR calculation of step {}, using mixing scheme {}'.format(self.ctx.loop_count, self.ctx.last_mixing_scheme)
        code = self.inputs.kkr
        remote = self.ctx.last_remote
        params = self.ctx.last_params
        options = {"max_wallclock_seconds": self.ctx.walltime_sec,
                   "resources": self.ctx.resources,
                   "queue_name" : self.ctx.queue}
        if self.ctx.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.ctx.custom_scheduler_commands
        inputs = get_inputs_kkr(code, remote, options, label, description, parameters=params, serial=(not self.ctx.use_mpi))

        # run the KKR calculation
        self.report('INFO: doing calculation')
        kkr_run = self.submit(KkrCalculation, **inputs)

        return ToContext(kkr=kkr_run, last_calc=kkr_run)


    def inspect_kkr(self):
        """
        check for convergence and store some of the results of the last calculation to context
        """
        self.report("INFO: inspecting kkr results step")

        self.ctx.calcs.append(self.ctx.last_calc)
        self.ctx.kkr_step_success = True

        # check calculation state
        if not self.ctx.last_calc.is_finished_ok:
            self.ctx.kkr_step_success = False
            self.report("ERROR: last calculation not finished correctly")

        self.report("INFO: kkr_step_success: {}".format(self.ctx.kkr_step_success))

        # extract convergence info about rms etc. (used to determine convergence behavior)
        try:
            self.report("INFO: trying to find output of last_calc: {}".format(self.ctx.last_calc))
            last_calc_output = self.ctx.last_calc.outputs.output_parameters.get_dict()
            found_last_calc_output = True
        except:
            found_last_calc_output = False
        self.report("INFO: found_last_calc_output: {}".format(found_last_calc_output))

        # try yo extract remote folder
        try:
            self.ctx.last_remote = self.ctx.kkr.outputs.remote_folder
        except:
            self.ctx.last_remote = None
            self.ctx.kkr_step_success = False

        self.report("INFO: last_remote: {}".format(self.ctx.last_remote))

        if self.ctx.kkr_step_success and found_last_calc_output:
            # check convergence
            self.ctx.kkr_converged = last_calc_output['convergence_group']['calculation_converged']
            # check rms
            self.ctx.rms.append(last_calc_output['convergence_group']['rms'])
            rms_all_iter_last_calc = list(last_calc_output['convergence_group']['rms_all_iterations'])
            #check charge neutrality
            self.ctx.neutr.append(last_calc_output['convergence_group']['charge_neutrality'])
            neutr_all_iter_last_calc = list(last_calc_output['convergence_group']['charge_neutrality_all_iterations'])

            # add lists of last iterations
            self.ctx.last_rms_all = rms_all_iter_last_calc
            self.ctx.last_neutr_all = neutr_all_iter_last_calc
            if self.ctx.kkr_step_success and self.convergence_on_track():
                self.ctx.rms_all_steps += rms_all_iter_last_calc
                self.ctx.neutr_all_steps += neutr_all_iter_last_calc
        else:
            self.ctx.kkr_converged = False

        self.report("INFO: kkr_converged: {}".format(self.ctx.kkr_converged))
        self.report("INFO: rms: {}".format(self.ctx.rms))
        self.report("INFO: last_rms_all: {}".format(self.ctx.last_rms_all))
        #self.report("INFO: rms_all_steps: {}".format(self.ctx.rms_all_steps))
        self.report("INFO: charge_neutrality: {}".format(self.ctx.neutr))
        self.report("INFO: last_neutr_all: {}".format(self.ctx.last_neutr_all))
        #self.report("INFO: neutr_all_steps: {}".format(self.ctx.neutr_all_steps))

        # turn off initial magnetization once one step was successful (update_kkr_params) used in
        if self.ctx.mag_init and self.ctx.kkr_step_success:
            self.ctx.mag_init_step_success = True

        # TODO: extract something else (maybe total energy, charge neutrality, magnetisation)?

        # store some statistics used to print table in the end of the report
        self.ctx.KKR_steps_stats['success'].append(self.ctx.kkr_step_success)
        try:
            isteps = self.ctx.last_calc.outputs.output_parameters.get_dict()['convergence_group']['number_of_iterations']
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

        self.ctx.KKR_steps_stats['isteps'].append(isteps)
        self.ctx.KKR_steps_stats['imix'].append(self.ctx.last_mixing_scheme)
        self.ctx.KKR_steps_stats['mixfac'].append(mixfac)
        self.ctx.KKR_steps_stats['qbound'].append(qbound)
        self.ctx.KKR_steps_stats['high_sett'].append(self.ctx.kkr_higher_accuracy)
        self.ctx.KKR_steps_stats['first_rms'].append(first_rms)
        self.ctx.KKR_steps_stats['last_rms'].append(last_rms)
        self.ctx.KKR_steps_stats['first_neutr'].append(first_neutr)
        self.ctx.KKR_steps_stats['last_neutr'].append(last_neutr)
        self.ctx.KKR_steps_stats['pk'].append(self.ctx.last_calc.pk)
        self.ctx.KKR_steps_stats['uuid'].append(self.ctx.last_calc.uuid)

        self.report("INFO: done inspecting kkr results step")


    def convergence_on_track(self):
        """
        Check if convergence behavior of the last calculation is on track (i.e. going down)
        """
        on_track = True
        threshold = 5. # used to check condition if at least one of charnge_neutrality, rms-error goes down fast enough

        # first check if previous calculation was stopped due to reaching the QBOUND limit
        try:
            calc_reached_qbound = self.ctx.last_calc.outputs.output_parameters.get_dict()['convergence_group']['calculation_converged']
        except AttributeError: # captures error when last_calc dies not have an output node
            calc_reached_qbound = False
        except KeyError: # captures
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
            r, n = last_rms/first_rms, last_neutr/first_neutr
            self.report("INFO convergence check: first/last rms {}, {}; first/last neutrality {}, {}".format(first_rms, last_rms, first_neutr, last_neutr))
            if r < 1 and n < 1:
                self.report("INFO convergence check: both rms and neutrality go down")
                on_track = True
            elif n > threshold or r > threshold:
                self.report("INFO convergence check: rms or neutrality goes up too fast, convergence is not expected")
                on_track = False
            elif n*r < 1:
                self.report("INFO convergence check: either rms goes up and neutrality goes down or vice versa")
                self.report("INFO convergence check: but product goes down fast enough")
                on_track = True
            elif len(self.ctx.last_rms_all) ==1:
                self.report("INFO convergence check: already converged after single iteration")
                on_track = True
            else:
                self.report("INFO convergence check: rms or neutrality do not shrink fast enough, convergence is not expected")
                on_track = False
        elif calc_reached_qbound:
            self.report("INFO convergence check: calculation reached QBOUND")
            on_track = True
        else:
            self.report("INFO convergence check: calculation unsuccessful")
            on_track = False

        self.report("INFO convergence check result: {}".format(on_track))

        return on_track



    def return_results(self):
        """
        return the results of the calculations
        This shoudl run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """

        self.report("INFO: entering return_results")

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
                self.ctx.warnings.append('cound not get pk of calc {}'.format(calc))


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
        self.report("INFO: collect outputnode_dict")
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['material'] = self.ctx.formula
        outputnode_dict['loop_count'] = self.ctx.loop_count
        outputnode_dict['warnings'] = self.ctx.warnings
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['last_params_nodeinfo'] = {'uuid':last_params_uuid, 'pk':last_params_pk}
        outputnode_dict['last_remote_nodeinfo'] = {'uuid':last_remote_uuid, 'pk':last_remote_pk}
        outputnode_dict['last_calc_nodeinfo'] = {'uuid':last_calc_uuid, 'pk':last_calc_pk}
        outputnode_dict['pks_all_calcs'] = all_pks
        outputnode_dict['errors'] = self.ctx.errors
        outputnode_dict['convergence_value'] = last_rms
        outputnode_dict['convergence_values_all_steps'] = array(self.ctx.rms_all_steps)
        outputnode_dict['convergence_values_last_step'] = array(self.ctx.last_rms_all)
        outputnode_dict['charge_neutrality'] = last_neutr
        outputnode_dict['charge_neutrality_all_steps'] = array(self.ctx.neutr_all_steps)
        outputnode_dict['charge_neutrality_last_step'] = array(self.ctx.last_neutr_all)
        outputnode_dict['dos_check_ok'] = self.ctx.dos_ok
        outputnode_dict['convergence_reached'] = self.ctx.kkr_converged
        outputnode_dict['voronoi_step_success'] = self.ctx.voro_step_success
        outputnode_dict['kkr_step_success'] = self.ctx.kkr_step_success
        outputnode_dict['used_higher_accuracy'] = self.ctx.kkr_higher_accuracy

        # report the status
        if self.ctx.successful:
            self.report('STATUS: Done, the convergence criteria are reached.\n'
                        'INFO: The charge density of the KKR calculation pk= {} '
                        'converged after {} KKR runs and {} iterations to {} \n'
                        ''.format(last_calc_pk, self.ctx.loop_count, self.ctx.loop_count, last_rms))
        else: # Termination ok, but not converged yet...
            if self.ctx.abort: # some error occured, donot use the output.
                self.report('STATUS/ERROR: I abort, see logs and '
                            'erros/warning/hints in output_kkr_scf_wc_para')
            else:
                self.report('STATUS/WARNING: Done, the maximum number of runs '
                            'was reached or something failed.\n INFO: The '
                            'charge density of the KKR calculation pk= '
                            'after {} KKR runs and {} iterations is {} "me/bohr^3"\n'
                            ''.format(self.ctx.loop_count, sum(self.ctx.KKR_steps_stats.get('isteps')), last_rms))

        # create results  node
        self.report("INFO: create results node") #: {}".format(outputnode_dict))
        outputnode_t = Dict(dict=outputnode_dict)
        outputnode_t.label = 'kkr_scf_wc_results'
        outputnode_t.description = 'Contains results of workflow (e.g. workflow version number, info about success of wf, lis tof warnings that occured during execution, ...)'


        # collect nodes in outputs dictionary
        # call helper function to create output nodes in correct AiiDA graph structure
        if last_calc_out is not None and last_RemoteData is not None and last_InputParameters is not None:
            if results_vorostart is not None and starting_dosdata_interpol is not None and final_dosdata_interpol is not None:
                outdict = create_scf_result_node(outpara=outputnode_t,
                                                 last_calc_out=last_calc_out,
                                                 last_RemoteData=last_RemoteData,
                                                 last_InputParameters=last_InputParameters,
                                                 final_dosdata_interpol=final_dosdata_interpol,
                                                 starting_dosdata_interpol=starting_dosdata_interpol,
                                                 results_vorostart=results_vorostart)
            elif results_vorostart is not None and starting_dosdata_interpol is not None:
                outdict = create_scf_result_node(outpara=outputnode_t,
                                                 last_calc_out=last_calc_out,
                                                 last_RemoteData=last_RemoteData,
                                                 last_InputParameters=last_InputParameters,
                                                 starting_dosdata_interpol=starting_dosdata_interpol,
                                                 results_vorostart=results_vorostart)
            elif results_vorostart is not None:
                outdict = create_scf_result_node(outpara=outputnode_t,
                                                 last_calc_out=last_calc_out,
                                                 last_RemoteData=last_RemoteData,
                                                 last_InputParameters=last_InputParameters,
                                                 results_vorostart=results_vorostart)
            elif final_dosdata_interpol is not None:
                outdict = create_scf_result_node(outpara=outputnode_t,
                                                 last_calc_out=last_calc_out,
                                                 last_RemoteData=last_RemoteData,
                                                 last_InputParameters=last_InputParameters,
                                                 final_dosdata_interpol=final_dosdata_interpol)
            else:
                outdict = create_scf_result_node(outpara=outputnode_t,
                                                 last_calc_out=last_calc_out,
                                                 last_RemoteData=last_RemoteData,
                                                 last_InputParameters=last_InputParameters)
        else:
            outdict = create_scf_result_node(outpara=outputnode_t)

        for link_name, node in outdict.items():
            #self.report("INFO: storing node {} {} with linkname {}".format(type(node), node, link_name))
            self.out(link_name, node)


        # print results table for overview
        # table layout:
        message = "INFO: overview of the result:\n\n"
        message += "|------|---------|--------|------|--------|-------------------|-----------------|-----------------|-------------\n"
        message += "| irun | success | isteps | imix | mixfac | accuracy settings |       rms       | abs(neutrality) | pk and uuid \n"
        message += "|      |         |        |      |        | qbound  | higher? | first  |  last  | first  |  last  |             \n"
        message += "|------|---------|--------|------|--------|---------|---------|--------|--------|--------|--------|-------------\n"
        #| %6i  | %9s     | %8i    | %6i  | %.2e   | %.3e    | %9s     | %.2e   |  %.2e  |  %.2e  |  %.2e  |
        KKR_steps_stats = self.ctx.KKR_steps_stats
        for irun in range(len(KKR_steps_stats.get('success'))):
            KKR_steps_stats.get('first_neutr')[irun] = abs(KKR_steps_stats.get('first_neutr')[irun])
            KKR_steps_stats.get('last_neutr')[irun] = abs(KKR_steps_stats.get('last_neutr')[irun])
            message += "|%6i|%9s|%8i|%6i|%.2e|%.3e|%9s|%.2e|%.2e|%.2e|%.2e|"%(irun+1,
                          KKR_steps_stats.get('success')[irun], KKR_steps_stats.get('isteps')[irun],
                          KKR_steps_stats.get('imix')[irun], KKR_steps_stats.get('mixfac')[irun],
                          KKR_steps_stats.get('qbound')[irun], KKR_steps_stats.get('high_sett')[irun],
                          KKR_steps_stats.get('first_rms')[irun], KKR_steps_stats.get('last_rms')[irun],
                          KKR_steps_stats.get('first_neutr')[irun], KKR_steps_stats.get('last_neutr')[irun])
            message += " {} | {}\n".format(KKR_steps_stats.get('pk')[irun], KKR_steps_stats.get('uuid')[irun])
            """
            message += "#|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".format(irun+1,
                          KKR_steps_stats.get('success')[irun], KKR_steps_stats.get('isteps')[irun],
                          KKR_steps_stats.get('imix')[irun], KKR_steps_stats.get('mixfac')[irun],
                          KKR_steps_stats.get('qbound')[irun], KKR_steps_stats.get('high_sett')[irun],
                          KKR_steps_stats.get('first_rms')[irun], KKR_steps_stats.get('last_rms')[irun],
                          KKR_steps_stats.get('first_neutr')[irun], KKR_steps_stats.get('last_neutr')[irun])
            """
        self.report(message)

        self.report("\nINFO: done with kkr_scf workflow!\n")


    def check_input_params(self, params, is_voronoi=False):
        """
        Checks input parameter consistency and aborts wf if check fails.
        """
        if params is None:
            return self.exit_codes.ERROR_NO_CALC_PARAMETERS_GIVEN
        else:
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
                return self.exit_codes.ERROR_CALC_PARAMTERS_INCONSISTENT

            # step 2: check if all mandatory keys are there
            missing_list = para_check.get_missing_keys(use_aiida=True)
            if missing_list != []:
                all_defaults = True
                for key in missing_list:
                    if key not in kkrparams.get_KKRcalc_parameter_defaults()[0]:
                        all_defaults = False
                if not all_defaults:
                    self.report('ERROR: calc_parameters given are not consistent! Missing mandatory keys: {}'.format(missing_list))
                    return self.exit_codes.ERROR_CALC_PARAMETERS_INCOMPLETE


    def get_dos(self):
        """
        call to dos sub workflow passing the appropriate input and submitting the calculation
        """
        if self.ctx.check_dos:
            self.report("INFO: Doing final DOS calculation")

            # fix emin/emax to include emin, ef of scf contour
            # remember: efermi, emin and emax are in internal units (Ry) but delta_e is in eV!
            eV2Ry = 1./get_Ry2eV()
            emin = self.ctx.dos_params['emin'] # from dos params
            emin_cont = self.ctx.last_calc.outputs.output_parameters.get_dict().get('energy_contour_group').get('emin') # from contour (sets limit of dos emin!)
            if emin_cont - self.ctx.delta_e*eV2Ry < emin:
                self.ctx.dos_params['emin'] = emin_cont - self.ctx.delta_e*eV2Ry
                self.report("INFO: emin ({} Ry) - delta_e ({} Ry) smaller than emin ({} Ry) of voronoi output. Setting automatically to {}Ry".format(emin_cont, self.ctx.delta_e*eV2Ry,  emin, emin_cont-self.ctx.delta_e*eV2Ry))
            self.ctx.efermi = get_ef_from_potfile(self.ctx.last_calc.outputs.retrieved.get_abs_path('out_potential'))
            emax = self.ctx.dos_params['emax']
            if emax < self.ctx.efermi + self.ctx.delta_e*eV2Ry:
                self.ctx.dos_params['emax'] = self.ctx.efermi + self.ctx.delta_e*eV2Ry
                self.report("INFO: self.ctx.efermi ({} Ry) + delta_e ({} Ry) larger than emax ({} Ry). Setting automatically to {}Ry".format(self.ctx.efermi, self.ctx.delta_e*eV2Ry, emax, self.ctx.efermi+self.ctx.delta_e*eV2Ry))

            # take subset of input and prepare parameter node for dos workflow
            wfdospara_dict = {'queue_name' : self.ctx.queue,
                              'resources': self.ctx.resources,
                              'max_wallclock_seconds' : self.ctx.walltime_sec,
                              'use_mpi' : self.ctx.use_mpi,
                              'custom_scheduler_commands' : self.ctx.custom_scheduler_commands,
                              'dos_params' : self.ctx.dos_params}
            wfdospara_node = Dict(dict=wfdospara_dict)
            wfdospara_node.label = 'DOS params'
            wfdospara_node.description = 'DOS parameter set for final DOS calculation of kkr_scf_wc'

            code = self.inputs.kkr
            remote = self.ctx.last_calc.outputs.remote_folder
            wf_label= ' final DOS calculation'
            wf_desc = ' subworkflow of a DOS calculation'
            future = self.submit(kkr_dos_wc, kkr=code, remote_data=remote,
                                 wf_parameters=wfdospara_node, label=wf_label,
                                 description=wf_desc, options=self.ctx.options_params_dict)

            return ToContext(doscal=future)


    def check_dos(self):
        """
        checks if dos of final potential is ok
        """
        # initialize dos_ok variable
        self.ctx.dos_ok = True

        # first check if dos should be checked or if the test is skipped
        if not self.ctx.check_dos:
            self.report("INFO: skipping DOS check")
            return

        self.report("INFO: checking DOS for consistency (EMIN position, negative DOS, etc.)")

        # check parser output
        doscal = self.ctx.doscal
        try:
            dos_outdict = doscal.outputs.results_wf.get_dict()
            if not dos_outdict['successful']:
                self.report("ERROR: DOS workflow unsuccessful")
                self.ctx.dos_ok = False
                return self.exit_codes.ERROR_DOS_RUN_UNSUCCESFUL

            if dos_outdict['list_of_errors'] != []:
                self.report("ERROR: DOS wf output contains errors: {}".format(dos_outdict['list_of_errors']))
                self.ctx.dos_ok = False
                return self.exit_codes.ERROR_DOS_RUN_UNSUCCESFUL
        except AttributeError:
            self.ctx.dos_ok = False
            return self.exit_codes.ERROR_DOS_RUN_UNSUCCESFUL

        # check for negative DOS
        try:
            dosdata = doscal.outputs.dos_data
            natom = self.ctx.last_calc.outputs.output_parameters.get_dict()['number_of_atoms_in_unit_cell']
            nspin = dos_outdict['nspin']

            ener = dosdata.get_x()[1] # shape= natom*nspin, nept
            totdos = dosdata.get_y()[0][1] # shape= natom*nspin, nept

            if len(ener) != nspin*natom:
                self.report("ERROR: DOS output shape does not fit nspin, natom information: len(energies)={}, natom={}, nspin={}".format(len(ener), natom, nspin))
                self.ctx.doscheck_ok = False
                return self.exit_codes.ERROR_DOS_RUN_UNSUCCESFUL

            # deal with snpin==1 or 2 cases and check negtive DOS
            for iatom in range(natom//nspin):
                for ispin in range(nspin):
                    x, y = ener[iatom*nspin+ispin], totdos[iatom*nspin+ispin]
                    if nspin == 2 and ispin == 0:
                        y = -y
                    if y.min() < 0:
                        self.report("INFO: negative DOS value found in (atom, spin)=({},{}) at iteration {}: {}".format(iatom, ispin, self.ctx.loop_count, y.min()))
                        self.ctx.dos_ok = False

            # check starting EMIN
            dosdata_interpol = doscal.outputs.dos_data_interpol

            ener = dosdata_interpol.get_x()[1] # shape= natom*nspin, nept
            totdos = dosdata_interpol.get_y()[0][1] # shape= natom*nspin, nept
            Ry2eV = get_Ry2eV()

            for iatom in range(natom//nspin):
                for ispin in range(nspin):
                    x, y = ener[iatom*nspin+ispin], totdos[iatom*nspin+ispin]
                    xrel = abs(x-(self.ctx.dos_params_dict['emin']-self.ctx.efermi)*Ry2eV)
                    mask_emin = where(xrel==xrel.min())
                    ymin = abs(y[mask_emin])
                    if ymin > self.ctx.threshold_dos_zero:
                        self.report("INFO: DOS at emin not zero! {}>{}".format(ymin,self.ctx.threshold_dos_zero))
                        self.ctx.dos_ok = False
        except AttributeError:
            self.ctx.dos_ok = False



@wf
def create_scf_result_node(**kwargs):
    """
    This is a pseudo wf, to create the right graph structure of AiiDA.
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
        if key == 'outpara': #  should always be there
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
        outputnode = outpara
        outputnode.label = 'workflow_Results'
        outputnode.description = ('Contains self-consistency results and '
                                  'information of an kkr_scf_wc run.')
        outdict['output_kkr_scf_wc_ParameterResults'] = outputnode

    if has_last_calc_out_dict:
        outputnode = last_calc_out_dict
        outputnode.label = 'last_calc_out'
        outputnode.description = ('Contains the Results Parameter node from the output '
                                   'of the last calculation done in the workflow.')
        outdict['last_calc_out'] = outputnode

    if has_last_RemoteData:
        outputnode = last_RemoteData_dict
        outputnode.label = 'last_RemoteData'
        outputnode.description = ('Contains a link to the latest remote data node '
                                   'where the output of the calculation can be accessed.')
        outdict['last_RemoteData'] = outputnode

    if has_last_InputParameters:
        outputnode = last_InputParameters_dict
        outputnode.label = 'last_InputParameters'
        outputnode.description = ('Contains the latest parameter data node '
                                   'where the input of the last calculation can be found.')
        outdict['last_InputParameters'] = outputnode

    if has_vorostart_output:
        outputnode = vorostart_output_dict
        outputnode.label = 'results_vorostart'
        outputnode.description = ('Contains the results parameter data node '
                                   'of the vorostart sub-workflow (sets up starting portentials).')
        outdict['results_vorostart'] = outputnode

    if has_starting_dos:
        outputnode = start_dosdata_interpol_dict
        outputnode.label = 'starting_dosdata_interpol'
        outputnode.description = ('Contains the interpolated DOS data note, computed '
                                   'from the starting portential.')
        outdict['starting_dosdata_interpol'] = outputnode

    if has_final_dos:
        outputnode = final_dosdata_interpol_dict
        outputnode.label = 'final_dosdata_interpol'
        outputnode.description = ('Contains the interpolated DOS data note, computed '
                                   'from the converged potential.')
        outdict['final_dosdata_interpol'] = outputnode

    return outdict
