#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for converging a kkr calculation and
some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory
from aiida.work.workchain import WorkChain, while_, if_, ToContext
from aiida.work.run import submit
from aiida.work import workfunction as wf
from aiida.work.process_registry import ProcessRegistry
from aiida.common.datastructures import calc_states
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import (test_and_get_codenode, get_inputs_kkr, 
                                                  get_parent_paranode, update_params_wf)
from aiida_kkr.workflows.voro_start import kkr_startpot_wc

__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.4"
__contributors__ = (u"Jens Broeder", u"Philipp Rüßmann")


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()
VoronoiProcess = VoronoiCalculation.process()

class kkr_scf_wc(WorkChain):
    """
    Workchain for converging a KKR calculation (SCF).

    It converges the charge density, potential?.
    Two paths are possible:

    (1) Start from a structure and run a voronoi calculation first,
    optional with calc_parameters
    (2) Start from an existing Voronoi or KKR calculation, with a remoteData

    :param wf_parameters: (ParameterData), Workchain Spezifications
    :param structure: (StructureData), Crystal structure
    :param calc_parameters: (ParameterData), Vornoi/Kkr Parameters
    :param remote_data: (RemoteData), from a KKR, or Vornoi calculation
    :param voronoi: (Code)
    :param kkr: (Code)

    :return output_kkr_scf_wc_para: (ParameterData), Information of workflow results 
        like Success, last result node, list with convergence behavior

    minimum input example:
    1. Code1, Code2, Structure, (Parameters), (wf_parameters)
    2. Code2, remote_data, (Parameters), (wf_parameters)

    maximum input example:
    1. Code1, Code2, Structure, Parameters
        wf_parameters: {
                        'queue' : String,
                        'resources' : dict({"num_machines": int, "num_mpiprocs_per_machine" : int})
                        'walltime' : int}
    2. Code2, (remote-data), wf_parameters as in 1.

    Hints:
    1. This workflow does not work with local codes!
    """

    _workflowversion = __version__
    _wf_default = {'kkr_runmax': 4,                           # Maximum number of kkr jobs/starts (defauld iterations per start)
                   'convergence_criterion' : 10**-6,          # Stop if charge denisty is converged below this value
                   'queue_name' : '',                         # Queue name to submit jobs too
                   'resources': {"num_machines": 1},          # resources to allowcate for the job
                   'walltime_sec' : 60*60,                    # walltime after which the job gets killed (gets parsed to KKR)
                   'use_mpi' : False,                         # execute KKR with mpi or without
                   'custom_scheduler_commands' : '',          # some additional scheduler commands 
                   'check_dos' : True,                        # check starting DOS for inconsistencies
                   'dos_params' : {'emax': 1,                 # DOS params: maximum of enery contour
                                   'emin': -1,                # DOS params: minimum of energy contour
                                   'kmesh': [50, 50, 50],     # DOS params: kmesh for dos run (typically higher than for scf contour)
                                   'nepts': 61,               # DOS params: number of energy points in DOS contour
                                   'tempr': 200},             # DOS params: temperature
                   'mag_init' : False,                        # initialize and converge magnetic calculation
                   'mixreduce': 0.5,                          # reduce mixing factor by this factor if calculaito fails due to too large mixing
                   'threshold_aggressive_mixing': 8*10**-3,   # threshold after which agressive mixing is used
                   'strmix': 0.01,                            # mixing factor of simple mixing
                   'brymix': 0.01,                            # mixing factor of aggressive mixing
                   'nsteps': 30,                              # number of iterations done per KKR calculation
                   'convergence_setting_coarse': {            # setting of the coarse preconvergence
                        'npol': 7, 
                        'n1': 3,
                        'n2': 11,
                        'n3': 3,
                        'tempr': 1400.0,
                        'kmesh': [10, 10, 10]},
                   'threshold_switch_high_accuracy': 1*10**-3,# threshold after which final conversion settings are used
                   'convergence_setting_fine': {              # setting of the final convergence (lower tempr, 48 epts, denser k-mesh)
                        'npol': 5, 
                        'n1': 7,
                        'n2': 29,
                        'n3': 7,
                        'tempr': 473.0,
                        'kmesh': [30, 30, 30]}
                   }
             
    # intended to guide user interactively in setting up a valid wf_params node
    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """
        print('Version of workflow: {}'.format(self._workflowversion))
        return self._wf_default
    
    
    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow. 
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_scf_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict=cls._wf_default))
        spec.input("structure", valid_type=StructureData, required=False)
        spec.input("calc_parameters", valid_type=ParameterData, required=False)
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
            # finalize calculation and create output nodes
            cls.return_results
        )


    def start(self):
        """
        init context and some parameters
        """
        self.report('INFO: started KKR convergence workflow version {}\n'
                    'INFO: Workchain node identifiers: {}'
                    ''.format(self._workflowversion, ProcessRegistry().current_calc_node))

        ####### init #######

        # internal para /control para
        self.ctx.loop_count = 0
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

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()

        if wf_dict == {}:
            wf_dict = self._wf_default
            self.report('INFO: using default wf parameter')

        # set values from input, or defaults
        self.ctx.use_mpi = wf_dict.get('use_mpi', self._wf_default['use_mpi'])
        self.ctx.max_number_runs = wf_dict.get('kkr_runmax', self._wf_default['kkr_runmax'])
        self.ctx.resources = wf_dict.get('resources', self._wf_default['resources'])
        self.ctx.walltime_sec = wf_dict.get('walltime_sec', self._wf_default['walltime_sec'])
        self.ctx.queue = wf_dict.get('queue_name', self._wf_default['queue_name'])
        self.ctx.custom_scheduler_commands = wf_dict.get('custom_scheduler_commands', self._wf_default['custom_scheduler_commands'])
        self.ctx.description_wf = self.inputs.get('_description', 'Workflow for '
                                                  'a KKR scf calculation starting '
                                                  'either from a structure with '
                                                  'automatic voronoi calculation '
                                                  'or a valid RemoteData node of '
                                                  'a previous calculation')
        self.ctx.label_wf = self.inputs.get('_label', 'kkr_scf_wc')
        self.ctx.strmix = wf_dict.get('strmix', self._wf_default['strmix'])
        self.ctx.brymix = wf_dict.get('brymix', self._wf_default['brymix'])
        self.ctx.check_dos = wf_dict.get('check_dos', self._wf_default['check_dos'])
        self.ctx.dos_params = wf_dict.get('dos_params', self._wf_default['dos_params'])
        self.ctx.convergence_criterion = wf_dict.get('convergence_criterion', self._wf_default['convergence_criterion'])
        self.ctx.convergence_setting_coarse = wf_dict.get('convergence_setting_coarse', self._wf_default['convergence_setting_coarse'])
        self.ctx.convergence_setting_fine = wf_dict.get('convergence_setting_fine', self._wf_default['convergence_setting_fine'])
        self.ctx.mag_init = wf_dict.get('mag_init', self._wf_default['mag_init'])
        self.ctx.mixreduce = wf_dict.get('mixreduce', self._wf_default['mixreduce'])
        self.ctx.nsteps = wf_dict.get('nsteps', self._wf_default['nsteps'])
        self.ctx.threshold_aggressive_mixing = wf_dict.get('threshold_aggressive_mixing', self._wf_default['threshold_aggressive_mixing'])
        self.ctx.threshold_switch_high_accuracy = wf_dict.get('threshold_switch_high_accuracy', self._wf_default['threshold_switch_high_accuracy'])
        
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
                    '\nAdditional parameter\n'
                    'check DOS between runs: {}\n'
                    'DOS parameters: {}\n'
                    'init magnetism in first step: {}\n'
                    'factor reduced mixing if failing calculation: {}\n'
                    ''.format(self.ctx.use_mpi, self.ctx.max_number_runs, 
                                self.ctx.resources, self.ctx.walltime_sec, 
                                self.ctx.queue, self.ctx.custom_scheduler_commands, 
                                self.ctx.description_wf, self.ctx.label_wf,
                                self.ctx.strmix, self.ctx.brymix, self.ctx.nsteps, 
                                self.ctx.convergence_criterion, 
                                self.ctx.threshold_aggressive_mixing, 
                                self.ctx.threshold_switch_high_accuracy,
                                self.ctx.convergence_setting_coarse, 
                                self.ctx.convergence_setting_fine, self.ctx.check_dos, 
                                self.ctx.dos_params, self.ctx.mag_init, self.ctx.mixreduce)
                    )

        # return para/vars
        self.ctx.successful = True
        self.ctx.rms = []
        self.ctx.warnings = []
        self.ctx.errors = []
        self.ctx.formula = ''

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
                error = 'ERROR: StructureData was provided, but no voronoi code was provided'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
        elif 'remote_data' in inputs:
            self.report('INFO: Found remote_data in input. Continue calculation without running voronoi step.')
            run_voronoi = False
        else:
            error = 'ERROR: No StructureData nor remote_data was provided as Input'
            self.ctx.errors.append(error)
            self.control_end_wc(error)

        if 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                error = ("The code you provided for voronoi does not "
                         "use the plugin kkr.voro")
                self.control_end_wc(error)

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                error = ("The code you provided for kkr does not "
                         "use the plugin kkr.kkr")
                self.control_end_wc(error)
                
        # set params and remote folder to input if voronoi step is skipped
        if not run_voronoi:
            self.ctx.last_remote = inputs.remote_data
            self.ctx.last_params = get_parent_paranode(inputs.remote_data)
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
            params = None # TODO: use defaults?
        
        # check consistency of input parameters before running calculation
        self.check_input_params(params, is_voronoi=True)
        
        # set parameters of voro_start sub workflow
        sub_wf_params_dict = kkr_startpot_wc.get_wf_defaults()
        label, description = "voro_start_default_params", "workflow parameters for voro_start"
        if 'wf_parameters' in self.inputs:
            wf_params_input = self.inputs.wf_parameters.get_dict()
            num_updated = 0
            for key in sub_wf_params_dict.keys():
                if key in wf_params_input.keys():
                    val = wf_params_input[key]
                    sub_wf_params_dict[key] = val
                    num_updated += 1
            if num_updated > 0:
                label = "voro_start_updated_params"
        sub_wf_params = ParameterData(dict=sub_wf_params_dict)
        sub_wf_params.label = label
        sub_wf_params.description = description
        
        
        self.report('INFO: run voronoi step')
        self.report('INFO: using calc_params ({}): {}'.format(params, params.get_dict()))
        self.report('INFO: using wf_parameters ({}): {}'.format(sub_wf_params, sub_wf_params.get_dict()))
        future = submit(kkr_startpot_wc, kkr=kkrcode, voronoi=voronoicode, 
                        calc_parameters=params, wf_parameters=sub_wf_params, 
                        structure=structure)

        return ToContext(voronoi=future, last_calc=future)
    
    
    def check_voronoi(self):
        """
        check output of kkr_startpot_wc workflow that creates starting potential, shapefun etc.
        """
        
        self.report("INFO: checking voronoi output")
        voro_step_ok = False
        
        # check some output
        kkrstartpot_results = self.ctx.voronoi.out.results_vorostart_wc.get_dict()
        if kkrstartpot_results['successful']:
            voro_step_ok = True
        # initialize last_remote and last_params (gets updated in loop of KKR calculations)
        self.ctx.last_params = self.ctx.voronoi.out.last_params_voronoi
        self.ctx.last_remote = self.ctx.voronoi.out.last_voronoi_remote
        
        # store result in context
        self.ctx.voro_step_success = voro_step_ok
        
        # abort calculation if something failed in voro_start step
        if not voro_step_ok:
            error = 'ERROR: kkr_startpot_wc step failed!'
            self.ctx.errors.append(error)
            self.control_end_wc(error)
        
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
            
        # check if previous calculation reached convergence criterion
        if self.ctx.kkr_converged:
            if not self.ctx.kkr_higher_accuracy:
                do_kkr_step = do_kkr_step & True 
            else:
                stopreason = 'KKR converged'
                do_kkr_step = False
        else:
            do_kkr_step = do_kkr_step & True
            # following check only if calculaiton not converged:
            # check if maximal number of iterations has been reached
            if self.ctx.loop_count <= self.ctx.max_number_runs:
                do_kkr_step = do_kkr_step & True
            else:
                error = 'INFO: maximal number of KKR restarts reached. Exiting now!'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
        
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
        if decrease_mixing_fac or switch_agressive_mixing or switch_higher_accuracy or initial_settings:
            
            label = ''
            description = ''
            
            # step 1: extract info from last input parameters and check consistency
            params = self.ctx.last_params
            input_dict = params.get_dict()
            para_check = kkrparams()
                    
            # step 1.1: try to fill keywords
            for key, val in input_dict.iteritems():
                para_check.set_value(key, val, silent=True)
                
            # step 1.2: check if all mandatory keys are there
            missing_list = para_check.get_missing_keys(use_aiida=True)
            if missing_list != []:
                error = 'ERROR: calc_parameters misses keys: {}'.format(missing_list)
                self.ctx.errors.append(error)
                self.control_end_wc(error)
                
            # step 2: change parameter (contained in new_params dictionary)
            last_mixing_scheme = para_check.get_value('IMIX')
            if last_mixing_scheme is None:
                last_mixing_scheme = 0
                
            strmixfac = self.ctx.strmix
            brymixfac = self.ctx.brymix
            nsteps = self.ctx.nsteps
            
            # add number of scf steps
            new_params = {'NSTEPS':nsteps}
            
            # step 2.1 fill new_params dict with values to be updated
            if decrease_mixing_fac:
                if last_mixing_scheme == 0:
                    self.report('(strmixfax, mixreduce)= ({}, {})'.format(strmixfac, self.ctx.mixreduce))
                    self.report('type(strmixfax, mixreduce)= {} {}'.format(type(strmixfac), type(self.ctx.mixreduce)))
                    strmixfac = strmixfac * self.ctx.mixreduce
                    self.ctx.strmix = strmixfac
                    label += 'decreased_mix_fac_str'
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
                label += '_switched_to_agressive_mixing'
                description += ' switched to agressive mixing scheme (IMIX={})'.format(last_mixing_scheme)
            
            # add mixing scheme
            new_params['IMIX'] = last_mixing_scheme
                
            if switch_higher_accuracy:
                self.ctx.kkr_higher_accuracy = True
                convergence_settings = self.ctx.convergence_setting_fine
                label += '_use_higher_accuracy'
                description += ' using higher accuracy settings goven in convergence_setting_fine'
            else:
                convergence_settings = self.ctx.convergence_setting_coarse
            
            # add convegence settings
            new_params['NPOL'] = convergence_settings['npol']
            new_params['NPT1'] = convergence_settings['n1']
            new_params['NPT2'] = convergence_settings['n2']
            new_params['NPT3'] = convergence_settings['n3']
            new_params['TEMPR'] = convergence_settings['tempr']
            new_params['BZDIVIDE'] = convergence_settings['kmesh']
            
            # step 2.2 update values
            try:
                for key, val in new_params.iteritems():
                    para_check.set_value(key, val, silent=True)
            except:
                error = 'ERROR: parameter update unsuccessful: some key, value pair not valid!'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
                
            # step 3: 
            self.report("INFO: update parameters to: {}".format(para_check.get_set_values()))
            updatenode = ParameterData(dict=para_check.get_dict())
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
        
        
        label = 'KKR calulation step {}'.format(self.ctx.loop_count)
        description = 'KKR calculation of step {}'.format(self.ctx.loop_count)         
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
        kkr_run = submit(KkrProcess, **inputs)

        return ToContext(kkr=kkr_run, last_calc=kkr_run)
    

    def inspect_kkr(self):
        """
        check for convergence and store some of the results of the last calculation to context
        """
        self.report("INFO: inspecting kkr results step")
        
        self.ctx.calcs.append(self.ctx.last_calc)
        self.ctx.kkr_step_success = True
        
        # try yo extract remote folder
        try:
            self.ctx.last_remote = self.ctx.kkr.out.remote_folder
        except:
            self.ctx.last_remote = None
            self.ctx.kkr_step_success = False
        
        self.report("INFO: last_remote: {}".format(self.ctx.last_remote))
        
        # check calculation state
        calc_state = self.ctx.last_calc.get_state()
        if calc_state != calc_states.FINISHED:
            self.ctx.kkr_step_success = False
            
        self.report("INFO: kkr_step_success: {}".format(self.ctx.kkr_step_success))
        
        # extract convergence info about rms etc. (used to determine convergence behavior)
        try:
            self.report("INFO: trying to find output of last_calc: {}".format(self.ctx.last_calc))
            last_calc_output = self.ctx.last_calc.out.output_parameters.get_dict()
            found_last_calc_output = True
        except:
            found_last_calc_output = False
        self.report("INFO: found_last_calc_output: {}".format(found_last_calc_output))
            
        if self.ctx.kkr_step_success and found_last_calc_output:
            # 
            self.ctx.kkr_converged = last_calc_output['convergence_group']['calculation_converged']
            self.ctx.rms.append(last_calc_output['convergence_group']['rms'])
            rms_all_iter_last_calc = list(last_calc_output['convergence_group']['rms_all_iterations'])
            
            self.ctx.last_rms_all = rms_all_iter_last_calc
            if self.ctx.kkr_step_success and self.convergence_on_track():
                self.ctx.rms_all_steps += rms_all_iter_last_calc
        else:
            self.ctx.kkr_converged = False
            
        self.report("INFO: kkr_converged: {}".format(self.ctx.kkr_converged))
        self.report("INFO: rms: {}".format(self.ctx.rms))
        self.report("INFO: last_rms_all: {}".format(self.ctx.last_rms_all))
        self.report("INFO: rms_all_steps: {}".format(self.ctx.rms_all_steps))
            
        # TODO: extract something else (maybe total energy, charge neutrality, magnetisation)?
            
        self.report("INFO: done inspecting kkr results step")
        
        
    def convergence_on_track(self):
        """
        Check if convergence behavior of the last calculation is on track (i.e. going down)
        """
        on_track = True
        
        if self.ctx.kkr_step_success:
            last_rms = self.ctx.last_rms_all[-1]
            if last_rms >= self.ctx.last_rms_all[0]:
                on_track = False
                
        else:
            on_track = False
            
        return on_track
        


    def return_results(self):
        """
        return the results of the calculations
        This shoudl run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """
        
        self.report("INFO: entering return_results")
        
        try:
            last_calc_uuid = self.ctx.last_calc.uuid
            last_calc_pk = self.ctx.last_calc.pk
        except AttributeError:
            last_calc_uuid = None
            last_calc_pk = None
        try:
            last_params_uuid = self.ctx.last_params.uuid
            last_params_pk = self.ctx.last_params.pk
        except AttributeError:
            last_params_uuid = None
            last_params_pk = None
        try:
            last_remote_uuid = self.ctx.last_remote.uuid
            last_remote_pk = self.ctx.last_remote.pk
        except AttributeError:
            last_remote_uuid = None
            last_remote_pk = None
        try: # if something failed, we still might be able to retrieve something
            last_calc_out = self.ctx.kkr.out['output_parameters']
            last_calc_out_dict = last_calc_out.get_dict()
        except:# KeyError:
            last_calc_out = None #self.ctx.last_calc.out.CALL.out['output_parameters']
            last_calc_out_dict = {}#last_calc_out.get_dict()
        #except AttributeError:
        #    last_calc_out = None
        #    last_calc_out_dict = {}
            
        try:
            last_rms = self.ctx.rms[-1]
        except:
            last_rms = None


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
        outputnode_dict['last_calc_nodeinfo']   = {'uuid':last_calc_uuid,   'pk':last_calc_pk}
        outputnode_dict['errors'] = self.ctx.errors
        outputnode_dict['convergence_value'] = last_rms
        outputnode_dict['convergence_values_all_steps'] = self.ctx.rms_all_steps
        outputnode_dict['convergence_values_last_step'] = self.ctx.last_rms_all
        outputnode_dict['dos_check_ok'] = self.ctx.dos_ok
        outputnode_dict['convergence_reached'] = self.ctx.kkr_converged
        outputnode_dict['voronoi_step_success'] = self.ctx.voro_step_success
        outputnode_dict['kkr_step_success'] = self.ctx.kkr_step_success
        outputnode_dict['used_higher_accuracy'] = self.ctx.kkr_higher_accuracy
        
        # maybe also store some information about the formula
        #also lognotes, which then can be parsed from subworkflow too workflow, list of calculations involved (pks, and uuids),
        #This node should contain everything you wish to plot, here iteration versus, total energy and distance.

        if self.ctx.successful:
            self.report('STATUS: Done, the convergence criteria are reached.\n'
                        'INFO: The charge density of the KKR calculation pk= '
                        'converged after {} KKR runs and {} iterations to {} '
                        '"me/bohr^3" \n'.format(self.ctx.loop_count,
                                       last_calc_out_dict.get('number_of_iterations_total', None),
                                       last_calc_out_dict.get('charge_density', None)))

        else: # Termination ok, but not converged yet...
            if self.ctx.abort: # some error occured, donot use the output.
                self.report('STATUS/ERROR: I abort, see logs and '
                            'erros/warning/hints in output_kkr_scf_wc_para')
            else:
                self.report('STATUS/WARNING: Done, the maximum number of runs '
                            'was reached or something failed.\n INFO: The '
                            'charge density of the KKR calculation pk= '
                            'after {} KKR runs and {} iterations is {} "me/bohr^3"\n'
                            ''.format(self.ctx.loop_count,
                            last_calc_out_dict.get('number_of_iterations_total', None),
                            last_calc_out_dict.get('charge_density', None)))
                

        self.report("INFO: create results node: {}".format(outputnode_dict))
        outputnode_t = ParameterData(dict=outputnode_dict)
        outputnode_t.label = 'kkr_scf_wc_results'
        outputnode_t.description = 'Contains results of workflow (e.g. workflow version number, info about success of wf, lis tof warnings that occured during execution, ...)'
        
         # this is unsafe so far, because last_calc_out could not exist...
        if last_calc_out:
            outdict = create_scf_result_node(outpara=outputnode_t, last_calc_out=last_calc_out, last_RemoteData=self.ctx.last_remote)
        else:
            outdict = create_scf_result_node(outpara=outputnode_t)

        if last_calc_out:
            outdict['last_kkr_calc_output'] = last_calc_out

        # TODO return other nodes? depending what calculation was run

        #outdict['output_scf_wc_para'] = outputnode
        for link_name, node in outdict.iteritems():
            self.report("INFO: storing node {} {} with linkname {}".format(type(node), node, link_name))
            self.out(link_name, node)
            
        self.report("INFO: done with kkr_scf workflow!")


    def control_end_wc(self, errormsg):
        """
        Controled way to shutdown the workchain. will initalize the output nodes
        """
        self.report('ERROR: shutting workchain down in a controlled way.')
        self.ctx.successful = False
        self.ctx.abort = True
        self.report(errormsg) # because return_results still fails somewhen
        self.return_results()
        #self.abort_nowait(errormsg)
        self.abort(errormsg)
        
        
    def check_input_params(self, params, is_voronoi=False):
        """
        Checks input parameter consistency and aborts wf if check fails.
        """
        if params is None:
            error = 'ERROR: calc_parameters not given as input but are needed!'
            self.ctx.errors.append(error)
            self.control_end_wc(error)
        else:
            input_dict = params.get_dict()
            if is_voronoi:
                para_check = kkrparams(params_type='voronoi')
            else:
                para_check = kkrparams()
                
            # step 1 try to fill keywords
            try:
                for key, val in input_dict.iteritems():
                    para_check.set_value(key, val, silent=True)
            except:
                error = 'ERROR: calc_parameters given are not consistent! Hint: did you give an unknown keyword?'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
                
            # step 2: check if all mandatory keys are there
            missing_list = para_check.get_missing_keys(use_aiida=True)
            if missing_list != []:
                error = 'ERROR: calc_parameters given are not consistent! Hint: are all mandatory keys set?'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
   
             

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
    for key, val in kwargs.iteritems():
        if key == 'outpara': #  should always be there
            outpara = val
            has_last_outpara = True
        elif key == 'last_calc_out':
            has_last_calc_out_dict = True
            last_calc_out_dict = val
        elif key =='last_RemoteData':
            last_RemoteData_dict = val
            has_last_RemoteData = True
            
    outdict = {}
    if has_last_outpara:
        outputnode = outpara.copy()
        outputnode.label = 'output_kkr_scf_wc_para'
        outputnode.description = ('Contains self-consistency results and '
                                  'information of an kkr_scf_wc run.')
        outdict['output_kkr_scf_wc_ParameterResults'] = outputnode
    if has_last_calc_out_dict:
        outputnode2 = last_calc_out_dict.copy()
        outputnode2.label = 'output_kkr_scf_wc_lastResults'
        outputnode2.description = ('Contains the Results Parameter node from the output '
                                   'of the last calculation done in the workflow.')
        outdict['output_kkr_scf_wc_lastResults'] = outputnode2
    if has_last_RemoteData:
        outputnode3 = last_RemoteData_dict.copy()
        outputnode3.label = 'output_kkr_scf_wc_lastRemoteData'
        outputnode3.description = ('Contains a link to the latest remote data node '
                                   'where the output of the calculation can be accessed.')
        outdict['output_kkr_scf_wc_lastRemoteData'] = outputnode3
    
    # copy, because we rather produce the same node twice then have a circle in the database for now...
    #output_para = args[0]
    #return {'output_eos_wc_para'}
    return outdict

