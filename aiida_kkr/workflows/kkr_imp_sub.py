# -*- coding: utf-8 -*-
"""
In this module you find the sub workflow for the kkrimp self consistency cycle  
and some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory, load_node
from aiida.work.workchain import WorkChain, ToContext, if_, while_
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode, update_params_wf, get_inputs_kkrimp, get_inputs_kkr
from aiida_kkr.calculations.kkrimp import KkrimpCalculation
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida.orm.calculation.job import JobCalculation
from aiida.common.datastructures import calc_states
from aiida.orm import WorkCalculation
from aiida.common.exceptions import InputValidationError
from numpy import array, where, ones


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = u"Fabian Bertoldo"

#TODO: work on return results function, so far output like in kkr_scf
#TODO: check get_inputs_kkr and modify for impurity calculation
#TODO: edit inspect_kkrimp function
#TODO: get rid of create_scf_result node and create output nodes differently
#TOTO: check if calculation parameters from previous calculation have to be 
#      loaded (in validate input, compare to kkr workflow)


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
SinglefileData = DataFactory('singlefile')
FolderData = DataFactory('folder')
KkrimpProcess = KkrimpCalculation.process()



class kkr_imp_sub_wc(WorkChain):
    """
    Workchain of a kkrimp sub self consistency calculation starting from the 
    host-impurity potential of the system.

    :param options_parameters: (ParameterData), Workchain specifications
    :param wf_parameters: (ParameterData), specifications for the calculation
    :param host_imp_startpot: (RemoteData), mandatory; host-impurity potential
    :param kkrimp: (Code), mandatory; KKRimp code converging the host-imp-potential

    :return result_kkr_imp_sub_wc: (ParameterData), Information of workflow results
                                   like success, last result node, list with 
                                   convergence behavior
    """
    
    _workflowversion = __version__
    _wf_label = 'kkr_imp_sub_wc'
    _wf_description = 'Workflow for a KKRimp self consistency calculation for converging a host-impurity potential'
       

    _options_default = {'queue_name' : '',                        # Queue name to submit jobs too
                        'resources': {"num_machines": 1},         # resources to allowcate for the job
                        'walltime_sec' : 60*60,                   # walltime after which the job gets killed (gets parsed to KKR)}
                        'custom_scheduler_commands' : '',         # some additional scheduler commands 
                        'use_mpi' : False}   
                        # execute KKR with mpi or without
    _wf_default = {'kkr_runmax': 5,                           # Maximum number of kkr jobs/starts (defauld iterations per start)
                   'convergence_criterion' : 3**-2,          # Stop if charge denisty is converged below this value
                   'mixreduce': 0.5,                          # reduce mixing factor by this factor if calculaito fails due to too large mixing
                   'threshold_aggressive_mixing': 3*10**-2,   # threshold after which agressive mixing is used
                   'strmix': 0.03,                            # mixing factor of simple mixing
                   'brymix': 0.05,                            # mixing factor of aggressive mixing
                   'nsteps': 100,                              # number of iterations done per KKR calculation
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
                   'fac_cls_increase' : 1.3,                  # factor by which the screening cluster is increased each iteration (up to num_rerun times)
                   'r_cls' : 1.3, # alat                      # default cluster radius, is increased iteratively
                   #'natom_in_cls_min' : 79                    # minimum number of atoms in screening cluster
                   }     

                        
    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create 
        set of wf_parameters.
        
        returns _wf_defaults
        """
    
        print('Version of workflow: {}'.format(self._workflowversion))
        return self._options_default, self._wf_default
        

        
    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """
    
        super(kkr_imp_sub_wc, cls).define(spec)
        
        # Define the inputs of the workflow
        spec.input("kkrimp", valid_type=Code, required=True)     
        spec.input("options_parameters", valid_type=ParameterData, required=False,
                       default=ParameterData(dict=cls._options_default))
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                       default=ParameterData(dict=cls._wf_default))
        spec.input("host_imp_startpot", valid_type=SinglefileData, required=True)
        spec.input("GF_remote_data", valid_type=RemoteData, required=True)
        
        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            cls.validate_input,
            while_(cls.condition)(
                cls.update_kkrimp_params,
                cls.run_kkrimp,
                cls.inspect_kkrimp),
            cls.return_results)
            
        # exit codes 
        spec.exit_code(121, 'ERROR_NO_HOST_IMP_POT', 
            message="ERROR: No host-impurity potential found in the inputs")
        spec.exit_code(122, 'ERROR_INVALID_INPUT_KKRIMP',
            message="ERROR: The code you provided for KKRimp does not "
                    "use the plugin kkr.kkrimp")
        spec.exit_code(123, 'ERROR_INVALID_HOST_IMP_POT',
            message="ERROR: Unable to extract parent paremeter node of "
                    "input remote folder")
        # probably not necessary
        spec.exit_code(124, 'ERROR_NO_CALC_PARAMS',
            message="ERROR: No calculation parameters provided")
        
        spec.exit_code(125, 'ERROR_SUB_FAILURE',
            message="ERROR: Last KKRcalc in SUBMISSIONFAILED state!\nstopping now")
        spec.exit_code(126, 'ERROR_MAX_STEPS_REACHED',
            message="ERROR: Maximal number of KKR restarts reached. Exiting now!")
        spec.exit_code(127, 'ERROR_SETTING_LAST_REMOTE',
            message="ERROR: Last_remote could not be set to a previous succesful calculation")
        spec.exit_code(127, 'ERROR_MISSING_PARAMS',
            message="ERROR: There are still missing calculation parameters")
        spec.exit_code(128, 'ERROR_PARAMETER_UPDATE',
            message="ERROR: Parameters could not be updated")
        spec.exit_code(129, 'ERROR_LAST_CALC_NOT_FINISHED',
            message="ERROR: Last calculation is not in finished state")
        
        # Define the outputs of the workflow
        spec.output('calculation_info', valid_type=ParameterData)
        

        
    def start(self):
        """
        init context and some parameters
        """
        self.report('INFO: started KKR impurity convergence workflow version {}'
                    ''.format(self._workflowversion))

        ####### init #######

        # internal para /control para
        self.ctx.loop_count = 0
        self.ctx.last_mixing_scheme = 0
        self.ctx.calcs = []
        self.ctx.abort = False
        # flags used internally to check whether the individual steps were successful
        self.ctx.kkr_converged = False
        self.ctx.kkr_step_success = False
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
        options_dict = self.inputs.options_parameters.get_dict()

        if options_dict == {}:
            options_dict = self._options_default
            self.report('INFO: using default options')
        
        if wf_dict == {}:
            wf_dict = self._wf_default
            self.report('INFO: using default wf parameter')

        # set option parameters from input, or defaults
        self.ctx.use_mpi = options_dict.get('use_mpi', self._options_default['use_mpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.walltime_sec = options_dict.get('walltime_sec', self._options_default['walltime_sec'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])
        
        # set workflow parameters from input, or defaults
        self.ctx.max_number_runs = wf_dict.get('kkr_runmax', self._wf_default['kkr_runmax'])
        self.ctx.description_wf = self.inputs.get('description', 'Workflow for '
                                                  'a KKR impurity calculation'
                                                  'starting from a host-impurity'
                                                  'potential')
        self.ctx.label_wf = self.inputs.get('label', 'kkr_imp_sub_wc')
        self.ctx.strmix = wf_dict.get('strmix', self._wf_default['strmix'])
        self.ctx.brymix = wf_dict.get('brymix', self._wf_default['brymix'])
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
        #self.ctx.delta_e = wf_dict.get('delta_e_min', self._wf_default['delta_e_min'])

        
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
                                self.ctx.mixreduce, self.ctx.mag_init,
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
        validate input and catch possible errors from the input
        """

        self.report('INFO: Validate input ...')
        
        inputs = self.inputs
        inputs_ok = True

        if not 'host_imp_startpot' in inputs:
            inputs_ok = False
            return self.exit_codes.ERROR_NO_HOST_IMP_POT

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkrimp', use_exceptions=True)
            except ValueError:
                inputs_ok = False
                return self.exit_codes.ERROR_INVALID_INPUT_KKRIMP        
                
        # set params and remote folder to input 
        self.ctx.last_remote = inputs.GF_remote_data
        
#        try: # first try parent of remote data output of a previous calc.
#            #parent_params = get_parent_paranode(inputs.host_imp_startpot)
#        except AttributeError:
#            inputs_ok = False
#            self.report('INFO: {}'.format(self.exit_codes.ERROR_INVALID_HOST_IMP_POT))
#            return self.exit_codes.ERROR_INVALID_HOST_IMP_POT
           
        if  'wf_parameters' in inputs:
            self.ctx.last_params = inputs.wf_parameters
        else:
            inputs_ok = False
            self.report('ERROR: {}'.format(self.exit_codes.ERROR_NO_CALC_PARAMS ))            
            return self.exit_codes.ERROR_NO_CALC_PARAMS 
            
        self.report('INFO: Validated input successfully: {}'.format(inputs_ok))   
        
        
        
    def condition(self):
        """
        check convergence condition
        """
        
        self.report("INFO: checking condition for kkrimp calculation")
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
                do_kkr_step = False
        else:
            do_kkr_step = do_kkr_step & True
            
        # check if previous calculation is in SUBMISSIONFAILED state
        if self.ctx.loop_count>1 and self.ctx.last_calc.get_state() == calc_states.SUBMISSIONFAILED:
            return self.exit_codes.ERROR_SUB_FAILURE

        # next check only needed if another iteration should be done after validating convergence etc. (previous checks)
        if do_kkr_step:
            # check if maximal number of iterations has been reached
            if self.ctx.loop_count <= self.ctx.max_number_runs:
                do_kkr_step = do_kkr_step & True
            else:
                return self.exit_codes.ERROR_MAX_STEPS_REACHED

        self.report("INFO: Done checking condition for kkr step (result={})".format(do_kkr_step))

        if not do_kkr_step:
            self.report("INFO: Stopreason={}".format(stopreason))

        return do_kkr_step
        
        
        
    def update_kkrimp_params(self):
        """
        update set of KKR parameters (check for reduced mixing, change of 
        mixing strategy, change of accuracy setting)
        """
        
        self.report("INFO: updating kkrimp param step")
        decrease_mixing_fac = False
        switch_agressive_mixing = False
        switch_higher_accuracy= False
        initial_settings = False

        # only do something other than simple mixing after first kkr run
        if self.ctx.loop_count != 1:
            # first determine if previous step was successful 
            # (otherwise try to find some rms value and decrease mixing to try again)
            if not self.ctx.kkr_step_success:
                decrease_mixing_fac = True
                self.report("INFO: last KKR calculation failed. Trying decreasing mixfac")

            convergence_on_track = self.convergence_on_track()

            # check if calculation was on its way to converge
            if not convergence_on_track:
                decrease_mixing_fac = True
                self.report("INFO: Last KKR did not converge. Trying decreasing mixfac")
                # reset last_remote to last successful calculation
                for icalc in range(len(self.ctx.calcs))[::-1]:
                    self.report("INFO: last calc success? {} {}".format(icalc, self.ctx.KKR_steps_stats['success'][icalc]))
                    if self.ctx.KKR_steps_stats['success'][icalc]:
                        self.ctx.last_remote = self.ctx.calcs[icalc].out.remote_folder
                        break # exit loop if last_remote was found successfully
                    else:
                        self.ctx.last_remote = None
                # if no previous calculation was succesful take voronoi output 
                # or remote data from input (depending on the inputs)
                self.report("INFO: Last_remote is None? {} {}".format(self.ctx.last_remote is None, 'structure' in self.inputs))
                if self.ctx.last_remote is None:
                    if 'structure' in self.inputs:
                        self.ctx.voronoi.out.last_voronoi_remote
                    else:
                        self.ctx.last_remote = self.inputs.remote_data
                # check if last_remote has finally been set and abort if this is not the case
                self.report("INFO: last_remote is still None? {}".format(self.ctx.last_remote is None))
                if self.ctx.last_remote is None:
                    error = 'ERROR: last_remote could not be set to a previous succesful calculation'
                    self.ctx.errors.append(error)
                    return self.exit_codes.ERROR_SETTING_LAST_REMOTE

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
            para_check = kkrparams(params_type='kkrimp')
            para_check.get_all_mandatory()
            self.report('INFO: use default kkrimp keywords: {}'.format(para_check.get_dict()))

            # step 1.1: try to fill keywords
            #for key, val in input_dict.iteritems():
            #    para_check.set_value(key, val, silent=True)

            # init new_params dict where updated params are collected
            new_params = {}
            
            # step 1.2: check if all mandatory keys are there and add defaults if missing
            missing_list = para_check.get_missing_keys(use_aiida=True)
            if missing_list != []:
                kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
                kkrdefaults_updated = []
                for key_default, val_default in kkrdefaults.items():
                    if key_default in missing_list:
                        new_params[key_default] = kkrdefaults.get(key_default)
                        kkrdefaults_updated.append(key_default)
                if len(kkrdefaults_updated)>0:
                    error = 'ERROR: Calc_parameters misses keys: {}'.format(missing_list)
                    self.ctx.errors.append(error)
                    print('ERROR: {}'.format(self.exit_codes.ERROR_MISSING_PARAMS))
                    return self.exit_codes.ERROR_MISSING_PARAMS
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
            new_params['SCFSTEPS'] = nsteps

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
            new_params['MIXFAC'] = strmixfac
            #new_params['BRYMIX'] = brymixfac

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

            # add convergence settings
            if self.ctx.loop_count == 1 or self.ctx.last_mixing_scheme == 0:
                new_params['QBOUND'] = self.ctx.threshold_aggressive_mixing
            else:
                new_params['QBOUND'] = self.ctx.convergence_criterion
            #new_params['NPOL'] = convergence_settings['npol']
            #new_params['NPT1'] = convergence_settings['n1']
            #new_params['NPT2'] = convergence_settings['n2']
            #new_params['NPT3'] = convergence_settings['n3']
            #new_params['TEMPR'] = convergence_settings['tempr']
            #new_params['BZDIVIDE'] = convergence_settings['kmesh']
            
            # initial magnetization
            if initial_settings and self.ctx.mag_init:
                if self.ctx.hfield <= 0:
                    self.report('\nWARNING: magnetization initialization chosen but hfield is zero. Automatically change back to default value (hfield={})\n'.format(self._wf_default['hfield']))
                    self.ctx.hfield = self._wf_default['hfield']
                xinipol = self.ctx.xinit
                if xinipol is None:
                    if 'structure' in self.inputs:
                        struc = self.inputs.structure
                    else:
                        struc, voro_parent = KkrCalculation.find_parent_structure(self.ctx.last_remote)
                    natom = len(struc.sites)
                    xinipol = ones(natom)
                #new_params['LINIPOL'] = True
                new_params['HFIELD'] = self.ctx.hfield
                new_params['XINIPOL'] = xinipol
            elif self.ctx.mag_init and self.ctx.mag_init_step_success: # turn off initialization after first (successful) iteration
                #new_params['LINIPOL'] = False
                new_params['HFIELD'] = [0.0, 0]
            elif not self.ctx.mag_init:
                self.report("INFO: mag_init is False. Overwrite 'HFIELD' to '0.0' and 'LINIPOL' to 'False'.")
                # reset mag init to avoid resinitializing
                #new_params['LINIPOL'] = False
                new_params['HFIELD'] = [0.0, 0]
                
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
                for key, val in new_params.iteritems():
                    para_check.set_value(key, val, silent=True)
            except:
                error = 'ERROR: parameter update unsuccessful: some key, value pair not valid!'
                self.ctx.errors.append(error)
                return self.exit_codes.ERROR_PARAMETER_UPDATE

            # step 3:
            self.report("INFO: update parameters to: {}".format(para_check.get_set_values()))
            
            #test
            self.ctx.last_params = ParameterData(dict={})
            
            updatenode = ParameterData(dict=para_check.get_dict())
            updatenode.label = label
            updatenode.description = description
            self.report('INFO: updatenode: {}'.format(updatenode.get_attrs()))
            self.report('INFO: old node: {}'.format(self.ctx.last_params.get_attrs()))
            
            paranode_new = updatenode#update_params_wf(self.ctx.last_params, updatenode)
            self.ctx.last_params = paranode_new
        else:
            self.report("INFO: reuse old settings")

        self.report("INFO: done updating kkr param step")


        
    def run_kkrimp(self):
        """
        submit a KKR impurity calculation
        """
        self.report("INFO: setting up kkrimp calculation step {}".format(self.ctx.loop_count))


        label = 'KKRimp calculation step {} (IMIX={})'.format(self.ctx.loop_count, self.ctx.last_mixing_scheme)
        description = 'KKRimp calculation of step {}, using mixing scheme {}'.format(self.ctx.loop_count, self.ctx.last_mixing_scheme)
        code = self.inputs.kkrimp
        remote = self.ctx.last_remote
        params = self.ctx.last_params
        print(type(params))
        host_GF = self.inputs.GF_remote_data
        imp_pot = self.inputs.host_imp_startpot
        
        options = {"max_wallclock_seconds": self.ctx.walltime_sec,
                   "resources": self.ctx.resources,
                   "queue_name" : self.ctx.queue}
        if self.ctx.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.ctx.custom_scheduler_commands
        inputs = get_inputs_kkrimp(code, options, label, description, params, not self.ctx.use_mpi, host_GF=host_GF, imp_pot=imp_pot)
        print(inputs)
        
        # run the KKR calculation
        self.report('INFO: doing calculation')
        kkrimp_run = self.submit(KkrimpProcess, **inputs)

        return ToContext(kkr=kkrimp_run, last_calc=kkrimp_run)
        
        
        
    def inspect_kkrimp(self):
        """
        check for convergence and store some of the results of the last calculation to context
        """
        self.report("INFO: inspecting kkrimp results step")

        self.ctx.calcs.append(self.ctx.last_calc)
        self.ctx.kkrimp_step_success = True

        # check calculation state
        calc_state = self.ctx.last_calc.get_state()
        if calc_state != calc_states.FINISHED:
            self.ctx.kkrimp_step_success = False
            return self.exit_codes.ERROR_LAST_CALC_NOT_FINISHED

        self.report("INFO: kkrimp_step_success: {}".format(self.ctx.kkrimp_step_success))

        # extract convergence info about rms etc. (used to determine convergence behavior)
        try:
            self.report("INFO: trying to find output of last_calc: {}".format(self.ctx.last_calc))
            last_calc_output = self.ctx.last_calc.out.output_parameters.get_dict()
            found_last_calc_output = True
        except:
            found_last_calc_output = False
        self.report("INFO: found_last_calc_output: {}".format(found_last_calc_output))
        
        # try yo extract remote folder
        try:
            self.ctx.last_remote = self.ctx.kkr.out.remote_folder
        except:
            self.ctx.last_remote = None
            self.ctx.kkrimp_step_success = False

        self.report("INFO: last_remote: {}".format(self.ctx.last_remote))

        if self.ctx.kkrimp_step_success and found_last_calc_output:
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
            if self.ctx.kkrimp_step_success and self.convergence_on_track():
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
        if self.ctx.mag_init and self.ctx.kkrimp_step_success:
            self.ctx.mag_init_step_success = True

        # TODO: extract something else (maybe total energy, charge neutrality, magnetisation)?

        # store some statistics used to print table in the end of the report
        self.ctx.KKR_steps_stats['success'].append(self.ctx.kkr_step_success)
        try:
            isteps = self.ctx.last_calc.out.output_parameters.get_dict()['convergence_group']['number_of_iterations']
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

        self.report("INFO: done inspecting kkrimp results step")

        
        
    def convergence_on_track(self):
        """
        Check if convergence behavior of the last calculation is on track (i.e. going down)
        """
        on_track = True
        threshold = 5. # used to check condition if at least one of charnge_neutrality, rms-error goes down fast enough

        # first check if previous calculation was stopped due to reaching the QBOUND limit
        try:
            calc_reached_qbound = self.ctx.last_calc.out.output_parameters.get_dict()['convergence_group']['calculation_converged']
        except AttributeError: # captures error when last_calc dies not have an output node
            calc_reached_qbound = False
        except KeyError: # captures
            calc_reached_qbound = False
           
        if self.ctx.kkrimp_step_success and not calc_reached_qbound:
            first_rms = self.ctx.last_rms_all[0]
            first_neutr = abs(self.ctx.last_neutr_all[0])
            last_rms = self.ctx.last_rms_all[-1]
            last_neutr = abs(self.ctx.last_neutr_all[-1])
            # use this trick to avoid division by zero
            if last_neutr == 0:
                last_neutr = 10**-16
            if last_rms == 0:
                last_rms = 10**-16
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
        Return the results of the calculations
        This should run through and produce output nodes even if everything failed,
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
            results_vorostart = self.ctx.voronoi.out.results_vorostart_wc
        except:
            results_vorostart = None
        try:
            starting_dosdata_interpol = self.ctx.voronoi.out.last_doscal_dosdata_interpol
        except:
            starting_dosdata_interpol = None

        try:
            final_dosdata_interpol = self.ctx.doscal.out.dos_data_interpol
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
        outputnode_t = ParameterData(dict=outputnode_dict)
        outputnode_t.label = 'kkr_scf_wc_results'
        outputnode_t.description = 'Contains results of workflow (e.g. workflow version number, info about success of wf, lis tof warnings that occured during execution, ...)'

        self.out('calculation_info', outputnode_t)

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
