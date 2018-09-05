# -*- coding: utf-8 -*-
"""
In this module you find the sub workflow for the kkrimp self consistency cycle  
and some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory, load_node
from aiida.work.workchain import WorkChain, ToContext, if_, while_
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode, update_params_wf, get_inputs_kkr
from aiida_kkr.calculations.kkr import KkrimpCalculation
from aiida.orm.calculation.job import JobCalculation
from aiida.common.datastructures import calc_states
from aiida.orm import WorkCalculation
from aiida.common.exceptions import InputValidationError


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = u"Fabian Bertoldo"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
FolderData = DataFactory('folder')
KkrimpProcess = KkrimpCalculation.process()



class kkr_imp_sub_wc(WorkChain):
    """
    Workchain of a kkrimp sub self consistency calculation starting from the 
    host-impurity potential of the system.

    :param options_parameters: (ParameterData), Workchain specifications
    :param wf_parameters: (ParameterData), specifications for the calculation
    :param host_imp_pot: (RemoteData), mandatory; host-impurity potential
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
                   'fac_cls_increase' : 1.3,                  # factor by which the screening cluster is increased each iteration (up to num_rerun times)
                   'r_cls' : 1.3, # alat                      # default cluster radius, is increased iteratively
                   'natom_in_cls_min' : 79                    # minimum number of atoms in screening cluster
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
        spec.input("host_imp_pot", valid_type=RemoteData, required=True)
        
        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            cls.validate_input,
            while_(cls.condition)(
                cls.update_kkrimp_params,
                cls.run_kkrimp,
                cls.inspect_kkrimp
            ),
            cls.return_results
        )
        
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
        spec.exit code(124, 'ERROR_NO_CALC_PARAMS',
            message="ERROR: No calculation parameters provided")
        
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
        self.ctx.delta_e = wf_dict.get('delta_e_min', self._wf_default['delta_e_min'])

        
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

        inputs = self.inputs
        inputs_ok = True

        if not 'host_imp_pot' in inputs:
            inputs_ok = False
            return self.exit_codes.ERROR_NO_HOST_IMP_POT

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkrimp', use_exceptions=True)
            except ValueError:
                inputs_ok = False
                return self.exit_codes.ERROR_INVALID_INPUT_KKRIMP

        # set params and remote folder to input 
        self.ctx.last_remote = inputs.remote_data
        try: # first try parent of remote data output of a previous calc.
            parent_params = get_parent_paranode(inputs.host_imp_pot)
        except AttributeError:
            inputs_ok = False
            return self.exit_codes.ERROR_INVALID_HOST_IMP_POT
        if  'calc_parameters' in inputs:
            self.ctx.last_params = inputs.calc_parameters
        else:
            inputs_ok = False
            return self.exit_codes.ERROR_NO_CALC_PARAMS
            
        self.report('Validated input successfully: {}'.format(inputs_ok))

        return inputs_ok
