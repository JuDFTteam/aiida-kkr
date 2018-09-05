# -*- coding: utf-8 -*-
"""
In this module you find the sub workflow for the kkrimp self consistency cycle  
and some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory, load_node
from aiida.work.workchain import WorkChain, ToContext, if_
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
        #spec.outline()
        
        # exit codes 
        #spec.exit_code()
        
        # Define the outputs of the workflow
        #spec.output()