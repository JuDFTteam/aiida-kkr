# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for writing out the kkr_flexfiles and
some helper methods to do so with AiiDA
"""
        

from aiida.orm import Code, DataFactory, load_node
from aiida.work.workchain import WorkChain, ToContext
from aiida.work.run import submit
from aiida.work import workfunction as wf
from aiida.work.process_registry import ProcessRegistry
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode, update_params_wf, get_inputs_kkr
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.orm.calculation.job import JobCalculation
from aiida.common.datastructures import calc_states
from aiida.orm import WorkCalculation
from aiida.common.exceptions import InputValidationError



__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.4"
__contributors__ = u"Fabian Bertoldo"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()


class kkr_flex_wc(WorkChain):
    """
    Workchain of a kkr_flex calculation with KKR starting from the RemoteData node 
    of a previous calculation (either Voronoi or KKR).

    :param wf_parameters: (ParameterData), Workchain specifications
    :param remote_data: (RemoteData), mandatory; from a converged KKR calculation
    :param kkr: (Code), mandatory; KKR code running the flexfile writeout
    :param imp_info: imp_info node specifying information of the impurities in the system

    :return result_kkr_flex_wc: (ParameterData), Information of workflow results
                                like success, last result node, list with convergence behavior
    """

    _workflowversion = __version__
    _wf_label = 'kkr_flex_wc'
    _wf_description = 'Workflow for a KKR flex calculation starting from RemoteData node of previous converged KKR calculation'
       

    _options_default = {'queue_name' : '',                        # Queue name to submit jobs too
                        'resources': {"num_machines": 1},         # resources to allowcate for the job
                        'walltime_sec' : 60*60,                   # walltime after which the job gets killed (gets parsed to KKR)}
                        'custom_scheduler_commands' : '',         # some additional scheduler commands 
                        'use_mpi' : False}                        # execute KKR with mpi or without

    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """
    
        print('Version of workflow: {}'.format(self._workflowversion))
        return self._options_default

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """
    
        # Take input of the workflow or use defaults defined above
        super(kkr_flex_wc, cls).define(spec)
        
        spec.input('code', valid_type=Code)     
        spec.input("options_parameters", valid_type=ParameterData, required=False,
                       default=ParameterData(dict=cls._options_default)
        spec.input("remote_data", valid_type=RemoteData, required=True)
        spec.input("kkr", valid_type=Code, required=True)
        spec.input("imp_info", valid_type=ParameterData, required=True)
    
        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            cls.validate_input,
            cls.set_params_flex,
            cls.get_flex, # calculate host GF and kkr-flexfiles
            cls.return_results
            )

        # ToDo: add implementation of exit codes (compare to aiida_quantumespresso)
        spec.exit_code(101, 'ERROR_INVALID_INPUT_IMP_INFO', 
            message="the 'imp_info' input ParameterData node could not be used")
        spec.exit_code(102, 'ERROR_INVALID_INPUT_KKR',
            message="the code you provided for kkr does not use the plugin kkr.kkr")


        # ToDo: add output of the workflow
        # Define the output of the workflow
        spec.output('remote_folder', valid_type=RemoteData)




    def start(self):
        """
        init context and some parameters
        """
    
        self.report('INFO: started KKRflex workflow version {}\n'
                    'INFO: Workchain node identifiers: {}'
                    ''.format(self._workflowversion, ProcessRegistry().current_calc_node))
    
        ####### init #######
        # internal para / control para
        self.ctx.abort = False
    
        # input both wf and options parameters
        options_dict = self.inputs.options_parameters.get_dict()
    
        if options_dict == {}:
            options_dict = self._options_default
            self.report('INFO: using default options parameters')

        # set values, or defaults
        # ToDo: arrange option assignment differently (look at scf.py from aiida-fleur)
        self.ctx.use_mpi = options_dict.get('use_mpi', self._options_default['use_mpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.walltime_sec = options_dict.get('walltime_sec', self._options_default['walltime_sec'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])

    
        # ToDo: exchange underscore in front of description and label?
        self.ctx.description_wf = self.inputs.get('_description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('_label', self._wf_label)
    

        self.report('INFO: use the following parameter:\n'
                    'use_mpi: {}\n'
                    'Resources: {}\n'
                    'Walltime (s): {}\n'
                    'queue name: {}\n'
                    'scheduler command: {}\n'
                    'description: {}\n'
                    'label: {}\n'.format(self.ctx.use_mpi, self.ctx.resources, self.ctx.walltime_sec, 
                                         self.ctx.queue, self.ctx.custom_scheduler_commands, 
                                         self.ctx.description_wf, self.ctx.label_wf))
    
        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''
    


    def validate_input(self):
        """
        Validate input
        """
        # ToDo: not modified enough yet! Taken from DOS workflow   

        inputs = self.inputs

        if not 'imp_info' in inputs:
            input_ok = False
            return self.exit_codes.ERROR_INVALID_INPUT_IMP_INFO
    
        if 'remote_data' in inputs:
            input_ok = True
        else:
            error = 'ERROR: No remote_data was provided as Input'
            self.ctx.errors.append(error)
            self.control_end_wc(error)
            input_ok = False
    
        # extract correct remote folder of last calculation if input remote_folder node
        # is not from KKRCalculation but kkr_scf_wc workflow
        input_remote = self.inputs.remote_data
        # check if input_remote has single KKRCalculation parent
        parents = input_remote.get_inputs(node_type=JobCalculation)
        nparents = len(parents)
        if nparents!=1:
            # extract parent workflow and get uuid of last calc from output node
            parent_workflow = input_remote.inp.last_RemoteData
            if not isinstance(parent_workflow, WorkCalculation):
                raise InputValidationError("Input remote_data node neither output of a KKR calculation nor of kkr_scf_wc workflow")
                parent_workflow_out = parent_workflow.out.output_kkr_scf_wc_ParameterResults
                uuid_last_calc = parent_workflow_out.get_dict().get('last_calc_nodeinfo').get('uuid')
                last_calc = load_node(uuid_last_calc)
                if not isinstance(last_calc, KkrCalculation):
                    raise InputValidationError("Extracted last_calc node not of type KkrCalculation: check remote_data input node")
                # overwrite remote_data node with extracted remote folder
                output_remote = last_calc.out.remote_folder
                self.inputs.remote_data = output_remote
            
            if 'kkr' in inputs:
                try:
                    test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
                except ValueError:
                    error = ("The code you provided for kkr does not "
                             "use the plugin kkr.kkr")
                    self.ctx.errors.append(error)
                    self.control_end_wc(error)
                    input_ok = False
            
            # set self.ctx.input_params_KKR
            self.ctx.input_params_KKR = get_parent_paranode(self.inputs.remote_data)
            
            return input_ok

            
            
    def set_params_flex(self):
        """
        Take input parameter node and change to input from wf_parameter and options 
        """
    
        params = self.ctx.input_params_KKR
        input_dict = params.get_dict()
        para_check = kkrparams()
    
        # step 1: try to fill keywords
        try:
            for key, val in input_dict.iteritems():
                para_check.set_value(key, val, silent=True)
        except:
            error = 'ERROR: calc_parameters given are not consistent! Hint: did you give an unknown keyword?'
            self.ctx.errors.append(error)
            self.control_end_wc(error)
    
        # step 2: check if all mandatory keys are there
        label = ''
        descr = ''
        missing_list = para_check.get_missing_keys(use_aiida=True)
        if missing_list != []:
            kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
            kkrdefaults_updated = []
            for key_default, val_default in kkrdefaults.items():
                if key_default in missing_list:
                    para_check.set_value(key_default, kkrdefaults.get(key_default), silent=True)
                    kkrdefaults_updated.append(key_default)
                    missing_list.remove(key_default)
            if len(missing_list)>0:
                error = 'ERROR: calc_parameters misses keys: {}'.format(missing_list)
                self.ctx.errors.append(error)
                self.control_end_wc(error)
            else:
                self.report('updated KKR parameter node with default values: {}'.format(kkrdefaults_updated))
                label = 'add_defaults_'
                descr = 'added missing default keys, '
    
        # add the RUNOPT = KKRFLEX to the params
        para_check.set_value(['RUNOPT'], ['KKRFLEX'], silent=True)

        #construct the final param node containing all of the params   
        updatenode = ParameterData(dict=para_check.get_dict())
        updatenode.label = label+'KKRparam_flex'
        updatenode.description = descr+'KKR parameter node extracted from parent parameters and wf_parameter and options input node.'
        paranode_flex = update_params_wf(self.ctx.input_params_KKR, updatenode)
        self.ctx.flex_kkrparams = paranode_flex

       

    def get_flex(self):
        """
        Submit a KKRFLEX calculation
        """

        label = 'KKRFLEX calc.'
        # flexdict = self.ctx.flex_params_dict
        description = 'KKRFLEX calculation to write out host GF using RUNOPT={}'.format()
        code = self.inputs.kkr
        remote = self.inputs.remote_data
        params = self.ctx.flex_kkrparams
        options = {"max_wallcloxk_seconds": self.ctx.walltime_sec,
                   "resources": self.ctx.resources,
                   "queue_name": self.ctx.queue}
        if self.ctx.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.ctx.custom_scheduler_commands
        inputs = get_inputs_kkr(code, remote, options, label, description, parameters=params, serial=(not self.ctx.use_mpi))

        # run the KKRFLEX calculation
        self.report('INFO: doing calculation')
        self.use_impurity_info(imp_info)
        flexrun = self.submit(KkrProcess, **inputs)

        return ToContext(flexrun=flexrun)


    def return_results(self):
        """
        Return the results of the KKRFLEX calculation.
        This should run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """ 
        # has to be worked on, ask Philipp

        # capture error of unsuccessful DOS run
        calc_state = self.ctx.flexrun.get_state()
        if calc_state != calc_states.FINISHED:
            self.ctx.successful = False
            error = ('ERROR: KKRFLEX calculation failed somehow it is '
                     'in state{}'.format(calc_state))
            self.ctx.errors.append(error)

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['use_mpi'] = self.ctx.use_mpi
        outputnode_dict['resources'] = self.ctx.resources
        outputnode_dict['walltime_sec'] = self.ctx.walltime_sec
        outputnode_dict['queue'] = self.ctx.queue
        outputnode_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        
        # ToDo: has to be modified
        # outputnode_dict['flex_params'] = self.ctx.flex_params_dict
        
        try:
            outputnode_dict['nspin'] = self.ctx.flex.res.nspin # TODO: how does that work?
        except:
            error = "ERROR: nspin not extracted"
            self.report(error)
            self.ctx.successful = False
            self.ctx.errors.append(error)
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['list_of_errors'] = self.ctx.errors
            
        outputnode = ParameterData(dict=outputnode_dict)
        outputnode.label = 'kkr_scf_wc_results'
        outputnode.description = ''
        outputnode.store()
        
        self.report("INFO: create GF writeout results nodes: outputnode={}".format(outputnode))
        try:
            self.report("INFO: create GF writeout results nodes. KKRFLEX calc retrieved node={}".format(self.ctx.flexrun.out.retrieved))
            has_flexrun = True
        except AttributeError as e:
            self.report("ERROR: no KKRFLEX calc retrieved node found")
            self.report("Caught AttributeError {}".format(e))
            has_flexrun = False

        # interpol dos file and store to XyData nodes
        if has_flexrun:
            # ask Philipp
            outdict = create_flex_result_node(outputnode, self.ctx.flexrun.out.retrieved)
        else:
            # ask Philipp how this function can be defined and where
            outdict = create_flex_result_node_minimal(outputnode)
    
        
        for link_name, node in outdict.iteritems():
            #self.report("INFO: storing node '{}' with link name '{}'".format(node, link_name))
            #self.report("INFO: node type: {}".format(type(node)))
            self.out(link_name, node)
            
        self.report("INFO: done with KKRFLEX GF writeout  workflow!\n")