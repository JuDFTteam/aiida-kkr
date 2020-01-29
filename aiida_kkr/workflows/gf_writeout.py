# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for writing out the kkr_flexfiles and
some helper methods to do so with AiiDA
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from aiida.orm import Code, load_node
from aiida.plugins import DataFactory
from aiida.engine import WorkChain, ToContext, if_
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode, update_params_wf, get_inputs_kkr
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode
from masci_tools.io.common_functions import get_Ry2eV
from aiida.orm import WorkChainNode
from aiida.common.exceptions import InputValidationError


__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.5.0"
__contributors__ = (u"Fabian Bertoldo", u"Philipp Ruessmann")

# ToDo: add more default values to wf_parameters

RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
Dict = DataFactory('dict')
FolderData = DataFactory('folder')

class kkr_flex_wc(WorkChain):
    """
    Workchain of a kkr_flex calculation to calculate the Green function with
    KKR starting from the RemoteData node of a previous calculation (either Voronoi or KKR).

    :param options: (Dict), Workchain specifications
    :param wf_parameters: (Dict), Workflow parameters that deviate from previous KKR RemoteData
    :param remote_data: (RemoteData), mandatory; from a converged KKR calculation
    :param kkr: (Code), mandatory; KKR code running the flexfile writeout
    :param impurity_info: Dict, mandatory: node specifying information
                          of the impurities in the system

    :return workflow_info: (Dict), Information of workflow results
                            like success, last result node, list with convergence behavior
    :return GF_host_remote: (RemoteData), host GF of the system
    """

    _workflowversion = __version__
    _wf_label = 'kkr_flex_wc'
    _wf_description = 'Workflow for a KKR flex calculation starting from RemoteData node of previous converged KKR calculation'


    _options_default = {'queue_name' : '',                        # Queue name to submit jobs too
                        'resources': {"num_machines": 1},         # resources to allowcate for the job
                        'max_wallclock_seconds' : 60*60,          # walltime after which the job gets killed (gets parsed to KKR)}
                        'custom_scheduler_commands' : '',         # some additional scheduler commands
                        'withmpi' : True}                         # execute KKR with mpi or without

    _wf_default = {'ef_shift': 0. ,                                  # set costum absolute E_F (in eV)
                   'dos_params': {'nepts': 61,                       # DOS params: number of points in contour
                                  'tempr': 200, # K                  # DOS params: temperature
                                  'emin': -1, # Ry                   # DOS params: start of energy contour
                                  'emax': 1,  # Ry                   # DOS params: end of energy contour
                                  'kmesh': [10, 10, 10]},            # DOS params: kmesh for DOS calculation (typically higher than in scf contour)
                   'dos_run': False
                   }

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
        Defines the outline of the workflow
        """

        # Take input of the workflow or use defaults defined above
        super(kkr_flex_wc, cls).define(spec)

        spec.input("kkr", valid_type=Code, required=True)
        spec.input("options", valid_type=Dict, required=False,
                       default=Dict(dict=cls._options_default))
        spec.input("wf_parameters", valid_type=Dict, required=False)
        spec.input("remote_data", valid_type=RemoteData, required=True)
        spec.input("impurity_info", valid_type=Dict, required=True)
        spec.input("params_kkr_overwrite", valid_type=Dict, required=False,
                   help="Set some input parameters of the KKR calculation.")

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            if_(cls.validate_input)(
                cls.set_params_flex,
                # TODO: encapsulate this in restarting mechanism (should be a base class of workflows that start calculations)
                # i.e. use base_restart_calc workchain as parent
                cls.get_flex), # calculate host GF and kkr-flexfiles
            cls.return_results)

        # ToDo: improve error codes
        spec.exit_code(101, 'ERROR_INVALID_INPUT_IMP_INFO',
            message="ERROR: the 'impurity_info' input Dict node could not be used")
        spec.exit_code(102, 'ERROR_INVALID_INPUT_KKR',
            message="ERROR: the code you provided for kkr does not use the plugin kkr.kkr")
        spec.exit_code(103, 'ERROR_INVALID_INPUT_REMOTE_DATA',
            message="ERROR: No remote_data was provided as Input")
        spec.exit_code(104, 'ERROR_INVALID_CALC_PARAMETERS',
            message='ERROR: calc_parameters given are not consistent! Hint: did you give an unknown keyword?')
        spec.exit_code(105, 'ERROR_CALC_PARAMETERS_INCOMPLETE',
            message='ERROR: calc_parameters misses keys')
        spec.exit_code(106, 'ERROR_KKR_CALCULATION_FAILED',
            message='ERROR: KKR calculation to write out kkrflex files unsuccessful')

        # specify the outputs
        #spec.output('remote_folder', valid_type=RemoteData)
        spec.output('workflow_info', valid_type=Dict)
        spec.output('GF_host_remote', valid_type=RemoteData)


    def start(self):
        """
        init context and some parameters
        """

        self.report('INFO: started KKR flex workflow version {}'
                    ''.format(self._workflowversion))

        ####### init #######
        # internal para / control para
        self.ctx.abort = False

        # input both wf and options parameters
        options_dict = self.inputs.options.get_dict()
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
            if wf_dict == {}:
                wf_dict = self._wf_default
                self.report('INFO: using default wf parameters')
        else:
            wf_dict = self._wf_default
            self.report('INFO: using default wf parameters')

        if options_dict == {}:
            options_dict = self._options_default
            self.report('INFO: using default options parameters')


        # set values, or defaults
        # ToDo: arrange option assignment differently (look at scf.py from aiida-fleur)
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get('max_wallclock_seconds', self._options_default['max_wallclock_seconds'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])

        self.ctx.ef_shift = wf_dict.get('ef_shift', self._wf_default['ef_shift'])
        self.ctx.dos_run = wf_dict.get('dos_run', self._wf_default['dos_run'])
        self.ctx.dos_params_dict = wf_dict.get('dos_params', self._wf_default['dos_params'])

        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)


        self.report('INFO: use the following parameter:\n'
                    'withmpi: {}\n'
                    'Resources: {}\n'
                    'Walltime (s): {}\n'
                    'queue name: {}\n'
                    'scheduler command: {}\n'
                    'description: {}\n'
                    'label: {}\n'.format(self.ctx.withmpi, self.ctx.resources, self.ctx.max_wallclock_seconds,
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

        inputs = self.inputs
        input_ok = True

        if not 'impurity_info' in inputs:
            input_ok = False
            return self.exit_codes.ERROR_INVALID_INPUT_IMP_INFO

        if 'remote_data' in inputs:
            input_ok = True
        else:
            input_ok = False
            return self.exit_codes.ERROR_INVALID_REMOTE_DATA

        # extract correct remote folder of last calculation if input remote_folder node
        # is not from KKRCalculation but kkr_scf_wc workflow
        input_remote = self.inputs.remote_data
        # check if input_remote has single KKRCalculation parent
        parents = input_remote.get_incoming(node_class=CalcJobNode)
        nparents = len(parents.all_link_labels())
        if nparents!=1:
            # extract parent workflow and get uuid of last calc from output node
            parent_workflow = input_remote.inputs.last_RemoteData
            if not isinstance(parent_workflow, WorkChainNode):
                raise InputValidationError("Input remote_data node neither output of a KKR calculation nor of kkr_scf_wc workflow")
                parent_workflow_out = parent_workflow.outputs.output_kkr_scf_wc_ParameterResults
                uuid_last_calc = parent_workflow_out.get_dict().get('last_calc_nodeinfo').get('uuid')
                last_calc = load_node(uuid_last_calc)
                if not isinstance(last_calc, KkrCalculation):
                    raise InputValidationError("Extracted last_calc node not of type KkrCalculation: check remote_data input node")
                # overwrite remote_data node with extracted remote folder
                output_remote = last_calc.outputs.remote_folder
                self.inputs.remote_data = output_remote

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                error = ("The code you provided for kkr does not "
                         "use the plugin kkr.kkr")
                self.ctx.errors.append(error)
                input_ok = False
                return self.exit_codes.ERROR_INVALID_INPUT_KKR

        # set self.ctx.input_params_KKR
        self.ctx.input_params_KKR = get_parent_paranode(self.inputs.remote_data)

        if input_ok:
            self.report('INFO: checking inputs successful')

        return input_ok



    def set_params_flex(self):
        """
        Take input parameter node and change to input from wf_parameter and options
        """

        self.report('INFO: setting parameters ...')

        params = self.ctx.input_params_KKR
        input_dict = params.get_dict()
        para_check = kkrparams()

        # step 1: try to fill keywords
        try:
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)
        except:
            return self.exit_codes.ERROR_INVALID_CALC_PARAMETERS

        # step 2: check if all mandatory keys are there
        label = ''
        descr = ''
        missing_list = para_check.get_missing_keys(use_aiida=True)
        if missing_list != []:
            kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
            kkrdefaults_updated = []
            for key_default, val_default in list(kkrdefaults.items()):
                if key_default in missing_list:
                    para_check.set_value(key_default, kkrdefaults.get(key_default), silent=True)
                    kkrdefaults_updated.append(key_default)
                    missing_list.remove(key_default)
            if len(missing_list)>0:
                self.report('ERROR: calc_parameters misses keys: {}'.format(missing_list))
                return self.exit_codes.ERROR_CALC_PARAMETERS_INCOMPLETE

            else:
                self.report('updated KKR parameter node with default values: {}'.format(kkrdefaults_updated))
                label = 'add_defaults_'
                descr = 'added missing default keys, '

        # updatedict contains all things that will be updated with update_params_wf later on
        updatedict = {}

        runopt = para_check.get_dict().get('RUNOPT', None)
        if runopt == None:
            runopt = []

        # overwrite some parameters of the KKR calculation by hand before setting mandatory keys
        if "params_kkr_overwrite" in self.inputs:
            for key, val in self.inputs.params_kkr_overwrite.get_dict().items():
                if key != runopt:
                    updatedict[key] = val
                else:
                    runopt = val
                self.report('INFO: overwriting KKR parameter: {} with {} from params_kkr_overwrite input node'.format(key, val))

        runopt = [i.strip() for i in runopt]
        if 'KKRFLEX' not in runopt:
            runopt.append('KKRFLEX')

        updatedict['RUNOPT'] = runopt

        self.report('INFO: RUNOPT set to: {}'.format(runopt))

        if 'wf_parameters' in self.inputs:
            if self.ctx.dos_run:
                for key, val in {'EMIN': self.ctx.dos_params_dict['emin'], 'EMAX': self.ctx.dos_params_dict['emax'], 
                                 'NPT2': self.ctx.dos_params_dict['nepts'], 'NPOL': 0, 'NPT1': 0, 'NPT3': 0, 
                                 'BZDIVIDE': self.ctx.dos_params_dict['kmesh'], 'IEMXD': self.ctx.dos_params_dict['nepts'],
                                 'TEMPR': self.ctx.dos_params_dict['tempr']}.items():
                    updatedict[key] = val
            elif self.ctx.ef_shift != 0:
                # extract old Fermi energy in Ry
                remote_data_parent = self.inputs.remote_data
                parent_calc = remote_data_parent.get_incoming(link_label_filter='remote_folder').first().node
                ef_old = parent_calc.outputs.output_parameters.get_dict().get('fermi_energy')
                # get Fermi energy shift in eV
                ef_shift = self.ctx.ef_shift #set new E_F in eV
                # calculate new Fermi energy in Ry
                ef_new = (ef_old + ef_shift/get_Ry2eV())
                self.report('INFO: ef_old + ef_shift = ef_new: {} eV + {} eV = {} eV'.format(ef_old*get_Ry2eV(), ef_shift, ef_new*get_Ry2eV()))
                updatedict['ef_set'] = ef_new

        #construct the final param node containing all of the params
        updatenode = Dict(dict=updatedict)
        updatenode.label = label+'KKRparam_flex'
        updatenode.description = descr+'KKR parameter node extracted from parent parameters and wf_parameter and options input node.'
        paranode_flex = update_params_wf(self.ctx.input_params_KKR, updatenode)
        self.ctx.flex_kkrparams = paranode_flex
        self.ctx.flex_runopt = runopt

        self.report('INFO: Updated params= {}'.format(paranode_flex.get_dict()))


    def get_flex(self):
        """
        Submit a KKRFLEX calculation
        """

        label = 'KKRFLEX calc.'
        description = 'KKRFLEX calculation to write out host GF'
        code = self.inputs.kkr
        remote = self.inputs.remote_data
        params = self.ctx.flex_kkrparams
        imp_info = self.inputs.impurity_info
        options = {"max_wallclock_seconds": self.ctx.max_wallclock_seconds,
                   "resources": self.ctx.resources,
                   "queue_name": self.ctx.queue}
        if self.ctx.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.ctx.custom_scheduler_commands
        inputs = get_inputs_kkr(code, remote, options, label, description, parameters=params, serial=(not self.ctx.withmpi), imp_info=imp_info)

        # run the KKRFLEX calculation
        self.report('INFO: doing calculation')
        flexrun = self.submit(KkrCalculation, **inputs)

        return ToContext(flexrun=flexrun)


    def return_results(self):
        """
        Return the results of the KKRFLEX calculation.
        This should run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """

        # capture error of unsuccessful flexrun
        if not self.ctx.flexrun.is_finished_ok:
            self.ctx.successful = False
            return self.exit_codes.ERROR_KKR_CALCULATION_FAILED

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['withmpi'] = self.ctx.withmpi
        outputnode_dict['resources'] = self.ctx.resources
        outputnode_dict['max_wallclock_seconds'] = self.ctx.max_wallclock_seconds
        outputnode_dict['queue_name'] = self.ctx.queue
        outputnode_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['pk_flexcalc'] = self.ctx.flexrun.pk
        outputnode_dict['list_of_errors'] = self.ctx.errors

        outputnode = Dict(dict=outputnode_dict)
        outputnode.label = 'kkr_flex_wc_results'
        outputnode.description = ''
        outputnode.store()

        # return the input remote_data folder as output node
        #self.out('remote_data', self.inputs.remote_data)
        # return Dict node containing information about previous calculation
        self.out('workflow_info', outputnode)
        # return retrieved data from kkrflex calculation
        self.out('GF_host_remote', self.ctx.flexrun.outputs.remote_folder)

        self.report('INFO: created GF writeout result nodes')

#        self.report("INFO: create GF writeout results nodes: outputnode={}".format(outputnode))
#        try:
#            self.report("INFO: create GF writeout results nodes. KKRFLEX calc retrieved node={}".format(self.ctx.flexrun.outputs.retrieved))
#            has_flexrun = True
#        except AttributeError as e:
#            self.report("ERROR: no KKRFLEX calc retrieved node found")
#            self.report("Caught AttributeError {}".format(e))
#            has_flexrun = False

        #for link_name, node in outdict.iteritems():
            #self.report("INFO: storing node '{}' with link name '{}'".format(node, link_name))
            #self.report("INFO: node type: {}".format(type(node)))
            #self.out(link_name, node)

        self.report("INFO: done with KKRFLEX GF writeout workflow!\n")
#        self.report("Successful run: {}".format(has_flexrun))
