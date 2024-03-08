# -*- coding: utf-8 -*-
# Workflow for simulating an STM around an impurity

from aiida.engine import WorkChain, ToContext, if_
from aiida.orm import Dict, RemoteData, Code, CalcJobNode, WorkChainNode, Float, Bool, XyData, SinglefileData, List
from aiida.orm import Group, load_group
from aiida_kkr.workflows import kkr_imp_wc, kkr_imp_dos_wc, kkr_dos_wc, kkr_flex_wc
from aiida_kkr.tools.find_parent import get_calc_from_remote
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode

__copyright__ = (u'Copyright (c), 2024, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.1'
__contributors__ = (u'Raffaele Aliberti, David Antognini Silva, Philipp Rüßmann')
_VERBOSE_ = True

class STM_wc(WorkChain):
    """
    Workchain for the Simulation of a (superconducting) STM tip around a magnetic impurity.
    The workchain uses the (converged) impurity calculation of a host system and combines them
    with vacuum sites in positions 

    inputs::
        :param options: (Dict), computer options
        :param tip_position: (Dict), specify the position of the STM tip
        :param impurity_info: (Dict), information about the impurity cluster
        :param host_calc: (RemoteData), information about the host structure of the sample 
        :param wf_parameters: (Dict), parameters that are used to run the workflow
        :param kkr: (Code), KKR host code for the writing out kkrflex files
        :param kkrimp: (Code), KKR impurity code for the normal state impurity scf and BdG impurity DOS calculation
        
     returns::
     
        :return workflow_info: (Dict), Information of workflow results
                            like success, last result node, list with convergence behavior
        :return STM_lmdos: (XYData), Returns the plot of the lmDOS of the calculation 
        :retrun STM_lmdos_interpol: (XYData), Returns the interpolated lmDOS of the calculation"""
    
    # TO DO: Ask how the stm tool works for the calculation of the impurity combined radius. 
    # TO DO: Ask how to run the code on the local machine for the testing.
    # TO DO: Find most efficient way to parallelize the calculations : matrix, dict, group (?)
    # TO DO: Add a parameter to include spherical cluster.
    # TO DO: Add BdG_Imp_dos setting for the workflow.
    # TO DO: Fix the KKRImpflex step, It loads random values
        
     
    _wf_version = __version__
    _wf_label = 'STM_wc'
    _wf_description = 'Workflow for simulating an STM measurement'

    _options_default = {
        'queue_name': '',  # Queue name to submit jobs too
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 48, 'num_cores_per_mpiproc': 1},  # resources to allocate for the job
        'max_wallclock_seconds': 3600*2,  # walltime after which the job gets killed (gets parsed to KKR)}
        'custom_scheduler_commands': '',  # some additional scheduler commands
        'withmpi': True
    }  # execute KKR with mpi or without

    _wf_default = {
        'jij_run': False,  # calculate Jij's energy resolved
        'lmdos': True,  # calculate also (l,m) or only l-resolved DOS, for the STM wf this is alwyas set on True as default
        'retrieve_kkrflex': False,  # retrieve kkrflex files to repository or leave on remote computer only
    }
    # add defaults of dos_params since they are passed onto that workflow
    _wf_default['dos_params'] = kkr_imp_dos_wc.get_wf_defaults()['dos_params']
    
    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(STM_wc, cls).define(spec)
        
        spec.input('kkr', valid_type=Code, required=False, 
             help='KKRhost code, needed if gf_dos_remote is not given.'
                  )
        
        spec.input('kkrimp', valid_type=Code, required=True, 
             help='KKRimp code, always needed.'
                  )
        
        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help='Computer options (resources, quene name, etc.).'
        )
        
        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
            help='Workflow parameter (see `kkr_dos_wc.get_wf_defaults()`).'
        )
        
        spec.input(
            'tip_position',
            valid_type=Dict,
            required=False,
            # Find a way to create an area of this position. 
            default=lambda:Dict({'ilayer': 0, 'nx': 0, 'ny': 0}),
            # In the previous line we have set to study the first layer
            # nx is the number of (symmetric) steps that we take in the x-direction starting from the impurity
            # ny is the number of (symmetric) steps that we take in the y-direction starting from the impurity
            # (0,0) correspond to calculate the DOS only on the impurity site
            help=
            'Position of the STM head'
        )
                
        spec.input(
            'imp_info',
            valid_type=Dict,
            required=True,
            help=
            'Information of the impurity like position in the unit cell, screening cluster, atom type.'
        )
        
        spec.input(
            'host_calc',
            valid_type=RemoteData,
            required=False,
            help=
            'The information about the clean host structure is required in order to continue the cluster'
            'Inside a bigger host structure with empty sites.'
        )
        
    
        spec.input(
            'host_remote',
            valid_type=RemoteData,
            required=True,
            help='Remote Data containing the remote folder from the outputs of the host calculation',
        )
        
        spec.input(
            'imp_potential_node',
            valid_type=SinglefileData,
            required=True,
            help=
            'Impurity potential node',
        )
        
        spec.input(
            'remote_data',
            valid_type=RemoteData,
            required=False,
            help=
            'Remote data from a converged kkr calculation, required for the gf writeout step',
        )
        
        spec.input(
            'kkrflex_files',
            valid_type=RemoteData,
            required=False,
            help=
            'with this input we can directly load the gf_dos files without calculating them')
        
        #spec.expose_input(
        #    kkr_flex_wc,
        #    namespace='gf_writeout',
        #    include=('wf_parameters')
        #)
        
        # Specify the possible outputs
        spec.output('tip_position', valid_type=Dict, required=True)
        
        #spec.output('combined_potential', valid_type=Dict, required=True)
                 
        #spec.output('combined_node',nvalid_type=SinglefileData,nrequired=True)
        
        spec.output('STM_data', valid_type=List, required=True)
        
        #spec.output('STM_lmdos_data_interpol', valid_type=XyData, required=True)
        
        #spec.output('workflow_info', valid_type=Dict)
        
        spec.output('kkrflex', valid_type=RemoteData)
            
        # Define all possible error messages
        
        spec.exit_code(
            100, 'ERROR_STM_POSITION_NOT_VALID', 'The position provided for the STM probe are incorrect'
        )
        spec.exit_code(
            101, 'ERROR_IMP_INFO_NOT_CORRECT', 'The node provided for the impurity info is not valid'
        )
            
            
        spec.outline(
            # For initializing workflow
            cls.start,
            # We first run the gf
            #cls.gf_writeout_run,
            # The gf is then used to evaluate the STM lmdos
            cls.STM_lmdos_run, 
            # Data aggregator, used to make the final result more user friendly
            cls.aggregate_results,
            
            cls.results
        )
        
            
    def combine_potentials(self,impurity_to_combine, da, db):
        from aiida_kkr.tools.tools_STM_scan import get_imp_info_add_position_cf
        """
        Here we want to combine the impurity potential and the host potential
        """
        tip_position = {}
        
        tip_position['ilayer'] = self.inputs.tip_position['ilayer']
        tip_position['da'] = da
        tip_position['db'] = db
        
        imp_info = self.inputs.imp_info #(impurity to combine)
        host_remote = self.inputs.host_remote
        
        combined_potential = get_imp_info_add_position_cf(host_remote, imp_info, tip_position)
        #self.ctx.combined_potential = combined_potential
        return combined_potential
            
    def combine_nodes(self, node_to_combine, da, db):
        from aiida_kkr.tools.tools_STM_scan import create_combined_potential_node_cf
        """
        Here we create a combined potential node from the host potential (no impurity)
        and from the impurity potential
        """
        tip_position = {}
    
        tip_position['ilayer'] = self.inputs.tip_position['ilayer']
        tip_position['da'] = da
        tip_position['db'] = db
        
        imp_potential_node = self.inputs.imp_potential_node # (node_to_combine)
        host_remote = self.inputs.host_remote
        
        combined_node = create_combined_potential_node_cf(tip_position, host_remote, imp_potential_node)
        #self.ctx.combined_node = combined_node
        return combined_node
    
    def start(self):
        """
        Initialise context and some parameters
        """
                
        self.report(f'INFO: started STM workflow version {self._wf_version}')
        if _VERBOSE_:
            self.report(f'inputs: {self.inputs}')
        
        # Input both wf and options parameters
        # We check if the inputs for the wc where given
        # Otherwise we assign the default values 
        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
            # empty dictionary evaluate as false to Python
            if options_dict == {}:
                options_dict = self._options_default
                self.report('INFO: using default options parameters')
                
        if 'wf_parameters' in self.inputs:
            wf_param_dict = self.inputs.wf_parameters.get_dict()
            # empty dictionary evaluate as false to Python
            if wf_param_dict == {}:
                options_dict = self._wf_default
                self.report('INFO: usign defalut wf parameters')
                
        # In this section we assign the computational resources to the builder
        
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

        # Set workflow parameters for the KKR imputrity calculations
        """This part is really important, this should always be set to True for an STM calculation"""
        self.ctx.lmdos = wf_param_dict.get('lmdos', self._wf_default['lmdos'])
        
        self.ctx.retrieve_kkrflex = wf_param_dict.get('retrieve_kkrflex', self._wf_default['retrieve_kkrflex'])

        self.ctx.dos_params_dict = wf_param_dict.get('dos_params', self._wf_default['dos_params'])
        
        # fill missing key, value pairs with defaults
        for k, v in self._wf_default['dos_params'].items():
            if k not in self.ctx.dos_params_dict.keys():
                self.ctx.dos_params_dict[k] = v

        # set workflow parameters for the KKR impurity calculation
        self.ctx.jij_run = wf_param_dict.get('jij_run', self._wf_default['jij_run'])

        # set workflow label and description
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)

        # whether or not to compute the GF writeout step
        self.ctx.skip_gfstep = False
            
        #if 'kkrflex_files' in self.inputs:
        #    self.ctx.skip_gfstep = True
        
        message = f"""
                        INFO: use the following parameter:
                        withmpi: {self.ctx.withmpi}
                        Resources: {self.ctx.resources}
                        Walltime (s): {self.ctx.max_wallclock_seconds}
                        queue name: {self.ctx.queue}
                        scheduler command: {self.ctx.custom_scheduler_commands}
                        description: {self.ctx.description_wf}
                        label: {self.ctx.label_wf}
                  """
        print(message)
        self.report(message)

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''
        
        
    #def validate_input(self):
    #    
    #    inputs = self.inputs
    #    inputs_ok = True
    #    gf_writeout_calc = None
    #    
    #    if 'imp_potential_node' in inputs:
    #        # check if input potential has incoming return link
    #        if len(inputs.imp_potential_node.get_incoming(link_type=LinkType.RETURN).all()) == 0:
    #            inputs_ok = self._imp_pot_not_from_wf(inputs_ok)
    #        else:
    #            gf_writeout_calc = self._imp_pot_from_wf()
    #
    #    if 'gf_dos_remote' in self.inputs:
    #        self.ctx.skip_gfstep = True
    #    else:
    #        inputs_ok = self._check_gf_writeout_inputs(inputs_ok, gf_writeout_calc)
    
    def gf_writeout_run(self):
        """This function would allow the workchian to evaluate the needed kkrflex files 
           And then use these for all the DOS calculation without the need for the kkrimpdos to evaluate
           The kkrflex files for each DOS step"""
            
        if not self.ctx.skip_gfstep:
            from aiida_kkr.workflows import kkr_flex_wc
            from aiida.engine import submit
            
            builder = kkr_flex_wc.get_builder()
            
            builder.kkr = self.inputs.kkr
            builder.options = self.ctx.options_params_dict
            builder.remote_data = self.inputs.remote_data
            builder.impurity_info = self.inputs.imp_info
            #builder.wf_parameters = self.ctx.kkr_params
            
            gf_run = self.submit(builder)
            
            message = f'INFO: running the gf_writeout step for impurity system in the STM workflow(pk: {gf_run.pk})'
            print(message)
            self.report(message)
            
            return ToContext(gf_node=gf_run)
        
    def impurity_cluster_evaluation(self, impurity_info, imp_potential_node, x, y):
        
        # Here we create an impurity cluster that has inside all the positions on which the STM will scan
        
        for da in np.arange(-x, x+1, 1):
            for db in np.arange(-y, y+1, 1):
                
                #tmp = self.combine_potentials(impurity_info, da, db)
                #impurity_info = tmp
                
                tmp = self.combine_nodes(imp_potential_node, da, db)
                imp_potential_node = tmp
                
        return  imp_potential_node
                     

    def STM_lmdos_run(self):
        """In this part of the worflow we want to simulate the lmdos which a STM is able to measure """
        
        # First we would like to distinguish between an impurity dos and a normal state calculation
         # Builder from tge     
        builder = kkr_imp_dos_wc.get_builder()
        
        # Code loading
        builder.kkrimp = self.inputs.kkrimp # needed for the kkr_imp_dos_wc
        builder.kkr = self.inputs.kkr # needed to evaluate the kkr_flex files in the DOS step
        
        # Builder options
        builder.options = self.ctx.options_params_dict
        
        #Testing kkrflex files
        #builder.gf_dos_remote = self.inputs.kkrflex_files
        #builder.kkrimp_remote = self.inputs.remote_data
        
        #if not self.ctx.skip_gfstep:
        #    # use computed gf_writeout
        #    if not self.ctx.gf_node.is_finished_ok:
        #        return self.exit_codes.ERROR_GF_WRITEOUT_UNSUCCESFUL  # pylint: disable=no-member
        #    else:
        #        builder.gf_dos_remote = self.ctx.gf_node.outputs.GF_host_remote
        #        kkrflex_out = self.ctx.gf_node.outputs.GF_host_remote
        #else:
        #    # use gf_writeout from input
        #    builder.gf_dos_remote = self.inputs.kkrflex_files
        #    kkrflex_out = self.inputs.kkrflex_files
        
        # load the information of the impurity
        #nspin = gf_writeout_calc.outputs.output_parameters.get_dict().get('nspin')
        #self.ctx.nspin = nspin
        #self.report(f'nspin: {nspin}')
        
        self.ctx.kkrimp_params_dict = Dict(dict={
            'nspin': 2, # take this information from the input 
            'nsteps': 1,
            'kkr_runmax': 1,
            'dos_run': True,
            'lmdos': self.ctx.lmdos,
            'jij_run': self.ctx.jij_run,
            'dos_params': self.ctx.dos_params_dict
        })
        
        
        # We want to set the energy to the Fermi level
        
        self.ctx.kkrimp_params_dict['dos_params']['emin'] = 0-0.005
        self.ctx.kkrimp_params_dict['dos_params']['emax'] = 0+0.005
            
        # Finally we overwrite the number of energy points to 1
        # This is because we want many epoints around the impurity position
            
        self.ctx.kkrimp_params_dict['dos_params']['nepts'] = 9 # Here 9 because of the interpolated files that aren't generated
        
        #builder.metadata.label = label_imp  # pylint: disable=no-member
        #builder.metadata.description = description_imp  # pylint: disable=no-member
        
        builder.wf_parameters = self.ctx.kkrimp_params_dict       
        # Host remote files that will be used for the actual plot step.
        builder.host_remote = self.inputs.host_remote
        
        # We now want to iterate over several in-plane positions. 
        # These are the number of vectors in which we want to move the STM tip.
        x = self.inputs.tip_position['nx']
        y = self.inputs.tip_position['ny']
        
        # Here we create the impurity cluster for the STM scanning tool
        imp_pot_sfd = self.impurity_cluster_evaluation(self.inputs.imp_info, self.inputs.imp_potential_node, x, y)
        
        Dict_result = {}
        # Indexing of the node number
        self.ctx.num_nodes = 0
        
        for da in np.arange(-x, x+1, 1):
                
            for db in np.arange(-y, y+1, 1):
                
                # Insertion of the combined potential at the position da, db
                builder.impurity_info = self.combine_potentials(self.inputs.imp_info, da, db)
                
                # Insertion of the combined potential node at the positon da, db
                builder.imp_pot_sfd = imp_pot_sfd
                
                # Insertion of the remote data for the gf step
                #builder.kkrimp_remote = self.inputs.remote_data
                                
                position = f'{da} {db}'
                builder.metadata.label = position
                #Here we store the information of the postion of the vector
                
                calc = self.submit(builder)
                message = f"""INFO: running DOS step for an STM measurement (pk: {calc.pk}) at position
                          (ilayer: {self.inputs.tip_position['ilayer']}, da: {da}, db: {db} )"""
                print(message)
                self.report(message)
                
                Dict_result[f'node_{self.ctx.num_nodes}'] = calc
                self.ctx.num_nodes += 1
                
                
                        
        """Next Gen features"""
        # LDA+U settings
        #if 'settings_LDAU' in self.inputs:
        #    self.report('Add settings_LDAU input node')
        #    builder.settings_LDAU = self.inputs.settings_LDAU

        #if 'params_overwrite' in self.inputs.BdG:
        #    builder.params_overwrite = self.inputs.BdG.params_overwrite
        #if 'initial_noco_angles' in self.inputs:
        #    builder.initial_noco_angles = self.inputs.initial_noco_angles
        #if 'rimpshift' in self.inputs:
        #    builder.rimpshift = self.inputs.rimpshift
        
        #lis = []
        #STM_lmdos = {}
        return ToContext(**Dict_result)
    
    def aggregate_results(self):
        
        # Function returns the data in a user-friendly way
        
        array = []
        
        for number in range(self.ctx.num_nodes):
            
            data = []
            # Here we select the node 
            attribute_name = f'node_{number}'
            current_node = getattr(self.ctx, attribute_name)
            
            # Now we want to extract the positions for the label 
            string = current_node.label            
            data.append(string)
                
            
            #for part in string:
            #    key, value = part.split(':')
            #    extracted_values[key] = value
            #    
            #da_value = extracted_values.get('x')
            #db_value = extracted_values.get('y')
            #
            #data.append(int(da_value))
            #data.append(int(db_value))
            
            # Here we append the lm dos
            data.append(current_node.outputs.dos_data_lm)
            
            array.append(data)
        
        
        
        
        #for k, v in self.ctx.node00.items():
        #    
        #    data = []
        #    
        #    string = k.strip('()')
        #    var = string.split(',')
        #    
        #    # Append the da position
        #    data.append(int(var[0]))
        #    # Append the db position
        #    data.append(int(var[1]))
        #
        #    # Access the da and db position of the lattice site
        #    data.append(v.outputs.dos_data_lm)
        #    
        #    lis.append(data)
        
        message = 'INFO: created output nodes for KKR STM workflow.'
        print(message)
        self.report(message)
        self.report(
            '\n'
            '|------------------------------------------------------------------------------------------------------------------|\n'
            '|-------------------------------------| Done with the STM workflow! |----------------------------------------------|\n'
            '|------------------------------------------------------------------------------------------------------------------|'
        )
        
        self.ctx.STM_data = orm.List(array)
        
        #last_calc_pk = self.ctx.kkrimp_dos.outputs.workflow_info.get_dict().get('last_calc_nodeinfo')['pk']
        #last_calc = load_node(last_calc_pk)
        #last_calc_output_params = last_calc.outputs.output_parameters
        #last_calc_info = self.ctx.kkrimp_dos.outputs.workflow_info
        #outputnode_dict = {}
        #outputnode_dict['impurity_info'] = self.ctx.imp_info.get_dict()
        #outputnode_dict['workflow_name'] = self.__class__.__name__
        #outputnode_dict['workflow_version'] = self._workflowversion
        #if not self.ctx.skip_gfstep:
        #    outputnode_dict['used_subworkflows'] = {'gf_writeout': self.ctx.gf_writeout.pk}
        #else:
        #    outputnode_dict['used_subworkflows'] = {}
        #outputnode_dict['used_subworkflows']['impurity_dos'] = self.ctx.kkrimp_dos.pk
        ## interpol dos file and store to XyData nodes
        #dos_extracted, dosXyDatas = self.extract_dos_data(last_calc)
        #message = f'INFO: extracted DOS data? {dos_extracted}'
        #print(message)
        #self.report(message)
        ## create results node and link rest of results
        #link_nodes = dosXyDatas.copy()
        #link_nodes['last_calc_output_parameters'] = last_calc_output_params
        #link_nodes['last_calc_remote'] = last_calc.outputs.remote_folder
        #outputnode_t = create_out_dict_node(Dict(dict=outputnode_dict), **link_nodes)
        #outputnode_t.label = 'kkr_imp_dos_wc_inform'
        #outputnode_t.description = 'Contains information for workflow'
        #if dos_extracted:
        #    self.out('dos_data', dosXyDatas['dos_data'])
        #    if 'dos_data_interpol' in dosXyDatas:
        #        self.out('dos_data_interpol', dosXyDatas['dos_data_interpol'])
        #    if self.ctx.lmdos:
        #        self.out('dos_data_lm', dosXyDatas['dos_data_lm'])
        #        if 'dos_data_interpol_lm' in dosXyDatas:
        #            self.out('dos_data_interpol_lm', dosXyDatas['dos_data_interpol_lm'])
        #message = f'INFO: workflow_info node: {outputnode_t.uuid}'
        #print(message)
        #self.report(message)
        #self.out('workflow_info', outputnode_t)
        #self.out('last_calc_output_parameters', last_calc_output_params)
        #self.out('last_calc_info', last_calc_info)
    
    def results(self):

        # Declaring the output
        self.out("STM_data", self.ctx.STM_data)
        #self.out("STM_lmdos_data_interpol", self.ctx.STM_lmdos.outputs.dos_data_interpol_lm)
        #self.out("workflow_info", self.ctx.STM_lmdos.outputs.workflow_info)
        self.out("tip_position", self.inputs.tip_position)
        #self.out("kkrflexfiles",kkrflex_out)
      
