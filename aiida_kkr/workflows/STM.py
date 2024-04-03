#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Workflow for STM swf_parameterson around a magnetic impurity

from aiida.engine import WorkChain, ToContext, if_, calcfunction
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
        :return STM_dos_data: (XYData), Returns the plot of the lmDOS of the calculation 
        :retrun STM_lmdos_data: (XYData), Returns the interpolated lmDOS of the calculation"""
    
    # TO DO: Add a parameter to include spherical cluster.
    # TO DO: Add BdG_Imp_dos setting for the workflow.    
    # TO DO: Add the initialize step.
    # TO DO: Add the 'path creation' step.
    # TO DO: Add the point symmetry analysis.
    # TO DO: Add check that between the ilayer and the actual number of layers in the structure.
    # TO DO: Add to the outputs the calculated imp_info and imp_potential.
     
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
    
        spec.expose_inputs(kkr_imp_dos_wc, namespace='NSHELD', include=('params_overwrite'))
    
        # Specify the possible outputs
        spec.output('tip_position', valid_type=Dict)
        
        spec.output('STM_dos_data', valid_type=XyData, required=True)
        
        spec.output('STM_dos_data_lmdos', valid_type=XyData, required=True)
        
        #spec.output('workflow_info', valid_type=Dict)
        
        spec.output('kkrflexfiles', valid_type=RemoteData)
        
        spec.output('combined_imp_info', valid_type=Dict)
        
        spec.output('combined_imp_potential', valid_type=SinglefileData)
            
        # Define all possible error messages
        
        spec.exit_code(
            100, 'ERROR_STM_POSITION_NOT_VALID', 'The position provided for the STM probe are incorrect'
        )
        spec.exit_code(
            101, 'ERROR_IMP_INFO_NOT_CORRECT', 'The node provided for the impurity info is not valid'
        )
        spec.exit_code(
            102, 'ERROR_NO_IMP_POT_SFD', 'No impurity node has been given in the intput'
        )
        spec.exit_code(
            103, 'ERROR_NO_IMPURITY_INFO', 'No impurity info has been given in the input'
        )
        spec.exit_code(
            104, 'ERROR_NO_DATA_FOR_THE_GF_STEP', """Neither the kkrflex files nor the KKR builder have been given. Please
                                                     provide already converged kkrflex files, or the kkr builder to evaluate them"""
        )
        spec.exit_code(
            201, 'ERROR_IMP_SUB_WORKFLOW_FAILURE', 'A step in the kkr_imp_dos workflow has failed' 
        )
            
            
        spec.outline(
            # For initializing workflow
            cls.start,
            # We first aggregate all the impurity data
            # The gf is then used to evaluate the STM lmdos
            #cls.gf_writeout_run,
            
            cls.STM_lmdos_run, 
            # Data aggregator, used to make the final result more user friendly
            # cls.finalize_results,
            
            cls.results
        )
        
            
    def combine_potentials(self, impurity_to_combine, da, db):
        from aiida_kkr.tools.tools_STM_scan import get_imp_info_add_position_cf
        """
        Here we want to combine the impurity information and the host information 
        """
       
        imp_info = self.inputs.imp_info #(impurity to combine)
        host_remote = self.inputs.host_remote

        # Since the objects in AiiDA are immutable we have to create a new dictionary and then convert  
        # it to the right AiiDA type
        tip_position = {}
        tip_position['ilayer'] = self.inputs.tip_position['ilayer']
        tip_position['da'] = da
        tip_position['db'] = db
        
        combined_imp_info = get_imp_info_add_position_cf(tip_position, host_remote, imp_info)
        # Add check to see if imp_cls is there
        if 'imp_cls' in impurity_to_combine:
            
            new_combined_imp_info = {}
            
            for key, value in impurity_to_combine.items():
                if key == 'Zimp':
                    new_combined_imp_info[key] = impurity_to_combine[key]
                    new_combined_imp_info[key].append(combined_imp_info[key][-1])
                else:
                    # Here we have lists of list that we need to confront
                    new_combined_imp_info[key] = impurity_to_combine[key]
                    set_tmp = [set(row) for row in impurity_to_combine[key]]
                    
                    new_combined_imp_info[key] += [row for row in combined_imp_info[key] if set(row) not in set_tmp]
                    
            new_combined_imp_info = orm.Dict(dict=new_combined_imp_info)

        else:
            
            new_combined_imp_info = combined_imp_info            
                        
        return new_combined_imp_info
            
    def combine_nodes(self, node_to_combine, da, db):
        from aiida_kkr.tools.tools_STM_scan import create_combined_potential_node_cf
        """
        Here we create a combined potential node from the host potential (no impurity)
        and from the impurity potential
        """
        #imp_potential_node = self.inputs.imp_potential_node # (node_to_combine).
        host_remote = self.inputs.host_remote # the remote host structure remains the same.
    
        # Since the objects in AiiDA are immutable we have to create a new dictionary and then convert  
        # it to the right AiiDA type
        
        tip_position = {}
        tip_position['ilayer'] = self.inputs.tip_position['ilayer'] # for now we require that the z position remains the same.
        tip_position['da'] = da
        tip_position['db'] = db

        combined_node = create_combined_potential_node_cf(tip_position, host_remote, node_to_combine)                
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
        
        
    def validate_input(self):
        
        inputs = self.inputs
        inputs_ok = True
        gf_writeout_calc = None
        
        if not 'imp_potential_node' in inputs:
            inputs_ok = False
            return self.exit_codes.ERROR_NO_IMP_POT_SFD
        
        if not 'imp_info' in inputs:
            inputs_ok = False
            return self.exit_codes.ERROR_NO_IMP_INFO
        
        if not 'kkrflex_files' and 'kkr' in inputs:
            inputs_ok = False
            return self.exit_codes.ERROR_NO_DATA_FOR_THE_GF_STEP
            
        
        if 'imp_potential_node' in inputs:
            # check if input potential has incoming return link
            if len(inputs.imp_potential_node.get_incoming(link_type=LinkType.RETURN).all()) == 0:
                inputs_ok = self._imp_pot_not_from_wf(inputs_ok)
            else:
                gf_writeout_calc = self._imp_pot_from_wf()
    
        if 'gf_dos_remote' in self.inputs:
            self.ctx.skip_gfstep = True
        else:
            inputs_ok = self._check_gf_writeout_inputs(inputs_ok, gf_writeout_calc)
    
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
        
    def impurity_cluster_evaluation(self):
        
        # Here we create an impurity cluster that has inside all the positions on which the STM will scan
        
        # We now want to iterate over several in-plane positions. 
        # These are the number of vectors in which we want to move the STM tip.
        x = self.inputs.tip_position['nx']
        y = self.inputs.tip_position['ny']
        
        impurity_info = self.inputs.imp_info # for the first step we combine the impurity info from the input
        imp_potential_node = self.inputs.imp_potential_node # for the first step we combine the impurity node from the input    
                
        # Aggregation of the impurity info
       
        for da in np.arange(x, x+1, 1):
            for db in np.arange(-y, y+1, 1):
        
                tmp_imp_info = self.combine_potentials(impurity_info, da, db)
                impurity_info = tmp_imp_info
                        
                # Aggregation the impurity nodes
                tmp_imp_pot = self.combine_nodes(imp_potential_node, da, db)
                imp_potential_node = tmp_imp_pot
                
        return impurity_info, imp_potential_node      

    def STM_lmdos_run(self):
        """In this part of the worflow we want to simulate the lmdos which a STM is able to measure """
        
        # First we would like to distinguish between an impurity dos and a normal state calculation     
        builder = kkr_imp_dos_wc.get_builder()
        
        # Code loading
        builder.kkrimp = self.inputs.kkrimp # needed for the kkr_imp_dos_wc
       
        # Builder options
        builder.options = self.ctx.options_params_dict
        
        # Check if the kkrflex files are already given in the outputs
        if 'kkrflex_files' in self.inputs:
            builder.gf_dos_remote = self.inputs.kkrflex_files
        else:
            builder.kkr = self.inputs.kkr # needed to evaluate the kkr_flex files in the DOS step
        
        self.ctx.kkrimp_params_dict = Dict(dict={
            'nsteps': 1,
            'kkr_runmax': 1,
            'dos_run': True,
            'lmdos': self.ctx.lmdos,
            'jij_run': self.ctx.jij_run,
            'dos_params': self.ctx.dos_params_dict
        })
        
        builder.params_kkr_overwrite = orm.Dict({'NSHELD': 1500})
        
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
    
        
        # Here we create the impurity cluster for the STM scanning tool
        impurity_info, imp_pot_sfd = self.impurity_cluster_evaluation()
        
        # impurity info for the workflow
        builder.impurity_info = impurity_info    
        builder.imp_pot_sfd = imp_pot_sfd
        
        x = self.inputs.tip_position['nx']
        y = self.inputs.tip_position['ny']
        
        calc = self.submit(builder)
        message = f"""INFO: running DOS step for an STM measurement (pk: {calc.pk}) at position
                          (ilayer: {self.inputs.tip_position['ilayer']}, da: {x}, db: {y} )"""
        
        print(message)
        self.report(message)
        
        return ToContext(STM_data = calc)
       
    
    def finalize_results(self):
        list_dos = []
        list_label = []
        
        for number in range(self.ctx.num_nodes):
            
            attribute_name = f'node_{number}'
            current_node = getattr(self.ctx, attribute_name)
            # Here we return the node where the calculation node itself
            list_dos.append(current_node.outputs.dos_data_lm)
            list_label.append(current_node.label)
        
        STM_data = aggregate_results(list_dos, list_label)
        return ToContext(STM_data = STM_data)
        
        
    def results(self):
        
        if not self.ctx.STM_data.is_finished_ok:
            
            message = 'ERROR: sub workflow for STM calculation failed'
            print(message)
            self.report(message)
            return self.exit_codes.ERROR_IMP_SUB_WORKFLOW_FAILURE
        
        else:

            # Declaring the output
            self.out("STM_dos_data", self.ctx.STM_data.outputs.dos_data)
            self.out("STM_dos_data_lmdos", self.ctx.STM_data.outputs.dos_data_lm)
            #self.out("workflow_info", self.ctx.STM_lmdos.outputs.workflow_info)
            self.out("tip_position", self.inputs.tip_position)
            try:
                self.out("kkrflexfiles",self.ctx.STM_data.outputs.gf_dos_remote)
            except:
                pass
        
        #self.out("combined_imp_info", impurity_info)
        #self.out("combined_imp_potential", imp_pot_sfd)
    
        
        message = 'INFO: created output nodes for KKR STM workflow.'
        print(message)
        self.report(message)
        self.report(
            '\n'
            '|------------------------------------------------------------------------------------------------------------------|\n'
            '|-----------------------------------------| Done with the STM workflow! |------------------------------------------|\n'
            '|------------------------------------------------------------------------------------------------------------------|'
        )