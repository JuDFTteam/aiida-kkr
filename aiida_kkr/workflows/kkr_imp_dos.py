#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a impurity DOS calculation and
some helper methods to do so with AiiDA
"""
from __future__ import print_function

from aiida.orm import Code, DataFactory, load_node
from aiida.work.workchain import if_, ToContext, WorkChain
from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc

__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.2"
__contributors__ = u"Fabian Bertoldo"

#TODO: improve workflow output node structure
#TODO: generalise search for imp_info and conv_host from startpot
    

ParameterData = DataFactory('parameter')  
RemoteData = DataFactory('remote') 
SinglefileData = DataFactory('singlefile') 


class kkr_imp_dos_wc(WorkChain):
    """
    Workchain of a DOS calculation for an impurity system starting from a
    converged impurity calculation or workflow
    
    :param options: (ParameterData), computer options
    :param wf_parameters: (ParameterData), specifications for the DOS
    :param kkr: (Code), mandatory: KKR code for gf_writeout step
    :param kkrimp: (Code), mandatory: KKRimp code for DOS calculation
    :param imp_host_pot: (SinglefileData), mandatory: impurity startpotential
    
    :return workflow_info: (ParameterData), Information on workflow results
    :return last_calc_output_parameters: (ParameterData), output parameters of 
                                         the last called calculation
    :return last_calc_info: (ParameterData), information of the last called calculation
    """
    
    _workflowversion = __version__
    _wf_label = 'kkr_imp_dos_wc'
    _wf_description = 'Workflow for a KKR impurity DOS calculation'
       

    _options_default = {'queue_name' : '',                        # Queue name to submit jobs too
                        'resources': {"num_machines": 1},         # resources to allocate for the job
                        'max_wallclock_seconds' : 60*60,          # walltime after which the job gets killed (gets parsed to KKR)}
                        'custom_scheduler_commands' : '',         # some additional scheduler commands 
                        'use_mpi' : False}                        # execute KKR with mpi or without
                        
    _wf_default = {'ef_shift': 0. ,                                  # set costum absolute E_F (in eV)
                   'dos_params': {'nepts': 61,                       # DOS params: number of points in contour
                                  'tempr': 200, # K                  # DOS params: temperature
                                  'emin': -1, # Ry                   # DOS params: start of energy contour
                                  'emax': 1,  # Ry                   # DOS params: end of energy contour
                                  'kmesh': [30, 30, 30]},            # DOS params: kmesh for DOS calculation (typically higher than in scf contour)
                   'non_spherical': 1,                      # use non-spherical parts of the potential (0 if you don't want that)
                   'born_iter': 2,                          # number of Born iterations for the non-spherical calculation
                   'init_pos' : None,                       # position in unit cell where magnetic field is applied [default (None) means apply to all]
                   'spinorbit' : False,                     # SOC calculation (True/False)
                   'newsol' : False }                       # new SOC solver is applied  
                   
                   
                   
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
        super(kkr_imp_dos_wc, cls).define(spec)
        
        spec.input("kkr", valid_type=Code, required=True)  
        spec.input("kkrimp", valid_type=Code, required=True)
        spec.input("host_imp_pot", valid_type=SinglefileData, required=True)        
        spec.input("options", valid_type=ParameterData, required=False,
                       default=ParameterData(dict=cls._options_default))
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                       default=ParameterData(dict=cls._wf_default))
    
        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,                  # start and initialise workflow
            if_(cls.validate_input)(    # validate the given input
                cls.run_gfstep,         # run GF step with DOS energy contour
                cls.run_imp_dos),       # run DOS for the impurity problem
            cls.return_results          # terminate workflow and return results
            )

        # Define possible exit codes for the workflow
        #TBD


        # specify the outputs
        spec.output('workflow_info', valid_type=ParameterData)
        spec.output('last_calc_output_parameters', valid_type=ParameterData)
        spec.output('last_calc_info', valid_type=ParameterData)
        
        
    def start(self):
        """
        Initialise context and some parameters
        """
    
        self.report('INFO: started KKR impurity DOS workflow version {}'
                    ''.format(self._workflowversion))
    
        # input both wf and options parameters
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
            if wf_dict == {}:
                wf_dict = self._wf_default
                self.report('INFO: using default wf parameters')
    
        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
            if options_dict == {}:
                options_dict = self._options_default
                self.report('INFO: using default options parameters')
            

        # set values, or defaults
        self.ctx.use_mpi = options_dict.get('use_mpi', self._options_default['use_mpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.walltime_sec = options_dict.get('max_wallclock_seconds', self._options_default['max_wallclock_seconds'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])
        self.ctx.options_params_dict = ParameterData(dict={'use_mpi': self.ctx.use_mpi, 'resources': self.ctx.resources, 
                                                           'max_wallclock_seconds': self.ctx.walltime_sec, 'queue_name': self.ctx.queue, 
                                                           'custom_scheduler_commands': self.ctx.custom_scheduler_commands})
        
        # set workflow parameters for the KKR imputrity calculations
        self.ctx.ef_shift = wf_dict.get('ef_shift', self._wf_default['ef_shift'])
        self.ctx.dos_params_dict = wf_dict.get('dos_params', self._wf_default['dos_params'])
        
        # set workflow parameters for the KKR impurity calculation
        self.ctx.nsteps = 1
        self.ctx.kkr_runmax = 1
        self.ctx.non_spherical = wf_dict.get('non_spherical', self._wf_default['non_spherical'])
        self.ctx.spinorbit = wf_dict.get('spinorbit', self._wf_default['spinorbit'])
        self.ctx.newsol = wf_dict.get('newsol', self._wf_default['newsol'])
        #self.ctx.kkrimp_params_dict = ParameterData(dict={#'nspin': self.ctx.nspin, 
        #                                                  'nsteps': self.ctx.nsteps, 
        #                                                  'kkr_runmax': self.ctx.kkr_runmax, 'non_spherical': self.ctx.non_spherical, 
        #                                                  'spinorbit': self.ctx.spinorbit, 'newsol': self.ctx.newsol,
        #                                                  'dosrun': True})

        # set workflow label and description
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)
    

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
        Validate input and catch possible errors
        """

        inputs = self.inputs
        inputs_ok = True

        if 'host_imp_pot' in inputs:
            if inputs.host_imp_pot.has_parents == False:
                self.report('WARNING: startpot has no parent and can not find '
                            'a converged host RemoteData node')
                inputs_ok = False
            else:
                self.report('INFO: get converged host RemoteData node and '
                            'impurity_info node from database')
                self.ctx.kkr_imp_wf = inputs.host_imp_pot.created_by.called_by
                self.report('INFO: found underlying kkr impurity workflow '
                            '(pk: {})'.format(self.ctx.kkr_imp_wf.pk))
                self.ctx.imp_info = self.ctx.kkr_imp_wf.inp.impurity_info
                self.report('INFO: found impurity_info node (pk: {})'.format(
                            self.ctx.imp_info.pk))
#                try:
#                    self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inp.gf_remote
#                except:
#                    self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inp.remote_converged_host
                try:
                    self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inp.remote_data_gf.inp.remote_folder.inp.parent_calc_folder.inp.remote_folder.out.remote_folder
                    self.report('INFO: imported converged_host_remote (pk: {}) and '
                                'impurity_info from database'.format(self.ctx.conv_host_remote.pk))
                except:
                    self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inp.gf_remote.inp.remote_folder.inp.parent_calc_folder.inp.remote_folder.out.remote_folder
                    self.report('INFO: imported converged_host_remote (pk: {}) and '
                                'impurity_info from database'.format(self.ctx.conv_host_remote.pk))
            
        self.report('INFO: validated input successfully: {}'.format(inputs_ok))  
        
        return inputs_ok
        
 
    def run_gfstep(self):
        """
        Start GF writeout step with DOS energy contour
        """
        
        options = self.ctx.options_params_dict
        kkrcode = self.inputs.kkr
        converged_host_remote = self.ctx.conv_host_remote
        imp_info = self.ctx.imp_info
        
        wf_params_gf = ParameterData(dict={'ef_shift':self.ctx.ef_shift, 'dos_run':True,
                                           'dos_params':self.ctx.dos_params_dict})
        label_gf = 'GF writeout for imp DOS'
        description_gf = 'GF writeout step with energy contour for impurity DOS'
        
        future = self.submit(kkr_flex_wc, label=label_gf, description=description_gf, 
                             kkr=kkrcode, options=options, 
                             remote_data=converged_host_remote, impurity_info=imp_info,
                             wf_parameters=wf_params_gf)
        
        self.report('INFO: running GF writeout (pid: {})'.format(future.pk))
        
        return ToContext(gf_writeout=future)
        
        
    def run_imp_dos(self):
        """
        Use previous GF step to calculate DOS for the impurity problem
        """

        options = self.ctx.options_params_dict
        kkrimpcode = self.inputs.kkrimp
        gf_writeout_wf = self.ctx.gf_writeout
        gf_writeout_calc = load_node(self.ctx.gf_writeout.out.workflow_info.get_attr('pk_flexcalc'))
        gf_writeout_remote = gf_writeout_wf.out.GF_host_remote
        impurity_pot = self.inputs.host_imp_pot
        imps = self.ctx.imp_info
        
        nspin = gf_writeout_calc.out.output_parameters.get_attr('nspin')
        self.report('nspin: {}'.format(nspin))
        self.ctx.kkrimp_params_dict = ParameterData(dict={'nspin': nspin, 
                                                          'nsteps': self.ctx.nsteps, 
                                                          'kkr_runmax': self.ctx.kkr_runmax, 'non_spherical': self.ctx.non_spherical, 
                                                          'spinorbit': self.ctx.spinorbit, 'newsol': self.ctx.newsol,
                                                          'dos_run': True})
        kkrimp_params = self.ctx.kkrimp_params_dict
        
        label_imp = 'KKRimp DOS (GF: {}, imp_pot: {}, Zimp: {}, ilayer_cent: {})'.format(
                    gf_writeout_wf.pk, impurity_pot.pk, imps.get_attr('Zimp'), imps.get_attr('ilayer_center'))
        description_imp = 'KKRimp DOS run (GF: {}, imp_pot: {}, Zimp: {}, ilayer_cent: {}, R_cut: {})'.format(
                    gf_writeout_wf.pk, impurity_pot.pk, imps.get_attr('Zimp'), imps.get_attr('ilayer_center'),
                    imps.get_attr('Rcut'))   
            
        future = self.submit(kkr_imp_sub_wc, label=label_imp, description=description_imp, 
                             kkrimp=kkrimpcode, options=options, 
                             wf_parameters=kkrimp_params, remote_data_gf=gf_writeout_remote, 
                             host_imp_startpot=impurity_pot)

        self.report('INFO: running DOS step for impurity system (pid: {})'.format(future.pk))
        
        return ToContext(kkrimp_dos=future)
        
    
    def return_results(self):
        """
        Return the results and create all of the output nodes
        """
        
        self.report('INFO: creating output nodes for the KKR imp DOS workflow ...')

        last_calc_pk = self.ctx.kkrimp_dos.out.workflow_info.get_attr('last_calc_nodeinfo')['pk']
        last_calc_output_params = load_node(last_calc_pk).out.output_parameters        
        last_calc_info = self.ctx.kkrimp_dos.out.workflow_info
        
        outputnode_dict = {}
        outputnode_dict['impurity_info'] = self.ctx.imp_info.get_attrs()
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['used_subworkflows'] = {'gf_writeout': self.ctx.gf_writeout.pk, 
                                                'impurity_dos': self.ctx.kkrimp_dos.pk}  
        outputnode_t = ParameterData(dict=outputnode_dict)
        outputnode_t.label = 'kkr_imp_dos_wc_inform'
        outputnode_t.description = 'Contains information for workflow'
        
        self.report('INFO: workflow_info node: {}'.format(outputnode_t.get_attrs()))

        self.out('workflow_info', outputnode_t)
        self.out('last_calc_output_parameters', last_calc_output_params)
        self.out('last_calc_info', last_calc_info)
        
        self.report('INFO: created output nodes for KKR imp DOS workflow.')      
        self.report('\n'
                    '|------------------------------------------------------------------------------------------------------------------|\n'
                    '|-------------------------------------| Done with the KKR imp DOS workflow! |--------------------------------------|\n'
                    '|------------------------------------------------------------------------------------------------------------------|')

            