#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory
from aiida.work.workchain import WorkChain, while_, if_, ToContext
from aiida.work.run import submit
from aiida.work import workfunction as wf
from aiida.work.process_registry import ProcessRegistry
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.workflows.dos import kkr_dos_wc
from aiida_kkr.tools.common_workfunctions import (test_and_get_codenode, 
                                                  update_params_wf, 
                                                  get_inputs_voronoi)
from aiida.common.datastructures import calc_states



__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.2"
__contributors__ = u"Philipp Rüßmann"


StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()
VoronoiProcess = VoronoiCalculation.process()

class kkr_startpot_wc(WorkChain):
    """
    Workchain  create starting potential for a KKR calculation by running 
    voronoi and getting the starting DOS for first checks on the validity of the input setting. 
    Starts from a structure together with a KKR parameter node.

    :param wf_parameters: (ParameterData), Workchain specifications
    :param structure: (StructureData), aiida structure node to begin 
        calculation from (needs to contain vacancies, if KKR needs empty spheres)
    :param kkr: (Code)
    :param voronoi: (Code)

    :return result_kkr_startpot_wc: (ParameterData), Information of workflow results 
        like Success, last result node, dos array data
    """

    _workflowversion = __version__
    _wf_default = {'queue_name' : '',                        # Queue name to submit jobs too
                   'resources': {"num_machines": 1},         # resources to allowcate for the job
                   'walltime_sec' : 60*60,                   # walltime after which the job gets killed (gets parsed to KKR)
                   'use_mpi' : False,                        # execute KKR with mpi or without
                   'custom_scheduler_commands' : '',         # some additional scheduler commands 
                   'dos_params' : {"nepts": 61,              # DOS params: number of points in contour
                                   "tempr": 200,             # DOS params: temperature
                                   "emin": -1,               # DOS params: start of energy contour
                                   "emax": 1,                # DOS params: end of energy contour
                                   "kmesh": [50, 50, 50]},   # DOS params: kmesh for DOS calculation (typically higher than in scf contour)
                   'num_rerun' : 3,                          # number of times voronoi+starting dos+checks is rerun to ensure non-negative DOS etc
                   'fac_cls_increase' : 1.3,                 # factor by which the screening cluster is increased each iteration (up to num_rerun times)
                   'r_cls' : 1.3,                            # default cluster radius, is increased iteratively
                   'natom_in_cls_min' : 79                   # minimum number of atoms in screening cluster
                }
                   
    _kkr_default_params = kkrparams.get_KKRcalc_parameter_defaults()

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
        super(kkr_startpot_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict=cls._wf_default))
        spec.input("structure", valid_type=StructureData, required=True)
        spec.input("kkr", valid_type=Code, required=True)
        spec.input("voronoi", valid_type=Code, required=True)
        spec.input("calc_parameters", valid_type=ParameterData, required=False)

        # Here the structure of the workflow is defined
        spec.outline(
            # initialize workflow and check inputs
            cls.start,
            # check if another iteration is done (in case of either voro_ok, doscheck_ok is False)
            while_(cls.do_iteration_check)(
                    # run voronoi calculation
                    cls.run_voronoi,
                    # check voronoi output (also sets ctx.voro_ok)
                    if_(cls.check_voronoi)(
                            # create starting DOS using dos sub-workflow
                            cls.get_dos,
                            # perform some checks and set ctx.doscheck_ok accordingly
                            cls.check_dos
                        )
                    ),
            # collect results and return
            cls.return_results
        )

       
    def start(self):
        """
        init context and some parameters
        """
        self.report('INFO: started VoroStart workflow version {}\n'
                    'INFO: Workchain node identifiers: {}'
                    ''.format(self._workflowversion, ProcessRegistry().current_calc_node))

        ####### init    #######

        # internal para /control para
        self.ctx.abort = False

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()

        #TODO: check for completeness
        if wf_dict == {}:
            wf_dict = self._wf_default
            self.report('INFO: using default wf parameter')

        # set values, or defaults
        self.ctx.use_mpi = wf_dict.get('use_mpi', self._wf_default['use_mpi'])
        self.ctx.resources = wf_dict.get('resources', self._wf_default['resources'])
        self.ctx.walltime_sec = wf_dict.get('walltime_sec', self._wf_default['walltime_sec'])
        self.ctx.queue = wf_dict.get('queue_name', self._wf_default['queue_name'])
        self.ctx.custom_scheduler_commands = wf_dict.get('custom_scheduler_commands', self._wf_default['custom_scheduler_commands'])
        
        self.ctx.dos_params_dict = wf_dict.get('dos_params', self._wf_default['dos_params'])
        
        self.ctx.description_wf = self.inputs.get('_description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('_label', self._wf_label)
        
        # iterative rerunning parameters
        self.ctx.Nrerun = self.inputs.get('num_rerun', self._wf_default['num_rerun'])
        self.ctx.fac_cls_increase = self.inputs.get('fac_cls_increase', self._wf_default['fac_cls_increase'])
        
        # initialize checking booleans
        self.ctx.is_starting_iter = True
        self.ctx.doscheck_ok = False
        self.ctx.voro_ok = False
        
        # some physical parameters that are reused
        self.ctx.r_cls = self.inputs.get('r_cls', self._wf_default['r_cls'])
        self.ctx.nclsmin = self.inputs.get('natom_in_cls_min', self._wf_default['natom_in_cls_min'])
        self.ctx.fac_clsincrease = self.inputs.get('fac_cls_increase', self._wf_default['fac_cls_increase'])
        
        
        #TODO add missing info
        # print the inputs
        self.report('INFO: use the following parameter:\n'
                    'use_mpi: {}\n'
                    'Resources: {}\n'
                    'Walltime (s): {}\n'
                    'queue name: {}\n'
                    'scheduler command: {}\n'
                    'description: {}\n'
                    'label: {}\n'
                    'dos_params: {}\n'.format(self.ctx.use_mpi, self.ctx.resources, self.ctx.walltime_sec, 
                                              self.ctx.queue, self.ctx.custom_scheduler_commands, 
                                              self.ctx.description_wf, self.ctx.label_wf, 
                                              self.ctx.dos_params_dict))

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''
        
        # get kkr and voronoi codes from input
        try:
            test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
        except ValueError:
            error = ("The code you provided for kkr does not "
                     "use the plugin kkr.kkr")
            self.ctx.errors.append(error)
            self.control_end_wc(error)
        try:
            test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
        except ValueError:
            error = ("The code you provided for voronoi does not "
                     "use the plugin kkr.voro")
            self.ctx.errors.append(error)
            self.control_end_wc(error)
            
            
       
    def run_voronoi(self):
        """
        run voronoi calculation with parameters from input
        """
        # incerement iteration counter
        self.ctx.iter += 1
        
        # increase some parameters
        if self.ctx.iter > 1:
            self.ctx.r_cls = self.ctx.r_cls * self.ctx.fac_clsincrease
            
        #
    
        structure = self.inputs.structure
        self.ctx.formula = structure.get_formula()
        label = 'voronoi calculation'
        description = '{} vornoi on {}'.format(self.ctx.description_wf, self.ctx.formula)

        voronoicode = self.inputs.voronoi
        
        # get valid KKR parameters
        if self.ctx.iter > 1:
            # take value from last run to continue
            params = self.ctx.last_params
            first_iter = False
        else:
            # used input or defaults in first iteration
            first_iter = True
            if 'calc_parameters' in self.inputs:
                params = self.inputs.calc_parameters
            else:
                kkrparams_default = kkrparams()
                para_version = self._kkr_default_params[0]
                for key, val in self._kkr_default_params[1].iteritems():
                    kkrparams_default.set_value(key, val, silent=True)
                # create ParameterData node
                params = ParameterData(dict={kkrparams_default})
                params.label = 'Defaults for KKR parameter node'
                params.description = 'defaults as defined in kkrparams of version {}'.format(para_version)
            #  set last_params accordingly (used below for provenance tracking)
            self.ctx.last_params = params
            
        # check if RCLUSTZ is set and use setting from wf_parameters instead (calls update_params_wf to keep track of provenance)
        updated_params = False
        update_list = []
        set_vals = params.get_dict().get_set_values()
        set_vals = [keyvalpair[0] for keyvalpair in set_vals]
        if 'RCLUSTZ' in set_vals:
            rcls_input = params.get_dict().get_value('RCLUSTZ')
            # set r_cls by default or from input in first iteration
            if self.ctx.r_cls < rcls_input and first_iter:
                self.ctx.r_cls = rcls_input
                updated_params = True
                update_list.append('RCLUSTZ')
            elif self.ctx.r_cls > rcls_input:
                # change rcls with iterations
                updated_params = True
                update_list.append('RCLUSTZ')
        else:
            updated_params = True
            update_list.append('RCLUSTZ')
        
        # store updated nodes
        if updated_params:
            updatenode = ParameterData(dict=params.get_dict())
            updatenode.label = 'updated params node'
            updatenode.description = 'changed values: {}'.format(update_list)
            # used workfunction for provenance tracking if parameters have been changed
            params = update_params_wf(self.ctx.last_params, updatenode)
            

        options = {"max_wallclock_seconds": self.ctx.walltime_sec,
                   "resources": self.ctx.resources,
                   "queue_name" : self.ctx.queue}

        inputs = get_inputs_voronoi(structure, voronoicode, options, label, description, params=params)
        self.report('INFO: run voronoi step')
        future = submit(VoronoiProcess, **inputs)
        
        
        # return remote_voro (passed to dos calculation as input)
        return ToContext(voro_calc=future, last_params = params)
    
       
    def check_voronoi(self):
        """
        check voronoi output. return True/False if voronoi output is ok/problematic
        if output is problematic try to increase some parameters (e.g. cluster radius) and rerun up tp N_rerun_max times
        initializes with returning True
        """
        
        #do some checks with the voronoi output (finally sets self.ctx.voro_ok)
        self.ctx.voro_ok = True
        
        # check calculation state (calculation must be completed)
        calc_state = self.ctx.voro_calc.get_state()
        if calc_state != calc_states.FINISHED:
            self.ctx.voro_ok = False
        
        #TODO check self.ctx.nclsmin condition
        voro_out_cls_info = res.cluster_info_group
        ncls = voro_out_cls_info.pop('number_of_clusters')
        
        if self.ctx.nclsmin > nclsmin_last_calc:
            self.ctx.voro_ok = False
        
        #TODO check radii condition
        voro_out_radii = self.ctx.voro_calc.res.radii
        
        #TODO other checks?
        
        # finally return result of check
        return self.ctx.voro_ok
    
    
    def do_iteration_check(self):
        """
        check if another iteration should be done
        """
        if self.ctx.is_starting_iter: 
            # initial iteration (at least one has to be done)
            # reset starting iter flag
            self.ctx.is_starting_iter = False
            return True
        elif self.ctx.iter >= self.ctx.Nrerun: 
            # check if maximal number of iterations is reached
            return False
        elif self.ctx.voro_ok and self.ctx.doscheck_ok:
            # if both steps succeed we are done
            return False
        else:
            return True
        
        
    def get_dos(self):
        """
        call to dos sub workflow passing the appropriate input and submitting the calculation
        """
        # take subset of input and prepare parameter node for dos workflow
        wfdospara_dict = {'queue_name' : self.ctx.queue, 
                          'resources': self.ctx.resources,
                          'walltime_sec' : self.ctx.walltime_sec,
                          'use_mpi' : self.ctx.use_mpi,
                          'custom_scheduler_commands' : self.ctx.custom_scheduler_commands,
                          'dos_params' : self.ctx.dos_params_dict}
        wfdospara_node = ParameterData(dict=wfdospara_dict)
        
        code = self.inputs.kkr
        remote = self.ctx.voro_calc.out.remote_folder
        future = submit(kkr_dos_wc, kkr=code, remote_data=remote, wf_parameters=wfdospara_node)
        
        return ToContext(doscal=future)
        
        
    def check_dos(self):
        """
        checks if dos of starting potential is ok
        """
        dos_ok = False
        #TODO implement checks for negative DOS, starting EMIN, semi-core-states
        
        # finally set the value in context (needed in do_iteration_check)
        if dos_ok:
            self.ctx.doscheck_ok = True
        else:
            self.ctx.doscheck_ok = False
        
        
    def control_end_wc(self, errormsg):
        """
        Controled way to shutdown the workchain. will initalize the output nodes
        """
        self.report('ERROR: shutting workchain down in a controlled way.')
        self.ctx.successful = False
        self.ctx.abort = True
        self.report(errormsg) # because return_results still fails somewhen
        self.return_results()
        self.abort(errormsg)
        
        
    def return_results(self):
        """
        return the results of the dos calculations
        This should run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['list_of_errors'] = self.ctx.errors
        outputnode_dict['use_mpi'] = self.ctx.use_mpi
        outputnode_dict['resources'] = self.ctx.resources
        outputnode_dict['walltime_sec'] = self.ctx.walltime_sec
        outputnode_dict['queue'] = self.ctx.queue
        outputnode_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        outputnode_dict['dos_params'] = self.ctx.dos_params_dict
        outputnode_dict['nspin'] = self.ctx.dosrun.res.nspin
        
        outputnode = ParameterData(dict=outputnode_dict)
        outputnode.label = 'kkr_scf_wc_results'
        outputnode.description = ''
        outputnode.store()
        
        self.report("INFO: create dos results nodes. outputnode={}, dos calc retrieved node={}".format(outputnode, self.ctx.dosrun.out.retrieved))
        self.report("INFO: type outputnode={}".format(type(outputnode)))
        
        # interpol dos file and store to XyData nodes
        outdict = create_dos_result_node(outputnode, self.ctx.dosrun.out.retrieved)
    
        
        for link_name, node in outdict.iteritems():
            self.report("INFO: storing node '{}' with link name '{}'".format(node, link_name))
            self.report("INFO: node type: {}".format(type(node)))
            self.out(link_name, node)
        