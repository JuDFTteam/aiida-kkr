# -*- coding: utf-8 -*-
"""
In this module you find the total workflow for a kkr impurity calculation  
and some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory, load_node
from aiida.work.workchain import WorkChain, ToContext
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, neworder_potential_wf
from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
import numpy as np

__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = u"Fabian Bertoldo"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
SinglefileData = DataFactory('singlefile')
FolderData = DataFactory('folder')



class kkr_imp_wc(WorkChain):
    """
    Workchain of a kkrimp calculation starting either from scratch (with a structure
    and impurity_info node), or with a converged host potential and impurity 
    startpotentials, ... to calculate the converged host-impurity potential of the system.

    :param options_parameters: (ParameterData), Workchain specifications
    :param wf_parameters: (ParameterData), specifications for the kkr impurity workflow
    :param voro_aux_parameters: (ParameterData), specification for the auxiliary voronoi calculation for the impurity
    :param kkrimpcode: (Code), mandatory: KKRimp code converging the host-imp-potential
    :param kkrcode: (Code), mandatory: KKR code for calculation the host potential
    :param vorocode: (Code), mandatory: Voronoi code to generate the impurity startpot
    :param GF_remote_data: (RemoteData): remote folder of a previous kkrflex calculation containing the flexfiles ...

    :return result_kkr_imp_wc: (ParameterData), Information of workflow results
    """
    
    
    _workflowversion = __version__
    _wf_label = 'kkr_imp_wc'
    _wf_description = 'Workflow for a KKRimp calculation'


    _options_default = {'queue_name' : '',                                          # Queue name to submit jobs too
                        'resources': {"num_machines": 1},                           # resources to allowcate for the job
                        'walltime_sec' : 60*60,                                     # walltime after which the job gets killed (gets parsed to KKR)}
                        'custom_scheduler_commands' : '',                           # some additional scheduler commands 
                        'use_mpi' : False}                                          # execute KKR with mpi or without
                        
    _wf_default = {'nspin':  1,                                                   # non-magnetic calculation, set nspin = 2 for magnetic case
                   'kkr_runmax': 3,                         # Maximum number of kkr jobs/starts (defauld iterations per start)
                   'threshold_aggressive_mixing': 5*10**-2, # threshold after which agressive mixing is used
                   'convergence_criterion' : 3*10**-2,      # Stop if charge denisty is converged below this value
                   'mixreduce': 0.5,                        # reduce mixing factor by this factor if calculaito fails due to too large mixing
                   'strmix': 0.03,                          # mixing factor of simple mixing                           
                   'aggressive_mix': 3,                     # type of aggressive mixing (3: broyden's 1st, 4: broyden's 2nd, 5: generalized anderson)
                   'aggrmix': 0.01,                         # mixing factor of aggressive mixing
                   'nsteps': 20,                            # number of iterations done per KKR calculation
                   'non_spherical': 1,                      # use non-spherical parts of the potential (0 if you don't want that)
                   'broyden_number': 20,                    # number of potentials to 'remember' for Broyden's mixing
                   'born_iter': 2,                          # number of Born iterations for the non-spherical calculation
                   'mag_init' : False,                      # initialize and converge magnetic calculation
                   'hfield' : [0.1, 10], # Ry               # external magnetic field used in initialization step
                   'init_pos' : None,                       # position in unit cell where magnetic field is applied [default (None) means apply to all]
                   'r_cls' : 1.3, # alat                    # default cluster radius, is increased iteratively
                   'calc_orbmom' : False,                   # defines of orbital moments will be calculated and written out
                   'spinorbit' : False,                     # SOC calculation (True/False)
                   'newsol' : False }                       # new SOC solver is applied                   
    
    _voro_aux_default = {'dos_params' : {'nepts': 61,                           # DOS params: number of points in contour
                                         'tempr': 200, # K                      # DOS params: temperature
                                         'emin': -1, # Ry                       # DOS params: start of energy contour
                                         'emax': 1,  # Ry                       # DOS params: end of energy contour
                                         'kmesh': [50, 50, 50]},                # DOS params: kmesh for DOS calculation (typically higher than in scf contour)
                        'num_rerun' : 4,                           # number of times voronoi+starting dos+checks is rerun to ensure non-negative DOS etc
                        'fac_cls_increase' : 1.3, # alat           # factor by which the screening cluster is increased each iteration (up to num_rerun times)
                        'r_cls' : 1.3,            # alat           # default cluster radius, is increased iteratively
                        'natom_in_cls_min' : 79,                   # minimum number of atoms in screening cluster
                        'delta_e_min' : 1., # eV                   # minimal distance in DOS contour to emin and emax in eV
                        'threshold_dos_zero' : 10**-3,             #states/eV  
                        'check_dos': True,                         # logical to determine if DOS is computed and checked
                        'delta_e_min_core_states' : 1.0,  # Ry     # minimal distance of start of energy contour to highest lying core state in Ry
                        'lmax': 3,
                        'gmax': 65.,
                        'rmax': 7.,
                        'rclustz': 2.5}
                                          
                                 
                                          
    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create 
        set of wf_parameters.
        
        returns _wf_defaults
        """
    
        print('Version of workflow: {}'.format(self._workflowversion))
        return self._options_default, self._wf_default, self._voro_aux_default                                          
                                                                  

        
    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """
    
        super(kkr_imp_wc, cls).define(spec)  
        
        # define the inputs of the workflow
        spec.input("vorocode", valid_type=Code, required=True) 
        spec.input("kkrcode", valid_type=Code, required=True)
        spec.input("kkrimpcode", valid_type=Code, required=True)
        spec.input("remote_converged_host", valid_type=RemoteData, required=True)
        spec.input("kkrflex_remote", valid_type=RemoteData, required=False)
        spec.input("impurity_info", valid_type=ParameterData, required=True)
        spec.input("options_parameters", valid_type=ParameterData, required=False)
        spec.input("voro_aux_parameters", valid_type=ParameterData, required=False)
        spec.input("wf_parameters", valid_type=ParameterData, required=False)
        
        
        # structure of the workflow
        spec.outline(
            cls.start,                                                          # initialize workflow
            cls.validate_input,                                                 # validate the input
            cls.run_gf_writeout,                                                # write out the host GF
            cls.run_voroaux,                                                    # calculate the auxiliary impurity potentials
            cls.construct_startpot,                                             # construct the host-impurity startpotential
            cls.run_kkrimp_scf,                                                 # run the kkrimp_sub workflow to converge the host-imp startpot
            cls.return_results)                                                 # check if the calculation was successful and return the result nodes
        
        
        # define the possible exit codes
        spec.exit_code(141, 'ERROR_INVALID_INPUT_CODE', 
            message="ERROR: one or more of the codes you provided do not "
                    "use the necessary plugins: kkr.voro, kkr.kkr, kkr.kkrimp")
        
    
        # define the outputs of the workflow
        spec.output('workflow_info', valid_type=ParameterData)
        spec.output('hostimp_remote_converged', valid_type=RemoteData)
        spec.output('hostimp_output_pot', valid_type=SinglefileData)
        
    
    
    
    def start(self):
        """
        Init context and some parameters
        """
        
        self.report('INFO: started KKR impurity workflow version {}'
                    ''.format(self._workflowversion))
        
     
        # get input parameters
        if 'options_parameters' in self.inputs:
            options_dict = self.inputs.options_parameters.get_dict() 
        else:
            options_dict = self._options_default
            self.report('INFO: using default options')
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
        else:
            wf_dict = self._wf_default
            self.report('INFO: using default workflow parameters')
        if 'voro_aux_parameters' in self.inputs:
            voro_aux_dict = self.inputs.voro_aux_parameters.get_dict()
        else:
            voro_aux_dict = self._voro_aux_default
            self.report('INFO: using default parameters auxiliary voronoi calculation')
        
 
        # set option parameters from input, or defaults
        self.ctx.use_mpi = options_dict.get('use_mpi', self._options_default['use_mpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.walltime_sec = options_dict.get('walltime_sec', self._options_default['walltime_sec'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])
        self.ctx.options_params_dict = ParameterData(dict={'use_mpi': self.ctx.use_mpi, 'resources': self.ctx.resources, 'walltime_sec': self.ctx.walltime_sec,
                                                           'queue_name': self.ctx.queue, 'custom_scheduler_commands': self.ctx.custom_scheduler_commands})

        # set label and description of the workflow
        self.ctx.description_wf = self.inputs.get('description', 'Workflow for a KKR impurity calculation starting from a host-impurity potential')
        self.ctx.label_wf = self.inputs.get('label', 'kkr_imp_sub_wc')            
            
        # set parameters for the auxiliary voronoi calculation
        self.ctx.voro_dos_params = voro_aux_dict.get('dos_params', self._voro_aux_default['dos_params'])
        self.ctx.voro_num_rerun = voro_aux_dict.get('num_rerun', self._voro_aux_default['num_rerun'])
        self.ctx.voro_fac_cls_increase = voro_aux_dict.get('fac_cls_increase', self._voro_aux_default['fac_cls_increase'])
        self.ctx.voro_r_cls = voro_aux_dict.get('r_cls', self._voro_aux_default['r_cls'])
        self.ctx.voro_natom_in_cls_min = voro_aux_dict.get('natom_in_cls_min', self._voro_aux_default['natom_in_cls_min'])
        self.ctx.voro_delta_e_min = voro_aux_dict.get('delta_e_min', self._voro_aux_default['delta_e_min'])
        self.ctx.voro_threshold_dos_zero = voro_aux_dict.get('threshold_dos_zero', self._voro_aux_default['threshold_dos_zero'])
        self.ctx.voro_check_dos = voro_aux_dict.get('check_dos', self._voro_aux_default['check_dos'])
        self.ctx.voro_delta_e_min_core_states = voro_aux_dict.get('delta_e_min_core_states', self._voro_aux_default['delta_e_min_core_states'])
        self.ctx.voro_lmax = voro_aux_dict.get('lmax', self._voro_aux_default['lmax'])
        self.ctx.voro_gmax = voro_aux_dict.get('gmax', self._voro_aux_default['gmax'])
        self.ctx.voro_rmax = voro_aux_dict.get('rmax', self._voro_aux_default['rmax'])
        self.ctx.voro_rclustz = voro_aux_dict.get('rclustz', self._voro_aux_default['rclustz'])
        # set up new parameter dict to pass to voronoi subworkflow later
        self.ctx.voro_params_dict = ParameterData(dict={'queue_name': self.ctx.queue, 'resources': self.ctx.resources, 'walltime_sec': self.ctx.walltime_sec, 
                                                        'use_mpi': self.ctx.use_mpi, 'custom_scheduler_commands': self.ctx.custom_scheduler_commands,
                                                        'dos_params': self.ctx.voro_dos_params, 'num_rerun': self.ctx.voro_num_rerun, 
                                                        'fac_cls_increase': self.ctx.voro_fac_cls_increase, 'r_cls': self.ctx.voro_r_cls,
                                                        'natom_in_cls_min': self.ctx.voro_natom_in_cls_min, 'delta_e_min': self.ctx.voro_delta_e_min,
                                                        'threshold_dos_zero': self.ctx.voro_threshold_dos_zero, 'check_dos': self.ctx.voro_check_dos,
                                                        'delta_e_min_core_states': self.ctx.voro_delta_e_min_core_states})
        
        # set workflow parameters for the KKR impurity calculation
        self.ctx.nspin = wf_dict.get('nspin', self._wf_default['nspin'])
        self.ctx.nsteps = wf_dict.get('nsteps', self._wf_default['nsteps'])
        self.ctx.kkr_runmax = wf_dict.get('kkr_runmax', self._wf_default['kkr_runmax'])
        self.ctx.threshold_aggressive_mixing = wf_dict.get('threshold_aggressive_mixing', self._wf_default['threshold_aggressive_mixing'])
        self.ctx.convergence_criterion = wf_dict.get('convergence_criterion', self._wf_default['convergence_criterion'])
        self.ctx.mixreduce = wf_dict.get('mixreduce', self._wf_default['mixreduce'])
        self.ctx.strmix = wf_dict.get('strmix', self._wf_default['strmix'])
        self.ctx.aggressive_mix = wf_dict.get('aggressive_mix', self._wf_default['aggressive_mix'])
        self.ctx.aggrmix = wf_dict.get('aggrmix', self._wf_default['aggrmix'])
        self.ctx.non_spherical = wf_dict.get('non_spherical', self._wf_default['non_spherical'])
        self.ctx.broyden_number = wf_dict.get('broyden_number', self._wf_default['broyden_number'])
        self.ctx.born_iter = wf_dict.get('born_iter', self._wf_default['born_iter'])
        self.ctx.mag_init = wf_dict.get('mag_init', self._wf_default['mag_init'])
        self.ctx.hfield = wf_dict.get('hfield', self._wf_default['hfield'])
        self.ctx.init_pos = wf_dict.get('init_pos', self._wf_default['init_pos'])
        self.ctx.r_cls = wf_dict.get('r_cls', self._wf_default['r_cls'])
        self.ctx.calc_orbmom = wf_dict.get('calc_orbmom', self._wf_default['calc_orbmom'])
        self.ctx.spinorbit = wf_dict.get('spinorbit', self._wf_default['spinorbit'])
        self.ctx.newsol = wf_dict.get('newsol', self._wf_default['newsol'])
        # set up new parameter dict to pass to kkrimp subworkflow later
        self.ctx.kkrimp_params_dict = ParameterData(dict={'nspin': self.ctx.nspin, 'nsteps': self.ctx.nsteps, 'kkr_runmax': self.ctx.kkr_runmax,
                                                          'threshold_aggressive_mixing': self.ctx.threshold_aggressive_mixing,
                                                          'convergence_criterion': self.ctx.convergence_criterion, 'mixreduce': self.ctx.mixreduce,
                                                          'strmix': self.ctx.strmix, 'aggressive_mix': self.ctx.aggressive_mix,
                                                          'aggrmix': self.ctx.aggrmix, 'non_spherical': self.ctx.non_spherical,
                                                          'broyden_number': self.ctx.broyden_number, 'born_iter': self.ctx.born_iter,
                                                          'mag_init': self.ctx.mag_init, 'hfield': self.ctx.hfield, 'init_pos': self.ctx.init_pos,
                                                          'r_cls': self.ctx.r_cls, 'calc_orbmom': self.ctx.calc_orbmom, 
                                                          'spinorbit': self.ctx.spinorbit, 'newsol': self.ctx.newsol})
        
        
        # report the chosen parameters to the user
        self.report('INFO: use the following parameter:\n'
                    '\nGeneral settings\n'
                    'use mpi: {}\n'
                    'resources: {}\n'
                    'walltime (s): {}\n'
                    'queue name: {}\n'
                    'scheduler command: {}\n'
                    'description: {}\n'
                    'label: {}\n'
                    'nspin: {}\n'
                    'parameters for the voroaux calculation: {}\n'
                    'parameters for the kkrimp scf: {}\n'
                    ''.format(self.ctx.use_mpi, self.ctx.resources, self.ctx.walltime_sec,
                              self.ctx.queue, self.ctx.custom_scheduler_commands,
                              self.ctx.description_wf, self.ctx.label_wf,
                              self.ctx.nspin, self.ctx.voro_params_dict.get_attrs(),
                              self.ctx.kkrimp_params_dict.get_attrs()))
            
            
            
    def validate_input(self):
        """
        Validate the input and catch possible errors from the input
        """

        inputs = self.inputs
        inputs_ok = True
        
        if 'kkrcode' and 'kkrimpcode' and 'vorocode' in inputs:
            try:
                test_and_get_codenode(inputs.kkrcode, 'kkr.kkr', use_exceptions=True)
                test_and_get_codenode(inputs.kkrimpcode, 'kkr.kkrimp', use_exceptions=True)
                test_and_get_codenode(inputs.vorocode, 'kkr.voro', use_exceptions=True)
            except ValueError:
                inputs_ok = False
                self.report(self.exit_codes.ERROR_INVALID_INPUT_CODE)
                return self.exit_codes.ERROR_INVALID_INPUT_CODE  
                
        if 'remote_converged_host' in inputs:
            self.report('INFO: found remote_data (pid: {}) in input'.format(inputs.remote_converged_host.pk))
            
        if 'impurity_info' in inputs:
            self.report('INFO: found impurity_info node in input: {}'.format(inputs.impurity_info.get_attrs()))
 
            
        self.report('INFO: validated input successfully: {}'.format(inputs_ok))
 
        
        
    def run_gf_writeout(self):
        """
        Run the gf_writeout workflow to calculate the host Green's function and the
        KKR flexfiles using the converged host remote folder and the impurity info node
        """
        
        # collect inputs
        kkrcode = self.inputs.kkrcode
        imp_info = self.inputs.impurity_info
        converged_host_remote = self.inputs.remote_converged_host
        options = self.ctx.options_params_dict
        
        # set label and description of the calc
        sub_label = 'GF writeout (conv. host pid: {}, imp_info pid: {})'.format(converged_host_remote.pk, imp_info.pk)
        sub_description = 'GF writeout sub workflow for kkrimp_wc using converged host remote data (pid: {}) and impurity_info node (pid: {})'.format(converged_host_remote.pk, imp_info.pk)
        
        self.report('INFO: run gf writeout workflow to get host GF and flexfiles')
        
        future = self.submit(kkr_flex_wc, label=sub_label, description=sub_description, kkr=kkrcode, options_parameters=options, 
                             remote_data=converged_host_remote, imp_info=imp_info)
        
        self.report('INFO: running gf writeout (pid: {})'.format(future.pk))
        
        return ToContext(gf_writeout=future, last_calc_gf=future)

      
    
    def run_voroaux(self):
        """
        Perform a voronoi calculation for every impurity charge using the structure
        from the converged KKR host calculation and return a list of auxiliary 
        voronoi potentials
        """
        
        # collect inputs
        vorocode = self.inputs.vorocode
        kkrcode = self.inputs.kkrcode
        imp_info = self.inputs.impurity_info
        voro_params = self.ctx.voro_params_dict
        converged_host_remote = self.inputs.remote_converged_host
        calc_params = ParameterData(dict=kkrparams(NSPIN=self.ctx.nspin, LMAX=self.ctx.voro_lmax, GMAX=self.ctx.voro_gmax, 
                                                   RMAX=self.ctx.voro_rmax, RCLUSTZ=self.ctx.voro_rclustz).get_dict())
        structure_host, voro_calc = VoronoiCalculation.find_parent_structure(converged_host_remote) 
        
        # for every impurity, generate a structure and launch the voronoi workflow
        # to get the auxiliary impurity startpotentials        
        self.ctx.voro_calcs = {}
#        for i in range(2):
        inter_struc = change_struc_imp_aux_wf(structure_host, imp_info)
        sub_label = 'voroaux calc for Zimp: {} in host-struc'.format(imp_info.get_attr('Zimp')[0])
        sub_description = 'Auxiliary voronoi calculation for an impurity with charge '
        sub_description += '{} in the host structure from pid: {}'.format(imp_info.get_attr('Zimp')[0], converged_host_remote.pk)
        
        self.report('INFO: run voroaux workflow for an impurity with Zimp= {}'.format(imp_info.get_attr('Zimp')[0]))
        
        future = self.submit(kkr_startpot_wc, label=sub_label, description=sub_description, structure=inter_struc,
                             voronoi=vorocode, kkr=kkrcode, wf_parameters=voro_params, calc_parameters=calc_params)
        
        tmp_calcname = 'voro_aux_{}'.format(1)    
        self.ctx.voro_calcs[tmp_calcname] = future    
        self.report('INFO: running voro aux (Zimp= {}, pid: {})'.format(imp_info.get_attr('Zimp')[0], future.pk))  
        
        return ToContext(last_voro_calc=future)   
        

                
    def construct_startpot(self):
        """
        Take the output of GF writeout and the converged host potential as well as the
        auxiliary startpotentials for the impurity to construct the startpotential for the
        KKR impurity sub workflow
        """
        
        nspin = self.ctx.nspin
        
        # collect all nodes necessary to construct the startpotential
        GF_host_calc_pk = self.ctx.gf_writeout.out.calculation_info.get_attr('pk_flexcalc')
        self.report('GF_host_calc_pk: {}'.format(GF_host_calc_pk))
        GF_host_calc = load_node(GF_host_calc_pk)
        converged_host_remote = self.inputs.remote_converged_host
        voro_calc_remote = self.ctx.last_voro_calc.out.last_voronoi_remote
        imp_info = self.inputs.impurity_info
        
        # prepare settings dict
        potname_converged = 'potential'
        potname_impvorostart = 'output.pot'
        potname_imp = 'potential_imp'
        if nspin < 2:
            replacelist_pot2 = [[0,0]]
        else:
            replacelist_pot2 = [[0,0],[1,1]]
        neworder_pot1 = [int(i) for i in np.loadtxt(GF_host_calc.out.retrieved.get_abs_path('scoef'), skiprows=1)[:,3]-1]    
        
        settings_label = 'startpot_KKRimp for imp_info node {}'.format(imp_info.pk)
        settings_description = 'starting potential for impurity info: {}'.format(imp_info)
        settings = ParameterData(dict={'pot1': potname_converged,  'out_pot': potname_imp, 'neworder': neworder_pot1,
                                       'pot2': potname_impvorostart, 'replace_newpos': replacelist_pot2, 'label': settings_label,
                                       'description': settings_description})
        startpot_kkrimp = neworder_potential_wf(settings_node=settings, parent_calc_folder=converged_host_remote,
                                                parent_calc_folder2=voro_calc_remote)
        
        # add starting potential for kkrimp calculation to context
        self.ctx.startpot_kkrimp = startpot_kkrimp
        
        

    def run_kkrimp_scf(self):
        """
        Uses both the previously generated host-impurity startpotential and the output from
        the GF writeout workflow as inputs to run the kkrimp_sub workflow in order to 
        converge the host-impurity potential
        """
        
        # collect all necessary input nodes
        kkrimpcode = self.inputs.kkrimpcode
        startpot = self.ctx.startpot_kkrimp
        gf_remote = self.ctx.gf_writeout.out.GF_host_remote
        kkrimp_params = self.ctx.kkrimp_params_dict
        options = self.ctx.options_params_dict
        
        # set label and description
        sub_label = 'kkrimp_sub scf wf (conv. host remote: {}, imp_info: {})'.format(self.inputs.remote_converged_host.pk, self.inputs.impurity_info.pk)
        sub_description = 'convergence of the host-impurity potential (pk: {}) using GF remote (pk: {})'.format(startpot.pk, gf_remote.pk)
        
        self.report('INFO: run kkrimp_sub_wf (startpot: {}, GF_remote: {})'.format(startpot.pk, gf_remote.pk))
        
        future = self.submit(kkr_imp_sub_wc, label=sub_label, description=sub_description, kkrimp=kkrimpcode, options_parameters=options, 
                             host_imp_startpot=startpot, GF_remote_data=gf_remote, wf_parameters=kkrimp_params)

        return ToContext(kkrimp_scf_sub=future)

        
        
    def return_results(self):
        """
        Return the results and create all of the output nodes
        """
        
        self.report('INFO: creating output nodes for the KKR impurity workflow ...')
        
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        #outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['used_subworkflows'] = {'gf_writeout': self.ctx.gf_writeout.pk, 'auxiliary_voronoi': self.ctx.last_voro_calc.pk,
                                                'kkr_imp_sub': self.ctx.kkrimp_scf_sub.pk}   
        outputnode_t = ParameterData(dict=outputnode_dict)
        outputnode_t.label = 'kkrimp_wc_inform'
        outputnode_t.description = 'Contains information for workflow'

        self.out('workflow_info', outputnode_t)
        #self.out('hostimp_output_pot', self.ctx.kkrimp_scf_sub.out.host_imp_pot)
        
        self.report('INFO: created output nodes successfully')        
        self.report('\n'
                    '|------------------------------------------------------------------------------------------------------------------|\n'
                    '|-------------------------------------| Done with the KKR impurity workflow! |-------------------------------------|\n'
                    '|------------------------------------------------------------------------------------------------------------------|')
            
        
        
def change_struc_imp_aux_wf(struc, imp_info): # Note: works for single imp at center only!
    from aiida.common.constants import elements as PeriodicTableElements
    _atomic_numbers = {data['symbol']: num for num, data in PeriodicTableElements.iteritems()}

    new_struc = StructureData(cell=struc.cell)
    isite = 0
    for site in struc.sites:
        sname = site.kind_name
        kind = struc.get_kind(sname)
        pos = site.position
        zatom = _atomic_numbers[kind.get_symbols_string()]
        if isite == imp_info.get_dict().get('ilayer_center'):
            zatom = imp_info.get_dict().get('Zimp')[0]
        symbol = PeriodicTableElements.get(zatom).get('symbol')
        new_struc.append_atom(position=pos, symbols=symbol)
        isite += 1

    return new_struc            
              