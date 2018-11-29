#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a EOS calculation and
some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory
from aiida.work.workchain import WorkChain, ToContext
from aiida.orm.data.base import Float
from aiida.work.workfunctions import workfunction as wf
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.common_workfunctions import update_params_wf
from aiida_kkr.tools.common_functions import get_Ry2eV
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
from masci_tools.io.kkr_params import kkrparams
from ase.eos import EquationOfState
from numpy import array, mean, std


__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = u"Philipp Rüßmann"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()
VoronoiProcess = VoronoiCalculation.process()


class kkr_eos_wc(WorkChain):
    """
    Workchain of an equation of states calculation with KKR.

    Layout of the workflow:
      1. determine V0, scale_range, etc. from input
      2. run voro_start for V0 and smallest volume
          2.1 get minimum for RMTCORE (needs to be fixed for all calculations to be able to compare total energies
      3. submit kkr_scf calculations for all volumes using RMTCORE setting determined in step 2
      4. collect results
    """

    _workflowversion = __version__
    _wf_label = 'kkr_eos_wc_{}' # replace with structure formula
    _wf_description = 'Equation of states workflow for {} using KKR' # replace with structure formula
    # workflow options (computer settings)
    _options_default = {'queue_name' : '',                # Queue name to submit jobs too
                        'resources': {"num_machines": 1}, # resources to allocate for the job
                        'walltime_sec' : 60*60,           # walltime in seconds after which the job gets killed (gets parsed to KKR)
                        'use_mpi' : True,                  # execute KKR with mpi or without
                        'custom_scheduler_commands' : '', # some additional scheduler commands (e.g. project numbers in job scripts, OpenMP settings, ...) 
                       }
    # workflow settings
    _wf_default = {'scale_range' : [0.94, 1.06],    # range around volume of starting structure which eos is computed
                   'nsteps': 7,                     # number of calculations around 
                   'ground_state_structure': True,  # create and return a structure which has the ground state volume determined by the fit used
                   'fitfunction': 'birchmurnaghan', # fitfunction used to determine ground state volume (see ase.eos.EquationOfState class for details)
                   'settings_kkr_startpot': kkr_startpot_wc.get_wf_defaults(), # settings for kkr_startpot behavior
                   'settings_kkr_scf': kkr_scf_wc.get_wf_defaults()            # settings for kkr_scf behavior
                   }

    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults, _options_default
        """
        print('Version of workflow: {}'.format(self._workflowversion))
        return self._wf_default, self._options_default


    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow. 
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_eos_wc, cls).define(spec)
        spec.input("options", valid_type=ParameterData, required=False,         # computer options
                   default=ParameterData(dict=cls._options_default))
        spec.input("wf_parameters", valid_type=ParameterData, required=False,   # workfunction settings
                   default=ParameterData(dict=cls._wf_default))
        spec.input("kkr", valid_type=Code, required=True)                       # KKRhost code
        spec.input("voronoi", valid_type=Code, required=True)                   # voronoi code
        spec.input("structure", valid_type=StructureData, required=True)        # starting structure node
        spec.input("calc_parameters", valid_type=ParameterData, required=False) # KKR input parameters (lmax etc.)

        # Here the structure of the workflow is defined
        spec.outline(
            # 1. initialize workflow and check input consistency
            cls.start,
            # 2. prepare structures
            cls.prepare_strucs,
            # 3. run voronoi calculation for smallest volume0
            cls.run_vorostart,
            # 4. check voronoi output and extract RMTCORE parameter
            cls.check_voro_out,
            # 5. submit KKR calculations for all steps
            cls.run_kkr_steps,
            # 6. collect output and fit results
            cls.collect_data_and_fit,
            # 7. collect and return output nodes
            cls.return_results
        )

        # ToDo: improve error codes
        spec.exit_code(221, 'ERROR_INVALID_INPUT', 
            message="ERROR: inputs invalid")


    def start(self):
        """
        initialize context and check input nodes
        """
        
        self.report('INFO: starting KKR eos workflow version {}'.format(self._workflowversion))

        # now extract information from input nodes
        try:
            self.ctx.wf_options = self.inputs.get('options').get_dict()
            self.ctx.wf_parameters = self.inputs.get('wf_parameters').get_dict()
            self.ctx.kkr = self.inputs.get('kkr')   #TODO: check if code is KKR code
            self.ctx.voro = self.inputs.get('voronoi') #TODO: check if code is voronoi code
            self.ctx.structure = self.inputs.get('structure')
            self.ctx.calc_parameters = self.inputs.get('calc_parameters') # optional, TODO: needs to be filled with defaults if not present
        except:
            # in case of failure, exit workflow here
            return self.exit_codes.ERROR_INVALID_INPUT

        # initialize some other things used to collect results etc.
        self.ctx.successful = True
        self.ctx.nsteps = self.ctx.wf_parameters.get('nsteps')
        self.ctx.scale_range = self.ctx.wf_parameters.get('scale_range')
        self.ctx.fitfunc_gs_out = self.ctx.wf_parameters.get('fitfunction') # fitfunction used to get ground state structure
        self.ctx.return_gs_struc = self.ctx.wf_parameters.get('ground_state_structure') # boolean, return output structure or not
        self.ctx.sub_wf_ids = {}             # for uuids of sub workflows, stored as (name: uuid) in dict
        self.ctx.scaled_structures = []      # filled in prepare_strucs

        # set scale_factors from scale_range and nsteps
        self.ctx.scale_factors = [] 
        for i in range(self.ctx.nsteps):
            scale_fac = self.ctx.scale_range[0]+i*(self.ctx.scale_range[1]-self.ctx.scale_range[0])/(self.ctx.nsteps-1)
            self.ctx.scale_factors.append(scale_fac)


    def prepare_strucs(self):
        """
        create new set of scaled structures using the 'rescale' workfunction (see end of the workflow)
        """
        for scale_fac in self.ctx.scale_factors:
            scaled_structure = rescale(self.ctx.structure, Float(scale_fac))
            self.ctx.scaled_structures.append(scaled_structure)


    def run_vorostart(self):
        """
        run vorostart workflow for smallest structure to determine rmtcore setting for all others
        """
        wfd = kkr_startpot_wc.get_wf_defaults()
        set_keys = []
        # first set options
        for key in self.ctx.wf_options.keys():
            wfd[key] = self.ctx.wf_options.get(key)
            set_keys.append(key)
        # then set ef_settings
        vorostart_settings = self.ctx.wf_parameters.get('settings_kkr_startpot')
        for key in vorostart_settings.keys():
            if key not in set_keys: # skip setting of options (done above already)
                wfd[key] = vorostart_settings[key]
        scaled_struc = self.ctx.scaled_structures[0]
        future = self.submit(kkr_startpot_wc, structure=scaled_struc, kkr=self.ctx.kkr, voronoi=self.ctx.voro, 
                             wf_parameters=ParameterData(dict=wfd), calc_parameters=self.ctx.calc_parameters)
        
        self.report('INFO: running kkr_startpot workflow (pk= {})'.format(future.pk))

        return ToContext(kkr_startpot=future)


    def check_voro_out(self):
        """
        check outout of vorostart workflow and create input for rest of calculations (rmtcore setting etc.)
        """
        self.report('INFO: checking voronoi output')
        # get output of kkr_startpot
        out_wc = self.ctx.kkr_startpot
        res = out_wc.out.results_vorostart_wc
        voro_params = out_wc.out.last_params_voronoi
        smallest_voro_remote = out_wc.out.last_voronoi_remote
        smallest_voro_results = out_wc.out.last_voronoi_results
        vorostart_success = res.get_dict()['successful']

        if vorostart_success:
            rmt = []
            radii = smallest_voro_results.get_dict()['radii_atoms_group']
            for rad_iatom in radii:
                if 'rmt0' in rad_iatom.keys():
                    rmt.append(rad_iatom['rmt0'])
            rmtcore_min = array(rmt) * smallest_voro_results.get_dict().get('alat') # needs to be mutiplied by alat in atomic units!
            self.report('INFO: extracted rmtcore_min')
        else:
            return self.error_code(222)

        # update parameter node with rmtcore setting
        voro_params_with_rmtcore = kkrparams(**voro_params.get_dict())
        voro_params_with_rmtcore.set_value('<RMTCORE>', rmtcore_min)
        voro_params_with_rmtcore_dict = voro_params_with_rmtcore.get_dict()
        voro_params_with_rmtcore = update_params_wf(voro_params, ParameterData(dict=voro_params_with_rmtcore_dict))
        self.report('INFO: updated kkr_parameters inlcuding RMTCORE setting')

        # store links to context
        self.ctx.params_kkr_run=voro_params_with_rmtcore
        self.ctx.smallest_voro_remote=smallest_voro_remote


    def run_kkr_steps(self):
        """
        submit KKR calculations for all structures, skip vorostart step for smallest structure
        """

        self.report('INFO: running kkr scf steps')
        # params for scf wfd
        wfd = kkr_scf_wc.get_wf_defaults()
        set_keys = []
        # first set options
        for key in self.ctx.wf_options.keys():
            wfd[key] = self.ctx.wf_options.get(key)
            set_keys.append(key)
        # then set ef_settings
        kkr_scf_settings = self.ctx.wf_parameters.get('settings_kkr_scf')
        for key in kkr_scf_settings.keys():
            if key not in set_keys: # skip setting of options (done above already)
                wfd[key] = kkr_scf_settings[key]

        # used to collect all submitted calculations
        calcs = {}
        
        # submit first calculation separately
        future = self.submit(kkr_scf_wc, kkr=self.ctx.kkr, remote_data=self.ctx.smallest_voro_remote, 
                             wf_parameters=ParameterData(dict=wfd), calc_parameters=self.ctx.params_kkr_run)
        scale_fac = self.ctx.scale_factors[0]
        calcs['kkr_{}_{}'.format(1, scale_fac)] = future

        # then also submit the rest of the calculations
        for i in range(len(self.ctx.scale_factors)-1):
            scale_fac = self.ctx.scale_factors[i+1]
            scaled_struc = self.ctx.scaled_structures[i+1]
            self.report('submit calc for scale fac= {} on {}'.format(scale_fac, scaled_struc.get_formula()))
            future = self.submit(kkr_scf_wc, structure=scaled_struc, kkr=self.ctx.kkr, voronoi=self.ctx.voro, 
                                 wf_parameters=ParameterData(dict=wfd), calc_parameters=self.ctx.params_kkr_run)
            calcs['kkr_{}_{}'.format(i+2, scale_fac)] = future

        self.report('INFO: submitted calculations: {}'.format(calcs))

        return ToContext(**calcs)


    def collect_data_and_fit(self):
        """
        collect output of KKR calculations and perform eos fitting to collect results
        """
        self.report('INFO: collect kkr results and fit data')
        calcs = self.ctx.kkr_scf_steps
        calc_pks = [cid.pk for cid in calcs]
        etot = []
        for iic in range(len(calc_pks)):
            pk = calc_pks[iic]
            n = load_node(pk)
            try:
                d_result = n.out.output_kkr_scf_wc_ParameterResults.get_dict()
                self.report(pk, d_result[u'successful'], d_result[u'convergence_value'])
                if d_result[u'successful']:
                    pk_last_calc = d_result['last_calc_nodeinfo']['pk']
                    n2 = load_node(pk_last_calc)
                    scale = self.ctx.scale_factors[iic]
                    ener = n2.out.output_parameters.get_dict()['total_energy_Ry']
                    rms = d_result[u'convergence_value']
                    scaled_struc = self.ctx.scaled_structures[iic]
                    v = scaled_struc.get_cell_volume()
                    etot.append([scale, ener, v, rms])
            except AttributeError:
                print pk, False, n.is_finished, n.is_finished_ok
               
        # collect calculation outcome
        etot = array(etot)
        scalings = etot[:,0]
        rms = etot[:,-1]
        # convert to eV and per atom units
        etot = etot/len(struc.sites) # per atom values
        etot[:,1] = etot[:,1] * get_Ry2eV() # convert energy from Ry to eV
        volumes, energies = etot[:,2], etot[:,1]-min(etot[:,1])
       
        # do multiple fits to data
        self.report('INFO: output of fits:')
        self.report('{:18} {:8} {:7} {:7}'.format('fitfunc', 'v0', 'e0', 'B'))
        self.report('-----------------------------------------')
        fitnames = ['sj', 'taylor', 'murnaghan', 'birch', 'birchmurnaghan', 'pouriertarantola', 'vinet', 'antonschmidt', 'p3']
        alldat = []
        fitdata = {}
        for fitfunc in fitnames:
            eos = EquationOfState(volumes, energies, eos=fitfunc)
            v0, e0, B = eos.fit()
            fitdata[fitfunc] = [v0, e0, B]
            alldat.append([v0, e0, B])
            self.report('{:16} {:8.3f} {:7.3f} {:7.3f}'.format(fitfunc, v0, e0, B))
        alldat = array(alldat)
        self.report('-----------------------------------------')
        self.report('{:16} {:8.3f} {:7.3f} {:7.3f}'.format('mean', mean(alldat[:,0]), mean(alldat[:,1]), mean(alldat[:,2])))
        self.report('{:16} {:8.3f} {:7.3f} {:7.3f}'.format('std', std(alldat[:,0]), std(alldat[:,1]), std(alldat[:,2])))

        # store results in context
        self.ctx.volumes=volumes
        self.ctx.energies=energies
        self.ctx.scalings=scalings
        self.ctx.fitdata=fitdata


    def return_results(self):
        """
        create output dictionary and run output node generation
        """
        self.report('INFO: create output node')
        outdict = {}
        outdict['successful'] = self.ctx.successful
        outdict['sub_workflow_uuids'] = self.ctx.sub_wf_ids
        outdict['nsteps_input'] = self.ctx.nsteps
        outdict['scale_range_input'] = self.ctx.scale_range
        outdict['scale_factors_all'] = self.ctx.scale_factors
        outdict['volumes'] = self.ctx.volumes
        outdict['energies'] = self.ctx.energies
        outdict['scalings'] = self.ctx.scalings
        outdict['rms'] = self.ctx.rms
        outdict['parameter_fits'] = self.ctx.fitdata
        if self.ctx.return_gs_struc:
            # final result: scaling factor for equilibium 
            scale_fac0 = v0/struc.get_cell_volume()*len(struc.sites)
            outdict['gs_scale_factor'] = scale_fac0
            outdict['gs_fitfunction'] = self.ctx.fitfunc_gs_out
            gs_structure = rescale(self.ctx.structure, Float(scale_fac0))
            gs_structure.label = 'ground_state_structure_{}'.format(gs_structure.get_formula())
            gs_structure.description = 'Ground state structure of {} after running eos workflow. Uses {} fit.'.format(gs_structure.get_formula(), self.ctx.fitfunc_gs_out)
            outdict['gs_structure'] = gs_structure

        # create output nodes with links
        for link_name, node in outdict.iteritems():
            self.out(link_name, node)


### Helper functions and workfunctions ###

    
# copied from aiida_fleur.tools.StructureData_util
def rescale_no_wf(structure, scale):
    """
    Rescales a crystal structure. DOES NOT keep the provanence in the database.

    :param structure, a StructureData node (pk, or uuid)
    :param scale, float scaling factor for the cell

    :returns: New StrcutureData node with rescalled structure, which is linked to input Structure
              and None if inp_structure was not a StructureData
    """

    the_ase = structure.get_ase()
    new_ase = the_ase.copy()
    new_ase.set_cell(the_ase.get_cell()*float(scale), scale_atoms=True)
    rescaled_structure = DataFactory('structure')(ase=new_ase)

    return rescaled_structure

# copied from aiida_fleur.tools.StructureData_util
@wf
def rescale(inp_structure, scale):
    """
    Rescales a crystal structure. Keeps the provanance in the database.

    :param inp_structure, a StructureData node (pk, or uuid)
    :param scale, float scaling factor for the cell

    :returns: New StrcutureData node with rescalled structure, which is linked to input Structure
              and None if inp_structure was not a StructureData
    """

    return rescale_no_wf(inp_structure, scale)
