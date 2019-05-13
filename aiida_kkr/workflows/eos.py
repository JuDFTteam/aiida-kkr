#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a EOS calculation and
some helper methods to do so with AiiDA
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from aiida.orm import Code, load_node
from aiida.plugins import DataFactory
from aiida.orm import Float, Bool
from aiida.engine import WorkChain, ToContext
from aiida.engine import calcfunction
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.common_workfunctions import update_params_wf
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
from masci_tools.io.kkr_params import kkrparams
from masci_tools.io.common_functions import get_Ry2eV
from ase.eos import EquationOfState
from numpy import array, mean, std, min, sort
from six.moves import range


__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7"
__contributors__ = u"Philipp Rüßmann"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
Dict = DataFactory('dict')

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
                        'max_wallclock_seconds' : 60*60,  # walltime in seconds after which the job gets killed (gets parsed to KKR)
                        'use_mpi' : True,                 # execute KKR with mpi or without
                        'custom_scheduler_commands' : ''  # some additional scheduler commands (e.g. project numbers in job scripts, OpenMP settings, ...)
                       }
    # workflow settings
    _wf_default = {'scale_range' : [0.94, 1.06],    # range around volume of starting structure which eos is computed
                   'nsteps': 7,                     # number of calculations around
                   'ground_state_structure': True,  # create and return a structure which has the ground state volume determined by the fit used
                   'use_primitive_structure': True, # use seekpath to get primitive structure after scaling to reduce computational time
                   'fitfunction': 'birchmurnaghan', # fitfunction used to determine ground state volume (see ase.eos.EquationOfState class for details)
                   'settings_kkr_startpot': kkr_startpot_wc.get_wf_defaults(silent=True), # settings for kkr_startpot behavior
                   'settings_kkr_scf': kkr_scf_wc.get_wf_defaults(silent=True)            # settings for kkr_scf behavior
                   }
    # change _wf_default of kkr_scf to deactivate DOS runs
    _wf_default['settings_kkr_scf']['check_dos'] = False

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults, _options_default
        """
        if not silent: print('Version of workflow: {}'.format(self._workflowversion))
        return self._wf_default, self._options_default


    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_eos_wc, cls).define(spec)
        spec.input("options", valid_type=Dict, required=False,         # computer options
                   default=Dict(dict=cls._options_default))
        spec.input("wf_parameters", valid_type=Dict, required=False,   # workfunction settings
                   default=Dict(dict=cls._wf_default))
        spec.input("kkr", valid_type=Code, required=True)                       # KKRhost code
        spec.input("voronoi", valid_type=Code, required=True)                   # voronoi code
        spec.input("structure", valid_type=StructureData, required=True)        # starting structure node
        spec.input("calc_parameters", valid_type=Dict, required=False) # KKR input parameters (lmax etc.)

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
        spec.exit_code(222, 'ERROR_NOT_ENOUGH_SUCCESSFUL_CALCS',
            message='ERROR: need at least 3 successful calculations')
        spec.exit_code(223, 'ERROR_NSTEPS_TOO_SMALL',
            message='ERROR: nsteps is smaller than 3, need at least three data points to do fitting')
        spec.exit_code(224, 'ERROR_INVALID_FITFUN',
            message='given fitfunction name not valid')
        spec.exit_code(225, 'ERROR_VOROSTART_NOT_SUCCESSFUL',
            message='ERROR: kkr_startpot was not successful. Check you inputs.')


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

        # add label and description if not given (contains structure name)
        #if self.label is None:
        self.ctx.label = self._wf_label.format(self.ctx.structure.get_formula())
        #if self.description is None:
        self.ctx.description = self._wf_description.format(self.ctx.structure.get_formula())
        if self.ctx.wf_parameters['settings_kkr_startpot'].get('_label', None) is None:
            self.ctx.wf_parameters['settings_kkr_startpot']['_label'] = self.ctx.label+'_kkr_startpot_{}'.format(self.ctx.structure.get_formula())
        if self.ctx.wf_parameters['settings_kkr_startpot'].get('_description', None) is None:
            self.ctx.wf_parameters['settings_kkr_startpot']['_description'] = self.ctx.description+' kkr_startpot step for {}'.format(self.ctx.structure.get_formula())
        if self.ctx.wf_parameters['settings_kkr_scf'].get('label', None) is None:
            self.ctx.wf_parameters['settings_kkr_scf']['label'] = self.ctx.label+'_kkr_scf_{}'.format(self.ctx.structure.get_formula())
        if self.ctx.wf_parameters['settings_kkr_scf'].get('description', None) is None:
            self.ctx.wf_parameters['settings_kkr_scf']['description'] = self.ctx.description+' kkr_scf step for {}'.format(self.ctx.structure.get_formula())

        # initialize some other things used to collect results etc.
        self.ctx.successful = True
        self.ctx.warnings = []
        self.ctx.rms_threshold = self.ctx.wf_parameters['settings_kkr_scf'].get('convergence_criterion', 10**-7)
        self.ctx.nsteps = self.ctx.wf_parameters.get('nsteps')
        self.ctx.scale_range = self.ctx.wf_parameters.get('scale_range')
        self.ctx.fitfunc_gs_out = self.ctx.wf_parameters.get('fitfunction') # fitfunction used to get ground state structure
        self.ctx.return_gs_struc = self.ctx.wf_parameters.get('ground_state_structure') # boolean, return output structure or not
        self.ctx.use_primitive_structure = self.ctx.wf_parameters.get('use_primitive_structure')
        self.ctx.scaled_structures = []      # filled in prepare_strucs
        self.ctx.fitnames = ['sj', 'taylor', 'murnaghan', 'birch', 'birchmurnaghan', 'pouriertarantola', 'vinet', 'antonschmidt', 'p3'] # list of allowed fits
        self.ctx.sub_wf_ids = {} # filled with workflow uuids

        # check input
        if self.ctx.nsteps<3:
            self.exit_codes.ERROR_NSTEPS_TOO_SMALL
        if self.ctx.fitfunc_gs_out not in self.ctx.fitnames:
            self.exit_codes.ERROR_INVALID_FITFUN

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
            if self.ctx.use_primitive_structure:
                scaled_structure = get_primitive_structure(scaled_structure, Bool(False))
            self.ctx.scaled_structures.append(scaled_structure)


    def run_vorostart(self):
        """
        run vorostart workflow for smallest structure to determine rmtcore setting for all others
        """
        wfd = kkr_startpot_wc.get_wf_defaults(silent=True)
        set_keys = []
        # first set options
        for key in list(self.ctx.wf_options.keys()):
            wfd[key] = self.ctx.wf_options.get(key)
            set_keys.append(key)
        # then set ef_settings
        vorostart_settings = self.ctx.wf_parameters.get('settings_kkr_startpot')
        for key in list(vorostart_settings.keys()):
            if key not in set_keys: # skip setting of options (done above already)
                wfd[key] = vorostart_settings[key]
        scaled_struc = self.ctx.scaled_structures[0]
        future = self.submit(kkr_startpot_wc, structure=scaled_struc, kkr=self.ctx.kkr,
                             voronoi=self.ctx.voro, wf_parameters=Dict(dict=wfd),
                             calc_parameters=self.ctx.calc_parameters,
                             options=Dict(dict=self.ctx.wf_options))

        self.report('INFO: running kkr_startpot workflow (pk= {})'.format(future.pk))
        self.ctx.sub_wf_ids['kkr_startpot_1'] = future.uuid

        return ToContext(kkr_startpot=future)


    def check_voro_out(self):
        """
        check outout of vorostart workflow and create input for rest of calculations (rmtcore setting etc.)
        """
        self.report('INFO: checking voronoi output')
        # get output of kkr_startpot
        out_wc = self.ctx.kkr_startpot
        try:
            res = out_wc.outputs.results_vorostart_wc
            voro_params = out_wc.outputs.last_params_voronoi
            smallest_voro_remote = out_wc.outputs.last_voronoi_remote
            smallest_voro_results = out_wc.outputs.last_voronoi_results
            vorostart_success = res.get_dict()['successful']
        except AttributeError:
            vorostart_success = False

        if vorostart_success:
            rmt = []
            radii = smallest_voro_results.get_dict()['radii_atoms_group']
            for rad_iatom in radii:
                if 'rmt0' in list(rad_iatom.keys()):
                    rmt.append(rad_iatom['rmt0'])
            rmtcore_min = array(rmt) * smallest_voro_results.get_dict().get('alat') # needs to be mutiplied by alat in atomic units!
            self.report('INFO: extracted rmtcore_min ({})'.format(rmtcore_min))
        else:
            return self.exit_codes.ERROR_VOROSTART_NOT_SUCCESSFUL

        # update parameter node with rmtcore setting
        voro_params_with_rmtcore = kkrparams(**voro_params.get_dict())
        voro_params_with_rmtcore.set_value('<RMTCORE>', rmtcore_min)
        voro_params_with_rmtcore_dict = voro_params_with_rmtcore.get_dict()
        voro_params_with_rmtcore = update_params_wf(voro_params, Dict(dict=voro_params_with_rmtcore_dict))
        self.report('INFO: updated kkr_parameters inlcuding RMTCORE setting (uuid={})'.format(voro_params_with_rmtcore.uuid))

        # store links to context
        self.ctx.params_kkr_run=voro_params_with_rmtcore
        self.ctx.smallest_voro_remote=smallest_voro_remote


    def run_kkr_steps(self):
        """
        submit KKR calculations for all structures, skip vorostart step for smallest structure
        """

        self.report('INFO: running kkr scf steps')
        # params for scf wfd
        wfd = kkr_scf_wc.get_wf_defaults(silent=True)
        set_keys = []
        # first set options
        for key in list(self.ctx.wf_options.keys()):
            wfd[key] = self.ctx.wf_options.get(key)
            set_keys.append(key)
        # then set ef_settings
        kkr_scf_settings = self.ctx.wf_parameters.get('settings_kkr_scf')
        for key in list(kkr_scf_settings.keys()):
            if key not in set_keys: # skip setting of options (done above already)
                wfd[key] = kkr_scf_settings[key]

        # used to collect all submitted calculations
        calcs = {}

        # submit first calculation separately
        self.report('submit calc for scale fac= {} on {}'.format(self.ctx.scale_factors[0], self.ctx.scaled_structures[0].get_formula()))
        future = self.submit(kkr_scf_wc, kkr=self.ctx.kkr, remote_data=self.ctx.smallest_voro_remote,
                             wf_parameters=Dict(dict=wfd), calc_parameters=self.ctx.params_kkr_run,
                             options=Dict(dict=self.ctx.wf_options))
        scale_fac = self.ctx.scale_factors[0]
        calcs['kkr_{}_{}'.format(1, scale_fac)] = future
        self.ctx.sub_wf_ids['kkr_scf_1'] = future.uuid

        # then also submit the rest of the calculations
        for i in range(len(self.ctx.scale_factors)-1):
            scale_fac = self.ctx.scale_factors[i+1]
            scaled_struc = self.ctx.scaled_structures[i+1]
            self.report('submit calc for scale fac= {} on {}'.format(scale_fac, scaled_struc.get_formula()))
            future = self.submit(kkr_scf_wc, structure=scaled_struc, kkr=self.ctx.kkr, voronoi=self.ctx.voro,
                                 wf_parameters=Dict(dict=wfd), calc_parameters=self.ctx.params_kkr_run,
                                 options=Dict(dict=self.ctx.wf_options))
            calcs['kkr_{}_{}'.format(i+2, scale_fac)] = future
            self.ctx.sub_wf_ids['kkr_scf_{}'.format(i+2)] = future.uuid

        # save uuids of calculations to context
        self.ctx.kkr_calc_uuids = []
        for name in sort(list(calcs.keys())): # sorting important to have correct assignment of scaling and structure info later on
            calc = calcs[name]
            self.ctx.kkr_calc_uuids.append(calc.uuid)

        self.report('INFO: submitted calculations: {}'.format(calcs))

        return ToContext(**calcs)


    def collect_data_and_fit(self):
        """
        collect output of KKR calculations and perform eos fitting to collect results
        """
        self.report('INFO: collect kkr results and fit data')
        calc_uuids = self.ctx.kkr_calc_uuids
        etot = []
        for iic in range(len(calc_uuids)):
            uuid = calc_uuids[iic]
            n = load_node(uuid)
            try:
                d_result = n.outputs.output_kkr_scf_wc_ParameterResults.get_dict()
                self.report('INFO: extracting output of calculation {}: successful={}, rms={}'.format(uuid, d_result[u'successful'], d_result[u'convergence_value']))
                if d_result[u'successful']:
                    pk_last_calc = d_result['last_calc_nodeinfo']['pk']
                    n2 = load_node(pk_last_calc)
                    scale = self.ctx.scale_factors[iic]
                    ener = n2.outputs.output_parameters.get_dict()['total_energy_Ry']
                    rms = d_result[u'convergence_value']
                    scaled_struc = self.ctx.scaled_structures[iic]
                    v = scaled_struc.get_cell_volume()
                    if rms<=self.ctx.rms_threshold: # only take those calculations which
                        etot.append([scale, ener, v, rms])
                    else:
                        warn = 'rms of calculation with uuid={} not low enough ({} > {})'.format(uuid, rms, self.ctx.rms_threshold)
                        self.report('WARNING: {}'.format(warn))
                        self.ctx.warnings.append(warn)
            except AttributeError:
                warn = 'calculation with uuid={} not successful'.format(uuid)
                self.report('WARNING: {}'.format(warn))
                self.ctx.warnings.append(warn)


        # collect calculation outcome
        etot = array(etot)
        self.report('INFO: collected data from calculations= {}'.format(etot))

        # check if at least 3 points were successful (otherwise fit does not work)
        if len(etot)<3:
            return self.exit_codes.ERROR_NOT_ENOUGH_SUCCESSFUL_CALCS

        scalings = etot[:,0]
        rms = etot[:,-1]
        # convert to eV and per atom units
        etot = etot/len(scaled_struc.sites) # per atom values
        etot[:,1] = etot[:,1] * get_Ry2eV() # convert energy from Ry to eV
        volumes, energies = etot[:,2], etot[:,1]

        # do multiple fits to data
        self.report('INFO: output of fits:')
        self.report('{:18} {:8} {:7} {:7}'.format('fitfunc', 'v0', 'e0', 'B'))
        self.report('-----------------------------------------')
        fitnames = self.ctx.fitnames
        alldat = []
        fitdata = {}
        for fitfunc in fitnames:
            try:
                eos = EquationOfState(volumes, energies, eos=fitfunc)
                v0, e0, B = eos.fit()
                fitdata[fitfunc] = [v0, e0, B]
                alldat.append([v0, e0, B])
                self.report('{:16} {:8.3f} {:7.3f} {:7.3f}'.format(fitfunc, v0, e0, B))
            except: # capture all errors and mark fit as unsuccessful
               self.ctx.warnings.append('fit unsuccessful for {} function'.format(fitfunc))
               if fitfunc == self.ctx.fitfunc_gs_out:
                   self.ctx.successful = False
        alldat = array(alldat)
        self.report('-----------------------------------------')
        self.report('{:16} {:8.3f} {:7.3f} {:7.3f}'.format('mean', mean(alldat[:,0]), mean(alldat[:,1]), mean(alldat[:,2])))
        self.report('{:16} {:8.3f} {:7.3f} {:7.3f}'.format('std', std(alldat[:,0]), std(alldat[:,1]), std(alldat[:,2])))

        # store results in context
        self.ctx.volumes=volumes
        self.ctx.energies=energies
        self.ctx.scalings=scalings
        self.ctx.rms = rms
        self.ctx.fitdata=fitdata
        self.ctx.fit_mean_values={'<v0>':mean(alldat[:,0]), '<e0>':mean(alldat[:,1]), '<B>':mean(alldat[:,2])}
        self.ctx.fit_std_values={'s_v0':std(alldat[:,0]), 's_e0':std(alldat[:,1]), 's_B':std(alldat[:,2])}


    def return_results(self):
        """
        create output dictionary and run output node generation
        """
        self.report('INFO: create output node')
        outdict = {}
        outdict['successful'] = self.ctx.successful
        outdict['warnings'] = self.ctx.warnings
        outdict['sub_workflow_uuids'] = self.ctx.sub_wf_ids
        outdict['nsteps_input'] = self.ctx.nsteps
        outdict['scale_range_input'] = self.ctx.scale_range
        outdict['scale_factors_all'] = self.ctx.scale_factors
        outdict['volumes'] = self.ctx.volumes
        outdict['energies'] = self.ctx.energies
        outdict['scalings'] = self.ctx.scalings
        outdict['rms'] = self.ctx.rms
        outdict['parameter_fits'] = self.ctx.fitdata
        outdict['fits_mean'] = self.ctx.fit_mean_values
        outdict['fits_std'] = self.ctx.fit_std_values
        outdict['formula'] = self.ctx.structure.get_formula()
        outdict['label'] = self.ctx.label
        if self.ctx.successful and self.ctx.return_gs_struc:
            # final result: scaling factor for equilibium
            v0, e0, B = self.ctx.fitdata.get(self.ctx.fitfunc_gs_out)
            scale_fac0 = v0/self.ctx.structure.get_cell_volume()*len(self.ctx.structure.sites)
            outdict['gs_scale_factor'] = scale_fac0
            outdict['gs_fitfunction'] = self.ctx.fitfunc_gs_out
            gs_structure = rescale(self.ctx.structure, Float(scale_fac0))
            if self.ctx.use_primitive_structure:
                tmpdict = get_primitive_structure(gs_structure, Bool(True))
                conv_structure, explicit_kpoints, parameters, gs_structure = tmpdict['conv_structure'], tmpdict['explicit_kpoints'], tmpdict['parameters'], tmpdict['primitive_structure']
                outdict['gs_kpoints_seekpath_params_uuid'] = parameters.uuid
            gs_structure.label = 'ground_state_structure_{}'.format(gs_structure.get_formula())
            gs_structure.description = 'Ground state structure of {} after running eos workflow. Uses {} fit.'.format(gs_structure.get_formula(), self.ctx.fitfunc_gs_out)
            outdict['gs_structure_uuid'] = gs_structure.uuid

        # create output nodes in dict with link names
        outnodes = {'eos_results': Dict(dict=outdict)}
        if self.ctx.successful and self.ctx.return_gs_struc:
            outnodes['gs_structure'] = gs_structure
            if self.ctx.use_primitive_structure:
                outnodes['explicit_kpoints'] = explicit_kpoints
                outnodes['get_explicit_kpoints_path_parameters'] = parameters
        # set out nodes and corresponding link names
        for link_name, node in outnodes.items():
            self.out(link_name, node)


### Helper functions and workfunctions ###

def rescale_no_wf(structure, scale):
    """
    Rescales a crystal structure. DOES NOT keep the provanence in the database.

    :param structure, a StructureData node (pk, or uuid)
    :param scale, float scaling factor for the cell

    :returns: New StrcutureData node with rescalled structure, which is linked to input Structure
              and None if inp_structure was not a StructureData

    copied and modified from aiida_fleur.tools.StructureData_util
    """

    the_ase = structure.get_ase()
    new_ase = the_ase.copy()
    new_ase.set_cell(the_ase.get_cell()*float(scale), scale_atoms=True)
    rescaled_structure = DataFactory('structure')(ase=new_ase)

    return rescaled_structure

@calcfunction
def rescale(inp_structure, scale):
    """
    Rescales a crystal structure. Keeps the provanance in the database.

    :param inp_structure, a StructureData node (pk, or uuid)
    :param scale, float scaling factor for the cell

    :returns: New StrcutureData node with rescalled structure, which is linked to input Structure
              and None if inp_structure was not a StructureData

    copied and modified from aiida_fleur.tools.StructureData_util
    """

    return rescale_no_wf(inp_structure, scale)


@calcfunction
def get_primitive_structure(structure, return_all):
    """
    calls get_explicit_kpoints_path which gives primitive structure
    auxiliary workfunction to keep provenance
    """
    from aiida.tools import get_explicit_kpoints_path
    output = get_explicit_kpoints_path(structure)
    conv_structure = output['conv_structure']
    explicit_kpoints = output['explicit_kpoints']
    parameters = output['parameters']
    primitive_structure = output['primitive_structure']
    if return_all:
        return {'conv_structure':conv_structure, 'explicit_kpoints':explicit_kpoints, 'parameters':parameters, 'primitive_structure':primitive_structure}
    else:
        return primitive_structure
