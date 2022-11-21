#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a EOS calculation and
some helper methods to do so with AiiDA
"""
import copy
import numpy as np
from masci_tools.io.kkr_params import kkrparams
from masci_tools.io.common_functions import get_Ry2eV
from ase.eos import EquationOfState
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistentAttributeError
from aiida.engine import WorkChain, ToContext
from aiida.engine import calcfunction
from aiida_kkr.tools.common_workfunctions import update_params_wf
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
from aiida_kkr.tools.save_output_nodes import create_out_dict_node

__copyright__ = (u'Copyright (c), 2018, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.9.2'
__contributors__ = u'Philipp Rüßmann'


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
    _wf_label = 'kkr_eos_wc_{}'  # replace with structure formula
    # replace with structure formula
    _wf_description = 'Equation of states workflow for {} using KKR'
    # workflow options (computer settings)
    _options_default = AttributeDict()
    # Queue name to submit jobs too
    _options_default.queue_name = ''
    _options_default.resources = AttributeDict()
    _options_default.resources.num_machines = 1
    # walltime in seconds after which the job gets killed (gets parsed to KKR)
    _options_default.max_wallclock_seconds = 60 * 60
    # execute KKR with mpi or without
    _options_default.withmpi = True
    # some additional scheduler commands
    # (e.g. project numbers in job scripts, OpenMP settings, ...)
    _options_default.custom_scheduler_commands = ''
    # workflow settings
    _wf_default = AttributeDict()
    # range around volume of starting structure which eos is computed
    _wf_default.scale_range = [0.94, 1.06]
    # number of calculations around
    _wf_default.nsteps = 7
    # create and return a structure which has the ground state volume
    # determined by the fit used
    _wf_default.ground_state_structure = True
    # use seekpath to get primitive structure after scaling to reduce
    # computational time
    _wf_default.use_primitive_structure = True
    # fitfunction used to determine ground state volume
    # (see ase.eos.EquationOfState class for details)
    _wf_default.fitfunction = 'birchmurnaghan'
    # settings for kkr_startpot behavior
    _wf_default.settings_kkr_startpot = kkr_startpot_wc.get_wf_defaults(silent=True)
    # settings for kkr_scf behavior
    _wf_default.settings_kkr_scf = kkr_scf_wc.get_wf_defaults(silent=True)[0]
    # change _wf_default of kkr_scf to deactivate DOS runs
    _wf_default.settings_kkr_scf['check_dos'] = False

    @classmethod
    def get_wf_defaults(cls, silent=False):
        """Print and return _wf_defaults dictionary.

        Can be used to easily create set of wf_parameters.
        returns _wf_defaults, _options_default
        """
        if not silent:
            print(f'Version of workflow: {cls._workflowversion}')
        return cls._wf_default, cls._options_default

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_eos_wc, cls).define(spec)
        spec.input(
            'options',
            valid_type=orm.Dict,
            required=False,  # computer options
            default=lambda: orm.Dict(dict=cls._options_default)
        )
        spec.input(
            'wf_parameters',
            valid_type=orm.Dict,
            required=False,  # workfunction settings
            default=lambda: orm.Dict(dict=cls._wf_default),
            help="""
            Workfunction specific parameters, controlling the behavior of the
            EOS workchain.
            """
        )
        spec.input(
            'kkr',
            valid_type=orm.Code,
            required=True,
            help="""
            Code entry for the KKRhost calculations.
            """
        )
        spec.input(
            'voronoi',
            valid_type=orm.Code,
            required=True,
            help="""
            Code entry for the Voronoi calculations.
            """
        )
        spec.input(
            'structure',
            valid_type=orm.StructureData,
            required=True,
            help="""
            Initial structure for which the EOS will be calculated.
            """
        )
        spec.input(
            'calc_parameters',
            valid_type=orm.Dict,
            required=False,
            help="""
            KKR input parameters. Everything (except structural factors) which
            would normally be written in the inputcard.
            """
        )

        # define output nodes
        spec.output(
            'eos_results',
            valid_type=orm.Dict,
            required=True,
        )
        spec.output(
            'gs_structure',
            valid_type=orm.StructureData,
            required=False,
        )
        spec.output(
            'explicit_kpoints',
            valid_type=orm.KpointsData,
            required=False,
        )
        spec.output(
            'get_explicit_kpoints_path_parameters',
            valid_type=orm.Dict,
            required=False,
        )

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
        spec.exit_code(
            221,
            'ERROR_INVALID_INPUT',
            message='ERROR: inputs invalid',
        )
        spec.exit_code(
            222,
            'ERROR_NOT_ENOUGH_SUCCESSFUL_CALCS',
            message='ERROR: need at least 3 successful calculations',
        )
        spec.exit_code(
            223,
            'ERROR_NSTEPS_TOO_SMALL',
            message='ERROR: nsteps is smaller than 3, need at least three data points to do fitting',
        )
        spec.exit_code(
            224,
            'ERROR_INVALID_FITFUN',
            message='given fitfunction name not valid',
        )
        spec.exit_code(
            225,
            'ERROR_VOROSTART_NOT_SUCCESSFUL',
            message='ERROR: kkr_startpot was not successful. Check you inputs.',
        )

    def start(self):
        """
        initialize context and check input nodes
        """

        self.report(f'INFO: starting KKR eos workflow version {self._workflowversion}')

        # now extract information from input nodes
        try:
            self.ctx.wf_options = self.inputs.get('options').get_dict()
            self.ctx.wf_parameters = self.inputs.get('wf_parameters').get_dict()
            # TODO: check if code is KKR code
            self.ctx.kkr = self.inputs.get('kkr')
            # TODO: check if code is voronoi code
            self.ctx.voro = self.inputs.get('voronoi')
            self.ctx.structure = self.inputs.get('structure')
            # optional, TODO: needs to be filled with defaults if not present
            self.ctx.calc_parameters = self.inputs.get('calc_parameters')
        except AttributeError:
            # in case of failure, exit workflow here
            return self.exit_codes.ERROR_INVALID_INPUT  # pylint: disable=no-member

        # add label and description if not given (contains structure name)
        # if self.label is None:
        self.ctx.label = self._wf_label.format(self.ctx.structure.get_formula())
        # if self.description is None:
        self.ctx.description = self._wf_description.format(self.ctx.structure.get_formula())
        if self.ctx.wf_parameters['settings_kkr_startpot'].get('_label', None) is None:
            self.ctx.wf_parameters['settings_kkr_startpot']['_label'] = (
                self.ctx.label + f'_kkr_startpot_{self.ctx.structure.get_formula()}'
            )
        if self.ctx.wf_parameters['settings_kkr_startpot'].get('_description', None) is None:
            self.ctx.wf_parameters['settings_kkr_startpot']['_description'] = (
                self.ctx.description + f' kkr_startpot step for {self.ctx.structure.get_formula()}'
            )
        if self.ctx.wf_parameters['settings_kkr_scf'].get('label', None) is None:
            self.ctx.wf_parameters['settings_kkr_scf']['label'] = (
                self.ctx.label + f'_kkr_scf_{self.ctx.structure.get_formula()}'
            )
        if self.ctx.wf_parameters['settings_kkr_scf'].get('description', None) is None:
            self.ctx.wf_parameters['settings_kkr_scf']['description'] = (
                self.ctx.description + f' kkr_scf step for {self.ctx.structure.get_formula()}'
            )

        # initialize some other things used to collect results etc.
        self.ctx.successful = True
        self.ctx.warnings = []
        self.ctx.rms_threshold = self.ctx.wf_parameters['settings_kkr_scf'].get(
            'convergence_criterion',
            10**-7,
        )
        self.ctx.nsteps = self.ctx.wf_parameters.get(
            'nsteps',
            self._wf_default.nsteps,
        )
        self.ctx.scale_range = self.ctx.wf_parameters.get(
            'scale_range',
            self._wf_default.scale_range,
        )
        # fitfunction used to get ground state structure
        self.ctx.fitfunc_gs_out = self.ctx.wf_parameters.get(
            'fitfunction',
            self._wf_default.fitfunction,
        )
        # boolean, return output structure or not
        self.ctx.return_gs_struc = self.ctx.wf_parameters.get(
            'ground_state_structure',
            self._wf_default.ground_state_structure,
        )
        self.ctx.use_primitive_structure = self.ctx.wf_parameters.get(
            'use_primitive_structure',
            self._wf_default.use_primitive_structure,
        )
        self.ctx.scaled_structures = []  # filled in prepare_strucs
        self.ctx.fitnames = [
            'sj',
            'taylor',
            'murnaghan',
            'birch',
            'birchmurnaghan',
            'pouriertarantola',
            'vinet',
            'antonschmidt',
            'p3',
        ]  # list of allowed fits
        self.ctx.sub_wf_ids = {}  # filled with workflow uuids

        # check input
        if self.ctx.nsteps < 3:
            return self.exit_codes.ERROR_NSTEPS_TOO_SMALL  # pylint: disable=no-member
        if self.ctx.fitfunc_gs_out not in self.ctx.fitnames:
            return self.exit_codes.ERROR_INVALID_FITFUN  # pylint: disable=no-member

        # set scale_factors from scale_range and nsteps
        self.ctx.scale_factors = np.linspace(
            start=self.ctx.scale_range[0],
            stop=self.ctx.scale_range[1],
            num=self.ctx.nsteps,
        )
        return None

    def prepare_strucs(self):
        """Create new set of scaled structures for the E-V curve.

        The structures are generated using the 'rescale' workfunction
        (see end of the workflow)
        """
        for scale_fac in self.ctx.scale_factors:
            scaled_structure = rescale(self.ctx.structure, orm.Float(scale_fac))
            if self.ctx.use_primitive_structure:
                scaled_structure = get_primitive_structure(
                    scaled_structure,
                    orm.Bool(False),
                )
            self.ctx.scaled_structures.append(scaled_structure)

    def run_vorostart(self):
        """Run vorostart workflow for smallest structure to determine rmtcore.

        One needs to run all the calculations with the same rmtcore to be able
        to compare the energies. This value is then set for all others.
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
            # skip setting of options (done above already)
            if key not in set_keys:
                wfd[key] = vorostart_settings[key]
        scaled_struc = self.ctx.scaled_structures[0]
        future = self.submit(
            kkr_startpot_wc,
            structure=scaled_struc,
            kkr=self.ctx.kkr,
            voronoi=self.ctx.voro,
            wf_parameters=orm.Dict(dict=wfd),
            calc_parameters=self.ctx.calc_parameters,
            options=orm.Dict(dict=self.ctx.wf_options)
        )

        self.report(f'INFO: running kkr_startpot workflow (pk= {future.pk})')
        self.ctx.sub_wf_ids['kkr_startpot_1'] = future.uuid

        return ToContext(kkr_startpot=future)

    def check_voro_out(self):
        """Check output of the vorostart workflow.

        The outputs are then used to create inputs for the next set of
        calculations (rmtcore setting etc.)
        """
        self.report('INFO: checking voronoi output')
        # get output of kkr_startpot
        out_wc = self.ctx.kkr_startpot
        res = out_wc.outputs.results_vorostart_wc
        voro_params = out_wc.outputs.last_params_voronoi
        smallest_voro_remote = out_wc.outputs.last_voronoi_remote
        smallest_voro_results = out_wc.outputs.last_voronoi_results
        vorostart_success = res.get_dict().get('successful', False)

        if vorostart_success:
            rmt = []
            radii = smallest_voro_results.get_dict()['radii_atoms_group']
            for rad_iatom in radii:
                if 'rmt0' in list(rad_iatom.keys()):
                    rmt.append(rad_iatom['rmt0'])
            # needs to be multiplied by alat in atomic units!
            rmtcore_min = np.array(rmt) * smallest_voro_results.get_dict().get('alat')
            self.report(f'INFO: extracted rmtcore_min ({rmtcore_min})')
        else:
            self.report(f'ERROR: kkr_startpot workflow not successful')
            return self.exit_codes.ERROR_VOROSTART_NOT_SUCCESSFUL  # pylint: disable=no-member

        # update parameter node with rmtcore setting
        voro_params_with_rmtcore = kkrparams(**voro_params.get_dict())
        voro_params_with_rmtcore.set_value('<RMTCORE>', rmtcore_min)
        voro_params_with_rmtcore_dict = voro_params_with_rmtcore.get_dict()
        voro_params_with_rmtcore = update_params_wf(
            voro_params,
            orm.Dict(dict=voro_params_with_rmtcore_dict),
        )
        self.report(f'INFO: updated kkr_parameters including RMTCORE setting (uuid={voro_params_with_rmtcore.uuid})')

        # store links to context
        self.ctx.params_kkr_run = voro_params_with_rmtcore
        self.ctx.smallest_voro_remote = smallest_voro_remote

        return None

    def run_kkr_steps(self):
        """Submit KKR calculations for all structures.

        This will skip the vorostart step for smallest structure.
        """

        self.report('INFO: running kkr scf steps')
        # params for scf wfd
        wfd = kkr_scf_wc.get_wf_defaults(silent=True)[0]
        set_keys = []
        # first set options
        for key in list(self.ctx.wf_options.keys()):
            wfd[key] = self.ctx.wf_options.get(key)
            set_keys.append(key)
        # then set ef_settings
        kkr_scf_settings = self.ctx.wf_parameters.get('settings_kkr_scf')
        for key in list(kkr_scf_settings.keys()):
            # skip setting of options (done above already)
            if key not in set_keys:
                wfd[key] = kkr_scf_settings[key]

        # used to collect all submitted calculations
        calcs = {}

        # submit first calculation separately
        self.report(
            f'submit calc for scale fac= {self.ctx.scale_factors[0]} on {self.ctx.scaled_structures[0].get_formula()}'
        )
        future = self.submit(
            kkr_scf_wc,
            kkr=self.ctx.kkr,
            remote_data=self.ctx.smallest_voro_remote,
            wf_parameters=orm.Dict(dict=wfd),
            calc_parameters=self.ctx.params_kkr_run,
            options=orm.Dict(dict=self.ctx.wf_options)
        )
        scale_fac = self.ctx.scale_factors[0]
        calcs[f'kkr_1_{scale_fac}'] = future
        self.ctx.sub_wf_ids['kkr_scf_1'] = future.uuid

        # then also submit the rest of the calculations
        for i in range(len(self.ctx.scale_factors) - 1):
            scale_fac = self.ctx.scale_factors[i + 1]
            scaled_struc = self.ctx.scaled_structures[i + 1]
            self.report(f'submit calc for scale fac= {scale_fac} on {scaled_struc.get_formula()}')
            future = self.submit(
                kkr_scf_wc,
                structure=scaled_struc,
                kkr=self.ctx.kkr,
                voronoi=self.ctx.voro,
                wf_parameters=orm.Dict(dict=wfd),
                calc_parameters=self.ctx.params_kkr_run,
                options=orm.Dict(dict=self.ctx.wf_options)
            )
            calcs[f'kkr_{i+2}_{scale_fac}'] = future
            self.ctx.sub_wf_ids[f'kkr_scf_{i+2}'] = future.uuid

        # save uuids of calculations to context
        self.ctx.kkr_calc_uuids = []
        # sorting important to have correct assignment of scaling and structure
        # info later on
        for name in np.sort(list(calcs.keys())):
            calc = calcs[name]
            self.ctx.kkr_calc_uuids.append(calc.uuid)

        self.report(f'INFO: submitted calculations: {calcs}')

        return ToContext(**calcs)

    def collect_data_and_fit(self):
        """Collect output of KKR calculations and perform eos fitting."""
        self.report('INFO: collect kkr results and fit data')
        calc_uuids = self.ctx.kkr_calc_uuids
        etot = []
        for iic, uuid in enumerate(calc_uuids):
            node = orm.load_node(uuid)
            try:
                d_result = node.outputs.output_kkr_scf_wc_ParameterResults.get_dict()
                self.report(
                    f'INFO: extracting output of calculation {uuid}: successful={d_result[u"successful"]}, rms={d_result[u"convergence_value"]}'
                )
                if d_result[u'successful']:
                    pk_last_calc = d_result['last_calc_nodeinfo']['pk']
                    node_2 = orm.load_node(pk_last_calc)
                    scale = self.ctx.scale_factors[iic]
                    ener = node_2.outputs.output_parameters.get_dict()['total_energy_Ry']
                    rms = d_result[u'convergence_value']
                    scaled_struc = self.ctx.scaled_structures[iic]
                    vol = scaled_struc.get_cell_volume()
                    alat = scaled_struc.cell_lengths[0]
                    # Only take those calculations which
                    if rms <= self.ctx.rms_threshold:
                        etot.append([scale, ener, vol, rms, alat])
                    else:
                        warn = f'rms of calculation with uuid={uuid} not low enough ({rms} > {self.ctx.rms_threshold})'

                        self.report(f'WARNING: {warn}')
                        self.ctx.warnings.append(warn)
            except (AttributeError, NotExistentAttributeError):
                warn = f'calculation with uuid={uuid} not successful'
                self.report(f'WARNING: {warn}')
                self.ctx.warnings.append(warn)

        # collect calculation outcome
        etot = np.array(etot)
        self.report(f'INFO: collected data from calculations= {etot}')

        # check if at least 3 points were successful (otherwise fit does not work)
        if len(etot) < 3:
            return self.exit_codes.ERROR_NOT_ENOUGH_SUCCESSFUL_CALCS  # pylint: disable=no-member

        scalings = etot[:, 0]
        rms = etot[:, -1]
        # convert to eV and per atom units
        etot = etot / len(scaled_struc.sites)  # per atom values
        etot[:, 1] = etot[:, 1] * get_Ry2eV()  # convert energy from Ry to eV
        volumes, energies = etot[:, 2], etot[:, 1]

        # do multiple fits to data
        self.report('INFO: output of fits:')
        self.report(f'{"fitfunc":18} {"v0":8} {"e0":7} {"B":7} {"alat":7}')
        self.report('-----------------------------------------')
        fitnames = self.ctx.fitnames
        alldat = []
        fitdata = {}
        for fitfunc in fitnames:
            try:
                eos = EquationOfState(volumes, energies, eos=fitfunc)
                v0, e0, B = eos.fit()
                py_func = scaled_struc.get_pymatgen()
                py_func.scale_lattice(v0)
                alat = py_func.lattice.a
                fitdata[fitfunc] = [v0, e0, B, alat]
                alldat.append([v0, e0, B, alat])
                self.report(f'{fitfunc:16} {v0:8.3f} {e0:7.3f} {B:7.3f} {alat:7.3f}')
            # capture all errors and mark fit as unsuccessful
            except (ValueError, RuntimeError, UserWarning, RuntimeWarning):
                self.ctx.warnings.append(f'fit unsuccessful for {fitfunc} function')
                if fitfunc == self.ctx.fitfunc_gs_out:
                    self.ctx.successful = False
        alldat = np.array(alldat)
        self.report('-----------------------------------------')
        self.report(
            f'{"mean":16} {np.mean(alldat[:, 0]):8.3f} {np.mean(alldat[:, 1]):7.3f} {np.mean(alldat[:, 2]):7.3f} {np.mean(alldat[:, 3]):7.3f}'
        )
        self.report(
            f'{"std":16} {np.std(alldat[:, 0]):8.3f} {np.std(alldat[:, 1]):7.3f} {np.std(alldat[:, 2]):7.3f} {np.std(alldat[:, 3]):7.3f}'
        )

        # store results in context
        self.ctx.volumes = volumes
        self.ctx.energies = energies
        self.ctx.scalings = scalings
        self.ctx.rms = rms
        self.ctx.fitdata = fitdata
        self.ctx.fit_mean_values = {
            '<v0>': np.mean(alldat[:, 0]),
            '<e0>': np.mean(alldat[:, 1]),
            '<B>': np.mean(alldat[:, 2]),
            '<alat>': np.mean(alldat[:, 3])
        }
        self.ctx.fit_std_values = {
            's_v0': np.std(alldat[:, 0]),
            's_e0': np.std(alldat[:, 1]),
            's_B': np.std(alldat[:, 2]),
            's_alat': np.std(alldat[:, 3])
        }

        return None

    def return_results(self):
        """Create output dictionary and run output node generation."""
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
            # final result: scaling factor for equilibrium
            v0, e0, B, alat = self.ctx.fitdata.get(self.ctx.fitfunc_gs_out)
            scale_fac0 = v0 / self.ctx.structure.get_cell_volume() * len(self.ctx.structure.sites)
            outdict['gs_scale_factor'] = scale_fac0
            outdict['gs_fitfunction'] = self.ctx.fitfunc_gs_out
            gs_structure = rescale(self.ctx.structure, orm.Float(scale_fac0))
            if self.ctx.use_primitive_structure:
                tmpdict = get_primitive_structure(gs_structure, orm.Bool(True))
                conv_structure, explicit_kpoints, parameters, gs_structure = tmpdict['conv_structure'], tmpdict[
                    'explicit_kpoints'], tmpdict['parameters'], tmpdict['primitive_structure']
                outdict['gs_kpoints_seekpath_params_uuid'] = parameters.uuid
            gs_structure.label = (f'ground_state_structure_{gs_structure.get_formula()}')
            gs_structure.description = (
                f'Ground state structure of {gs_structure.get_formula()} after running eos workflow. Uses {self.ctx.fitfunc_gs_out} fit.'
            )
            outdict['gs_structure_uuid'] = gs_structure.uuid

        # create output nodes in dict with link names
        outnodes = {}
        if self.ctx.successful and self.ctx.return_gs_struc:
            outnodes['gs_structure'] = gs_structure
            if self.ctx.use_primitive_structure:
                outnodes['explicit_kpoints'] = explicit_kpoints
                outnodes['get_explicit_kpoints_path_parameters'] = parameters

        # create results node with calcfunction for data provenance
        link_nodes = outnodes.copy()
        for wf_label, sub_wf_uuid in self.ctx.sub_wf_ids.items():
            if 'kkr_scf' in wf_label:
                link_nodes[wf_label] = orm.load_node(sub_wf_uuid).outputs.output_kkr_scf_wc_ParameterResults
            else:
                link_nodes[wf_label] = orm.load_node(sub_wf_uuid).outputs.results_vorostart_wc
        outnodes['eos_results'] = create_out_dict_node(orm.Dict(dict=outdict), **link_nodes)

        # set out nodes and corresponding link names
        for link_name, node in outnodes.items():
            self.out(link_name, node)


### Helper functions and workfunctions ###


def rescale_no_wf(structure, scale) -> orm.StructureData:
    """
    Rescales a crystal structure. DOES NOT keep the provenance in the database.

    :param structure, a StructureData node (pk, or uuid)
    :param scale, float scaling factor for the cell

    :returns: New StructureData node with rescalled structure, which is linked to input Structure
              and None if inp_structure was not a StructureData

    copied and modified from aiida_fleur.tools.StructureData_util
    """

    scaled_volume = structure.get_cell_volume() * scale.value
    scaled_structure = copy.deepcopy(structure.get_pymatgen())
    scaled_structure.scale_lattice(scaled_volume)

    return orm.StructureData(pymatgen=scaled_structure)


@calcfunction
def rescale(inp_structure, scale):
    """
    Rescales a crystal structure. Keeps the provenance in the database.

    :param inp_structure, a StructureData node (pk, or uuid)
    :param scale, float scaling factor for the cell

    :returns: New StructureData node with rescalled structure, which is linked to input Structure
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
        return {
            'conv_structure': conv_structure,
            'explicit_kpoints': explicit_kpoints,
            'parameters': parameters,
            'primitive_structure': primitive_structure
        }

    return primitive_structure
