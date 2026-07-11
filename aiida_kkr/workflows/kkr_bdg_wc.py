#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
WorkChain for automated Bogoliubov-de Gennes (BdG) calculations for host materials
using the AiiDA-KKR plugin.

Sequence:
    1. Normal KKR SCF      (optional, bypass with `remote_data_normal`)
    2. Semi-circle contour restart calculation
    3. BdG initialization  (one-shot raw KkrCalculation, NSTEPS=1)
    4. Final BdG SCF       (full self-consistency of the anomalous density)
"""

from aiida import orm
from aiida.engine import WorkChain, ToContext, if_, calcfunction
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_inputs_kkr
from aiida_kkr.tools.save_output_nodes import create_out_dict_node
from aiida_kkr.calculations.kkr import KkrCalculation
from masci_tools.io.kkr_params import kkrparams

__copyright__ = (u'Copyright (c), 2026, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.2.3'
# Changelog:
#   0.2.3 — cell-aware BZDIVIDE: the semi-circle step now inherits the parent (normal-SCF) k-mesh
#           instead of a hardcoded [100,100,100]. A too-dense mesh overruns the compiled KPOIBZ once
#           BdG lowers the symmetry (250001 irreducible k-points -> abort). Generalises the 0.2.1
#           RCLUSTZ inheritance to both cell-tuned keys. See update_params_semi_circle.
#   0.2.2 — NSPIN-aware spin coupling: the BdG steps now drop DECOUPLE_SPIN_CHEBY whenever
#           NSPIN=2 (BdG pairing couples spin up/down, so decoupled Chebyshev spin channels are
#           invalid and the KKRhost BdG solver hard-errors). NSPIN=1 is unaffected. See
#           update_params_bdg.
#   0.2.1 — semi-circle RCLUSTZ is now cell-aware: inherited from calc_parameters/normal-SCF
#           params instead of a hardcoded 3.5 default (explicit semi_circle_settings['RCLUSTZ']
#           still overrides). Prevents oversized screening clusters exceeding NACLSD on
#           multi-atom cells. See update_params_semi_circle. Also: rank/energy-grid sanity
#           WARNING in validate_inputs (MPI ranks must divide the semi-circle NPT1 grid).
__contributors__ = u'Philipp Rüßmann, Mohammad Hemmati'


# ── Helper calcfunctions (maintain provenance) ────────────────────────────────

@calcfunction
def update_params_semi_circle(params_node, semi_circle_settings):
    _unregistered_keys = ['NPT1']
    # Strip SCF-specific keys inherited from kkr_scf_wc convergence loop.
    # The tutorial never sets these — KKR uses its own internal defaults.
    _scf_inherited_keys = ['STRMIX', 'BRYMIX', 'QBOUND', 'HFIELD', 'LINIPOL']
    clean_dict = {k: v for k, v in params_node.get_dict().items()
                  if v is not None and k not in _scf_inherited_keys}
    para = kkrparams(**clean_dict)
    for k in ['TEMPR', 'NPT2', 'NPT3', 'NPOL']:
        try:
            para.remove_value(k)
        except KeyError:
            pass
    settings_dict = semi_circle_settings.get_dict()
    # Cell-aware defaults: RCLUSTZ (screening-cluster radius) and BZDIVIDE (k-mesh) tuned for
    # 1-atom bulk are wrong for multi-atom cells. A too-large RCLUSTZ overruns the compiled
    # NACLSD (RCLUSTZ=3.5 -> 225-atom cluster, abort [302]); a too-dense BZDIVIDE overruns the
    # compiled KPOIBZ once BdG lowers the symmetry ([100,100,100] -> 250001 irreducible k-points
    # at the BdG step, abort [302]). If the caller did not explicitly set them, inherit the parent
    # (normal-SCF / calc_parameters) cell-appropriate values.
    parent_dict = params_node.get_dict()
    for _key in ('RCLUSTZ', 'BZDIVIDE'):
        if settings_dict.get(_key) is None:
            settings_dict.pop(_key, None)
            _pv = parent_dict.get(_key)
            if _pv is not None:
                settings_dict[_key] = _pv
    unregistered = {k: settings_dict.pop(k) for k in _unregistered_keys if k in settings_dict}
    para.set_multiple_values(**settings_dict)
    result_dict = para.get_dict()
    result_dict.update(unregistered)
    return orm.Dict(dict=result_dict)

@calcfunction
def update_params_bdg(params_node, bdg_settings):
    """
    Update KKR parameters for a BdG step.

    Uses a dict-level merge (not kkrparams reconstruction) to preserve ALL
    keys from the parent params — including RUNOPT, bracket keys
    (<USE_SEMI_CIRCLE_CONTOUR>, <DECOUPLE_SPIN_CHEBY>, …), and list values
    (BZDIVIDE) — which kkrparams(**parent_dict) would silently drop.

    Strategy:
      1. Start from the full parent dict (preserves everything)
      2. Run bdg_settings through an EMPTY kkrparams to normalize aliases
      3. Merge normalized BdG keys on top
    """
    # Step 1: preserve every key from the parent
    result_dict = {k: v for k, v in params_node.get_dict().items() if v is not None}

    # Step 2: normalize BdG aliases through an empty kkrparams instance
    para = kkrparams()
    para.set_multiple_values(**bdg_settings.get_dict())
    bdg_normalized = {k: v for k, v in para.get_dict().items() if v is not None}

    # Step 3: BdG settings take precedence over parent
    result_dict.update(bdg_normalized)

    # Step 4 — NSPIN-aware spin coupling: BdG pairing couples spin up/down, so decoupled
    # Chebyshev spin channels are physically invalid for NSPIN=2 (the KKRhost BdG solver
    # hard-errors: "BdG formalism for nspin=2 works only with coupled spin channels"). Drop the
    # decouple option for spin-polarised BdG so KKR uses its coupled default; NSPIN=1 is
    # unaffected (single channel — decoupling is a harmless no-op, as in the 1-atom Nb probe).
    nspin = result_dict.get('NSPIN', result_dict.get('<NSPIN>'))
    if nspin == 2:
        for _k in ('<DECOUPLE_SPIN_CHEBY>', 'DECOUPLE_SPIN_CHEBY', 'decouple_spin_cheby'):
            result_dict.pop(_k, None)

    return orm.Dict(dict=result_dict)
# ── WorkChain ─────────────────────────────────────────────────────────────────

class kkr_bdg_wc(WorkChain):
    """
    Workchain for performing a Bogoliubov-de Gennes (BdG) host calculation.

    :param kkr:                  (Code) KKRhost code for normal-state & semi-circle steps.
    :param kkr_bdg:              (Code) KKRhost BdG code for the superconducting solver.
    :param voronoi:              (Code, optional) Voronoi code (only needed from scratch).
    :param structure:            (StructureData, optional) Crystal structure.
    :param calc_parameters:      (Dict) Initial KKR parameters (LMAX, RMAX, …).
    :param remote_data_normal:   (RemoteData, optional) Bypass normal SCF.
    :param remote_data_semi_circle: (RemoteData, optional) Bypass normal + semi-circle SCF.
    :param semi_circle_settings: (Dict) Parameters injected for the semi-circle step.
    :param bdg_init_settings:    (Dict) Parameters for the one-shot BdG initialisation.
    :param bdg_settings:         (Dict) Parameters for the final BdG SCF.
    :param options:              (Dict) Computer options (queue, wallclock, …).
    :param wf_parameters:        (Dict) Workflow parameters forwarded to sub-workflows.

    :return output_kkr_bdg_wc_ParameterResults: (Dict) Summary of the workflow.
    :return last_RemoteData:     (RemoteData) Remote folder of the final BdG calculation.
    :return last_InputParameters: (Dict) Input parameters of the final BdG calculation.
    """

    _workflowversion = __version__
    _wf_label = 'kkr_bdg_wc'
    _wf_description = (
        'Workflow for a BdG KKR calculation starting either from a structure '
        '(with automatic Voronoi + normal SCF) or from a converged RemoteData node.'
    )

    # ── Default workflow settings ──────────────────────────────────────────────
    _wf_default = {
        # semi-circle contour defaults
        'semi_circle': {
            'USE_SEMI_CIRCLE_CONTOUR': True,
            'IM_E_CIRC_MIN': 5e-5,
            'NPT1': 32,
            'MAX_NUM_KMESH': 4,
            # RCLUSTZ and BZDIVIDE intentionally omitted: the semi-circle step now INHERITS these
            # cell-appropriate values from the parent (calc_parameters/normal SCF) — see
            # update_params_semi_circle. Pass them in semi_circle_settings to override.
            'NSTEPS': 200,
            'IMIX': 4,
            'DISABLE_CHARGE_NEUTRALITY': True,
            'RUNOPT': ['NEWSOSOL'],
            'R_LOG': 0.6,
            'NPAN_EQ': 7,
            'NPAN_LOG': 18,
            'NCHEB': 12,
            'DECOUPLE_SPIN_CHEBY': True,
        },
        # BdG one-shot init defaults
        'bdg_init': {
            'NSTEPS':              1,
            'use_BdG':             True,
            'Delta_BdG':           5e-4,
            'use_e_symm_BdG':      True,
            'at_scale_BdG':        [1.0],           # must match number of atoms
            'RUNOPT':              ['NEWSOSOL'],
            'DECOUPLE_SPIN_CHEBY': True,
        },
        # BdG SCF defaults
        'bdg_scf': {
            'NSTEPS': 200,
            'use_BdG': True,
            'lambda_BdG': 0.025,
            'mixfac_BdG': 0.3,
            'memlen_Broyden_BdG': 20,
            'Ninit_Broyden_BdG': 5,
            'NSIMPLEMIXFIRST': 25,
            'IMIX': 4,
            'DISABLE_CHARGE_NEUTRALITY': True,
            'use_e_symm_BdG': True,
            'FORCE_BZ_SYMM': False,
        },
    }

    _options_default = {
        'queue_name': '',
        'resources': {'num_machines': 1},
        'max_wallclock_seconds': 60 * 60 * 4,
        'withmpi': True,
        'custom_scheduler_commands': '',
    }

    @classmethod
    def get_wf_defaults(cls, silent=False):
        """
        Print and return the default workflow settings.

        :returns: (_wf_default dict, _options_default dict)
        """
        if not silent:
            print(f'Version of workflow: {cls._workflowversion}')
        return cls._wf_default.copy(), cls._options_default.copy()

    # ── Spec ──────────────────────────────────────────────────────────────────

    @classmethod
    def define(cls, spec):
        super(kkr_bdg_wc, cls).define(spec)

        # Codes
        spec.input('kkr', valid_type=orm.Code, required=True,
                   help='KKRhost code for the normal-state and semi-circle contour steps.')
        spec.input('kkr_bdg', valid_type=orm.Code, required=True,
                   help='KKRhost BdG code for the superconducting solver.')
        spec.input('voronoi', valid_type=orm.Code, required=False,
                   help='Voronoi code (required only when starting from scratch).')

        # Optional bypass inputs
        spec.input('remote_data_normal', valid_type=orm.RemoteData, required=False,
                   help='Converged RemoteData of the normal-state SCF (skips normal SCF step).')
        spec.input('remote_data_semi_circle', valid_type=orm.RemoteData, required=False,
                   help='Converged RemoteData of the semi-circle SCF (skips normal + semi-circle steps).')

        # Structure and parameters
        spec.input('structure', valid_type=orm.StructureData, required=False,
                   help='Crystal structure (required only when starting from scratch).')
        spec.input('calc_parameters', valid_type=orm.Dict, required=True,
                   help='Initial KKR parameters (LMAX, NSPIN, RMAX, GMAX, …).')

        # Settings dicts
        spec.input('semi_circle_settings', valid_type=orm.Dict, required=False,
                   default=lambda: orm.Dict(dict=cls._wf_default['semi_circle']),
                   help='Parameters injected for the semi-circle contour step.')
        spec.input('bdg_init_settings', valid_type=orm.Dict, required=False,
                   default=lambda: orm.Dict(dict=cls._wf_default['bdg_init']),
                   help='Parameters for the one-shot BdG initialisation (NSTEPS=1).')
        spec.input('bdg_settings', valid_type=orm.Dict, required=False,
                   default=lambda: orm.Dict(dict=cls._wf_default['bdg_scf']),
                   help='Parameters for the final BdG SCF convergence loop.')

        # Computer/workflow options
        spec.input('options', valid_type=orm.Dict, required=False,
                   default=lambda: orm.Dict(dict=cls._options_default),
                   help='Computer options (queue, wallclock, resources, …).')
        spec.input('wf_parameters', valid_type=orm.Dict, required=False,
                   help='Workflow parameters forwarded to kkr_scf_wc sub-workflows.')

        # Outputs
        spec.output('output_kkr_bdg_wc_ParameterResults', valid_type=orm.Dict, required=True,
                    help='Summary dictionary of the BdG workflow.')
        spec.output('last_RemoteData', valid_type=orm.RemoteData, required=True,
                    help='Remote folder of the final converged BdG calculation.')
        spec.output('last_InputParameters', valid_type=orm.Dict, required=True,
                    help='Input parameters used in the final BdG calculation.')

        # Outline
        spec.outline(
            cls.start,
            cls.validate_inputs,
            if_(cls.should_run_normal_scf)(
                cls.run_normal_scf,
                cls.check_normal_scf,
            ),
            if_(cls.should_run_semi_circle)(
                cls.run_semi_circle_scf,
                cls.check_semi_circle_scf,
            ),
            cls.run_bdg_init,
            cls.check_bdg_init,
            cls.run_bdg_scf,
            cls.check_bdg_scf,
            cls.results,
        )

        # Exit codes
        spec.exit_code(301, 'ERROR_NORMAL_SCF_FAILED',
                       message='The normal-state KKR SCF step failed.')
        spec.exit_code(302, 'ERROR_SEMI_CIRCLE_SCF_FAILED',
                       message='The semi-circle contour SCF step failed.')
        spec.exit_code(303, 'ERROR_BDG_INIT_FAILED',
                       message='The BdG one-shot initialisation failed.')
        spec.exit_code(304, 'ERROR_BDG_SCF_FAILED',
                       message='The final BdG SCF step failed.')

    # ── Workflow steps ────────────────────────────────────────────────────────

    def start(self):
        """Initialise context variables and parse compute options."""
        self.report(f'INFO: Started kkr_bdg_wc version {self._workflowversion}')
        self.ctx.current_remote = None
        self.ctx.current_params = self.inputs.calc_parameters
        self.ctx.successful = True
        self.ctx.errors = []

        # Parse compute options into ctx (mirrors pattern from dos.py / bs.py)
        options_dict = self.inputs.options.get_dict() if 'options' in self.inputs else {}
        if not options_dict:
            options_dict = self._options_default
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get('max_wallclock_seconds', self._options_default['max_wallclock_seconds'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)

        self.report(
            f'INFO: use the following settings:\n'
            f'withmpi: {self.ctx.withmpi}\n'
            f'Resources: {self.ctx.resources}\n'
            f'Walltime (s): {self.ctx.max_wallclock_seconds}\n'
            f'queue name: {self.ctx.queue}\n'
            f'description: {self.ctx.description_wf}\n'
            f'label: {self.ctx.label_wf}\n'
        )

    def validate_inputs(self):
        test_and_get_codenode(self.inputs.kkr, 'kkr.kkr', use_exceptions=True)
        test_and_get_codenode(self.inputs.kkr_bdg, 'kkr.kkr', use_exceptions=True)

        # Rank / energy-grid sanity (added 2026-07-07).
        # The semi-circle & BdG steps parallelise over NPT1 energy points, while the normal SCF
        # (kkr_scf_wc) uses its own (smaller) fine contour. MPI ranks must divide EVERY step's
        # energy-point count or KKR aborts "No rest ranks allowed". We can only see NPT1 here, so
        # warn (do not hard-fail, since the normal-step contour is not visible to this workchain).
        _res = self.ctx.resources or {}
        _nranks = _res.get('tot_num_mpiprocs') or (
            (_res.get('num_machines', 1) or 1) * (_res.get('num_mpiprocs_per_machine', 1) or 1))
        _npt1 = (self.inputs.semi_circle_settings.get_dict().get('NPT1')
                 if 'semi_circle_settings' in self.inputs else None)
        if _npt1 and _nranks and (_nranks > _npt1 or _npt1 % _nranks != 0):
            self.report(
                f'WARNING: MPI ranks ({_nranks}) do not divide the semi-circle NPT1 ({_npt1}); the '
                f'BdG steps may abort with "No rest ranks allowed". The normal SCF uses kkr_scf_wc\'s '
                f'own (smaller) contour, so ranks must be a COMMON divisor of BOTH grids — pick a '
                f'small divisor (e.g. 8).')

        if 'remote_data_semi_circle' in self.inputs:
            self.ctx.current_remote = self.inputs.remote_data_semi_circle
            self.report('INFO: Bypassing normal and semi-circle SCF, starting from remote_data_semi_circle.')
            # ← NEW: load full params from the parent semi-circle KkrCalculation
            try:
                parent_calc = self.inputs.remote_data_semi_circle.creator
                self.ctx.current_params = parent_calc.inputs.parameters
                self.report('INFO: Loaded calc_parameters from parent semi-circle KKR calculation.')
            except Exception as e:
                self.report(f'WARNING: Could not load params from parent semi-circle calc ({e}), '
                            f'falling back to input calc_parameters.')
    
        elif 'remote_data_normal' in self.inputs:
            self.ctx.current_remote = self.inputs.remote_data_normal
            self.report('INFO: Bypassing normal SCF, starting from remote_data_normal.')
            try:
                parent_calc = self.inputs.remote_data_normal.creator
                self.ctx.current_params = parent_calc.inputs.parameters
                self.report('INFO: Loaded calc_parameters from parent KKR calculation.')
            except Exception as e:
                self.report(f'WARNING: Could not load params from parent calc ({e}), '
                            f'falling back to input calc_parameters.')
    
        else:
            if 'structure' not in self.inputs:
                self.report('ERROR: `structure` is required when starting from scratch.')
                return self.exit_codes.ERROR_NORMAL_SCF_FAILED

    def should_run_normal_scf(self):
        """Return True if no starting remote data was provided."""
        return self.ctx.current_remote is None

    def run_normal_scf(self):
        """Submit the normal-state KKR SCF."""
        self.report('INFO: Submitting normal-state KKR SCF.')
        builder = kkr_scf_wc.get_builder()
        builder.kkr = self.inputs.kkr
        builder.calc_parameters = self.ctx.current_params
        builder.structure = self.inputs.structure
        if 'voronoi' in self.inputs:
            builder.voronoi = self.inputs.voronoi
        if 'options' in self.inputs:
            builder.options = self.inputs.options
        if 'wf_parameters' in self.inputs:
            builder.wf_parameters = self.inputs.wf_parameters
        return ToContext(normal_scf=self.submit(builder))

    def check_normal_scf(self):
        """Check normal SCF result and update context."""
        if not self.ctx.normal_scf.is_finished_ok:
            self.report('ERROR: normal_scf failed.')
            return self.exit_codes.ERROR_NORMAL_SCF_FAILED
        self.ctx.current_remote = self.ctx.normal_scf.outputs.last_RemoteData
        self.ctx.current_params = self.ctx.normal_scf.outputs.last_InputParameters

    def should_run_semi_circle(self):
        """Return True if a semi-circle SCF is still needed."""
        return 'remote_data_semi_circle' not in self.inputs

    def run_semi_circle_scf(self):
        """
        Submit the semi-circle contour SCF as a raw KkrCalculation.

        Using a raw KkrCalculation (not kkr_scf_wc) is essential here because
        kkr_scf_wc internally overrides the energy contour parameters (NPT1,
        NPT2, NPT3, NPOL, TEMPR) with its own convergence_setting_fine values,
        silently destroying the semi-circle contour settings (e.g. NPT1=32
        becomes 7, causing 'too many ranks' MPI errors with 32 processes).
        """
        self.report('INFO: Submitting semi-circle contour KKR SCF (raw KkrCalculation).')
        new_params = update_params_semi_circle(self.ctx.current_params,
                                               self.inputs.semi_circle_settings)
        options_dict = self.inputs.options.get_dict() if 'options' in self.inputs else {}
        inputs = get_inputs_kkr(
            code=self.inputs.kkr,
            remote=self.ctx.current_remote,
            options=options_dict,
            label='semi_circle_scf',
            description='Semi-circle contour KKR SCF (NEWSOSOL + BdG code)',
            parameters=new_params,
        )
        return ToContext(semi_circle_scf=self.submit(KkrCalculation, **inputs))

    def check_semi_circle_scf(self):
        """Check semi-circle SCF result and extract remote folder and parameters."""
        if not self.ctx.semi_circle_scf.is_finished_ok:
            self.report('ERROR: semi_circle_scf failed.')
            return self.exit_codes.ERROR_SEMI_CIRCLE_SCF_FAILED
        # Read outputs directly from the KkrCalculation (same pattern as check_bdg_init)
        self.ctx.current_remote = self.ctx.semi_circle_scf.outputs.remote_folder
        self.ctx.current_params = self.ctx.semi_circle_scf.inputs.parameters

    def run_bdg_init(self):
        """
        Submit the one-shot BdG initialisation using a raw KkrCalculation.

        Using a raw KkrCalculation (not kkr_scf_wc) ensures that NSTEPS=1 is
        respected exactly — kkr_scf_wc would silently override NSTEPS via its
        own convergence loop logic, causing the anomalous density to oscillate
        and collapse to zero.
        """
        self.report('INFO: Submitting BdG one-shot initialisation (NSTEPS=1).')
        new_params = update_params_bdg(self.ctx.current_params, self.inputs.bdg_init_settings)
        options_dict = self.inputs.options.get_dict() if 'options' in self.inputs else {}
        inputs = get_inputs_kkr(
            code=self.inputs.kkr_bdg,
            remote=self.ctx.current_remote,
            options=options_dict,
            label='BdG_init_1step',
            description='BdG one-shot initialisation (NSTEPS=1)',
            parameters=new_params,
        )
        return ToContext(bdg_init=self.submit(KkrCalculation, **inputs))

    def check_bdg_init(self):
        """Check BdG init and extract remote folder and parameters."""
        if not self.ctx.bdg_init.is_finished_ok:
            self.report('ERROR: bdg_init failed.')
            return self.exit_codes.ERROR_BDG_INIT_FAILED
        self.ctx.current_remote = self.ctx.bdg_init.outputs.remote_folder
        self.ctx.current_params = self.ctx.bdg_init.inputs.parameters

    def run_bdg_scf(self):
        """
        Submit the final BdG SCF as a raw KkrCalculation.
    
        kkr_scf_wc must NOT be used here: it internally resets NPT1, BZDIVIDE,
        and the energy contour to convergence_setting_fine values, destroying the
        semi-circle contour and causing 'too many ranks' MPI errors.
        The tutorial pattern is: raw KkrCalculation with NSTEPS=200.
        """
        self.report('INFO: Submitting final BdG SCF (raw KkrCalculation, NSTEPS=200).')
        new_params = update_params_bdg(self.ctx.current_params, self.inputs.bdg_settings)
        options_dict = self.inputs.options.get_dict() if 'options' in self.inputs else {}
        inputs = get_inputs_kkr(
            code=self.inputs.kkr_bdg,
            remote=self.ctx.current_remote,
            options=options_dict,
            label='BdG_scf',
            description='BdG SCF convergence (NSTEPS=200)',
            parameters=new_params,
        )
        return ToContext(bdg_scf=self.submit(KkrCalculation, **inputs))
    
    
    def check_bdg_scf(self):
        """Check final BdG SCF result."""
        if not self.ctx.bdg_scf.is_finished_ok:
            self.report('ERROR: bdg_scf failed.')
            self.ctx.successful = False
            self.ctx.errors.append('Final BdG SCF step failed.')
            return self.exit_codes.ERROR_BDG_SCF_FAILED
        self.ctx.current_remote = self.ctx.bdg_scf.outputs.remote_folder
        self.ctx.current_params = self.ctx.bdg_scf.inputs.parameters

    def results(self):
        """Collect outputs and publish them with full data provenance."""
        last_remote = self.ctx.bdg_scf.outputs.remote_folder
        last_params = self.ctx.bdg_scf.inputs.parameters
        last_output = self.ctx.bdg_scf.outputs.output_parameters
    
        outputnode_dict = {
            'workflow_name': self.__class__.__name__,
            'workflow_version': self._workflowversion,
            'successful': self.ctx.successful,
            'list_of_errors': self.ctx.errors,
            'withmpi': self.ctx.withmpi,
            'resources': self.ctx.resources,
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'queue_name': self.ctx.queue,
            'last_calc_pk': self.ctx.bdg_scf.pk,
        }
        outputnode = orm.Dict(dict=outputnode_dict)
        outputnode.label = 'kkr_bdg_wc_results'
        outputnode.description = 'Summary of the kkr_bdg_wc workflow.'
    
        result_node = create_out_dict_node(outputnode, last_output_parameters=last_output)
    
        self.out('output_kkr_bdg_wc_ParameterResults', result_node)
        self.out('last_RemoteData', last_remote)
        self.out('last_InputParameters', last_params)
        self.report('INFO: kkr_bdg_wc finished successfully.')