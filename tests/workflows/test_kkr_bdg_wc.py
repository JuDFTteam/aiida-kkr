#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Local (offline) tests for the kkr_bdg_wc WorkChain.

Tests the spec, calcfunctions, default settings, and outline branching logic
without any daemon, real KKR calculations, or cluster access.

Tests are divided into:
  - DB-free tests: spec validation, defaults, outline logic (no aiida_profile needed)
  - DB tests: calcfunction logic (need aiida_profile for Dict node storage)
"""

import pytest
from unittest.mock import MagicMock
from plumpy.utils import AttributesFrozendict


# ═══════════════════════════════════════════════════════════════════════════════
# Group 1 — Spec validation (no DB needed — just class introspection)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBdgSpec:
    """Verify that kkr_bdg_wc spec is correctly defined."""

    def test_spec_required_inputs(self):
        """kkr, kkr_bdg, calc_parameters must be required; others optional."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        spec = kkr_bdg_wc.spec()

        # Required inputs
        for name in ('kkr', 'kkr_bdg', 'calc_parameters'):
            port = spec.inputs[name]
            assert port.required, f"'{name}' should be required"

        # Optional inputs
        for name in ('voronoi', 'structure', 'remote_data_normal',
                      'remote_data_semi_circle', 'semi_circle_settings',
                      'bdg_init_settings', 'bdg_settings', 'options',
                      'wf_parameters'):
            port = spec.inputs[name]
            assert not port.required, f"'{name}' should be optional"

    def test_spec_input_types(self):
        """All inputs should have the correct valid_type."""
        from aiida import orm
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        spec = kkr_bdg_wc.spec()

        def _check_valid_type(port, expected_type):
            """Check that a port's valid_type includes the expected type.
            Optional ports may have valid_type as a tuple including NoneType."""
            vt = port.valid_type
            if isinstance(vt, tuple):
                assert expected_type in vt, (
                    f"Expected {expected_type} in valid_type tuple {vt}"
                )
            else:
                assert issubclass(vt, expected_type), (
                    f"Expected subclass of {expected_type}, got {vt}"
                )

        # Code inputs (required)
        _check_valid_type(spec.inputs['kkr'], orm.Code)
        _check_valid_type(spec.inputs['kkr_bdg'], orm.Code)

        # Dict inputs
        _check_valid_type(spec.inputs['calc_parameters'], orm.Dict)

        # Optional inputs (valid_type may be tuple with NoneType)
        _check_valid_type(spec.inputs['structure'], orm.StructureData)
        _check_valid_type(spec.inputs['remote_data_normal'], orm.RemoteData)
        _check_valid_type(spec.inputs['remote_data_semi_circle'], orm.RemoteData)

    def test_spec_outputs(self):
        """All three declared outputs must exist with correct valid_type."""
        from aiida import orm
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        spec = kkr_bdg_wc.spec()

        assert 'output_kkr_bdg_wc_ParameterResults' in spec.outputs
        assert spec.outputs['output_kkr_bdg_wc_ParameterResults'].valid_type is orm.Dict

        assert 'last_RemoteData' in spec.outputs
        assert spec.outputs['last_RemoteData'].valid_type is orm.RemoteData

        assert 'last_InputParameters' in spec.outputs
        assert spec.outputs['last_InputParameters'].valid_type is orm.Dict

    def test_spec_exit_codes(self):
        """Exit codes 301–304 must be defined with correct labels."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        spec = kkr_bdg_wc.spec()
        exit_codes = spec.exit_codes

        # Exit codes are keyed by name (str), so look them up by name
        expected = {
            'ERROR_NORMAL_SCF_FAILED': 301,
            'ERROR_SEMI_CIRCLE_SCF_FAILED': 302,
            'ERROR_BDG_INIT_FAILED': 303,
            'ERROR_BDG_SCF_FAILED': 304,
        }
        for label, expected_status in expected.items():
            assert label in exit_codes, f"Exit code '{label}' not defined"
            assert exit_codes[label].status == expected_status, (
                f"Exit code '{label}' has status {exit_codes[label].status}, "
                f"expected {expected_status}"
            )

    def test_spec_outline_has_all_steps(self):
        """The outline should reference all expected step methods."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        # Verify the class has all the step methods we expect
        expected_steps = [
            'start', 'validate_inputs',
            'should_run_normal_scf', 'run_normal_scf', 'check_normal_scf',
            'should_run_semi_circle', 'run_semi_circle_scf', 'check_semi_circle_scf',
            'run_bdg_init', 'check_bdg_init',
            'run_bdg_scf', 'check_bdg_scf',
            'results',
        ]
        for step in expected_steps:
            assert hasattr(kkr_bdg_wc, step), (
                f"WorkChain missing step method '{step}'"
            )
            assert callable(getattr(kkr_bdg_wc, step)), (
                f"'{step}' should be callable"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Group 2 — Default settings (no DB needed)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBdgDefaults:
    """Verify the default workflow settings are correct."""

    def test_get_wf_defaults_returns_two_dicts(self):
        """get_wf_defaults() should return a 2-tuple of dicts."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        result = kkr_bdg_wc.get_wf_defaults(silent=True)
        assert isinstance(result, tuple) and len(result) == 2
        wf_defaults, options_defaults = result
        assert isinstance(wf_defaults, dict)
        assert isinstance(options_defaults, dict)

    def test_wf_defaults_structure(self):
        """wf_defaults should have three sub-dicts: semi_circle, bdg_init, bdg_scf."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        wf_defaults, _ = kkr_bdg_wc.get_wf_defaults(silent=True)

        assert 'semi_circle' in wf_defaults
        assert 'bdg_init' in wf_defaults
        assert 'bdg_scf' in wf_defaults

    def test_semi_circle_defaults(self):
        """Semi-circle defaults should have the key physics settings."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        sc = kkr_bdg_wc.get_wf_defaults(silent=True)[0]['semi_circle']

        assert sc['USE_SEMI_CIRCLE_CONTOUR'] is True
        assert sc['NPT1'] == 32
        assert sc['BZDIVIDE'] == [100, 100, 100]
        assert sc['NSTEPS'] == 200
        assert sc['DECOUPLE_SPIN_CHEBY'] is True
        assert 'NEWSOSOL' in sc['RUNOPT']

    def test_bdg_init_defaults(self):
        """BdG init defaults should have NSTEPS=1 and use_BdG=True."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        bdg_init = kkr_bdg_wc.get_wf_defaults(silent=True)[0]['bdg_init']

        assert bdg_init['NSTEPS'] == 1
        assert bdg_init['use_BdG'] is True
        assert bdg_init['Delta_BdG'] == 5e-4

    def test_bdg_scf_defaults(self):
        """BdG SCF defaults should have NSTEPS=200 and mixing settings."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        bdg_scf = kkr_bdg_wc.get_wf_defaults(silent=True)[0]['bdg_scf']

        assert bdg_scf['NSTEPS'] == 200
        assert bdg_scf['use_BdG'] is True
        assert bdg_scf['lambda_BdG'] == 0.025
        assert bdg_scf['mixfac_BdG'] == 0.3

    def test_options_defaults(self):
        """Options defaults should have standard compute keys."""
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        _, options = kkr_bdg_wc.get_wf_defaults(silent=True)

        assert options['resources'] == {'num_machines': 1}
        assert options['withmpi'] is True
        assert options['max_wallclock_seconds'] == 60 * 60 * 4


# ═══════════════════════════════════════════════════════════════════════════════
# Group 3 — Outline & step logic (mock-based, no DB needed)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBdgOutlineLogic:
    """Test the branching/conditional logic in the workchain steps."""

    def _make_wc_instance(self):
        """
        Create a bare workchain instance for testing step methods.

        Uses __new__ to skip __init__, then sets up the backing stores
        for the `ctx` property (_context) and `inputs` property
        (_parsed_inputs) plus other mocked attributes.
        """
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        wc = kkr_bdg_wc.__new__(kkr_bdg_wc)
        # Set backing stores for ctx and inputs properties
        wc._context = AttributesFrozendict()
        wc._parsed_inputs = AttributesFrozendict()  # empty → 'x' in wc.inputs is False
        # Mock report and provide exit_codes
        wc.report = MagicMock()
        wc.exit_codes = kkr_bdg_wc.spec().exit_codes
        return wc

    def test_should_run_normal_scf_true_when_no_remote(self):
        """should_run_normal_scf() returns True when ctx.current_remote is None."""
        wc = self._make_wc_instance()
        wc.ctx.current_remote = None
        assert wc.should_run_normal_scf() is True

    def test_should_run_normal_scf_false_when_remote_provided(self):
        """should_run_normal_scf() returns False when ctx.current_remote is set."""
        wc = self._make_wc_instance()
        wc.ctx.current_remote = MagicMock()  # simulate a RemoteData
        assert wc.should_run_normal_scf() is False

    def test_should_run_semi_circle_true_when_no_bypass(self):
        """
        should_run_semi_circle() returns True when 'remote_data_semi_circle'
        is NOT in inputs.
        """
        wc = self._make_wc_instance()
        # MagicMock(spec=[]) has no attributes →
        # 'remote_data_semi_circle' in wc.inputs → False
        assert wc.should_run_semi_circle() is True

    def test_should_run_semi_circle_false_when_bypass(self):
        """
        should_run_semi_circle() returns False when 'remote_data_semi_circle'
        IS in inputs.
        """
        wc = self._make_wc_instance()
        # Set inputs backing store with remote_data_semi_circle present
        wc._parsed_inputs = AttributesFrozendict(
            {'remote_data_semi_circle': MagicMock()}
        )
        assert wc.should_run_semi_circle() is False

    def test_validate_inputs_missing_structure_returns_error(self):
        """
        When starting from scratch (no remote_data_*) and no structure is
        provided, validate_inputs() should return ERROR_NORMAL_SCF_FAILED (301).
        """
        wc = self._make_wc_instance()

        # Simulate: codes exist with correct plugin, but no remote/structure
        mock_code = MagicMock()
        mock_code.get_input_plugin_name.return_value = 'kkr.kkr'

        # Set inputs with kkr/kkr_bdg but no remote_data or structure
        wc._parsed_inputs = AttributesFrozendict(
            {'kkr': mock_code, 'kkr_bdg': mock_code}
        )

        result = wc.validate_inputs()
        assert result is not None, "validate_inputs should return an exit code"
        assert result.status == 301, (
            f"Expected exit code 301 (ERROR_NORMAL_SCF_FAILED), got {result.status}"
        )

    def test_check_normal_scf_returns_error_on_failure(self):
        """check_normal_scf should return ERROR_NORMAL_SCF_FAILED when sub-wc fails."""
        wc = self._make_wc_instance()
        mock_scf = MagicMock()
        mock_scf.is_finished_ok = False
        wc.ctx.normal_scf = mock_scf

        result = wc.check_normal_scf()
        assert result is not None
        assert result.status == 301

    def test_check_normal_scf_updates_context_on_success(self):
        """check_normal_scf should update ctx on success."""
        wc = self._make_wc_instance()
        mock_scf = MagicMock()
        mock_scf.is_finished_ok = True
        mock_remote = MagicMock(name='remote')
        mock_params = MagicMock(name='params')
        mock_scf.outputs.last_RemoteData = mock_remote
        mock_scf.outputs.last_InputParameters = mock_params
        wc.ctx.normal_scf = mock_scf

        result = wc.check_normal_scf()
        assert result is None  # no exit code on success
        assert wc.ctx.current_remote is mock_remote
        assert wc.ctx.current_params is mock_params

    def test_check_semi_circle_scf_returns_error_on_failure(self):
        """check_semi_circle_scf should return ERROR_SEMI_CIRCLE_SCF_FAILED."""
        wc = self._make_wc_instance()
        wc.ctx.semi_circle_scf = MagicMock()
        wc.ctx.semi_circle_scf.is_finished_ok = False

        result = wc.check_semi_circle_scf()
        assert result is not None
        assert result.status == 302

    def test_check_bdg_init_returns_error_on_failure(self):
        """check_bdg_init should return ERROR_BDG_INIT_FAILED."""
        wc = self._make_wc_instance()
        wc.ctx.bdg_init = MagicMock()
        wc.ctx.bdg_init.is_finished_ok = False

        result = wc.check_bdg_init()
        assert result is not None
        assert result.status == 303

    def test_check_bdg_scf_returns_error_on_failure(self):
        """check_bdg_scf should return ERROR_BDG_SCF_FAILED and track error."""
        wc = self._make_wc_instance()
        wc.ctx.successful = True
        wc.ctx.errors = []
        wc.ctx.bdg_scf = MagicMock()
        wc.ctx.bdg_scf.is_finished_ok = False

        result = wc.check_bdg_scf()
        assert result is not None
        assert result.status == 304
        assert wc.ctx.successful is False
        assert len(wc.ctx.errors) == 1

    def test_check_bdg_scf_updates_context_on_success(self):
        """check_bdg_scf should update ctx.current_remote and params on success."""
        wc = self._make_wc_instance()
        wc.ctx.successful = True
        wc.ctx.errors = []
        mock_bdg = MagicMock()
        mock_bdg.is_finished_ok = True
        mock_remote = MagicMock(name='remote')
        mock_params = MagicMock(name='params')
        mock_bdg.outputs.remote_folder = mock_remote
        mock_bdg.inputs.parameters = mock_params
        wc.ctx.bdg_scf = mock_bdg

        result = wc.check_bdg_scf()
        assert result is None
        assert wc.ctx.current_remote is mock_remote
        assert wc.ctx.current_params is mock_params


# ═══════════════════════════════════════════════════════════════════════════════
# Group 4 — Calcfunction logic (needs AiiDA profile for Dict node storage)
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.usefixtures('aiida_profile')
class TestBdgCalcfunctions:
    """
    Test the @calcfunction helpers used by kkr_bdg_wc.

    These tests require an AiiDA profile (via aiida_profile fixture) because
    @calcfunction needs to store Dict nodes in the database.
    """

    def test_update_params_semi_circle_strips_scf_keys(self):
        """
        update_params_semi_circle should:
        - Strip SCF-inherited keys (STRMIX, BRYMIX, QBOUND, HFIELD, LINIPOL)
        - Remove TEMPR, NPT2, NPT3, NPOL
        - Set semi-circle settings
        - Preserve unregistered keys like NPT1
        """
        from aiida import orm
        from aiida_kkr.workflows.kkr_bdg_wc import update_params_semi_circle
        from masci_tools.io.kkr_params import kkrparams

        # Build a parent parameter set mimicking post-SCF output
        parent_dict = kkrparams(LMAX=2, NSPIN=1, RMAX=10.0, GMAX=100.0).get_dict()
        # Add SCF-inherited keys that should be stripped
        parent_dict['STRMIX'] = 0.03
        parent_dict['BRYMIX'] = 0.05
        parent_dict['QBOUND'] = 1e-7
        parent_dict['TEMPR'] = 800.0
        parent_dict['NPT2'] = 50

        params_node = orm.Dict(dict=parent_dict)

        semi_settings = orm.Dict(dict={
            'USE_SEMI_CIRCLE_CONTOUR': True,
            'IM_E_CIRC_MIN': 5e-5,
            'NPT1': 32,
            'NSTEPS': 200,
            'IMIX': 4,
        })

        result = update_params_semi_circle(params_node, semi_settings)
        result_dict = result.get_dict()

        # SCF-inherited keys should be gone
        for key in ('STRMIX', 'BRYMIX', 'QBOUND'):
            assert key not in result_dict or result_dict[key] is None, (
                f"SCF key '{key}' should have been stripped"
            )

        # TEMPR should be removed
        assert result_dict.get('TEMPR') is None, "TEMPR should have been removed"

        # Semi-circle settings should be present
        assert result_dict.get('NSTEPS') == 200
        assert result_dict.get('IMIX') == 4

        # NPT1 is an unregistered key, should be preserved
        assert result_dict.get('NPT1') == 32, "Unregistered key NPT1 should be preserved"

    def test_update_params_bdg_preserves_parent_keys(self):
        """
        update_params_bdg should:
        - Preserve all parent keys including bracket keys
        - Merge BdG settings on top (BdG takes precedence)
        """
        from aiida import orm
        from aiida_kkr.workflows.kkr_bdg_wc import update_params_bdg

        parent_dict = {
            'LMAX': 2,
            'NSPIN': 1,
            'NSTEPS': 200,
            '<USE_SEMI_CIRCLE_CONTOUR>': True,
            '<DECOUPLE_SPIN_CHEBY>': True,
            'BZDIVIDE': [100, 100, 100],
            'RUNOPT': ['NEWSOSOL'],
        }
        params_node = orm.Dict(dict=parent_dict)

        bdg_settings = orm.Dict(dict={
            'NSTEPS': 1,
            'use_BdG': True,
            'Delta_BdG': 5e-4,
        })

        result = update_params_bdg(params_node, bdg_settings)
        result_dict = result.get_dict()

        # BdG setting should override NSTEPS
        assert result_dict['NSTEPS'] == 1, "NSTEPS should be overridden to 1"

        # Parent bracket keys should be preserved
        assert result_dict.get('<USE_SEMI_CIRCLE_CONTOUR>') is True, (
            "Bracket key <USE_SEMI_CIRCLE_CONTOUR> should be preserved"
        )

        # BZDIVIDE and RUNOPT should survive (dict-level merge)
        assert result_dict.get('BZDIVIDE') == [100, 100, 100], (
            "BZDIVIDE should be preserved"
        )
        assert result_dict.get('RUNOPT') == ['NEWSOSOL'], (
            "RUNOPT should be preserved"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Group 5 — Entry point (needs installed package)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBdgEntryPoint:
    """Verify the kkr.bdg entry point is registered correctly."""

    def test_bdg_entry_point(self):
        """
        WorkflowFactory('kkr.bdg') must resolve to kkr_bdg_wc.

        NOTE: This test only passes when aiida-kkr is pip-installed (not just
        on PYTHONPATH), because entry points require package registration.
        If it fails with MissingEntryPointError, run: pip install -e .
        """
        from aiida_kkr.workflows.kkr_bdg_wc import kkr_bdg_wc

        try:
            from aiida.plugins import WorkflowFactory
            wf = WorkflowFactory('kkr.bdg')
            assert wf is kkr_bdg_wc
        except Exception as e:
            if 'MissingEntryPointError' in type(e).__name__:
                pytest.skip(
                    "Entry point 'kkr.bdg' not found — aiida-kkr may need "
                    "reinstall: pip install -e ."
                )
            raise
