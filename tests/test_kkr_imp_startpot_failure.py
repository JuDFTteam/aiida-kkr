#!/usr/bin/env python
"""Failure routing in kkr_imp_wc when a sub-workflow did not finish (no database needed).

Regression guard: construct_startpot used to record an exit code in ctx.exit_code and then
fall through into code that assumes the sub-workflow succeeded, so the workchain excepted
with `ValueError: max() arg is an empty sequence` (or NotExistentAttributeError for the GF
writeout) instead of returning the exit code it had already prepared.
"""

from types import SimpleNamespace

from aiida_kkr.workflows.kkr_imp import kkr_imp_wc


def _outgoing_stub(nodes):
    """Stand-in for the result of Node.get_outgoing(...)."""
    return SimpleNamespace(all=lambda: nodes)


def _workchain_stub(voro_finished_ok=True, gf_finished_ok=True, do_gf_calc=True, exit_code=None):
    """Minimal stand-in for a kkr_imp_wc instance at the construct_startpot step."""
    last_voro_calc = SimpleNamespace(
        is_finished_ok=voro_finished_ok,
        pk=1234,
        get_outgoing=lambda **kwargs: _outgoing_stub([]),
    )
    gf_writeout = SimpleNamespace(is_finished_ok=gf_finished_ok, pk=5678)
    return SimpleNamespace(
        ctx=SimpleNamespace(
            exit_code=exit_code,
            do_gf_calc=do_gf_calc,
            last_voro_calc=last_voro_calc,
            gf_writeout=gf_writeout,
        ),
        inputs={},  # only ever probed with `'startpot' in self.inputs`
        exit_codes=kkr_imp_wc.exit_codes,
        report=lambda message: None,
    )


def test_construct_startpot_returns_145_when_voronoi_failed():
    """A failed kkr_startpot_wc ends the workchain with 145 instead of raising.

    What this asserts is the early return, not any particular downstream exception. The
    observed production failure was `ValueError: max() arg is an empty sequence` at the
    `max([i.node.pk for i in all_nodes])` line, and the stub's get_outgoing returns no
    calculations to match that state -- but with do_gf_calc set, an unguarded run now dies
    one branch earlier, on the GF writeout outputs. Either way it must not get that far.
    """
    stub = _workchain_stub(voro_finished_ok=False)
    result = kkr_imp_wc.construct_startpot(stub)

    assert result.status == 145
    assert stub.ctx.exit_code.status == 145


def test_construct_startpot_returns_146_when_gf_writeout_failed():
    """A failed kkr_flex_wc ends the workchain with 146 instead of raising.

    The stub's gf_writeout carries no outputs attribute, so reaching for
    `.outputs.workflow_info` raises AttributeError here, standing in for the
    NotExistentAttributeError the real workchain raised on a missing workflow_info output.
    """
    stub = _workchain_stub(gf_finished_ok=False)
    result = kkr_imp_wc.construct_startpot(stub)

    assert result.status == 146
    assert stub.ctx.exit_code.status == 146


def test_bail_on_error_returns_a_stored_exit_code():
    """An exit code set by input validation stops the workchain before any work is done."""
    stub = _workchain_stub(exit_code=kkr_imp_wc.exit_codes.ERROR_MISSING_REMOTE)

    assert kkr_imp_wc.bail_on_error(stub).status == 143


def test_bail_on_error_is_transparent_on_the_happy_path():
    """With no exit code recorded, the step returns None and the outline carries on."""
    assert kkr_imp_wc.bail_on_error(_workchain_stub()) is None


def test_has_starting_potential_input_ignores_a_stored_exit_code():
    """The predicate decides only whether a startpotential must be built.

    It used to return True when ctx.exit_code was set, which *entered* the if_ block and ran
    the voronoi step on inputs already known to be invalid; bail_on_error now handles that
    case before the predicate is reached.
    """
    stub = _workchain_stub(exit_code=kkr_imp_wc.exit_codes.ERROR_MISSING_REMOTE)

    assert kkr_imp_wc.has_starting_potential_input(stub) is True
    assert stub.ctx.create_startpot is True


# run test manually
if __name__ == '__main__':
    test_construct_startpot_returns_145_when_voronoi_failed()
    test_construct_startpot_returns_146_when_gf_writeout_failed()
    test_bail_on_error_returns_a_stored_exit_code()
    test_bail_on_error_is_transparent_on_the_happy_path()
    test_has_starting_potential_input_ignores_a_stored_exit_code()
    print('ok')
