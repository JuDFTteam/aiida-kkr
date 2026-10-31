#!/usr/bin/env python
# coding: utf-8
"""
Run options that a workflow does not set itself are passed on to the calculation
"""

from types import SimpleNamespace
from aiida.common import AttributeDict
from aiida_kkr.tools.common_workfunctions import get_extra_options

OPTIONS = {
    'withmpi': True,
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,
        'num_cores_per_mpiproc': 2
    },
    'max_wallclock_seconds': 3600,
    'queue_name': 'q',
    'custom_scheduler_commands': '',
    'environment_variables': {
        'OMP_NUM_THREADS': '2'
    },
    'prepend_text': 'module load intel',
    'account': 'proj',
}
EXTRA = {'environment_variables': {'OMP_NUM_THREADS': '2'}, 'prepend_text': 'module load intel', 'account': 'proj'}


def test_get_extra_options():
    assert get_extra_options(OPTIONS) == EXTRA
    assert get_extra_options(OPTIONS, ('prepend_text',)) == {
        k: v for k, v in EXTRA.items() if k != 'prepend_text'
    }
    # the five options every workflow sets itself give nothing extra, so default runs keep their options (and hashes)
    assert get_extra_options({k: OPTIONS[k] for k in list(OPTIONS)[:5]}) == {}


def test_gf_writeout_passes_extra_options(monkeypatch):
    """environment_variables etc. reach the KKRFLEX calculation"""
    import aiida_kkr.workflows.gf_writeout as gf_writeout

    captured = {}

    def fake_get_inputs_kkr(code, remote, options, *args, **kwargs):
        captured.update(options)
        return AttributeDict()

    monkeypatch.setattr(gf_writeout, 'get_inputs_kkr', fake_get_inputs_kkr)
    fake_self = SimpleNamespace(
        inputs=SimpleNamespace(kkr='code', remote_data='remote', impurity_info='imp_info'),
        ctx=AttributeDict(
            dict(
                flex_kkrparams='params',
                withmpi=OPTIONS['withmpi'],
                resources=OPTIONS['resources'],
                max_wallclock_seconds=OPTIONS['max_wallclock_seconds'],
                queue=OPTIONS['queue_name'],
                custom_scheduler_commands=OPTIONS['custom_scheduler_commands'],
                extra_options=get_extra_options(OPTIONS),
                retrieve_kkrflex=True,
            )
        ),
        report=print,
        submit=lambda *args,
        **kwargs: 'flexrun',
    )
    gf_writeout.kkr_flex_wc.get_flex(fake_self)
    assert captured == {
        'max_wallclock_seconds': 3600,
        'resources': OPTIONS['resources'],
        'queue_name': 'q',
        **EXTRA,
    }
