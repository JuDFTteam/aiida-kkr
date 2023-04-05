#!/usr/bin/env python
# coding: utf-8
"""
Tests for Jij workflow
"""

import pytest
from ..dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
from aiida.manage.tests.pytest_fixtures import (
    aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile, clear_database, clear_database_after_test,
    clear_database_before_test
)
from ..conftest import kkrhost_local_code, data_dir, import_with_migration


@pytest.mark.timeout(240, method='thread')
def test_jij(clear_database_before_test, kkrhost_local_code, run_with_cache, ndarrays_regression):
    """
    Jij test with SOC
    """

    from aiida.orm import load_node, Dict
    from aiida_kkr.workflows.jijs import kkr_jij_wc

    # for runing in local computer
    options = Dict({
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    })

    # import calculation which is used as parent calculation
    import_with_migration('files/jij_input.aiida')
    parent = load_node('d3e62972-2f38-470b-8bbf-300c21847a6a')

    # now set up process builder
    builder = kkr_jij_wc.get_builder()
    wfd = kkr_jij_wc.get_wf_defaults()
    wfd['jijrad_ang'] = 5.0
    builder.wf_parameters = Dict(wfd)
    builder.remote_data = parent.outputs.remote_folder
    builder.kkr = kkrhost_local_code
    builder.options = options

    # make test run faster
    builder.params_kkr_overwrite = Dict({
        'LMAX': 2,
        'BZDIVIDE': [10, 10, 10],
        'NPT1': 3,
        'NPT2': 10,
        'NPT3': 3,
    })

    # run the calculation using cached data is available
    out, _ = run_with_cache(builder, data_dir=data_dir)

    # check results
    print('check results', out)
    n = out['results_wf']
    assert n.get_dict().get('successful')
    # check Jij data arrays
    check_dict = {
        'Jij_expanded': out['jij_data'].get_array('Jij_expanded'),
        'positions_expanded': out['jij_data'].get_array('positions_expanded')
    }
    print(check_dict['Jij_expanded'][:10])
    ndarrays_regression.check(check_dict)


@pytest.mark.timeout(240, method='thread')
def test_jij_soc(clear_database_before_test, kkrhost_local_code, run_with_cache, ndarrays_regression):
    """
    Jij test with SOC
    """

    from aiida.orm import load_node, Dict
    from aiida_kkr.workflows.jijs import kkr_jij_wc

    # for runing in local computer
    options = Dict({
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 12 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    })

    # import calculation which is used as parent calculation
    import_with_migration('files/jij_input.aiida')
    parent = load_node('d3e62972-2f38-470b-8bbf-300c21847a6a')

    # now set up process builder
    builder = kkr_jij_wc.get_builder()
    wfd = kkr_jij_wc.get_wf_defaults()
    wfd['jijrad_ang'] = 5.0
    builder.wf_parameters = Dict(wfd)
    builder.remote_data = parent.outputs.remote_folder
    builder.kkr = kkrhost_local_code
    builder.options = options

    # make test run faster
    builder.params_kkr_overwrite = Dict({
        'RUNOPT': ['NEWSOSOL'],
        'LMAX': 2,
        'BZDIVIDE': [10, 10, 10],
        'NPT1': 3,
        'NPT2': 10,
        'NPT3': 3,
        'NPAN_LOG': 5,
        'NPAN_EQ': 12,
        'NCHEB': 12,
        'R_LOG': 0.6,
        '<USE_CHEBYCHEV_SOLVER>': True,
        '<SET_CHEBY_NOSOC>': True
    })

    # run the calculation using cached data is available
    out, _ = run_with_cache(builder, data_dir=data_dir)

    # check results
    print('check results', out)
    n = out['results_wf']
    assert n.get_dict().get('successful')
    # check Jij data arrays
    check_dict = {
        'Jij_expanded': out['jij_data'].get_array('Jij_expanded'),
        'positions_expanded': out['jij_data'].get_array('positions_expanded')
    }
    print(check_dict['Jij_expanded'][:10])
    ndarrays_regression.check(check_dict)
