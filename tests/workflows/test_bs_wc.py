#!/usr/bin/env python
# coding: utf-8
"""
Tests for the bandstructure workflow
"""

import pytest
from ..dbsetup import *
from aiida_test_cache.archive_cache import enable_archive_cache
from aiida.engine import run_get_node
from ..conftest import voronoi_local_code, kkrhost_local_code, data_dir
from ..conftest import import_with_migration


@pytest.mark.timeout(240, method='thread')
def test_bs_wc_Cu(clear_database_before_test, kkrhost_local_code, enable_archive_cache, ndarrays_regression):
    """
    minimal bandstructure calculation for Cu bulk
    """

    from aiida import get_version
    from aiida.orm import Code, load_node, Dict, StructureData, Computer
    from aiida.plugins import DataFactory
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.bs import kkr_bs_wc
    import numpy as np

    print(f'AiiDA version: {get_version()}')
    Dict = DataFactory('core.dict')
    StructureData = DataFactory('core.structure')

    # create workflow parameters
    wfbs = kkr_bs_wc.get_wf_defaults()
    wfbs['nepts'] = 12
    wfbs['emax'] = 5
    wfbs['emin'] = -10
    wfbs['RCLUSTZ'] = 2.3
    wfbs['tempr'] = 50.0
    params_bs = Dict(wfbs)

    label = 'bs calc Cu bulk'
    descr = 'testing bs workflow for Cu bulk'

    # import calculation which is used as parent calculation
    import_with_migration('files/db_dump_bs/db_dump_kkrcalc_bs.tar.gz')
    kkr_calc_remote = load_node('d5782162-8393-4212-9340-c8ee8b725474').outputs.remote_folder

    # now set up process builder
    builder = kkr_bs_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkr = kkrhost_local_code
    builder.wf_parameters = params_bs
    builder.options = Dict({
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    })
    builder.remote_data = kkr_calc_remote

    # run the calculation using cached data is available
    with enable_archive_cache(data_dir / 'bs_wc_Cu.aiida'):
        out, node = run_get_node(builder)

    # check results
    print('check results', out)
    n = out['results_wf']
    n = n.get_dict()
    assert n.get('successful')
    assert n.get('list_of_errors') == []
    # check band structure data arrays
    check_dict = {
        'Kpts': out['BS_Data'].get_array('Kpts'),
        'energy_points': out['BS_Data'].get_array('energy_points'),
        'BlochSpectralFunction': out['BS_Data'].get_array('BlochSpectralFunction'),
    }
    print(check_dict)
    ndarrays_regression.check(check_dict)
