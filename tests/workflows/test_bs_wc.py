#!/usr/bin/env python
# coding: utf-8
"""
Tests for the bandstructure workflow
"""

from __future__ import absolute_import
from __future__ import print_function
import pytest
from ..dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
from ..conftest import voronoi_local_code, kkrhost_local_code, data_dir
from ..conftest import import_with_migration
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test
from six.moves import range


@pytest.mark.timeout(240, method='thread')
def test_bs_wc_Cu(clear_database_before_test, kkrhost_local_code, run_with_cache, ndarrays_regression):
    """
    minimal bandstructure calualtion for Cu bulk
    """

    from aiida import get_version
    from aiida.orm import Code, load_node, Dict, StructureData, Computer
    from aiida.plugins import DataFactory
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.bs import kkr_bs_wc
    import numpy as np

    print(f'AiiDA version: {get_version()}')
    Dict = DataFactory('dict')
    StructureData = DataFactory('structure')

    # create workflow parameters
    wfbs = kkr_bs_wc.get_wf_defaults()
    wfbs['nepts'] = 12
    wfbs['emax'] = 5
    wfbs['emin'] = -10
    wfbs['RCLUSTZ'] = 2.3
    wfbs['tempr'] = 50.0
    params_bs = Dict(dict=wfbs)

    # for runing in local computer
    options2 = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    options = Dict(dict=options2)

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
    builder.options = options
    builder.remote_data = kkr_calc_remote

    # run the calculation using cached data is available
    out, _ = run_with_cache(builder, data_dir=data_dir)

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
    ndarrays_regression.check(check_dict)
