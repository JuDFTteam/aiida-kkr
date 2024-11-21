#!/usr/bin/env python

import pytest
from aiida_test_cache.archive_cache import enable_archive_cache
from aiida.engine import run_get_node
from ..dbsetup import *
from ..conftest import voronoi_local_code, kkrhost_local_code, data_dir, import_with_migration


@pytest.mark.timeout(240, method='thread')
def test_kkrflex_writeout_wc(clear_database_before_test, kkrhost_local_code, enable_archive_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData, RemoteData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
    from numpy import array
    import os

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd = kkr_flex_wc.get_wf_defaults()
    print(wfd)
    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'custom_scheduler_commands': '',
        'withmpi': False
    }
    options = Dict(options)

    imp_info = Dict({'Rcut': 2.5533, 'ilayer_center': 0, 'Zimp': [29.]})

    label = 'GF_writeout Cu bulk'
    descr = 'GF_writeout workflow for Cu bulk'

    import_with_migration('files/db_dump_kkrcalc.tar.gz')
    kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

    # create process builder to set parameters
    builder = kkr_flex_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkr = kkrhost_local_code
    builder.options = options
    builder.remote_data = kkr_calc_remote
    builder.impurity_info = imp_info

    # now run calculation
    with enable_archive_cache(data_dir / 'kkrflex_writeout_wc.aiida'):
        out, node = run_get_node(builder)
    print(out)

    n = out['workflow_info']
    n = n.get_dict()

    assert n.get('successful')
    assert n.get('list_of_errors') == []

    d = out['GF_host_remote']
    assert isinstance(d, RemoteData)

    kkrflex_calc = load_node(n.get('pk_flexcalc'))
    kkrflex_retrieved = kkrflex_calc.outputs.retrieved
    for name in 'tmat green atominfo intercell_cmoms intercell_ref'.split():
        assert 'kkrflex_' + name in kkrflex_retrieved.list_object_names()


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    test_kkrflex_writeout_wc()
