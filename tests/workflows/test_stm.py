#!/usr/bin/env python

import pytest
from ..dbsetup import *
from aiida.engine import run_get_node
from ..conftest import voronoi_local_code, kkrhost_local_code, test_dir, data_dir, import_with_migration


@pytest.mark.timeout(900, method='thread')
def test_stm_wc(
    clear_database_before_test, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, enable_archive_cache,
    ndarrays_regression
):
    """
    STM workflow test for a simple Cu host in host impurity adding a minimal scanning region
    """
    from aiida.orm import Code, load_node, Dict, StructureData, load_group
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows import kkr_STM_wc
    from numpy import array

    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    options = Dict(options)

    # import parent calculation (converged host system)
    group_pk = import_with_migration('data_dir/kkrimp_full_wc.aiida')
    kkr_imp_scf = [n for n in load_group(group_pk).nodes if n.label == 'kkrimp_scf full Cu host_in_host'][0]

    # create process builder to set parameters
    builder = kkr_STM_wc.get_builder()
    builder.metadata.label = 'stm test'
    builder.kkrimp = kkrimp_local_code
    builder.voronoi = voronoi_local_code
    builder.kkr = kkrhost_local_code
    builder.options = options

    builder.host_remote = kkr_imp_scf.inputs.remote_data_host
    builder.imp_info = kkr_imp_scf.inputs.impurity_info
    builder.imp_potential_node = kkr_imp_scf.outputs.converged_potential

    builder.tip_position = Dict({'ilayer': 0, 'nx': 2, 'ny': 2})

    builder.wf_parameters = Dict({
        'jij_run': False,
        'lmdos': True,
        'retrieve_kkrflex': True,
        'dos_params': {
            'nepts': 7,
            'tempr': 200.0,
            'emin': -1.0,
            'emax': 1.0,
            'kmesh': [5, 5, 5]
        }
    })

    # now run calculation
    with enable_archive_cache(data_dir / 'stm.aiida'):
        out, node = run_get_node(builder)
    print(out)
    print(list(node.called))

    # check outcome
    assert 'STM_dos_data_lmdos' in out

    # check dos data
    check_dict = {
        'x': out['STM_dos_data_lmdos'].get_x()[1],
        'y': out['STM_dos_data_lmdos'].get_y()[5][1],
    }
    print(check_dict)
    ndarrays_regression.check(check_dict)

    # print('x', out['STM_dos_data_lmdos'].get_x()[1])
    # print('y', out['STM_dos_data_lmdos'].get_y()[5][1])
