#!/usr/bin/env python

import pytest
from ..dbsetup import *
from ..conftest import kkrimp_local_code, data_dir, import_with_migration
from aiida.engine import run_get_node


@pytest.mark.timeout(1200, method='thread')
def test_kkrimp_sub_wc(clear_database_before_test, kkrimp_local_code, enable_archive_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
    from numpy import array
    from aiida_kkr.tools import neworder_potential_wf
    from numpy import loadtxt
    from aiida.manage.caching import enable_caching

    wfd = kkr_imp_sub_wc.get_wf_defaults()

    wfd['nsteps'] = 20
    wfd['strmix'] = 0.05
    # deactivate final cleanup to be able to use caching
    wfd['do_final_cleanup'] = False

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

    # import previous GF writeout
    import_with_migration('files/db_dump_kkrflex_create.tar.gz')
    GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

    # now create a SingleFileData node containing the impurity starting potential
    with GF_host_calc.outputs.retrieved.open('scoef') as _f:
        neworder_pot1 = [int(i) for i in loadtxt(_f, skiprows=1)[:, 3] - 1]
    settings_dict = {'pot1': 'out_potential', 'out_pot': 'potential_imp', 'neworder': neworder_pot1}
    settings = Dict(settings_dict)
    with enable_caching():
        startpot_imp_sfd = neworder_potential_wf(
            settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder
        )

    label = 'kkrimp_scf Cu host_in_host'
    descr = 'kkrimp_scf workflow for Cu bulk'

    # create process builder to set parameters
    builder = kkr_imp_sub_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkrimp = kkrimp_local_code
    builder.options = options
    builder.remote_data = GF_host_calc.outputs.remote_folder
    builder.wf_parameters = Dict(wfd)
    builder.host_imp_startpot = startpot_imp_sfd

    print('builder', builder)

    # now run calculation
    with enable_archive_cache(data_dir / 'kkrimp_sub_wc.aiida'):
        out, node = run_get_node(builder)
    print('out', out)
    print('node', node)
    print(node.process_status)

    n = out['workflow_info']
    n = n.get_dict()

    assert n.get('successful')
    assert n.get('convergence_reached')


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    test_kkrimp_sub_wc()
