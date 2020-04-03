#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
if __name__!='__main__':
    import pytest
    from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
    from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test
    from ..conftest import kkrimp_local_code, kkrhost_local_code, test_dir, data_dir
from aiida.orm import load_node, Dict
from aiida_kkr.workflows import combine_imps_wc
from aiida.tools.importexport import import_data

def test_combine_imps(clear_database_before_test, kkrhost_local_code, kkrimp_local_code, run_with_cache, nopytest=False):
    """
    test for combine_imps_wc (place two imps next to each other)
    """
    # import single imp calculations
    imported_nodes = import_data(test_dir/'data_dir/kkr_imp_wc-nodes-2ae76fce3a2cdcded6341f368ef0f00b.tar.gz', silent=True)['Node']
    for _, pk in imported_nodes['new']+imported_nodes['existing']:
        node = load_node(pk)
        if node.label=='kkrimp_scf full Cu host_in_host':
            imp1 = node
    imp1_out = imp1.outputs.workflow_info
    imp2_out = imp1_out # use the same impurity and create a dimer

    #set up combine_imps_wc workflow
    builder = combine_imps_wc.get_builder()
    builder.impurity1_output_node = imp1_out
    builder.impurity2_output_node = imp2_out
    builder.offset_imp2 = Dict(dict={'index': 1})
    builder.metadata.label = 'test_combine_imps'
    # set codes
    builder.host_gf.kkr = kkrhost_local_code
    builder.scf.kkrimp = kkrimp_local_code
    # set computer options
    options = {'queue_name' : '', 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
    builder.host_gf.options = Dict(dict=options)
    builder.scf.options = builder.host_gf.options

    # now submit
    print(builder, type(builder))
    if not nopytest:
        out, node = run_with_cache(builder, data_dir=data_dir)
    else:
        from aiida.engine import run_get_node
        out, node = run_get_node(builder)
    print((out, node))

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert results['successful']
    assert results['convergence_reached']



# run manual:
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    from aiida.orm import Code
    import pathlib

    test_dir = pathlib.Path('.')
    #define codes
    try:
        # on mac
        kkrhost_local_code = Code.get_from_string('kkrhost_intel19@localhost')
        kkrimp_local_code = Code.get_from_string('kkrimp_intel19@localhost')
    except:
        # on iff desktop
        kkrhost_local_code = Code.get_from_string('kkrhost@localhost')
        kkrimp_local_code = Code.get_from_string('kkrimp@localhost')

    # run test
    test_combine_imps(None, kkrhost_local_code, kkrimp_local_code, None, nopytest=True)
