#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
if __name__ != '__main__':
    import pytest
    from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
    from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test
    from ..conftest import kkrimp_local_code, kkrhost_local_code, test_dir, data_dir
from aiida.orm import load_node, Dict
from aiida_kkr.workflows import combine_imps_wc
from ..conftest import import_with_migration


def write_graph(node, label=''):
    #if create_graph_file:
    from aiida.tools.visualization import Graph
    builder = node.get_builder_restart()
    name = str(builder.process_class).split('.')[-1].strip("'>")
    name += label
    graph = Graph(engine='dot', node_id_type='uuid')
    graph.recurse_ancestors(
        node, depth=None, link_types=(), annotate_links='both', include_process_outputs=True, print_func=None
    )
    graph.recurse_descendants(
        node, depth=None, link_types=(), annotate_links='both', include_process_inputs=True, print_func=None
    )
    output_file_name = graph.graphviz.render(str(data_dir) + '/' + name, format='pdf', view=False, cleanup=True)
    print('wrote graph to', output_file_name)


def get_single_imp_inputs():
    # import single imp calculations
    imported_nodes = import_with_migration(
        test_dir / 'data_dir/kkr_imp_wc-nodes-1e7804d6388fea2ca1e492d2f1a148c8.tar.gz'
    )['Node']
    for _, pk in imported_nodes['new'] + imported_nodes['existing']:
        node = load_node(pk)
        if node.label == 'kkrimp_scf full Cu host_in_host':
            imp1 = node
    imp1_out = imp1.outputs.workflow_info
    imp2_out = imp1_out  # use the same impurity and create a dimer

    return imp1_out, imp2_out


def get_builder_basic(label, kkrhost_local_code, kkrimp_local_code):

    imp1_out, imp2_out = get_single_imp_inputs()

    #set up combine_imps_wc workflow
    builder = combine_imps_wc.get_builder()
    builder.metadata.label = label
    builder.impurity1_output_node = imp1_out
    builder.impurity2_output_node = imp2_out
    builder.offset_imp2 = Dict(dict={'index': 1})
    # set code
    builder.host_gf.kkr = kkrhost_local_code  # should not be required if gf_host_remote is given, seems to be a problem of aiida-testing's run_with_cache
    builder.scf.kkrimp = kkrimp_local_code
    # set computer options
    options = {
        'queue_name': '',
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    builder.scf.options = Dict(dict=options)
    builder.scf.wf_parameters = Dict(dict={'do_final_cleanup': False})  # this is needed to allow for caching
    builder.host_gf.options = builder.scf.options

    return builder


def test_combine_imps(
    clear_database_before_test, kkrhost_local_code, kkrimp_local_code, run_with_cache, nopytest=False
):
    """
    test for combine_imps_wc (place two imps next to each other)
    """

    builder = get_builder_basic('test_combine_imps', kkrhost_local_code, kkrimp_local_code)

    # now submit
    print(builder, type(builder))
    if not nopytest:
        out, node = run_with_cache(builder, data_dir=data_dir)
    else:
        from aiida.engine import run_get_node
        out, node = run_get_node(builder)
    print((out, node))
    write_graph(node)

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert results['successful']
    assert results['convergence_reached']
    assert 'remote_data_gf' in out  # make sure GF writeout step was done


def test_combine_imps_params_kkr_overwrite(
    clear_database_before_test, kkrhost_local_code, kkrimp_local_code, run_with_cache, nopytest=False
):
    """
    test for combine_imps_wc overwriting the k-mesh with hte params_kkr_overwrite input to the gf writeout step
    """

    builder = get_builder_basic('test_combine_imps_params_kkr_overwrite', kkrhost_local_code, kkrimp_local_code)
    # increase k-mesh for GF writeout step with params_kkr_overwrite
    builder.host_gf.params_kkr_overwrite = Dict(dict={'BZDIVIDE': [20, 20, 20]})

    # now submit
    print(builder, type(builder))
    if not nopytest:
        out, node = run_with_cache(builder, data_dir=data_dir)
    else:
        from aiida.engine import run_get_node
        out, node = run_get_node(builder)
    print((out, node))
    write_graph(node, '_params_kkr_overwrite')

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert results['successful']
    assert results['convergence_reached']
    assert 'remote_data_gf' in out  # make sure GF writeout step was done


def test_combine_imps_reuse_gf(
    clear_database_before_test, kkrhost_local_code, kkrimp_local_code, run_with_cache, nopytest=False
):
    """
    test for combine_imps_wc reusing the host gf from a previous calculation
    """

    # import previous combine_imps workflow and reuse the host GF
    imported_nodes = import_with_migration(
        test_dir / 'data_dir/combine_imps_wc-nodes-76b8a1118abf371e7948b7ef15bf6e55.tar.gz'
    )['Node']
    for _, pk in imported_nodes['new'] + imported_nodes['existing']:
        node = load_node(pk)
        if node.label == 'test_combine_imps':
            imp_combine_imported = node
    host_gf_remote = imp_combine_imported.outputs.remote_data_gf

    builder = get_builder_basic('test_combine_imps_reuse_gf', kkrhost_local_code, kkrimp_local_code)
    # set pre-calculated host_gf
    builder.gf_host_remote = host_gf_remote

    # now submit
    print(builder, type(builder))
    if not nopytest:
        out, node = run_with_cache(builder, data_dir=data_dir)
    else:
        from aiida.engine import run_get_node
        out, node = run_get_node(builder)
    print((out, node))
    write_graph(node, '_reuse_gf')

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert results['successful']
    assert results['convergence_reached']
    assert 'remote_data_gf' not in out  # make sure GF writeout step was skipped


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
    #test_combine_imps(None, kkrhost_local_code, kkrimp_local_code, None, nopytest=True)
    test_combine_imps_reuse_gf(None, kkrhost_local_code, kkrimp_local_code, None, nopytest=True)
