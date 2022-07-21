#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from ..dbsetup import *
if __name__ != '__main__':
    from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
    from ..conftest import voronoi_local_code, kkrhost_local_code, test_dir, data_dir, import_with_migration
    from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile
    from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test


@pytest.mark.timeout(600, method='thread')
def test_kkrimp_full_wc(
    clear_database_before_test, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, run_with_cache
):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
    from numpy import array

    # settings for workflow
    _, wfd, voro_aux_settings = kkr_imp_wc.get_wf_defaults()

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
    options = Dict(dict=options)
    voro_aux_settings['check_dos'] = False
    voro_aux_settings['natom_in_cls_min'] = 50
    voro_aux_settings['rclustz'] = 1.5

    voro_aux_settings = Dict(dict=voro_aux_settings)
    wf_inputs = Dict(dict=wfd)

    imp_info = Dict(dict={'Rcut': 2.5533, 'ilayer_center': 0, 'Zimp': [30.]})

    # import parent calculation (converged host system)
    import_with_migration('files/db_dump_kkrflex_create.tar.gz')
    GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

    # give workflow label and description
    label = 'kkrimp_scf full Cu host_in_host'
    descr = 'kkrimp_scf full workflow for Cu bulk inlcuding GF writeout and vorostart for starting potential'

    # create process builder to set parameters
    builder = kkr_imp_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkrimp = kkrimp_local_code
    builder.voronoi = voronoi_local_code
    builder.kkr = kkrhost_local_code
    builder.options = options
    builder.voro_aux_parameters = voro_aux_settings
    builder.wf_parameters = wf_inputs
    builder.impurity_info = imp_info
    builder.remote_data_host = GF_host_calc.outputs.remote_folder

    # now run calculation
    out, node = run_with_cache(builder, data_dir=data_dir)
    print(out)
    print(list(node.called))

    # check outcome
    n = out['workflow_info']
    n = n.get_dict()
    for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
        assert sub in list(n.get('used_subworkflows').keys())

    kkrimp_sub = load_node(n['used_subworkflows']['kkr_imp_sub'])
    assert kkrimp_sub.outputs.workflow_info.get_dict().get('successful')


def test_kkrimp_full_Ag_Cu_onsite(
    clear_database_before_test, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, run_with_cache
):
    """
    Simple Ag_Cu (bulk) noSOC, FP, lmax2 example  where impurity cluster contains only the impurity atom
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
    from numpy import array

    # settings for workflow
    options, wfd, voro_aux_settings = kkr_imp_wc.get_wf_defaults()

    # workflow behavior
    wfd['nsteps'] = 50
    wfd['strmix'] = 0.05
    wfd['do_final_cleanup'] = False
    wfd['convergence_criterion'] = 10**-4
    # computer settings
    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    options = Dict(dict=options)
    # voronoi settings for impurity startpot
    voro_aux_settings['check_dos'] = False
    voro_aux_settings['natom_in_cls_min'] = 50
    voro_aux_settings['rclustz'] = 1.5

    # make cluster radius small so that only the impurity is inside
    imp_info = Dict(dict={'Rcut': 3.5, 'ilayer_center': 0, 'Zimp': [47.]})

    # import parent calculation (converged host system)
    imported_nodes = import_with_migration('data_dir/kkr_scf_wc-nodes-db396f0dabbf666d9a247b3dca766421.tar.gz')['Node']
    for _, pk in imported_nodes['new'] + imported_nodes['existing']:
        node = load_node(pk)
        if node.label == 'KKR-scf for Cu bulk':
            kkr_scf_wc = node
    kkr_converged = load_node(kkr_scf_wc.outputs.output_kkr_scf_wc_ParameterResults['last_calc_nodeinfo']['uuid'])
    kkrhost_calc_remote = kkr_converged.outputs.remote_folder

    # give workflow label and description
    label = 'kkrimp_scf full Cu host_in_host'
    descr = 'kkrimp_scf full workflow for Cu bulk inlcuding GF writeout and vorostart for starting potential'

    # create process builder to set parameters
    builder = kkr_imp_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkrimp = kkrimp_local_code
    builder.voronoi = voronoi_local_code
    builder.kkr = kkrhost_local_code
    builder.options = options
    builder.voro_aux_parameters = Dict(dict=voro_aux_settings)
    builder.wf_parameters = Dict(dict=wfd)
    builder.impurity_info = imp_info
    builder.remote_data_host = kkrhost_calc_remote

    # now run calculation
    out, node = run_with_cache(builder, data_dir=data_dir)
    print(out)

    # check outcome
    n = out['workflow_info']
    n = n.get_dict()
    for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
        assert sub in list(n.get('used_subworkflows').keys())

    kkrimp_sub = load_node(n['used_subworkflows']['kkr_imp_sub'])
    assert kkrimp_sub.outputs.workflow_info.get_dict().get('successful')


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    from aiida.engine import run_get_node
    from aiida.orm import Code

    #define codes
    try:
        # on mac
        kkrhost_local_code = Code.get_from_string('kkrhost_intel19@localhost')
        kkrimp_local_code = Code.get_from_string('kkrimp_intel19@localhost')
    except:
        # on iff desktop
        kkrhost_local_code = Code.get_from_string('kkrhost@localhost')
        kkrimp_local_code = Code.get_from_string('kkrimp@localhost')

    voronoi_local_code = Code.get_from_string('voronoi@localhost')

    #test_kkrimp_full_wc()
    test_kkrimp_full_Ag_Cu_onsite(None, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, run_get_node)
