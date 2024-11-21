#!/usr/bin/env python

import pytest
from ..dbsetup import *
if __name__ != '__main__':
    from aiida.engine import run_get_node
    from ..conftest import voronoi_local_code, kkrhost_local_code, test_dir, data_dir, import_with_migration


@pytest.mark.timeout(900, method='thread')
def test_kkrimp_full_wc(
    clear_database_before_test, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, enable_archive_cache
):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from aiida.orm import Code, load_node, Dict, StructureData, load_group
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
    from numpy import array

    # settings for workflow
    _, wfd, voro_aux_settings = kkr_imp_wc.get_wf_defaults()

    wfd['nsteps'] = 5
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
    voro_aux_settings['check_dos'] = False
    voro_aux_settings['natom_in_cls_min'] = 50
    voro_aux_settings['rclustz'] = 1.5

    voro_aux_settings = Dict(voro_aux_settings)
    wf_inputs = Dict(wfd)

    imp_info = Dict({'Rcut': 2.5533, 'ilayer_center': 0, 'Zimp': [30.]})

    # import parent calculation (converged host system)
    group_pk = import_with_migration('files/db_dump_kkrflex_create.tar.gz')
    gf_writeout_workflow = [i for i in load_group(group_pk).nodes if i.label == 'GF_writeout Cu bulk'][0]

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
    builder.remote_data_host = gf_writeout_workflow.outputs.GF_host_remote
    builder.scf.params_overwrite = Dict({'TOL_ALAT_CHECK': 1e-8})

    # now run calculation
    with enable_archive_cache(data_dir / 'kkrimp_full_wc.aiida'):
        out, node = run_get_node(builder)
    print(out)
    print(list(node.called))

    # check outcome
    n = out['workflow_info']
    n = n.get_dict()
    print(n)
    for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
        assert sub in list(n.get('used_subworkflows').keys())
    # check outcome (calculation cannot be converged, but should be on the way)
    assert not n['converged']
    assert n['gf_wc_success']
    assert n['voro_wc_success']
    assert not n['kkrimp_wc_success']
    assert n['number_of_rms_steps'] == 25
    rms = n['convergence_values_all_steps']
    assert rms[-1] < rms[0]
    assert rms[-1] < 1.40


# @pytest.mark.timeout(900, method='thread')
# def test_kkrimp_full_Ag_Cu_onsite(
#     clear_database_before_test, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, enable_archive_cache
# ):
#     """
#     Simple Ag_Cu (bulk) noSOC, FP, lmax2 example  where impurity cluster contains only the impurity atom
#     """
#     from aiida.orm import Code, load_node, Dict, StructureData, load_group
#     from aiida.orm.querybuilder import QueryBuilder
#     from masci_tools.io.kkr_params import kkrparams
#     from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
#     from numpy import array

#     # settings for workflow
#     options, wfd, voro_aux_settings = kkr_imp_wc.get_wf_defaults()

#     # workflow behavior
#     wfd['nsteps'] = 10
#     wfd['strmix'] = 0.05
#     wfd['do_final_cleanup'] = False
#     wfd['convergence_criterion'] = 10**-4
#     # computer settings
#     options = {
#         'queue_name': queuename,
#         'resources': {
#             'num_machines': 1
#         },
#         'max_wallclock_seconds': 5 * 60,
#         'withmpi': False,
#         'custom_scheduler_commands': ''
#     }
#     options = Dict(options)
#     # voronoi settings for impurity startpot
#     voro_aux_settings['check_dos'] = False
#     voro_aux_settings['natom_in_cls_min'] = 50
#     voro_aux_settings['rclustz'] = 1.5

#     # make cluster radius small so that only the impurity is inside
#     imp_info = Dict({'Rcut': 3.5, 'ilayer_center': 0, 'Zimp': [47.]})

#     # import parent calculation (converged host system)
#     group_pk = import_with_migration('data_dir/kkr_scf_wc-nodes-31a2e00e231215133475de79d47f7c0b.tar.gz')
#     for node in load_group(group_pk).nodes:
#         if node.label == 'KKR-scf for Cu bulk':
#             kkr_scf_wc = node
#     kkr_converged = load_node(kkr_scf_wc.outputs.output_kkr_scf_wc_ParameterResults['last_calc_nodeinfo']['uuid'])
#     kkrhost_calc_remote = kkr_converged.outputs.remote_folder

#     # give workflow label and description
#     label = 'kkrimp_scf full Cu host_in_host'
#     descr = 'kkrimp_scf full workflow for Cu bulk inlcuding GF writeout and vorostart for starting potential'

#     # create process builder to set parameters
#     builder = kkr_imp_wc.get_builder()
#     builder.metadata.description = descr
#     builder.metadata.label = label
#     builder.kkrimp = kkrimp_local_code
#     builder.voronoi = voronoi_local_code
#     builder.kkr = kkrhost_local_code
#     builder.options = options
#     builder.voro_aux_parameters = Dict(voro_aux_settings)
#     builder.wf_parameters = Dict(wfd)
#     builder.impurity_info = imp_info
#     builder.remote_data_host = kkrhost_calc_remote
#     builder.scf.params_overwrite = Dict({'TOL_ALAT_CHECK': 1e-8})

#     # now run calculation
#     with enable_archive_cache(data_dir / 'kkrimp_full_Ag_Cu_onsite.aiida'):
#         out, node = run_get_node(builder)
#     print(out)

#     # check outcome
#     n = out['workflow_info']
#     n = n.get_dict()
#     print(n)
#     for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
#         assert sub in list(n.get('used_subworkflows').keys())

#     kkrimp_sub = load_node(n['used_subworkflows']['kkr_imp_sub'])
#     assert kkrimp_sub.outputs.workflow_info.get_dict().get('successful')

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
