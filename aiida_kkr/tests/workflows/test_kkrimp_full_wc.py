#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache
from ..conftest import voronoi_local_code, kkrhost_local_code
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile


@pytest.mark.timeout(600, method='thread')
def test_kkrimp_full_wc(aiida_profile, voronoi_local_code, kkrhost_local_code, kkrimp_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
    from numpy import array

    # settings for workflow
    options, wfd, voro_aux_settings =kkr_imp_wc.get_wf_defaults()

    wfd['nsteps'] = 20
    wfd['strmix'] = 0.05
    wfd['do_final_cleanup'] = False
    options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
    options = Dict(dict=options)
    voro_aux_settings['check_dos'] = False
    voro_aux_settings['dos_params']['kmesh'] = [10,10,10]
    voro_aux_settings['dos_params']['nepts'] = 10
    voro_aux_settings['natom_in_cls_min'] = 50
    voro_aux_settings['rclustz'] = 1.5

    voro_aux_settings = Dict(dict=voro_aux_settings)
    wf_inputs = Dict(dict=wfd)

    imp_info = Dict(dict={'Rcut':2.5533, 'ilayer_center': 0, 'Zimp':[30.]})

    # import parent calculation (converged host system)
    from aiida.tools.importexport import import_data
    import_data('files/db_dump_kkrflex_create.tar.gz')
    kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

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
    builder.remote_data_host = kkr_calc_remote

    # now run calculation
    out, node = run_with_cache(builder)
    print(out)

    # check outcome
    n = out['workflow_info']
    n = n.get_dict()
    for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
        assert sub in list(n.get('used_subworkflows').keys())

    kkrimp_sub = load_node(n['used_subworkflows']['kkr_imp_sub'])
    assert kkrimp_sub.outputs.workflow_info.get_dict().get('successful')


#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   test_kkrimp_full_wc()
