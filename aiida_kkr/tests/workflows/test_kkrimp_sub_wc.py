#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache
from ..conftest import kkrimp_local_code
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile



@pytest.mark.timeout(600, method='thread')
def test_kkrimp_sub_wc(aiida_profile, kkrimp_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
    from numpy import array

    wfd =kkr_imp_sub_wc.get_wf_defaults()

    wfd['nsteps'] = 20
    wfd['strmix'] = 0.05
    # deactivate final cleanup to be able to use caching
    wfd['do_final_cleanup'] = False

    options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
    options = Dict(dict=options)

    # import previous GF writeout
    from aiida.tools.importexport import import_data
    import_data('files/db_dump_kkrflex_create.tar.gz')
    GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

    # now create a SingleFileData node containing the impurity starting potential
    from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
    from numpy import loadtxt
    neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.outputs.retrieved.open('scoef'), skiprows=1)[:,3]-1]
    settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
    settings = Dict(dict=settings_dict)
    from aiida.manage.caching import enable_caching
    with enable_caching(): # should enable caching globally in this python interpreter 
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder)

    label = 'kkrimp_scf Cu host_in_host'
    descr = 'kkrimp_scf workflow for Cu bulk'

    # create process builder to set parameters
    builder = kkr_imp_sub_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkrimp = kkrimp_local_code
    builder.options = options
    builder.remote_data = GF_host_calc.outputs.remote_folder
    builder.wf_parameters = Dict(dict=wfd)
    builder.host_imp_startpot = startpot_imp_sfd

    # now run calculation
    out, node = run_with_cache(builder)
    print(out)
    print(node)
    print(node.process_status)

    n = out['workflow_info']
    n = n.get_dict()

    assert n.get('successful')
    assert n.get('convergence_reached')


#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   test_kkrimp_sub_wc()
