#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache
from ..conftest import voronoi_local_code, kkrhost_local_code
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile


@pytest.mark.timeout(500, method='thread')
def test_vorostart_wc_Cu(aiida_profile, voronoi_local_code, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.voro_start import kkr_startpot_wc
    from numpy import array

    """
    # import data from previous run to use caching
    from aiida.tools.importexport import import_data
    import_data('files/export_kkr_startpot.tar.gz', extras_mode_existing='ncu', extras_mode_new='import')


    # need to rehash after import, otherwise cashing does not work
    from aiida.orm import Data, ProcessNode, QueryBuilder
    entry_point = (Data, ProcessNode)
    qb = QueryBuilder()
    qb.append(ProcessNode, tag='node') # query for ProcessNodes
    to_hash = qb.all()
    num_nodes = qb.count()
    print(num_nodes, to_hash)
    for node in to_hash:
        node[0].rehash()
    #"""


    # prepare computer and code (needed so that
    #prepare_code(voro_codename, codelocation, computername, workdir)
    #if kkr_codename=='kkrhost':
    #prepare_code(kkr_codename, codelocation, computername, workdir)

    # Then set up the structure
    alat = 6.83 # in a_Bohr
    abohr = 0.52917721067 # conversion factor to Angstroem units
    bravais = array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])# bravais vectors
    a = 0.5*alat*abohr
    Cu = StructureData(cell=[[a, a, 0.0], [a, 0.0, a], [0.0, a, a]])
    Cu.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')

    Cu.store()
    print(Cu)

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd = kkr_startpot_wc.get_wf_defaults()
    wfd['check_dos'] = True 
    wfd['natom_in_cls_min'] = 20
    wfd['num_rerun'] = 2
    options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
    params_vorostart = Dict(dict=wfd)

    # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
    #VoroCode = Code.get_from_string(voro_codename+'@'+computername)

    # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a Dict object for the use in aiida
    ParaNode = Dict(dict=kkrparams(LMAX=2, NSPIN=1, RCLUSTZ=1.9).get_dict())

    # create process builder to set parameters
    builder = kkr_startpot_wc.get_builder()
    builder.calc_parameters = ParaNode
    builder.metadata.description = 'voronoi startpot workflow for Cu bulk'
    builder.metadata.label = 'startpot for Cu bulk'
    #builder.voronoi = VoroCode
    builder.voronoi = voronoi_local_code
    builder.structure = Cu
    builder.wf_parameters = params_vorostart
    builder.options = Dict(dict=options)
    #KKRCode = Code.get_from_string(kkr_codename+'@'+computername)
    #builder.kkr = KKRCode
    builder.kkr = kkrhost_local_code

    """
    # now run calculation
    from aiida.engine import run_get_node #run
    from aiida.manage.caching import enable_caching
    from aiida.manage.caching import get_use_cache
    with enable_caching(): # should enable caching globally in this python interpreter 
        out, node = run_get_node(builder)
    print(out)
    from aiida.orm import Node
    for _,node in builder.items():
        if isinstance(node,Node):
            if not node.is_stored:
                node.store()
    """


    out, node = run_with_cache(builder)

    # check output
    n = out['results_vorostart_wc']
    n = n.get_dict()
    assert n.get('successful')
    assert n.get('last_voro_ok')
    assert n.get('list_of_errors') == []
    assert abs(n.get('starting_fermi_energy') - 0.409241) < 10**-14

    #kkrcalc = [i for i in node.called_descendants if i.process_label=='KkrCalculation'][0]
    #vorocalc = [i for i in node.called_descendants if i.process_label=='VoronoiCalculation'][0]
    #print('hashes of computed voro and kkr calcs', vorocalc.get_hash(), kkrcalc.get_hash())
    #print('was cached (voro/kkr)?', vorocalc.get_cache_source(), kkrcalc.get_cache_source())
    #print('caching info calcul.d voro', vorocalc._get_objects_to_hash())
    #print('caching info calcul.d kkr ', kkrcalc._get_objects_to_hash())

    # export data to reuse it later
    #from aiida.tools.importexport import export
    #export([node], outfile='export_data.aiida.tar.gz', overwrite=True) # add export of extras automatically

#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   #Test = Test_vorostart_workflow()
   #Test.test_vorostart_wc_Cu()
   test_vorostart_wc_Cu()
