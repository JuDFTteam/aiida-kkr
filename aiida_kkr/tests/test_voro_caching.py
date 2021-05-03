#!/usr/bin/env python

import pytest
from aiida_kkr.tests.dbsetup import *
from aiida.manage.tests.pytest_fixtures import aiida_profile, clear_database, clear_database_after_test

# tests
@pytest.mark.timeout(120, method='thread')
def test_vorocalc_caching(aiida_profile, voronoi_local_code):
    """
    simple test to see if caching works
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from numpy import array
    from aiida_kkr.calculations import VoronoiCalculation

    # Then set up the structure
    a = 0.5*6.83*0.52917721067 # alat in Ang.
    Cu = StructureData(cell=[[a, a, 0.0], [a, 0.0, a], [0.0, a, a]])
    Cu.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')
    Cu.store()

    # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a Dict object for the use in aiida
    ParaNode = Dict(dict=kkrparams(LMAX=2, NSPIN=1, RCLUSTZ=1.9).get_dict())

    # create process builder to set parameters
    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.parameters = ParaNode
    builder.structure = Cu
    builder.metadata.options.resources = {'num_machines': 1}
    builder.metadata.options.withmpi = False
    print(builder)

    # now run 1st calculation without caching
    from aiida.engine import run_get_node
    out1, node1 = run_get_node(builder)
    print(out1, node1)

    #run now with caching
    from aiida.manage.caching import enable_caching
    with enable_caching(): # should enable caching globally in this python interpreter 
        out2, node2 = run_get_node(builder)
    print(out2, node2)

    #run now another time with caching
    with enable_caching():
        out3, node3 = run_get_node(builder)
    print(out3, node3)

    # check if caching works
    no1 = node1.outputs.retrieved.list_object_names()
    no2 = node2.outputs.retrieved.list_object_names()
    no3 = node3.outputs.retrieved.list_object_names()
    
    # write stderr and stdout
    print('std.out')
    with node1.outputs.retrieved.open('_scheduler-stdout.txt') as f:
        print(f.readlines())
    print('std.err')
    with node1.outputs.retrieved.open('_scheduler-stderr.txt') as f:
        print(f.readlines())

    print(node1.get_cache_source(), no1)
    assert len(no1)==9
    assert node1.get_cache_source() is None
    print(node2.get_cache_source(), no2)
    assert len(no2)==9
    assert node2.get_cache_source() is not None
    print(node3.get_cache_source(), no3)
    assert len(no3)==9
    assert node3.get_cache_source() is not None

#run test manually
if __name__=='__main__':
    from aiida import load_profile
    load_profile()
    test_vorocalc_caching()
