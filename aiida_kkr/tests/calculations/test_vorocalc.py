#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object
from aiida_kkr.tests.dbsetup import *
from ..conftest import voronoi_local_code
from aiida_testing.export_cache._fixtures import run_with_cache
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile
import pytest


kkr_codename = 'kkrhost'

#TODO
# implement missing tests:
# * test_vca_structure
# * test_overwrite_alat_input
# * test_voronoi_after_kkr
# * test_overwrite_potential

def wait_for_it(calc, maxwait=300, dT=10):
    """
    helper function used to wait until calculation reaches FINISHED state
    wait for maximally <maxwait> seconds and check the calculation's state every <dT> seconds
    """
    from time import sleep
    nsteps = maxwait/dT
    print('waiting for calculation to finish (maximally wait for {} seconds)'.format(maxwait))
    istep = 0
    calcstate = u'UNKNOWN'
    while istep < nsteps:
        print('checking status')
        sleep(dT)
        calcstate = calc.get_state()
        istep += 1
        if calcstate == u'FINISHED' or calcstate == u'FAILED':
            break

    if calcstate == u'FINISHED':
        print('calculation reached FINISHED state')
    elif calcstate == u'FAILED':
        print('calculation in FAILED state')
    else:
        print('maximum waiting time exhausted')


# tests
def test_startpot_Cu_simple(aiida_profile, voronoi_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example
    """
    from aiida.orm import Code, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # create StructureData instance for Cu
    alat = 3.61 # lattice constant in Angstroem
    bravais = [[0.5*alat, 0.5*alat, 0], [0.5*alat, 0, 0.5*alat], [0, 0.5*alat, 0.5*alat]] # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')

    # create Dict input node using kkrparams class from masci-tools
    params = kkrparams(params_type='voronoi')
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    ParaNode = Dict(dict=params.get_dict())
    
    options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.metadata.options = options
    builder.parameters = ParaNode
    builder.structure = Cu
    builder._hash_ignored_inputs = ['code']
    builder.metadata.dry_run = True
    from aiida.engine import run
    run(builder)
    #out, node = run_with_cache(builder)
    #print(out, node)
    #print(node.get_cache_source())
    #print(node.get_hash())
    #print(voronoi_local_code.get_hash())
    #print(voronoi_local_code._get_objects_to_hash())
    #print(voronoi_local_code._hash_ignored_attributes)


def test_vca_structure(aiida_profile, voronoi_local_code):
    """
    test for vca_structure behaviour
    """
    pass

def test_overwrite_alat_input(aiida_profile, voronoi_local_code):
    """
    test using 'use_alat_input' keyword in input parameters
    """
    pass

def test_voronoi_after_kkr(aiida_profile, voronoi_local_code):
    """
    test voronoi run from parent kkr calculation (e.g. to update to a higher lmax value)
    """
    pass

def test_overwrite_potential(aiida_profile, voronoi_local_code):
    """
    test providing overwirte_potential input node which overwrites the starting potentai with the given input
    """
    pass


#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   Test = Test_voronoi_calculation()
   Test.test_startpot_Cu_simple()
