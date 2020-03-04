#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache
from ..conftest import voronoi_local_code, kkrhost_local_code
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile



@pytest.mark.timeout(500, method='thread')
def test_eos_wc_Cu_simple(aiida_profile, voronoi_local_code, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.eos import kkr_eos_wc
    from pprint import pprint
    from numpy import array

    # create structure
    alat = 6.83 # in a_Bohr
    abohr = 0.52917721067 # conversion factor to Angstroem units
    # bravais vectors
    bravais = array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])

    a = 0.5*alat*abohr
    Cu = StructureData(cell=[[a, a, 0.0], [a, 0.0, a], [0.0, a, a]])
    Cu.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')

    Cu.store()

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd, options = kkr_eos_wc.get_wf_defaults()
    wfd['nsteps'] = 3
    wfd['settings_kkr_scf']['convergence_criterion'] = 10**-4
    wfd['settings_kkr_scf']['convergence_setting_fine'] = wfd['settings_kkr_scf']['convergence_setting_coarse']
    wfd['settings_kkr_scf']['nsteps'] = 80
    wfd['settings_kkr_scf']['num_rerun'] = 2
    wfd['settings_kkr_scf']['natom_in_cls_min'] = 20
    wfd['settings_kkr_startpot']['natom_in_cls_min'] = 20
    wfd['settings_kkr_startpot']['num_rerun'] = 2
    wfd['fitfunction'] = 'sj' # for only three points only sj fit works

    KKReos_wf_parameters = Dict(dict=wfd)
    options['queue_name'] = queuename
    options['max_wallclock_seconds'] = 5*60
    options['withmpi'] = False
    options = Dict(dict=options)

    # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a Dict object for the use in aiida
    ParaNode = Dict(dict=kkrparams(LMAX=2, RMAX=7, GMAX=65, NSPIN=1, RCLUSTZ=1.9).get_dict())

    label = 'KKR-eos for Cu bulk'
    descr = 'KKR equation of states for Cu bulk'

    # create process builder to set parameters
    builder = kkr_eos_wc.get_builder()
    builder.calc_parameters = ParaNode
    builder.voronoi = voronoi_local_code
    builder.kkr = kkrhost_local_code
    builder.structure = Cu
    builder.wf_parameters = KKReos_wf_parameters
    builder.options = options
    builder.metadata.label = label
    builder.metadata.description = descr

    # now run calculation
    out, node = run_with_cache(builder)

    # load node of workflow
    print(out)
    n = out['eos_results']

    # get output dictionary
    out = n.get_dict()
    print('\n\noutput dictionary:\n-------------------------------------------------')
    pprint(out)

    # finally check some output
    print('\n\ncheck values ...\n-------------------------------------------------')

    print('successful', out['successful'])
    assert out['successful']

    print('rms', out['rms'])
    assert max(out['rms'])<10**-4

    print('gs_scale_factor', out['gs_scale_factor'])
    assert abs(out['gs_scale_factor']-1.0707334700693) < 10**-7

    print('\ndone with checks\n')

#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   test_eos_wc_Cu_simple()
