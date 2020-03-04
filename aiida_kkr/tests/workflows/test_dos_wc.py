#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache
from ..conftest import voronoi_local_code, kkrhost_local_code
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile


# change kkr_condename for testing (on mac)
kkr_codename = 'kkrhost_intel19'


@pytest.mark.timeout(240, method='thread')
def test_dos_wc_Cu(aiida_profile, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida import get_version
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm import Computer
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.dos import kkr_dos_wc
    from numpy import array

    print('AiiDA version: {}'.format(get_version()))

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
    wfd = kkr_dos_wc.get_wf_defaults()
    wfd['dos_params']['kmesh'] = [10, 10, 10]
    wfd['dos_params']['nepts'] = 10
    params_dos = Dict(dict=wfd)

    options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
    options = Dict(dict=options)

    label = 'dos Cu bulk'
    descr = 'DOS workflow for Cu bulk'

    from aiida.tools.importexport import import_data
    import_data('files/db_dump_kkrcalc.tar.gz')
    kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

    # create process builder to set parameters
    builder = kkr_dos_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkr = kkrhost_local_code
    builder.wf_parameters = params_dos
    builder.options = options
    builder.remote_data = kkr_calc_remote

    # now run calculation
    out, node = run_with_cache(builder)

    # check outcome
    n = out['results_wf']
    n = n.get_dict()
    assert n.get('successful')
    assert n.get('list_of_errors') == []
    assert n.get('dos_params').get('nepts') == 10

    d = out['dos_data']
    x = d.get_x()
    y = d.get_y()

    assert sum(abs(x[1][0] - array([-19.24321191, -16.2197246 , -13.1962373 , -10.17274986, -7.14926255,  -4.12577525,  -1.10228794,   1.9211995 , 4.94468681,   7.96817411]))) < 10**-7
    assert sum(abs(y[0][1][0] - array([  9.86819781e-04,   1.40981029e-03,   2.27894713e-03,         4.79231363e-03,   3.59368494e-02,   2.32929524e+00,         3.06973485e-01,   4.17629157e-01,   3.04021941e-01,         1.24897739e-01]))) < 10**-8

#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   test_dos_wc_Cu()
