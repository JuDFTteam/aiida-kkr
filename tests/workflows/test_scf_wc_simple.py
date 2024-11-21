#!/usr/bin/env python

import pytest
from aiida.engine import run_get_node
from ..dbsetup import *
from ..conftest import voronoi_local_code, kkrhost_local_code, test_dir, data_dir


# helper function
def print_clean_inouts(node):
    from pprint import pprint
    inputs = node.get_inputs_dict()
    outputs = node.get_outputs_dict()
    for key in outputs.keys():
        try:
            int(key.split('_')[-1])
            has_id = True
        except:
            has_id = False
        if 'CALL' in key or 'CREATE' in key or has_id:
            outputs.pop(key)
    print('inputs:')
    pprint(inputs)
    print('\noutputs:')
    pprint(outputs)


@pytest.mark.timeout(300, method='thread')
def test_scf_wc_Cu_simple(clear_database_before_test, voronoi_local_code, kkrhost_local_code, enable_archive_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
    from pprint import pprint
    from numpy import array

    # create structure
    alat = 6.83  # in a_Bohr
    abohr = 0.52917721067  # conversion factor to Angstroem units
    a = 0.5 * alat * abohr
    Cu = StructureData(cell=[[a, a, 0.0], [a, 0.0, a], [0.0, a, a]])
    Cu.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')

    Cu.store()
    print(Cu)

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd = kkr_scf_wc.get_wf_defaults()[0]

    wfd['convergence_criterion'] = 10**-4
    wfd['check_dos'] = False
    wfd['kkr_runmax'] = 5
    wfd['nsteps'] = 50

    wfd['num_rerun'] = 2
    wfd['natom_in_cls_min'] = 20

    KKRscf_wf_parameters = Dict(wfd)

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

    # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a Dict object for the use in aiida
    ParaNode = Dict(kkrparams(LMAX=2, RMAX=7, GMAX=65, NSPIN=1, RCLUSTZ=1.9).get_dict())

    label = 'KKR-scf for Cu bulk'
    descr = 'KKR self-consistency workflow for Cu bulk'

    # create process builder to set parameters
    builder = kkr_scf_wc.get_builder()
    builder.calc_parameters = ParaNode
    builder.voronoi = voronoi_local_code
    builder.kkr = kkrhost_local_code
    builder.structure = Cu
    builder.wf_parameters = KKRscf_wf_parameters
    builder.options = options
    builder.metadata.label = label
    builder.metadata.description = descr

    # now run calculation
    with enable_archive_cache(data_dir / 'scf_wc_Cu_simple.aiida'):
        out, node = run_get_node(builder)

    # load node of workflow
    print(out)
    n = out['output_kkr_scf_wc_ParameterResults']

    # get output dictionary
    out = n.get_dict()
    print('\n\noutput dictionary:\n-------------------------------------------------')
    pprint(out)

    # finally check some output
    print('\n\ncheck values ...\n-------------------------------------------------')

    print('voronoi_step_success', out['voronoi_step_success'])
    assert out['voronoi_step_success']

    print('kkr_step_success', out['kkr_step_success'])
    assert out['kkr_step_success']

    print('successful', out['successful'])
    assert out['successful']

    print('error', out['errors'])
    assert out['errors'] == []

    print('warning', out['warnings'])
    assert out['warnings'] == []

    print('convergence_reached', out['convergence_reached'])
    assert out['convergence_reached']

    print('convergence_value', out['convergence_value'])
    assert out['convergence_value'] < 10**-4

    print('charge_neutrality', abs(out['charge_neutrality']))
    assert abs(out['charge_neutrality']) < 8 * 10**-4

    print('used_higher_accuracy', out['used_higher_accuracy'])
    assert out['used_higher_accuracy']

    print('\ndone with checks\n')


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    test_scf_wc_Cu_simple()
