#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
if __name__!='__main__':
    import pytest
    from ..conftest import kkrimp_local_code, kkrhost_local_code
    from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test

from aiida.orm import load_node, Dict
from aiida.engine import run_get_node
from aiida_kkr.workflows import combine_imps_wc
from aiida.tools.importexport import import_data

def test_combine_imps(clear_database_before_test, kkrhost_local_code, kkrimp_local_code):
    """
    test for combine_imps_wc (place two imps next to each other)
    """
    # import single imp calculations
    #import_data('data_dir/kkr_imp_wc-nodes-a2780cfccc03ac4373b9a1169ad605c2.tar.gz', silent=True)['Node']
    imp1_out = load_node('ead0a84a-57a4-4315-9165-ecd8c4307e75')
    imp2_out = imp1_out # use tha same impurity and create a dimer

    #set up combine_imps_wc workflow
    builder = combine_imps_wc.get_builder()
    builder.impurity1_output_node = imp1_out
    builder.impurity2_output_node = imp2_out
    builder.offset_imp2 = Dict(dict={'index': 1})
    builder.metadata.label = 'test_combine_imps'
    # set codes
    builder.host_gf.kkr = kkrhost_local_code
    builder.scf.kkrimp = kkrimp_local_code
    # set computer options
    options = {'queue_name' : '', 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
    builder.host_gf.options = Dict(dict=options)
    builder.scf.options = builder.host_gf.options

    # now submit
    print(builder)
    out, node = run_get_node(builder)
    print((out, node))

    # check outcome
    #assert  ...

def test_create_combined_imo_info():
    from aiida.orm import load_node
    host_structure1, impinfo1 = load_node('e51ee6a1-bd27-4901-9612-7bac256bf117'), load_node('1a5c2ebc-05db-4b76-bb0c-ebe2cbf3f67d')
    from aiida_kkr.tools.combine_imps import create_combined_imp_info_cf
    o = create_combined_imp_info_cf(host_structure1, impinfo1, impinfo1, Dict(dict={'index': 1}))
    o['imp_info_combined'].get_dict()
    o['imp_info_combined']['imp_cls']
    imp_cls_out = o['imp_info_combined']['imp_cls']
    return imp_cls_out


# run manual:
if __name__=='__main__':
    from aiida import load_profile
    load_profile()
    from aiida.orm import Code
    
    #define codes
    kkrhost_local_code = Code.get_from_string('kkrhost_intel19@localhost')
    kkrimp_local_code = Code.get_from_string('kkrimp_intel19@localhost')

    # run test
    test_combine_imps('dummy', kkrhost_local_code, kkrimp_local_code)
    test_create_combined_imo_info()