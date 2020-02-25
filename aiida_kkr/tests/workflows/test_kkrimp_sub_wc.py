#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *

# tests
class Test_kkrimp_scf_workflow():
    """
    Tests for the kkrimp_scf workflow
    """

    @pytest.mark.usefixtures("fresh_aiida_env")
    @pytest.mark.timeout(600, method='thread')
    def test_kkrimp_sub_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node, Dict, StructureData
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
        from numpy import array


        # import data from previous run to use caching
        from aiida.tools.importexport import import_data
        import_data('files/export_kkrimp_sub.tar.gz', extras_mode_existing='ncu', extras_mode_new='import')

        # need to rehash after import, otherwise cashing does not work
        from aiida.orm import Data, ProcessNode, QueryBuilder
        entry_point = (Data, ProcessNode)
        qb = QueryBuilder()
        qb.append(entry_point, tag='node') # query for ProcessNodes
        to_hash = qb.all()
        num_nodes = qb.count()
        print(num_nodes, to_hash)
        for node in to_hash:
            node[0].rehash()


        # prepare computer and code (needed so that
        if kkrimp_codename=='kkrimp':
            prepare_code(kkrimp_codename, codelocation, computername, workdir)


        wfd =kkr_imp_sub_wc.get_wf_defaults()

        wfd['nsteps'] = 20
        wfd['strmix'] = 0.05

        options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
        options = Dict(dict=options)

        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)

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
        builder.kkrimp = KKRimpCode
        builder.options = options
        builder.remote_data = GF_host_calc.outputs.remote_folder
        builder.wf_parameters = Dict(dict=wfd)
        builder.host_imp_startpot = startpot_imp_sfd

        # now run calculation
        from aiida.engine import run_get_node #run
        from aiida.manage.caching import enable_caching
        from aiida.manage.caching import get_use_cache
        with enable_caching(): # should enable caching globally in this python interpreter 
            out, node = run_get_node(builder)
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
   Test = Test_kkrimp_scf_workflow()
   Test.test_kkrimp_sub_wc()
