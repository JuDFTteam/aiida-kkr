#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *

# tests
class Test_kkrimp_dos_workflow():
    """
    Tests for the kkrimp_scf workflow
    """

    @pytest.mark.timeout(300, method='thread')
    @pytest.mark.usefixtures("fresh_aiida_env")
    def test_dos_startpot_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node, Dict, StructureData
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp_dos import kkr_imp_dos_wc
        from numpy import array


        #"""
        # import data from previous run to use caching
        from aiida.tools.importexport import import_data
        import_data('files/export_kkrimp_dos.tar.gz', extras_mode_existing='ncu', extras_mode_new='import')
        #import_data('export_kkrimp_dos.tar.gz', extras_mode_existing='ncu', extras_mode_new='import')

        # import previous GF writeout, need to do this before rehashing
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrflex_create.tar.gz')

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
        #"""

        GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

        # prepare computer and code (needed so that
        if kkrimp_codename=='kkrimp':
            prepare_code(kkrimp_codename, codelocation, computername, workdir)
        if kkr_codename=='kkrhost':
            prepare_code(kkr_codename, codelocation, computername, workdir)


        wfd =kkr_imp_dos_wc.get_wf_defaults()
        wfd['clean_impcalc_retrieved'] = False # deactivate cleaning of unused data to regain cachability

        options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
        options = Dict(dict=options)

        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)
        KKRCode = Code.get_from_string(kkr_codename+'@'+computername)


        # now create a SingleFileData node containing the impurity starting potential
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from numpy import loadtxt
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.outputs.retrieved.open('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(dict=settings_dict)

        from aiida.manage.caching import enable_caching
        with enable_caching(): # should enable caching globally in this python interpreter 
            startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder)

        label = 'kkrimp_dos Cu host_in_host'
        descr = 'kkrimp_dos workflow for Cu bulk'

        imp_info = GF_host_calc.inputs.impurity_info.get_dict()
        imp_info ['Rcut'] = 2.5533
        print(imp_info)

        # create process builder to set parameters
        builder = kkr_imp_dos_wc.get_builder()
        builder.metadata.description = descr
        builder.metadata.label = label
        builder.options = options
        builder.kkr = KKRCode
        builder.kkrimp = KKRimpCode
        builder.imp_pot_sfd = startpot_imp_sfd
        builder.wf_parameters = Dict(dict=wfd)
        builder.impurity_info = Dict(dict=imp_info)
        builder.host_remote = GF_host_calc.outputs.remote_folder

        # now run calculation
        from aiida.engine import run_get_node #run
        from aiida.manage.caching import enable_caching
        from aiida.manage.caching import get_use_cache
        with enable_caching(): # should enable caching globally in this python interpreter
            print(builder)
            out, node = run_get_node(builder)
        print(node)
        print(out)

        assert 'last_calc_info' in list(out.keys())
        assert 'last_calc_output_parameters' in list(out.keys())
        assert 'workflow_info' in list(out.keys())
        assert 'dos_data' in list(out.keys())
        assert 'dos_data_interpol' in list(out.keys())
        assert len(out['dos_data_interpol'].get_y()) == 5
        assert len(out['dos_data_interpol'].get_y()[0]) == 3
        assert len(out['dos_data_interpol'].get_y()[0][0]) == 20


        # create export file
        #from aiida.tools.importexport import export
        #export([startpot_imp_sfd, node], outfile='export_kkrimp_dos.tar.gz', overwrite=True, silent=False)

    """
    @pytest.mark.timeout(300, method='thread')
    def test_dos_reuse_gf_writeout(self):
        pass


    @pytest.mark.timeout(300, method='thread')
    def test_dos_from_kkrimp_sub(self):
        pass


    @pytest.mark.timeout(300, method='thread')
    def test_dos_from_kkrimp_full(self):
        pass
    """

#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   Test = Test_kkrimp_dos_workflow()
   Test.test_dos_startpot_wc()
   #Test.test_dos_reuse_gf_writeout()
   #Test.test_dos_from_kkrimp_sub()
   #Test.test_dos_from_kkrimp_full()
