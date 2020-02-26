#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *

# tests
@pytest.mark.usefixtures("fresh_aiida_env")
class Test_kkrimp_full_workflow():
    """
    Tests for the full kkrimp_scf workflow with GF writeout and voroaux steps
    """

    @pytest.mark.timeout(600, method='thread')
    def test_kkrimp_full_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node, Dict, StructureData
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
        from numpy import array


        #"""
        # import data from previous run to use caching
        from aiida.tools.importexport import import_data
        import_data('files/export_kkrimp_full.tar.gz', extras_mode_existing='ncu', extras_mode_new='import')

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


        # prepare computer and code (needed so that
        prepare_code(voro_codename, codelocation, computername, workdir)
        if kkrimp_codename=='kkrimp':
            prepare_code(kkrimp_codename, codelocation, computername, workdir)
        if kkr_codename=='kkrhost':
            prepare_code(kkr_codename, codelocation, computername, workdir)

        options, wfd, voro_aux_settings =kkr_imp_wc.get_wf_defaults()

        wfd['nsteps'] = 20
        wfd['strmix'] = 0.05
        wfd['do_final_cleanup'] = False
        options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
        options = Dict(dict=options)
        voro_aux_settings['check_dos'] = False
        voro_aux_settings['dos_params']['kmesh'] = [10,10,10]
        voro_aux_settings['dos_params']['nepts'] = 10
        voro_aux_settings['natom_in_cls_min'] = 50
        voro_aux_settings['rclustz'] = 1.5

        voro_aux_settings = Dict(dict=voro_aux_settings)
        wf_inputs = Dict(dict=wfd)

        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRhostCode = Code.get_from_string(kkr_codename+'@'+computername)
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)
        VoroCode = Code.get_from_string(voro_codename+'@'+computername)

        imp_info = Dict(dict={'Rcut':2.5533, 'ilayer_center': 0, 'Zimp':[30.]})

        kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

        label = 'kkrimp_scf full Cu host_in_host'
        descr = 'kkrimp_scf full workflow for Cu bulk inlcuding GF writeout and vorostart for starting potential'

        # create process builder to set parameters
        builder = kkr_imp_wc.get_builder()
        builder.metadata.description = descr
        builder.metadata.label = label
        builder.kkrimp = KKRimpCode
        builder.voronoi = VoroCode
        builder.kkr = KKRhostCode
        builder.options = options
        builder.voro_aux_parameters = voro_aux_settings
        builder.wf_parameters = wf_inputs
        builder.impurity_info = imp_info
        builder.remote_data_host = kkr_calc_remote

        # now run calculation
        from aiida.engine import run_get_node #run
        from aiida.manage.caching import enable_caching
        from aiida.manage.caching import get_use_cache
        with enable_caching(): # should enable caching globally in this python interpreter 
            out, node = run_get_node(builder)
        print(out)

        # check outcome
        n = out['workflow_info']
        n = n.get_dict()
        for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
            assert sub in list(n.get('used_subworkflows').keys())

        kkrimp_sub = load_node(n['used_subworkflows']['kkr_imp_sub'])
        assert kkrimp_sub.outputs.workflow_info.get_dict().get('successful')

        # create export file
        #from aiida.tools.importexport import export
        #export([node], outfile='export_kkrimp_full.tar.gz', overwrite=True, silent=False)


#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   Test = Test_kkrimp_full_workflow()
   Test.test_kkrimp_full_wc()
