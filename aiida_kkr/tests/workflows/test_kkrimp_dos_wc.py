#!/usr/bin/env python

from __future__ import absolute_import
import pytest
from aiida_kkr.tests.dbsetup import *

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_dos_workflow():
    """
    Tests for the kkrimp_scf workflow
    """

    @pytest.mark.timeout(300, method='thread')
    def test_dos_startpot_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp_dos import kkr_imp_dos_wc
        from numpy import array

        Dict = DataFactory('dict')
        StructureData = DataFactory('structure')

        # prepare computer and code (needed so that
        prepare_code(kkrimp_codename, codelocation, computername, workdir)


        wfd =kkr_imp_dos_wc.get_wf_defaults()

        options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'use_mpi' : False, 'custom_scheduler_commands' : ''}
        options = Dict(dict=options)

        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)
        KKRCode = Code.get_from_string(kkr_codename+'@'+computername)

        # import previous GF writeout
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('de9b5093-25e7-407e-939e-9282c4431343')

        # now create a SingleFileData node containing the impurity starting potential
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from numpy import loadtxt
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.outputs.retrieved.open('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder)

        label = 'kkrimp_dos Cu host_in_host'
        descr = 'kkrimp_dos workflow for Cu bulk'

        # create process builder to set parameters
        builder = kkr_imp_dos_wc.get_builder()
        builder.metadata.description = descr
        builder.metadata.label = label
        builder.options = options
        builder.kkr = KKRCode
        builder.kkrimp = KKRimpCode
        builder.host_imp_pot = startpot_imp_sfd
        builder.wf_parameters = Dict(dict=wfd)
        builder.impurity_info = GF_host_calc.inputs.impurity_info
        builder.host_remote = GF_host_calc.outputs.remote_folder

        # now run calculation
        from aiida.engine import run
        print(builder)
        out = run(builder)

        assert 'last_calc_info' in out.keys()
        assert 'last_calc_output_parameters' in out.keys()
        assert 'workflow_info' in out.keys()
        assert 'dos_data' in out.keys()
        assert 'dos_data_interpol' in out.keys()
        assert len(out['dos_data_interpol'].get_y()) == 5
        assert len(out['dos_data_interpol'].get_y()[0]) == 3
        assert len(out['dos_data_interpol'].get_y()[0][0]) == 20


    @pytest.mark.timeout(300, method='thread')
    def test_dos_from_kkrimp_sub(self):
        pass


    @pytest.mark.timeout(300, method='thread')
    def test_dos_from_kkrimp_full(self):
        pass

#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkrimp_dos_workflow()
   Test.test_dos_startpot_wc()
   Test.test_dos_from_kkrimp_sub()
   Test.test_dos_from_kkrimp_full()

