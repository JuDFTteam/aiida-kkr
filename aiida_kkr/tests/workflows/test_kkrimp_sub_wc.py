#!/usr/bin/env python

from __future__ import absolute_import
import pytest
from aiida_kkr.tests.dbsetup import *

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_scf_workflow():
    """
    Tests for the kkrimp_scf workflow
    """

    @pytest.mark.timeout(300, method='thread')
    def test_kkrimp_sub_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
        from numpy import array

        Dict = DataFactory('dict')
        StructureData = DataFactory('structure')

        # prepare computer and code (needed so that
        prepare_code(kkrimp_codename, codelocation, computername, workdir)


        wfd =kkr_imp_sub_wc.get_wf_defaults()

        wfd['nsteps'] = 20
        wfd['strmix'] = 0.05

        options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'use_mpi' : False, 'custom_scheduler_commands' : ''}
        options = Dict(dict=options)

        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)

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

        label = 'kkrimp_scf Cu host_in_host'
        descr = 'kkrimp_scf workflow for Cu bulk'

        # create process builder to set parameters
        builder = kkr_imp_sub_wc.get_builder()
        builder.description = descr
        builder.label = label
        builder.kkrimp = KKRimpCode
        builder.options = options
        builder.remote_data = GF_host_calc.outputs.remote_folder
        builder.wf_parameters = Dict(dict=wfd)
        builder.host_imp_startpot = startpot_imp_sfd

        # now run calculation
        from aiida.engine import run
        print(builder)
        out = run(builder)

        n = out['workflow_info']
        n = n.get_dict()

        assert n.get('successful')
        assert n.get('convergence_reached')

#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkrimp_scf_workflow()
   Test.test_kkrimp_sub_wc()
