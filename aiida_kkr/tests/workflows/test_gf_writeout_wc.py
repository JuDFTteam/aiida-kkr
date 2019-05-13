#!/usr/bin/env python

from __future__ import absolute_import
import pytest
from aiida_kkr.tests.dbsetup import *

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_gf_writeout_workflow():
    """
    Tests for the kkr_startpot workflow
    """

    @pytest.mark.timeout(240, method='thread')
    def test_kkrflex_writeout_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
        from numpy import array
        import os

        Dict = DataFactory('dict')
        StructureData = DataFactory('structure')

        # prepare computer and code (needed so that
        prepare_code(kkr_codename, codelocation, computername, workdir)

        # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
        wfd =kkr_flex_wc.get_wf_defaults()
        options = {'queue_name' : queuename, 'resources': {"num_machines": 1},'max_wallclock_seconds' : 5*60, 'custom_scheduler_commands' : '', 'use_mpi' : False}
        options = Dict(dict=options)

        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRCode = Code.get_from_string(kkr_codename+'@'+computername)

        imp_info = Dict(dict={'Rcut':1.01, 'ilayer_center': 0, 'Zimp':[79.]})

        label = 'GF_writeout Cu bulk'
        descr = 'GF_writeout workflow for Cu bulk'

        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')
        kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

        # create process builder to set parameters
        builder = kkr_flex_wc.get_builder()
        builder.description = descr
        builder.label = label
        builder.kkr = KKRCode
        builder.options = options
        builder.remote_data = kkr_calc_remote
        builder.impurity_info = imp_info

        # now run calculation
        from aiida.engine import run
        out = run(builder)

        n = out['workflow_info']
        n = n.get_dict()

        assert n.get('successful')
        assert n.get('list_of_errors') == []

        d = out['GF_host_remote']
        assert isinstance(d, DataFactory('remote'))

        kkrflex_calc = load_node(n.get('pk_flexcalc'))
        kkrflex_retrieved = kkrflex_calc.outputs.retrieved
        for name in 'tmat green atominfo intercell_cmoms intercell_ref'.split():
            assert 'kkrflex_'+name in kkrflex_retrieved.list_object_names()


#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_gf_writeout_workflow()
   Test.test_kkrflex_writeout_wc()
