#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from dbsetup import *

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_vorostart_workflow():
    """
    Tests for the kkr_startpot workflow
    """
    # make sure running the workflow exists after at most 5 minutes
    import timeout_decorator
    @timeout_decorator.timeout(300, use_signals=False)
    def run_timeout(self, builder):
        from aiida.engine.launch import run
        out = run(builder)
        return out
    
    def test_vorostart_wc_Cu(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow
        """
        from aiida.plugins import Code, load_node, DataFactory
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.voro_start import kkr_startpot_wc
        from numpy import array
       
        ParameterData = DataFactory('parameter')
        StructureData = DataFactory('structure')

        # prepare computer and code (needed so that 
        prepare_code(voro_codename, codelocation, computername, workdir)

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
        wfd = kkr_startpot_wc.get_wf_defaults()
        wfd['check_dos'] = False
        wfd['natom_in_cls_min'] = 20
        wfd['num_rerun'] = 2
        options = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'use_mpi' : False, 'custom_scheduler_commands' : ''}
        params_vorostart = Dict(dict=wfd)
        options = Dict(dict=options)
       
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        VoroCode = Code.get_from_string(voro_codename+'@'+computername)
       
        # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a ParameterData object for the use in aiida
        ParaNode = Dict(dict=kkrparams(LMAX=2, NSPIN=1, RCLUSTZ=1.9).get_dict())
       
        # create process builder to set parameters
        builder = kkr_startpot_wc.get_builder()
        builder.calc_parameters = ParaNode
        builder.description = 'voronoi startpot workflow for Cu bulk'
        builder.label = 'startpot for Cu bulk'
        builder.voronoi = VoroCode
        builder.structure = Cu
        builder.wf_parameters = params_vorostart
        builder.options = options
       
        # now run calculation
        out = self.run_timeout(builder)
        print(out)
       
        # check output
        n = out['results_vorostart_wc']
        n = n.get_dict()
        assert n.get('successful')
        assert n.get('last_voro_ok')
        assert n.get('list_of_errors') == []
        assert abs(n.get('starting_fermi_energy') - 0.409241) < 10**-14
 
#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_vorostart_workflow()
   Test.test_vorostart_wc_Cu()
