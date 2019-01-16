#!/usr/bin/env python

import pytest

# some global settings

voro_codename = 'voronoi'
computername = 'localhost'
queuename = ''
workdir = '/temp/ruess/aiida_run_iff734/'
codelocation = '/Users/ruess/sourcecodes/aiida/codes_localhost'

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_vorostart_workflow():
    """
    Tests for the kkr_startpot workflow
    """
    
    def test_vorostart_wc_Cu(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow
        """
        from aiida.orm import Code, load_node, DataFactory
        from aiida.orm.computers import Computer
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.voro_start import kkr_startpot_wc
        from numpy import array
       
        ParameterData = DataFactory('parameter')
        StructureData = DataFactory('structure')

        # create or read computer and code
        # first check if computer exists already in database
        qb = QueryBuilder()
        qb.append(Computer, tag='computer')
        all_computers = qb.get_results_dict()
        computer_found_in_db = False
        if len(all_computers)>0:
            for icomp in range(len(all_computers)):
                c = all_computers[icomp].get('computer').get('*')
                if c.get_name() == computername:
                    computer_found_in_db = True
                    comp = Computer.from_backend_entity(c)
        # if it is not there create a new one
        if not computer_found_in_db:
            comp = Computer(computername, 'test computer', transport_type='local', scheduler_type='direct', workdir='/temp/ruess/aiida_run_iff734/')
            comp.set_default_mpiprocs_per_machine(4)
            comp.store()
            print 'computer stored now cofigure'
            comp.configure()
        else:
            print 'found computer in database'

        # then get code from database or create a new code
        from aiida.common.exceptions import NotExistent
        try:
            code = Code.get_from_string(voro_codename+'@'+computername)
        except NotExistent as exception:
            code = Code()
            code.label = voro_codename
            code.description = ''
            code.set_remote_computer_exec((comp, '/Users/ruess/sourcecodes/aiida/codes_localhost/voronoi.exe'))
            code.set_input_plugin_name('kkr.voro')
            code.set_prepend_text('ln -s /Users/ruess/sourcecodes/aiida/codes_localhost/ElementDataBase .')
            code.store()

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
        wfd['queue_name'] = queuename
        wfd['resources']['num_machines'] = 1
        params_vorostart = ParameterData(dict=wfd)
       
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        VoroCode = Code.get_from_string(voro_codename+'@'+computername)
       
        # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a ParameterData object for the use in aiida
        ParaNode = ParameterData(dict=kkrparams(LMAX=2, NSPIN=1, RCLUSTZ=1.9).get_dict())
       
        # create process builder to set parameters
        builder = kkr_startpot_wc.get_builder()
        builder.calc_parameters = ParaNode
        builder.description = 'voronoi startpot workflow for Cu bulk'
        builder.label = 'startpot for Cu bulk'
        builder.voronoi = VoroCode
        builder.structure = Cu
        builder.wf_parameters = params_vorostart
       
        # now run calculation
        from aiida.work.launch import run, submit
        out = run(builder)
       
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
