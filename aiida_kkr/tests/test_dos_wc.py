#!/usr/bin/env python

import pytest

# some global settings

kkr_codename = 'KKRhost'
computername = 'localhost'
queuename = ''

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_dos_workflow():
    """
    Tests for the kkr_startpot workflow
    """
    
    def test_dos_wc_Cu(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow
        """
        from aiida.orm import Code, load_node, DataFactory
        from aiida.orm.computers import Computer
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.dos import kkr_dos_wc
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
                if c.get_hostname() == computername:
                    computer_found_in_db = True
                    comp = c
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
            code = Code.get_from_string(kkr_codename+'@'+computername)
        except NotExistent as exception:
            code = Code()
            code.label = kkr_codename
            code.description = ''
            code.set_remote_computer_exec((comp, '/Users/ruess/sourcecodes/aiida/codes_localhost/kkr.x'))
            code.set_input_plugin_name('kkr.kkr')
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
        wfd = kkr_dos_wc.get_wf_defaults()
        wfd['dos_params']['kmesh'] = [10, 10, 10]
        wfd['dos_params']['nepts'] = 10
        wfd['queue_name'] = queuename
       
        params_dos = ParameterData(dict=wfd)
       
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRCode = Code.get_from_string(kkr_codename+'@'+computername)
       
        label = 'dos Cu bulk'
        descr = 'DOS workflow for Cu bulk'
       
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')
        kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').out.remote_folder
       
        # create process builder to set parameters
        builder = kkr_dos_wc.get_builder()
        builder.description = descr
        builder.label = label
        builder.kkr = KKRCode
        builder.wf_parameters = params_dos
        builder.remote_data = kkr_calc_remote
       
        # now run calculation
        from aiida.work.launch import run, submit
        out = run(builder)
       
        # check outcome
        n = out['results_wf']
        n = n.get_dict()
        assert n.get('successful')
        assert n.get('list_of_errors') == []
        assert n.get('dos_params').get('nepts') == 10
       
        d = out['dos_data']
        x = d.get_x()
        y = d.get_y()
       
        assert sum(abs(x[1][0] - array([-19.24321191, -16.2197246 , -13.1962373 , -10.17274986, -7.14926255,  -4.12577525,  -1.10228794,   1.9211995 , 4.94468681,   7.96817411]))) < 10**-7
        assert sum(abs(y[0][1][0] - array([  9.86819781e-04,   1.40981029e-03,   2.27894713e-03,         4.79231363e-03,   3.59368494e-02,   2.32929524e+00,         3.06973485e-01,   4.17629157e-01,   3.04021941e-01,         1.24897739e-01]))) < 10**-8

#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_dos_workflow()
   Test.test_dos_wc_Cu()
