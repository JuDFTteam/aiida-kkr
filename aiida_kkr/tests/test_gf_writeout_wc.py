#!/usr/bin/env python

import pytest

# some global settings

kkr_codename = 'KKRhost'
computername = 'localhost'
queuename = ''
workdir = '/temp/ruess/aiida_run_iff734/'
codelocation = '/Users/ruess/sourcecodes/aiida/codes_localhost'

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_gf_writeout_workflow():
    """
    Tests for the kkr_startpot workflow
    """
    
    def test_kkrflex_writeout_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow
        """
        from aiida.orm import Code, load_node, DataFactory
        from aiida.orm.computers import Computer
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
        from numpy import array
        import os
       
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
            code = Code.get_from_string(kkr_codename+'@'+computername)
        except NotExistent as exception:
            code = Code()
            code.label = kkr_codename
            code.description = ''
            code.set_remote_computer_exec((comp, '/Users/ruess/sourcecodes/aiida/codes_localhost/kkr.x'))
            code.set_input_plugin_name('kkr.kkr')
            code.store()

       
        # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
        wfd =kkr_flex_wc.get_wf_defaults()
        wfd['queue_name'] = queuename
        wfd['use_mpi'] = True
        options = ParameterData(dict=wfd)
       
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRCode = Code.get_from_string(kkr_codename+'@'+computername)
       
        imp_info = ParameterData(dict={'Rcut':1.01, 'ilayer_center': 0, 'Zimp':[79.]})
       
        label = 'GF_writeout Cu bulk'
        descr = 'GF_writeout workflow for Cu bulk'
       
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')
        kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').out.remote_folder
       
        # create process builder to set parameters
        builder = kkr_flex_wc.get_builder()
        builder.description = descr
        builder.label = label
        builder.kkr = KKRCode
        builder.options_parameters = options
        builder.remote_data = kkr_calc_remote
        builder.imp_info = imp_info
       
        # now run calculation
        from aiida.work.launch import run, submit
        out = run(builder)
       
        n = out['calculation_info']
        n = n.get_dict()
       
        assert n.get('successful')
        assert n.get('list_of_errors') == []
       
        d = out['GF_host_remote']
        assert isinstance(d, DataFactory('remote'))
       
        kkrflex_retrieved = load_node(n.get('pk_flexcalc'))
        kkrflex_retrieved = kkrflex_retrieved.out.retrieved
        kkrflex_path = kkrflex_retrieved.get_abs_path('')
        for name in 'tmat green atominfo intercell_cmoms intercell_ref'.split():
            assert 'kkrflex_'+name in os.listdir(kkrflex_path)

#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_gf_writeout_workflow()
   Test.test_kkrflex_writeout_wc()
