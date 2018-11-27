#!/usr/bin/env python

import pytest

# some global settings

kkrimp_codename = 'KKRimp'
computername = 'localhost'
queuename = ''

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_scf_workflow():
    """
    Tests for the kkrimp_scf workflow
    """
    
    def test_kkrimp_sub_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node, DataFactory
        from aiida.orm.computers import Computer
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
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
                    print comp
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
            code = Code.get_from_string(kkrimp_codename+'@'+computername)
        except NotExistent as exception:
            code = Code()
            code.label = kkrimp_codename
            code.description = ''
            code.set_remote_computer_exec((comp, '/Users/ruess/sourcecodes/aiida/codes_localhost/kkrflex.exe'))
            code.set_input_plugin_name('kkr.kkrimp')
            code.store()


        options, wfd =kkr_imp_sub_wc.get_wf_defaults()
       
        wfd['nsteps'] = 20
        wfd['strmix'] = 0.05
        options['queue_name'] = queuename
        options['use_mpi'] = True
       
        options = ParameterData(dict=options)
       
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)
       
        # import previous GF writeout
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('de9b5093-25e7-407e-939e-9282c4431343')
       
        # now create a SingleFileData node containing the impurity starting potential
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from numpy import loadtxt
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.out.retrieved.get_abs_path('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = ParameterData(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.out.remote_folder)
       
        label = 'kkrimp_scf Cu host_in_host'
        descr = 'kkrimp_scf workflow for Cu bulk'
       
        # create process builder to set parameters
        builder = kkr_imp_sub_wc.get_builder()
        builder.description = descr
        builder.label = label
        builder.kkrimp = KKRimpCode
        builder.options_parameters = options
        builder.GF_remote_data = GF_host_calc.out.remote_folder
        builder.wf_parameters = ParameterData(dict=wfd)
        builder.host_imp_startpot = startpot_imp_sfd
       
        # now run calculation
        from aiida.work.launch import run, submit
        out = run(builder)
       
        n = out['calculation_info']
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
