#!/usr/bin/env python

import pytest

# some global settings

kkr_codename = 'KKRhost'
kkrimp_codename = 'KKRimp'
voro_codename = 'voronoi'
computername = 'localhost'
queuename = ''
workdir = '/temp/ruess/aiida_run_iff734/'
codelocation = '/Users/ruess/sourcecodes/aiida/codes_localhost'

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_full_workflow():
    """
    Tests for the full kkrimp_scf workflow with GF writeout and voroaux steps
    """
    
    def test_kkrimp_full_wc(self):
        """
        simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
        """
        from aiida.orm import Code, load_node, DataFactory
        from aiida.orm.computers import Computer
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.kkr_imp import kkr_imp_wc
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
       
        # then get code from database or create a new code
        from aiida.common.exceptions import NotExistent
        try:
            code = Code.get_from_string(voro_codename+'@'+computername)
        except NotExistent as exception:
            code = Code()
            code.label = voro_codename
            code.description = ''
            code.set_remote_computer_exec((comp, '/Users/ruess/sourcecodes/aiida/codes_localhost/voronoi.exe'))
            code.set_prepend_text('ln -s /Users/ruess/sourcecodes/aiida/codes_localhost/ElementDataBase .')
            code.set_input_plugin_name('kkr.voro')
            code.store()


        options, wfd, voro_aux_settings =kkr_imp_wc.get_wf_defaults()
       
        wfd['nsteps'] = 20
        wfd['strmix'] = 0.05
        options['queue_name'] = queuename
        options['use_mpi'] = True
        voro_aux_settings['check_dos'] = False
        voro_aux_settings['dos_params']['kmesh'] = [10,10,10]
        voro_aux_settings['dos_params']['nepts'] = 10
        voro_aux_settings['natom_in_cls_min'] = 50
        voro_aux_settings['rclustz'] = 1.5
       
        options = ParameterData(dict=options)
        voro_aux_settings = ParameterData(dict=voro_aux_settings)
        wf_inputs = ParameterData(dict=wfd)
       
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRhostCode = Code.get_from_string(kkr_codename+'@'+computername)
        KKRimpCode = Code.get_from_string(kkrimp_codename+'@'+computername)
        VoroCode = Code.get_from_string(voro_codename+'@'+computername)
       
        imp_info = ParameterData(dict={'Rcut':1.01, 'ilayer_center': 0, 'Zimp':[30.]}) 

        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')
        kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').out.remote_folder
       
        label = 'kkrimp_scf full Cu host_in_host'
        descr = 'kkrimp_scf full workflow for Cu bulk inlcuding GF writeout and vorostart for starting potential'
       
        # create process builder to set parameters
        builder = kkr_imp_wc.get_builder()
        builder.description = descr
        builder.label = label
        builder.kkrimpcode = KKRimpCode
        builder.vorocode = VoroCode
        builder.kkrcode = KKRhostCode
        builder.options_parameters = options
        builder.voro_aux_parameters = voro_aux_settings
        builder.wf_parameters = wf_inputs
        builder.impurity_info = imp_info
        builder.remote_converged_host = kkr_calc_remote
       
        # now run calculation
        from aiida.work.launch import run, submit
        out = run(builder)
       
        # check outcome
        n = out['workflow_info']
        n = n.get_dict()
        for sub in 'auxiliary_voronoi gf_writeout kkr_imp_sub'.split():
            assert sub in n.get('used_subworkflows').keys()
       
        kkrimp_sub = load_node(n['used_subworkflows']['kkr_imp_sub'])
        assert kkrimp_sub.out.calculation_info.get_attr('successful')


#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkrimp_full_workflow()
   Test.test_kkrimp_full_wc()
