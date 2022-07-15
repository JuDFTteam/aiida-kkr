# -*- coding: utf-8 -*-
'''
Module to test all CLI commands to launch for calcjob and workflows of aiida-kkr.

Comment: while 'launch process' is mocked to do nothing these tests are still quite slow
but execute large parts of the workchain code base.
'''
import os
from aiida.orm import Dict
from masci_tools.io.kkr_params import kkrparams

THISFILE_PATH = os.path.dirname(os.path.abspath(__file__))


def test_launch_voro_base(run_cli_process_launch_command, fixture_code):
    """Test invoking the voro launch command with only required inputs."""
    from aiida_kkr.cmdline.launch import launch_voro

    code = fixture_code('kkr.voro').store()
    params = kkrparams(params_type='voronoi')
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    para_node = Dict(dict=params.get_dict()).store()

    options = ['--voro', code.uuid, '--parameters', para_node.uuid]
    run_cli_process_launch_command(launch_voro, options=options)


def test_launch_bs(run_cli_process_launch_command, fixture_code, generate_remote_data):
    """A test for invoking the bs launch commmand with required inputs."""
    from aiida_kkr.cmdline.launch import launch_bs
    from aiida.engine.processes.calcjobs.tasks import PreSubmitException
    code = fixture_code('kkr.kkr').store()
    params = kkrparams(params_type='kkr')
    params.set_multiple_values(
        MIN=-10,
        MAX=5,
        NPT2=12,
        RCLUSTZ=2.3,
        TEMPR=50,
    )

    param_node = Dict(dict=params.get_dict()).store()
    path = os.path.abspath(os.path.join(THISFILE_PATH, '../../files/bd_dump_bs/parent_kkr'))
    remote = generate_remote_data(code.computer, path).store()
    options = ['--kkr', code.uuid, '--parameters', param_node.uuid, '--parent-folder', remote.uuid]
    run_cli_process_launch_command(launch_bs, options=options, raises=PreSubmitException)


def test_launch_kkr_base(run_cli_process_launch_command, fixture_code, generate_remote_data):
    """Test invoking the kkr launch command with only required inputs."""
    from aiida_kkr.cmdline.launch import launch_kkr
    from aiida.engine.processes.calcjobs.tasks import PreSubmitException

    code = fixture_code('kkr.kkr').store()
    params = kkrparams(params_type='voronoi')
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    para_node = Dict(dict=params.get_dict()).store()
    path = os.path.abspath(os.path.join(THISFILE_PATH, '../../files/voronoi'))
    remote = generate_remote_data(code.computer, path).store()

    # TODO generate real full Voronoi parent kkr gets at the inputs...
    #
    options = ['--kkr', code.uuid, '--parameters', para_node.uuid, '-P', remote.uuid]
    # This will raise, because parent is only a folder.
    run_cli_process_launch_command(launch_kkr, options=options, raises=PreSubmitException)


'''
def test_launch_kkrimp_base(run_cli_process_launch_command, fixture_code, generate_remote_data):
    """Test invoking the kkrimp launch command with only required inputs."""
    from aiida_kkr.cmdline.launch import launch_kkrimp
    from aiida.engine.processes.calcjobs.tasks import PreSubmitException

    code = fixture_code('kkr.imp').store()
    params = kkrparams(params_type='voronoi')
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    para_node = Dict(dict=params.get_dict()).store()
    path = os.path.abspath(os.path.join(THISFILE_PATH, '../../files/voronoi'))
    remote = generate_remote_data(code.computer, path).store()

    # TODO generate real full Voronoi parent kkr gets at the inputs...
    #
    options = ['--kkrimp', code.uuid, '--parameters', para_node.uuid, '-P', remote.uuid]
    # This will raise, because parent is only a folder.
    run_cli_process_launch_command(launch_kkr, options=options, raises=PreSubmitException)
'''
