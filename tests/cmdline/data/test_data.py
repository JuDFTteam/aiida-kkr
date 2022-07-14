'''
Module to test the cmds from the data cmd group from the commandline
'''
'''
def test_data_list(run_cli_command):
    """Test data list command"""
    from aiida_kkr.cmdline.data import cmd_list_kkrstructures

    # empty
    run_cli_command(cmd_list_kkrstructures)

    # generate kkrstructure and test again
'''
