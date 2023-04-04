# -*- coding: utf-8 -*-
'''
Module to test the plot cmd from the commandline
'''


def test_cmd_plot_base(run_cli_process_launch_command):
    from aiida_kkr.cmdline.visualization import cmd_plot
    from aiida.orm import Dict
    # test Dict node
    para_node = Dict({'test1': 1})
    options = [para_node.uuid, '--o', 'silent=False']
    run_cli_process_launch_command(cmd_plot, options=options, raises=ValueError)
