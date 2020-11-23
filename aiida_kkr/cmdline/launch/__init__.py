# -*- coding: utf-8 -*-
'''
Module with CLI commands for calcjob types of aiida-kkr.
'''
import click
from .launch import launch_voro
from .launch import launch_kkr
from .launch import launch_kkr_imp
#from .launch import launch_dos
#from .launch import launch_eos
#from .launch import launch_gf_writeout
#from .launch import launch_kkr_imp_dos
#from .launch import launch_kkr_imp_wc
#from .launch import launch_kkr_imp_dos
#from .launch import launch_kkr_imp_sub
#from .launch import launch_kkr_scf
#from .launch import launch_vor_start


@click.group('launch')
def cmd_launch():
    """Commands to launch workflows and calcjobs of aiida-kkr."""


# we do it like this and not in short from to avoid cyclic imports
# and get the full bash completion working
cmd_launch.add_command(launch_voro)
cmd_launch.add_command(launch_kkr)
cmd_launch.add_command(launch_kkr_imp)
#cmd_launch.add_command(launch_dos)
#cmd_launch.add_command(launch_eos)
#cmd_launch.add_command(launch_gf_writeout)
#cmd_launch.add_command(launch_kkr_imp_dos)
#cmd_launch.add_command(launch_kkr_imp_wc)
#cmd_launch.add_command(launch_kkr_imp_dos)
#cmd_launch.add_command(launch_kkr_imp_sub)
#cmd_launch.add_command(launch_kkr_scf)
#cmd_launch.add_command(launch_vor_start)
