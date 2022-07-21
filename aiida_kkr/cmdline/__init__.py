# -*- coding: utf-8 -*-
'''
Module for the command line interface of AiiDA-KKR
'''
import click
import click_completion

from aiida.cmdline.params import options, types
from .launch import cmd_launch
from .data import cmd_data
from .visualization import cmd_plot

# Activate the completion of parameter types provided by the click_completion package
# for bash: eval "$(_AIIDA_KKR_COMPLETE=source_bash aiida-kkr)"
click_completion.init()

# Instead of using entrypoints and directly injecting verdi commands into aiida-core
# we created our own separete CLI because verdi will prob change and become
# less material science specific


@click.group('aiida-kkr', context_settings={'help_option_names': ['-h', '--help']})
@options.PROFILE(type=types.ProfileParamType(load_profile=True))
def cmd_root(profile):  # pylint: disable=unused-argument
    """CLI for the `aiida-kkr` plugin."""


# To avoid circular imports all commands are not yet connected to the root
# but they have to be here because of bash completion

cmd_root.add_command(cmd_launch)
cmd_root.add_command(cmd_data)
cmd_root.add_command(cmd_plot)
