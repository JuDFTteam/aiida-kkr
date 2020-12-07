# -*- coding: utf-8 -*-
"""Tests if all available CLI commands can print help."""
from click import Context, Group
from aiida_kkr.cmdline import cmd_root


def test_base_all_commands():
    """Test that all commands in ``cmd_root`` are reachable and can print the help message.

    This doesn't guarantee that the command works but at least that it can be successfully called and there are no
    import errors or other basic problems.
    """

    def recursively_print_help(ctx):
        assert isinstance(ctx.get_help(), str)

        if isinstance(ctx.command, Group):
            for subcommand in ctx.command.commands.values():
                ctx.command = subcommand
                recursively_print_help(ctx)

    ctx = Context(cmd_root)
    recursively_print_help(ctx)
