# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c), Forschungszentrum Jülich GmbH, IAS-1/PGI-1, Germany.         #
#                All rights reserved.                                         #
# This file is part of the AiiDA-KKR package.                                 #
#                                                                             #
# For further information on the license, see the LICENSE.txt file            #
###############################################################################
# pylint: disable=cyclic-import
"""Command line utilities to create and inspect `Dict` nodes with FLAPW parameters."""
import click
from aiida.cmdline.params import options
from aiida.cmdline.utils import decorators, echo
from . import cmd_data
from aiida.orm import Dict
from aiida_kkr.tools import kkrparams


@cmd_data.group('parameter')
def cmd_parameter():
    """Commands to create and inspect `Dict` nodes containing FLAPW parameters ('calc_parameters')."""


@cmd_parameter.command('import')
@click.argument('filename', type=click.Path(exists=True))
@click.option('--show/--no-show', default=True, show_default=True, help='Print the contents from the extracted dict.')
@click.option('--verbose', '-v', is_flag=True, default=False, show_default=True, help='Print verbose output.')
@options.DRY_RUN()
@decorators.with_dbenv()
def cmd_param_import(filename, show, verbose, dry_run):
    """
    Import kkr params Dict from a KKRhost input file.

    FILENAME is the name/path of the inputcard file to use.
    """

    kkr_para = kkrparams()
    if verbose:
        print(f'start reading from file: {filename}')
    kkr_para.read_keywords_from_inputcard(filename)
    if verbose:
        print(f'read kkr params: {kkr_para.get_set_values()}')
    missing = kkr_para.get_missing_keys(use_aiida=True)
    if len(missing) > 0:
        # TODO replace by click echo
        echo.echo_warning(
            f"""inputcard did not contain all necessary keys to run a kkr calculation!
The following keys are missing:
{missing}"""
        )

    struc_keys = ['BRAVAIS', '<RBASIS>', '<ZATOM>', 'ALATBASIS', 'NAEZ', 'NATYP']
    no_struc_kkrparams = {k: v for k, v in kkr_para.get_set_values() if k not in struc_keys}
    if verbose:
        print(f'keys after taking out the ones used in the structure node:: {no_struc_kkrparams}')

    if len(no_struc_kkrparams.keys()) == 0:
        echo.echo_critical('failed to extract kkr params')
    else:
        param_node = Dict(no_struc_kkrparams)
        if dry_run:
            echo.echo_success('parsed kkr params from input file')
        else:
            param_node.store()
            message = f'stored kkr params Dict node <{param_node.pk}>'
            echo.echo_success(message)

    if show:
        echo.echo_dictionary(param_node.get_dict())
