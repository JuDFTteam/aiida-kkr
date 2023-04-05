# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c), Forschungszentrum Jülich GmbH, IAS-1/PGI-1, Germany.         #
#                All rights reserved.                                         #
# This file is part of the AiiDA-KKR package.                                 #
#                                                                             #
# For further information on the license, see the LICENSE.txt file            #
###############################################################################
"""Command line utilities to create and inspect `StructureData` nodes."""
import click
from aiida.cmdline.params import options
from aiida.cmdline.utils import decorators, echo
from aiida.orm import Dict
from aiida_kkr.tools import kkrparams
from aiida_kkr.tools.common_workfunctions import structure_from_params
from . import cmd_data


@click.group('structure')
def cmd_structure():
    """Commands to create and inspect `StructureData` nodes."""


cmd_data.add_command(cmd_structure)


@cmd_structure.command('import')
@click.argument('filename', type=click.Path(exists=True))
@click.option(
    '--kkrpara/--no-kkrpara', default=False, show_default=True, help='Store also the kkrparams Dict in the db.'
)
@click.option('--verbose', '-v', is_flag=True, default=False, show_default=True, help='Print verbose output.')
@options.DRY_RUN()
@decorators.with_dbenv()
def cmd_import(filename, kkrpara, verbose, dry_run):
    """
    Import a `StructureData` from a KKRhost input file.

    FILENAME is the name/path of the inputcard file to use.

    If you want to import a structure from another file type you can use
    'verdi data structure import -ase <filename>' instead.
    """

    kkr_para = kkrparams()
    if verbose:
        print(f'start reading from file: {filename}')
    kkr_para.read_keywords_from_inputcard(filename)
    if verbose:
        print(f'read kkr params: {kkr_para.get_set_values()}')

    # extract kkr params
    missing = kkr_para.get_missing_keys(use_aiida=True)
    if kkrpara:
        if len(missing) > 0:
            # TODO replace by click echo
            echo.echo_warning(
                f"""inputcard did not contain all necessary keys to run a kkr calculation!
The following keys are missing:
{missing}"""
            )

        struc_keys = ['BRAVAIS', '<RBASIS>', '<ZATOM>', 'ALATBASIS', 'NAEZ', 'NATYP']
        no_struc_kkrparams = {k: v for k, v in kkr_para.get_set_values() if k not in struc_keys}
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

    # extract structure
    success, struc_node = structure_from_params(kkr_para)

    if success:
        formula = struc_node.get_formula()
        if dry_run:
            echo.echo_success(f'parsed structure from kkr params from input file with formula {formula}')
        else:
            struc_node.store()
            message = f'parsed and stored StructureData <{struc_node.pk}> with formula {formula}'
            echo.echo_success(message)
    else:
        echo.echo_critical('Error in getting the structure, check the input file.')
