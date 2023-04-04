# -*- coding: utf-8 -*-
'''
This module contains click option types specific to aiida-kkr
'''
import click
from aiida.cmdline.params import types
from aiida.cmdline.utils import echo
from aiida.common.exceptions import NotExistent
from aiida.plugins import DataFactory


class StructureNodeOrFileParamType(click.ParamType):
    """
    The ParamType for identifying a structure by node or to extract it from a given file

    Pro: It is convenient
    Con: If users only use other formats to launch their workflows it will create many
    more structures in the database.
    """

    name = 'StructureFile'

    def convert(self, value, param, ctx):
        is_path = False
        # Alternative one could check if int or uuid
        # aiida allows also for shorten uuids

        try:
            structure = types.DataParamType(sub_classes=('aiida.data:core.structure',)).convert(value, param, ctx)
        except (NotExistent, click.exceptions.BadParameter) as er:
            echo.echo(
                f'Tried to load node, could not fine one for {value}. '
                'I will further check if it is a filepath.'
            )
            is_path = True

        if is_path:
            StructureData = DataFactory('core.structure')
            # If it is a path to a file try to convert the structure
            pathtype = click.Path(exists=True, dir_okay=False, resolve_path=True)
            filename = pathtype.convert(value, param, ctx)
            try:
                import ase.io
            except ImportError:
                echo.echo_critical('You have not installed the package ase. \nYou can install it with: pip install ase')

            try:
                asecell = ase.io.read(filename)
                structure = StructureData(ase=asecell)
            except ValueError as err:
                echo.echo_critical(str(err))
            # do not store structure, since this option is for calculation and workflow
            # input, which will store the structure anyway.
        return structure
