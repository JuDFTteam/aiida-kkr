# -*- coding: utf-8 -*-
"""
Module with CLI commands for various data types.
"""
import click


@click.group('data')
def cmd_data():
    """Commands to create and inspect data nodes."""


'''
@cmd_data.command('list')
def cmd_list_kkrstructures():
    """
    Display all KkrstructureData nodes
    """
    from aiida import load_profile
    load_profile()

    from aiida.orm.querybuilder import QueryBuilder
    from aiida.plugins import DataFactory
    KkrStructure = DataFactory('kkr.kkrstructure')

    qb = QueryBuilder()
    qb.append(KkrStructure)
    results = qb.all()

    s = ""
    for result in results:
        obj = result[0]
        s += "{}, pk: {}\n".format(str(obj), obj.pk)
    click.echo(s)
'''
