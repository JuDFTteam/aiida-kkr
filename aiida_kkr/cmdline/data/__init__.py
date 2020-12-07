# -*- coding: utf-8 -*-
"""
Module with CLI commands for various data types.
"""
import click


@click.group('data')
def cmd_data():
    """Commands to create and inspect data nodes."""

#TODO
# aiida-kkr data create kkrparams: create a Dict node using the kkrparams class to have help and verification
# aiida-kkr data create options: create options node for workflows (withmpi, resources etc.)
# aiida-kkr data create impinfo: create imputiy info node which has help and verification
# aiida-kkr data create kpoints: create k-points for bandstructure claculations
# aiida-kkr data create noco_angles: create k-points for bandstructure claculations
#also corresponding list and show commands needed
# for list and show we should have additionally imp pot, maybe something else 

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
