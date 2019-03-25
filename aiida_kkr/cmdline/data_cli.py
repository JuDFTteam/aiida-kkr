from __future__ import absolute_import
from builtins import str
import click
import sys


@click.group()
def cli():
    """Command line interface for template plugin"""
    pass


@cli.command()
def list():  # pylint: disable=redefined-builtin
    """
    Display all KkrstructureData nodes
    """
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()

    from aiida.orm.querybuilder import QueryBuilder
    from aiida.orm import DataFactory
    KKrStructure = DataFactory('kkr.kkrstructure')

    qb = QueryBuilder()
    qb.append(KkrStructure)
    results = qb.all()

    s = ""
    for result in results:
        obj = result[0]
        s += "{}, pk: {}\n".format(str(obj), obj.pk)
    sys.stdout.write(s)
