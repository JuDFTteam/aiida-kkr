from __future__ import print_function

from builtins import range
import os


# some global settings

kkr_codename = 'kkrhost'
kkrimp_codename = 'kkrimp'
voro_codename = 'voronoi'
computername = 'localhost'
queuename = ''

cwd = os.getcwd()
workdir = cwd+'/jukkr/work'
codelocation = cwd+'/jukkr/'



def prepare_computer(computername, workdir):
    """Create new computer in db or read computer from db if it already exists."""
    from aiida.orm import Computer
    from aiida.orm.backend import construct_backend
    from aiida.orm.querybuilder import QueryBuilder
    # 
    # first check if computer exists already in database
    qb = QueryBuilder()
    qb.append(Computer, tag='computer')
    all_computers = qb.get_results_dict()
    computer_found_in_db = False
    if len(all_computers)>0:
        for icomp in range(len(all_computers)):
            c = all_computers[icomp].get('computer').get('*')
            if c.get_name() == computername:
                computer_found_in_db = True
                comp = c 
    # if it is not there create a new one
    if not computer_found_in_db:
        #comp = Computer(computername, 'test computer', transport_type='local', scheduler_type='direct', workdir=workdir)
        b = construct_backend()
        comp = b.computers.create(computername, 'test computer', transport_type='local', scheduler_type='direct', workdir=workdir)
        comp.set_default_mpiprocs_per_machine(4)
        comp.store()
        print('computer stored now cofigure')
        comp.configure()
    else:
        print('found computer in database')
    # return computer
    return comp


def prepare_code(codename, codelocation, computername, workdir):
    """."""
    # first create or read computer
    comp = prepare_computer(computername, workdir)
    # now decide which code to add
    if codename == 'kkrhost':
        execname = 'kkr.x'
        pluginname = 'kkr.kkr'
    elif codename == 'voronoi':
        execname = 'voronoi.exe'
        pluginname = 'kkr.voro'
    elif codename == 'kkrimp':
        execname = 'kkrflex.exe'
        pluginname = 'kkr.kkrimp'
    else:
        raise ValueError('unknown codename')
    # then get code from database or create a new code
    from aiida.orm import Code
    from aiida.common.exceptions import NotExistent
    try:
        code = Code.get_from_string(codename+'@'+computername)
    except NotExistent as exception:
        code = Code()
        code.label = codename
        code.description = ''
        code.set_remote_computer_exec((comp, codelocation+execname))
        code.set_input_plugin_name(pluginname)
        if codename == 'voronoi':
            code.set_prepend_text('ln -s '+codelocation+'ElementDataBase .')
        code.store()

