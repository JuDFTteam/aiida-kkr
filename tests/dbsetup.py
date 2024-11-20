import os

# some global settings

kkr_codename = 'kkrhost'
kkrimp_codename = 'kkrimp'
voro_codename = 'voronoi'
computername = 'localhost'
queuename = ''

cwd = os.getcwd()
workdir = cwd + '/jukkr/work'
codelocation = cwd + '/jukkr/'

# change kkr_condename for testing (on mac)
# also used for caching
kkr_codename = 'kkrhost_intel19'
kkrimp_codename = 'kkrimp_intel19'


def prepare_computer(computername, workdir):
    """Create new computer in db or read computer from db if it already exists."""
    from aiida.orm import Computer
    from aiida.orm.querybuilder import QueryBuilder

    # first check if computer exists already in database
    qb = QueryBuilder()
    qb.append(Computer, tag='computer')
    all_computers = qb.dict()
    computer_found_in_db = False
    if len(all_computers) > 0:
        for icomp in range(len(all_computers)):
            c = all_computers[icomp].get('computer').get('*')
            if c.get_name() == computername:
                computer_found_in_db = True
                comp = c
    # if it is not there create a new one
    if not computer_found_in_db:
        comp = Computer(computername, 'test computer', transport_type='local', scheduler_type='direct', workdir=workdir)
        comp.set_default_mpiprocs_per_machine(4)
        comp.store()
        print('computer stored now cofigure')
    else:
        print('found computer in database')
    # configure for cases where computer was imported
    comp.configure()
    comp.set_minimum_job_poll_interval(0.)
    # return computer
    return comp


def prepare_code(codename, codelocation, computername, workdir):
    """Prepare a code, either create entry in AiiDA DB or load it from DB."""
    # first create or read computer
    comp = prepare_computer(computername, workdir)
    # now decide which code to add
    if 'kkrhost' in codename:
        execname = 'kkr.x'
        pluginname = 'kkr.kkr'
    elif 'voronoi' in codename:
        execname = 'voronoi.exe'
        pluginname = 'kkr.voro'
    elif 'kkrimp' in codename:
        execname = 'kkrflex.exe'
        pluginname = 'kkr.kkrimp'
    else:
        raise ValueError('unknown codename')
    # then get code from database or create a new code
    from aiida.orm import Code
    from aiida.common.exceptions import NotExistent
    try:
        code = Code.get_from_string(codename + '@' + computername)
    except NotExistent as exception:
        code = Code()
        code.label = codename
        code.description = ''
        code.set_remote_computer_exec((comp, codelocation + execname))
        code.set_input_plugin_name(pluginname)
        if 'voronoi' in codename:
            code.set_prepend_text('ln -s ' + codelocation + 'ElementDataBase .')
        code.store()
