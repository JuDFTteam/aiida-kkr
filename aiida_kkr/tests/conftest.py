"""
Here we define the fixtures for the tests
"""

from __future__ import absolute_import
import pytest
from aiida.common.hashing import make_hash


import tempfile
import shutil
#from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory



@pytest.fixture(scope='function')
def temp_dir_session():
    """Get a temporary directory.

    E.g. to use as the working directory of an AiiDA computer.

    :return: The path to the directory
    :rtype: str

    """
    try:
        dirpath = tempfile.mkdtemp()
        yield dirpath
    finally:
        # after the test function has completed, remove the directory again
        shutil.rmtree(dirpath)



@pytest.fixture(scope='function')
def aiida_localhost_session(temp_dir_session):  # pylint: disable=redefined-outer-name
    """Get an AiiDA computer for localhost.

    Usage::

      def test_1(aiida_localhost):
          name = aiida_localhost.get_name()
          # proceed to set up code or use 'aiida_local_code_factory' instead


    :return: The computer node
    :rtype: :py:class:`aiida.orm.Computer`
    """
    from aiida.orm import Computer
    from aiida.common.exceptions import NotExistent

    name = 'localhost-test'

    try:
        computer = Computer.objects.get(name=name)
    except NotExistent:
        computer = Computer(
            name=name,
            description='localhost computer set up by test manager',
            hostname=name,
            workdir=temp_dir_session,
            transport_type='local',
            scheduler_type='direct'
        )
        computer.store()
        computer.set_minimum_job_poll_interval(0.)
        computer.set_default_mpiprocs_per_machine(1)
        computer.configure()

    return computer



@pytest.fixture(scope='function')
def aiida_local_code_factory_session(aiida_localhost_session):  # pylint: disable=redefined-outer-name
    """Get an AiiDA code on localhost.

    Searches in the PATH for a given executable and creates an AiiDA code with provided entry point.

    Usage::

      def test_1(aiida_local_code_factory):
          code = aiida_local_code_factory('pw.x', 'quantumespresso.pw')
          # use code for testing ...

    :return: A function get_code(executable, entry_point) that returns the Code node.
    :rtype: object
    """

    def get_code(entry_point, executable, computer=aiida_localhost_session, prepend_text=None):
        """Get local code.
        Sets up code for given entry point on given computer.

        :param entry_point: Entry point of calculation plugin
        :param executable: name of executable; will be searched for in local system PATH.
        :param computer: (local) AiiDA computer
        :return: The code node
        :rtype: :py:class:`aiida.orm.Code`
        """
        from aiida.orm import Code

        codes = Code.objects.find(filters={'label': executable})  # pylint: disable=no-member
        if codes:
            return codes[0]

        executable_path = shutil.which(executable)

        if not executable_path:
            raise ValueError('The executable "{}" was not found in the $PATH.'.format(executable))

        code = Code(
            input_plugin_name=entry_point,
            remote_computer_exec=[computer, executable_path]
        )
        code.label = executable
        if prepend_text is not None:
            code.set_prepend_text(prepend_text)
        return code.store()

    return get_code


@pytest.fixture(scope='function')
def voronoi_local_code(aiida_local_code_factory_session):
    import os, pathlib
    from aiida.tools.importexport import import_data, export
    from aiida.orm import ProcessNode, QueryBuilder, Code, load_node

    cwd = pathlib.Path(os.getcwd())                  # Might be not the best idea.
    data_dir = (cwd / 'data_dir')                    # TODO: get from config?
    full_import_path = str(data_dir)+'/voro_code.tar.gz'
    # check if exported code exists and load it, otherwise create new code (will have different has due to different working directory)
    if pathlib.Path(full_import_path).exists():
        import_data(full_import_path)
        codes = Code.objects.find(filters={'label': 'voronoi.exe'})  # pylint: disable=no-member
        voro_code = codes[0]
        voro_code.computer.configure()

    else:
        # get voro code using aiida_local_code_factory fixture
        _voro_exe_path = os.path.abspath('jukkr')
        os.environ['PATH']+=':'+_voro_exe_path
        voro_code = aiida_local_code_factory_session('kkr.voro', 'voronoi.exe', prepend_text="ln -s "+_voro_exe_path+"/ElementDataBase .")
        
        #export for later use
        export([voro_code], outfile=full_import_path, overwrite=True) # add export of extras automatically

    return voro_code

@pytest.fixture(scope='function')
def kkrhost_local_code(aiida_local_code_factory_session):
    import os, pathlib
    from aiida.tools.importexport import import_data, export
    from aiida.orm import ProcessNode, QueryBuilder, Code, load_node

    cwd = pathlib.Path(os.getcwd())                  # Might be not the best idea.
    data_dir = (cwd / 'data_dir')                    # TODO: get from config?
    full_import_path = str(data_dir)+'/kkrhost_code.tar.gz'
    # check if exported code exists and load it, otherwise create new code (will have different has due to different working directory)
    if pathlib.Path(full_import_path).exists():
        import_data(full_import_path)
        codes = Code.objects.find(filters={'label': 'kkr.x'})  # pylint: disable=no-member
        kkrhost_code = codes[0]
        kkrhost_code.computer.configure()

    else:
        # get voro code using aiida_local_code_factory fixture
        _kkr_exe_path = os.path.abspath('jukkr/build_new_kkrhost/')
        os.environ['PATH']+=':'+_kkr_exe_path
        kkrhost_code = aiida_local_code_factory_session('kkr.kkr', 'kkr.x', prepend_text="""
        ulimit -s hard
        export OMP_STACKSIZE=2G
        source compiler-select intel
        """)
        
        #export for later use
        export([kkrhost_code], outfile=full_import_path, overwrite=True) # add export of extras automatically

    return kkrhost_code

    

