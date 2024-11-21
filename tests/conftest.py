"""
Here we define the fixtures for the tests
"""

import pytest
import tempfile
import shutil
import pathlib
from aiida import __version__ as aiida_core_version
from aiida.orm import RemoteData, CalcJobNode
from aiida.common.hashing import make_hash
# from aiida.manage.tests.pytest_fixtures import aiida_profile, temp_dir
# from aiida.tools.pytest_fixtures import *
import aiida_kkr

pytest_plugins = [
    'aiida.manage.tests.pytest_fixtures',
    # 'aiida_testing.mock_code',
    # 'aiida_testing.export_cache',
]

# pytest_plugins = 'aiida.tools.pytest_fixtures'
# pytest_plugins = 'aiida.tools.pytest_fixtures'

# test settings:
# paths where the tests are located and where the test input data is stored
test_dir = pathlib.Path(aiida_kkr.__file__).parent.parent / 'tests'
data_dir = (test_dir / 'data_dir')  # TODO: get from config?

# fixtures

# @pytest.fixture(scope='function', autouse=True)
# def clear_database_auto(clear_database):
#     """Automatically clear database in between tests."""
#     pass


# need fixed aiida_localhost to have set_default_mpiprocs_per_machine set to 1
@pytest.fixture(scope='function')
def aiida_localhost_serial(tmp_path):  # pylint: disable=redefined-outer-name
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
        computer = Computer.objects.get(label=name)
    except NotExistent:
        if int(aiida_core_version.split('.')[0]) < 2:
            # old name before aiida-core v2.0
            transport_type = 'local'
            scheduler_type = 'direct'
        else:
            # new names in aiida-core >= 2.0
            transport_type = 'core.local'
            scheduler_type = 'core.direct'
        computer = Computer(
            label=name,
            description='localhost computer set up by test manager',
            hostname=name,
            workdir=str(tmp_path),
            transport_type=transport_type,
            scheduler_type=scheduler_type
        )
        computer.store()
        computer.set_minimum_job_poll_interval(0.)
        computer.set_default_mpiprocs_per_machine(1)
        computer.configure()

    return computer


@pytest.fixture(scope='function')
def aiida_local_code_factory_prepend(aiida_localhost_serial):  # pylint: disable=redefined-outer-name
    """Get an AiiDA code on localhost.

    Searches in the PATH for a given executable and creates an AiiDA code with provided entry point.

    Usage::

      def test_1(aiida_local_code_factory):
          code = aiida_local_code_factory('pw.x', 'quantumespresso.pw')
          # use code for testing ...

    :return: A function get_code(executable, entry_point) that returns the Code node.
    :rtype: object
    """

    def get_code(entry_point, executable, computer=aiida_localhost_serial, prepend_text=None):
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
            raise ValueError(f'The executable "{executable}" was not found in the $PATH.')

        code = Code(input_plugin_name=entry_point, remote_computer_exec=[computer, executable_path])
        code.label = executable
        if prepend_text is not None:
            code.set_prepend_text(prepend_text)
        return code.store()

    return get_code


@pytest.fixture(scope='function')
def reuse_local_code(aiida_local_code_factory_prepend):

    def _get_code(executable, exec_relpath, entrypoint, prepend_text=None, use_export_file=True):
        import os
        from aiida.orm import ProcessNode, QueryBuilder, Code, load_node

        try:
            # this is for aiida-core < 2.0
            from aiida.tools.importexport import import_data, export
            aiida_core_ver_old = True
        except ImportError:
            # This is the case for aiida >= 2.0.0
            from aiida.tools.archive import import_archive
            import_kwargs = {'group': None}
            aiida_core_ver_old = False

        full_import_path = str(data_dir) + '/' + executable + '.tar.gz'

        # check if exported code exists and load it, otherwise create new code (will have different has due to different working directory)
        if use_export_file and pathlib.Path(full_import_path).exists():
            if aiida_core_ver_old:
                import_data(full_import_path, silent=True)
            else:
                import_archive(archive_path, **import_kwargs)
            codes = Code.objects.find(filters={'label': executable})  # pylint: disable=no-member
            code = codes[0]
            code.computer.configure()

        else:
            # make sure code is found in PATH
            _exe_path = (test_dir / pathlib.Path(exec_relpath)).absolute()
            print(_exe_path)
            os.environ['PATH'] += ':' + str(_exe_path)
            # get code using aiida_local_code_factory fixture
            code = aiida_local_code_factory_prepend(entrypoint, executable, prepend_text=prepend_text)

            if use_export_file and aiida_core_ver_old:
                #export for later reuse
                export([code], outfile=full_import_path, overwrite=True)  # add export of extras automatically

        return code

    return _get_code


@pytest.fixture(scope='function')
def voronoi_local_code_import(reuse_local_code):
    """
    Create or load KKRhost code
    """
    import os
    executable = 'voronoi.exe'  # name of the Voronoi executable
    exec_rel_path = 'jukkr/'  # location where it is found
    entrypoint = 'kkr.voro'  # entrypoint
    # prepend text to be added before execution
    prepend_text = f"""
ulimit -s unlimited
ln -s {os.path.abspath(exec_rel_path)}/ElementDataBase ."""
    voro_code = reuse_local_code(executable, exec_rel_path, entrypoint, prepend_text)

    return voro_code


@pytest.fixture(scope='function')
def voronoi_local_code(reuse_local_code):
    """
    Create or load KKRhost code
    """
    import os
    executable = 'voronoi.exe'  # name of the Voronoi executable
    exec_rel_path = 'jukkr/'  # location where it is found
    entrypoint = 'kkr.voro'  # entrypoint
    # prepend text to be added before execution
    prepend_text = f"""
ulimit -s unlimited
ln -s {os.path.abspath(exec_rel_path)}/ElementDataBase .
"""
    voro_code = reuse_local_code(executable, exec_rel_path, entrypoint, prepend_text, use_export_file=False)

    return voro_code


@pytest.fixture(scope='function')
def kkrhost_local_code(reuse_local_code):
    """
    Create or load KKRhost code
    """
    executable = 'kkr.x'  # name of the KKRhost executable
    exec_rel_path = 'jukkr/'  # location where it is found
    entrypoint = 'kkr.kkr'  # entrypoint
    # prepend text to be added before execution
    prepend_text = """
ulimit -s unlimited
export OMP_STACKSIZE=2G
"""
    kkrhost_code = reuse_local_code(executable, exec_rel_path, entrypoint, prepend_text, use_export_file=False)

    return kkrhost_code


@pytest.fixture(scope='function')
def kkrimp_local_code(reuse_local_code):
    """
    Create or load KKRimp code
    """
    executable = 'kkrflex.exe'  # name of the KKRimp executable
    exec_rel_path = 'jukkr/'  # location where it is found
    entrypoint = 'kkr.kkrimp'  # entrypoint
    # prepend text to be added before execution
    prepend_text = """
ulimit -s unlimited
export OMP_STACKSIZE=2G
"""
    kkrimp_code = reuse_local_code(executable, exec_rel_path, entrypoint, prepend_text, use_export_file=False)

    return kkrimp_code


@pytest.fixture
def generate_remote_data():
    """Return a `RemoteData` node."""

    def _generate_remote_data(computer, remote_path, entry_point_name=None):
        """Return a `RemoteData` node which points to some dir."""
        from aiida.common.links import LinkType
        from aiida.plugins.entry_point import format_entry_point_string

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        remote = RemoteData(remote_path=remote_path)
        remote.computer = computer

        if entry_point_name is not None:
            creator = CalcJobNode(computer=computer, process_type=entry_point)
            creator.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
            remote.add_incoming(creator, link_type=LinkType.CREATE, link_label='remote_folder')
            creator.store()

        return remote

    return _generate_remote_data


def import_with_migration(archive_path):
    """Import aiida export file and try migration if version is incompatible"""
    try:
        # this is the case for aiida-core<=1.6
        from aiida.tools.importexport import (
            detect_archive_type, EXPORT_VERSION, import_data, IncompatibleArchiveVersionError
        )
        from aiida.tools.importexport.archive.migrators import get_migrator
        from aiida.common.folders import SandboxFolder
        import_kwargs = dict(extras_mode_existing='nnl', silent=True)
        try:
            imported_nodes = import_data(archive_path, **import_kwargs)
        except IncompatibleArchiveVersionError as exception:
            print('incompatible version detected for import file, trying migration')
            with SandboxFolder() as temp_folder:
                try:
                    migrator = get_migrator(detect_archive_type(archive_path))(archive_path)
                    archive_path = migrator.migrate(
                        EXPORT_VERSION, None, out_compression='none', work_dir=temp_folder.abspath
                    )
                except Exception as exception:
                    print('an exception occurred while migrating the archive', exception)

                print('proceeding with import of migrated archive')
                try:
                    imported_nodes = import_data(archive_path, **import_kwargs)
                except Exception as exception:
                    print('an exception occurred while trying to import the migrated archive', exception)

    except ImportError:
        # This is the case for aiida >= 2.0.0
        from click import echo
        from aiida.tools.archive import import_archive, get_format
        from aiida.common.exceptions import IncompatibleStorageSchema

        import_kwargs = {'group': None}

        try:
            imported_nodes = import_archive(archive_path, **import_kwargs)
        except IncompatibleStorageSchema:
            echo(f'incompatible version detected for {archive_path}, trying migration')
            archive_format = get_format()
            version = archive_format.latest_version
            archive_format.migrate(archive_path, archive_path, version, force=True, compression=6)
            imported_nodes = import_archive(archive_path, **import_kwargs)

    return imported_nodes
