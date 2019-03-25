"""
Here we define the fixtures for the tests
"""

from __future__ import absolute_import
import pytest
from aiida.manage.fixtures import fixture_manager

@pytest.fixture(scope='session')
def aiida_env():
    with fixture_manager() as manager:
        yield manager

@pytest.fixture()
def fresh_aiida_env(aiida_env):
    yield
    aiida_env.reset_db()


# for computers and codes
@pytest.fixture(scope='session')
def computers_and_codes(aiida_env):
    pass

# for previous data
@pytest.fixture(scope='session')
def import_data(aiida_env):
    from aiida.orm.importexport import import_data
    for db_export_file in ['db_dump_kkrcalc.tar.gz', 'db_dump_kkrflex_create.tar.gz', 'db_dump_vorocalc.tar.gz']:
        import_data('files/'+db_export_file)

