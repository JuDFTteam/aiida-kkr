import pytest
#from aiida.utils.fixtures import fixture_manager

@pytest.fixture(scope='session')
def aiida_profile():
    from aiida.utils.fixtures import fixture_manager
    
    with fixture_manager() as fixture_mgr:
        fixture_mgr.create_profile()
        yield fixture_mg
        

#@pytest.fixture(scope='function')
#def test_data(aiida_profile):
#    # load my test data
#    yield
#    aiida_profile.reset_db()

##def test_my_stuff(test_data):
#   # run a test
#   print('test_my_stuf works')
