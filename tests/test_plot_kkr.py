#!/usr/bin/env python

from builtins import object
import pytest
import matplotlib

matplotlib.use('Agg')  # needs to be done before pyplot is imported
from matplotlib.pyplot import gcf, title
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test
from .conftest import import_with_migration


@pytest.mark.usefixtures('aiida_profile', 'clear_database')
class Test_plot_kkr(object):
    """
    Tests for the plot_kkr tool
    """

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='kkr.png', tolerance=8)
    def test_plot_kkr_calc(self):
        basic_test('1792144e-746c-4575-a1e1-40125a26778c', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='voro.png')
    def test_plot_voro_calc(self):
        basic_test('cd88cad8-16a0-4eb4-b3fb-8887a857c376', strucplot=False)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='dos.png', remove_text=True)
    def test_plot_dos(self):
        basic_test('d96d2db2-8276-4fec-836e-789e8710c487', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='qdos.png', remove_text=True)
    def test_plot_qdos(self):
        import_with_migration('files/db_dump_kkrcalc_qdos.tar.gz')
        basic_test('a0d0d29f-7b22-4ca4-ba55-6b97569d94af', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='vorostart.png', tolerance=5)
    def test_plot_vorostart_wc(self):
        basic_test('0c52eff5-3c5a-4623-a278-8febab698d30', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='scf.png', tolerance=8)
    def test_plot_scf_wc(self):
        basic_test('224e0f00-6f81-4e63-a142-d45e83ec33e8', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='scf_grouped.png', tolerance=8)
    def test_plot_scf_wc_grouped(self):
        basic_test([
            'bf644c1a-10a8-4b35-96e2-e8d1a3a97d9f', '1d04ac8b-9b84-4cbf-94ce-8761eef2d05c',
            'e3d206ef-4ffc-40aa-9fbb-79e93fab56a0'
        ],
                   strucplot=False,
                   nolegend=True,
                   noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='eos.png', remove_text=True)
    def test_plot_eos_wc(self):
        basic_test('1419fe6f-cdff-4fb8-b1c1-d8d2c06e33f4', strucplot=False, nolegend=True, noshow=True)
        return gcf()


# helper functions
def basic_test(node_id, **kwargs):
    from aiida_kkr.tools.plot_kkr import plot_kkr
    from aiida.orm import load_node
    from aiida.common.exceptions import NotExistent
    # import database if not there
    try:
        if type(node_id) == list:
            load_node(node_id[0])
        else:
            load_node(node_id)
    except NotExistent:
        print('Node not yet in database. Import test database')
        import_with_migration('files/export_kkr_eos.tar.gz')
        import_with_migration('files/export_kkr_startpot.tar.gz')
        import_with_migration('files/export_kkr_dos.tar.gz')
        import_with_migration('files/export_kkr_scf.tar.gz')

    # now clear old figure and do plotting
    gcf().clear()
    plot_kkr(node_id, **kwargs)


if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    t = Test_plot_kkr()
    t.test_plot_kkr_calc()
