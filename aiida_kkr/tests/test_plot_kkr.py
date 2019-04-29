#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import object
import pytest
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import gcf, title

@pytest.mark.usefixtures("aiida_env")
class Test_plot_kkr(object):
    """
    Tests for the plot_kkr tool
    """

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='kkr.png')
    def test_plot_kkr_calc(self):
        basic_test('d507d133-faec-4b31-857e-b0e6e7e99a18', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='voro.png')
    def test_plot_voro_calc(self):
        basic_test('086eb074-3275-4e80-9c14-811058c641ff', strucplot=False)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='dos.png')
    def test_plot_dos(self):
        basic_test('b4745834-93ff-4ac5-88a4-337a1de57168', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='qdos.png', remove_text=True)
    def test_plot_qdos(self):
        # load necessary files from db_dump files
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrcalc_qdos.tar.gz')
        basic_test('a0d0d29f-7b22-4ca4-ba55-6b97569d94af', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='vorostart.png')
    def test_plot_vorostart_wc(self):
        basic_test('165e9636-07d1-4827-a0e0-6dc2d5b5ec5a', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='scf.png')
    def test_plot_scf_wc(self):
        basic_test('d3f45122-f3f8-4726-a817-be7c98f9447e', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='scf_grouped.png')
    def test_plot_scf_wc_grouped(self):
        basic_test(['23d7f92d-b952-4267-95fa-9ebba354ec85', 'c1e70d2b-5142-49d1-8a19-4a8a3b56d9cc', '99f00f1c-21a0-41ae-aaac-4921d16f9ee2'], strucplot=False, nolegend=True, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='eos.png')
    def test_plot_eos_wc(self):
        basic_test('58a3e5e4-aeba-400d-a1f4-2f4756bd6f9e', strucplot=False, nolegend=True, noshow=True)
        return gcf()


# helper functions
def basic_test(node_id, **kwargs):
    from aiida_kkr.tools.plot_kkr import plot_kkr
    from aiida.orm import load_node
    from aiida.common.exceptions import NotExistent
    # import database if note there
    try:
        if type(node_id)==list:
            load_node(node_id[0])
        else:
            load_node(node_id)
    except NotExistent:
        print('Node not yet in database. Import test database')
        from aiida.orm.importexport import import_data
        import_data('files/export_eos_workflow.tar.gz')
    # now clear old figure and do plotting
    gcf().clear()
    plot_kkr(node_id, **kwargs)

if __name__=='__main__':
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()
    t = Test_plot_kkr()
    t.test_plot_kkr_calc()
