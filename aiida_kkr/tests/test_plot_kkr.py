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
        basic_test('65d578a5-7227-4413-a606-472ae2c597f6', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='voro.png')
    def test_plot_voro_calc(self):
        basic_test('40d6b054-e522-423c-9ae0-2ff765ca2a51', strucplot=False)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='dos.png')
    def test_plot_dos(self):
        basic_test('daeb58d4-c2de-4309-8ee7-7c298fa591a8', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='qdos.png', remove_text=True)
    def test_plot_qdos(self):
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrcalc_qdos.tar.gz')
        basic_test('a0d0d29f-7b22-4ca4-ba55-6b97569d94af', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='vorostart.png')
    def test_plot_vorostart_wc(self):
        basic_test('2a523b90-45e6-47d6-a263-c1c3c88d88f7', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='scf.png')
    def test_plot_scf_wc(self):
        basic_test('bbd9d4d4-30e4-43e5-bfd6-ab9822a1df9b', strucplot=False, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='scf_grouped.png')
    def test_plot_scf_wc_grouped(self):
        basic_test(['b84742fb-6dcc-4667-b995-3faf0c841901', 'bbd9d4d4-30e4-43e5-bfd6-ab9822a1df9b', '3166f9b7-bfcd-4ca4-bd8b-19d31cf15f5c'], strucplot=False, nolegend=True, noshow=True)
        return gcf()

    @pytest.mark.mpl_image_compare(baseline_dir='files/baseline_images/', filename='eos.png')
    def test_plot_eos_wc(self):
        basic_test('edc9c875-2632-437c-b23c-864134675e22', strucplot=False, nolegend=True, noshow=True)
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
        from aiida.tools.importexport import import_data
        import_data('files/export_eos_workflow.tar.gz')
        import_data('files/db_dump_vorostart.tar.gz')
    # now clear old figure and do plotting
    gcf().clear()
    plot_kkr(node_id, **kwargs)

if __name__=='__main__':
    from aiida import load_profile
    load_profile()
    t = Test_plot_kkr()
    t.test_plot_kkr_calc()
