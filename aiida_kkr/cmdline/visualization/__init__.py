# -*- coding: utf-8 -*-
"""
Module with CLI commands for various visualizations of data types.
"""
import click
from aiida.cmdline.utils import decorators
from aiida.cmdline.params import arguments

@click.command('plot')
@arguments.NODES('nodes')
@click.option('-f', 'filename', type=click.File('r'), default=None)
@click.option('--show/--no-show', default=False, show_default=True, help='')
@click.option('--only/--no-only', default=False, show_default=True, help='')
@click.option('--silent/--no-silent', default=False, show_default=True, help='print information about input node including inputs and outputs.')
@click.option('--strucplot/--no-strucplot', default=True, show_default=True, help='plot structure using ase’s view function')
@click.option('--interpol/--no-interpol', default=True, show_default=True, help='use interpolated data for DOS plots')
@click.option('--all_atoms/--no-all_atoms', default=False, show_default=True, help='plot all atoms in DOS plots (default: False, i.e. plot total DOS only)')
@click.option('--l_channels/--no-l_channels', default=True, show_default=True, help='plot l-channels in addition to total DOS')
@click.option('--sum_spins/--no-sum_spins', default=False, show_default=True, help='sum up both spin channels or plot both?')
@click.option('--logscale/--no-logscale', default=True, show_default=True, help='plot rms and charge neutrality curves on a log-scale')
@click.option('--switch_xy/--no-switch_xy', default=False, show_default=True, help='')
@click.option('--iatom', multiple=True, default=[], show_default=True, help='list of atom indices which are supposed to be plotted')
def cmd_plot(nodes, filename, noshow, only, silent, strucplot, interpol, all_atoms, l_channels, sum_spins, logscale, switch_xy, iatom):
    """
    Invoke the plot_kkr command on given nodes and kwargs

    Parsing additional keyword arguments onto the plotting function is not supported, for example,
    to change the markers used in a DOS plot to crosses via `marker='x'`
    """

    # TODO: make iatom input different from multiple
    from aiida_kkr.tools.plot_kkr import plot_kkr

    nodes = list(nodes)
    nodesf = []
    if filename is not None:
        nodesf = filename.read().split()
        filename.close()
    nodes = nodes + nodesf

    if iatom is None:
        iatom = []
    kwargs = {'silent': silent, 
              'strucplot': strucplot, 
              'interpol' : interpol,
              'all_atoms': all_atoms,
              'l_channels': l_channels,
              'sum_spins': sum_spins, 
              'logscale': logscale, 
              'switch_xy': switch_xy, 
              'iatom': iatom}

    plot_kkr(nodes, noshow=noshow, only=only, **kwargs)
