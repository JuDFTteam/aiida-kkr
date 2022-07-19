#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np


def fermi_plot(ax1, energies, yvalues, label, usefermi, efermi, pingap, scatter, color, linestyle):
    """ Plot single entry, w.r.t. to Fermi energy or not. Possible to pin gap (tricky, though,untested!!!) """
    if pingap:
        egap = energies[np.argmin(yvalues)]
        energies = energies - egap
    if usefermi:
        ax1.plot(energies - efermi, yvalues, label=f'{label}', c=color, linestyle=linestyle)
        if scatter:
            ax1.scatter(energies - efermi, yvalues, c=color)
    else:
        ax1.plot(energies, yvalues, label=f'{label}', c=color, linestyle=linestyle)
        if scatter:
            ax1.scatter(energies, yvalues, c=color)


def plot_DOS_atoms(
    output_dict,
    ytitle,
    usefermi=True,
    pingap=False,
    allatomsDOS=False,
    allatomsDOSlabel='total DOS, all atoms',
    f=None,
    ax1=None,
    scatter=True,
    scale=1.0,
    colorlist=[],
    showlegend=True,
    linestyle='solid',
    manual_fermi=False
):
    """plots DOS relative to Fermi level from Dict node"

    :param output_dict: Dict node to start from
    :param ytitle: plot key in output_dict, e. g. `s`,`p`,`d`,`total_DOS`, `non-spherical`
    :param usefermi:bool, True means all energies are plotted relative to the Fermi level
    :param pingap: bool, True means that the bandgap is pinned and all energies and DOS values are given relative to this
    :param allatomsDOS: bool, True for getting the total DOS of all atoms' added up
    :param f: matplotlib figure object, can be passed to add to existing figure and axis ax1
    :param ax1: matplotlib axis object, can be passed to add to existing figure and axis ax1
    :param scatter: bool, True means that besides the lineplot, also the scatter plot is shown (useful to highlight the finite set of energy points)
    :param  scale: float, multiply y-values with factor scale (helpful for plotting different calcs together)
    :param colorlist: list of colors readable by matplotlib, to specify the plotting colors (follows order of input dict)
    :param showlegend: bool
    :param linestyle: string
    :param manual_fermi: float, specify Fermi energy to be used
    """
    if manual_fermi == False:
        try:
            parent_calc = output_dict.get_incoming().all_nodes()[0].inputs.parent_folder.get_incoming().all_nodes()[0]
            efermi = parent_calc.outputs.output_parameters.get_dict()['fermi_energy_in_ryd'][-1]
        except:
            print('WARNING: Reading Efermi failed!')
            efermi = 0.0
            usefermi = False
            #efermi=output_dict.get_dict()['fermi_energy_in_ryd'][-1]
    else:
        efermi = manual_fermi

    DOS_dict = output_dict.get_dict()['DOS']
    if DOS_dict['atom 1']['spin_directions'] == 2:
        energies = DOS_dict['atom 1']['spin_up']['energy_in_ryd']
    else:
        energies = DOS_dict['atom 1']['energy_in_ryd']
    energies = np.array(energies)

    set_fontsize = 11
    set_labelpad = 8

    #Check if new subplot needs to be created
    if f is None and ax1 is None:
        f, ax1 = plt.subplots(1, 1)
    elif f is None and ax1 is not None:
        print('ERROR: No fig argument passed, but axis')
    elif ax1 is None and f is not None:
        print('ERROR: No axis argument passed, but fig')

    if pingap:
        usefermi = False
        energylabel = r'$\Delta E_{\text{min}}$'

    if usefermi:
        energylabel = 'E - E_F'
    else:
        energylabel = 'E'
    ax1.set_xlabel(f'{energylabel} (Ryd)', fontsize=set_fontsize, labelpad=set_labelpad)
    # eV + ryd axis
    ryd = 13.605684958731  # 1 Ryd in eV
    ryd2eV = lambda x: x * ryd
    eV2ryd = lambda x: x / ryd
    ax1.tick_params(axis='x', labelsize=set_fontsize)
    asecax = ax1.secondary_xaxis('top', functions=(ryd2eV, eV2ryd))
    asecax.tick_params(axis='x', labelsize=set_fontsize)
    asecax.set_xlabel(f'{energylabel} (eV)', fontsize=set_fontsize, labelpad=set_labelpad)

    totalyvalues = np.zeros(np.shape(energies))

    counter = 0
    for key in DOS_dict:
        if len(colorlist) == 0:
            color = None
        else:
            color = colorlist[counter]
        print(key, DOS_dict[key]['element'])
        label_atom = f"{key}, {DOS_dict[key]['element']}"
        #try:
        if DOS_dict[key]['spin_directions'] == 2:
            yvalues = scale * np.abs(np.array(DOS_dict[key]['spin_up'][ytitle])) + np.abs(
                np.array((DOS_dict[key]['spin_down'][ytitle]))
            )
        else:
            yvalues = scale * np.array(DOS_dict[key][ytitle])
        if allatomsDOS:
            totalyvalues += yvalues
        else:
            fermi_plot(ax1, energies, yvalues, label_atom, usefermi, efermi, pingap, scatter, color, linestyle)
        counter += 1

    if allatomsDOS:
        fermi_plot(ax1, energies, totalyvalues, allatomsDOSlabel, usefermi, efermi, pingap, scatter, color, linestyle)
        #except: pass

    ax1.axvline(0, color='black', linewidth=0.5)
    ax1.axhline(0, color='black', linewidth=0.5)
    ax1.set_ylabel(f'DOS:  {ytitle} (e/Ryd)', fontsize=set_fontsize, labelpad=set_labelpad)
    ax1.tick_params(axis='y', labelsize=set_fontsize)
    if showlegend:
        ax1.legend(loc=3, fontsize='x-small')  #set_fontsize)

    return f, ax1
    f.set_size_inches(8, 4)
    f.savefig(f'DOSplot_{ytitle}.png', dpi=100)
