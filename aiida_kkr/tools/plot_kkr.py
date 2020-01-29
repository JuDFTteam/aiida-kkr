# -*- coding: utf-8 -*-
"""
contains plot_kkr class for node visualization
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object, str
from six.moves import range
__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.5.0"
__contributors__ = ("Philipp Rüßmann")


class plot_kkr(object):
    """
    Class grouping all functionality to plot typical nodes (calculations, workflows, ...) of the aiida-kkr plugin.

    :param nodes: node identifier which is to be visualized

    optional arguments:

    :param silent: print information about input node including inputs and outputs (default: False)
    :type silent: bool
    :param strucplot: plot structure using ase’s view function (default: True)
    :type strucplot: bool
    :param interpol: use interpolated data for DOS plots (default: True)
    :type interpol: bool
    :param all_atoms: plot all atoms in DOS plots (default: False, i.e. plot total DOS only)
    :type all_atoms: bool
    :param l_channels: plot l-channels in addition to total DOS (default: True, i.e. plot all l-channels)
    :type l_channels: bool
    :param sum_spins: sum up both spin channels or plot both? (default: False, i.e. plot both spin channels)
    :type sum_spins: bool
    :param logscale: plot rms and charge neutrality curves on a log-scale (default: True)
    :type locscale: bool
    :param switch_xy:  (default: False)
    :type switch_xy: bool
    :param iatom: list of atom indices which are supposed to be plotted (default: [], i.e. show all atoms)
    :type iatom: list

    additional keyword arguments are passed onto the plotting function which allows, for example,
    to change the markers used in a DOS plot to crosses via `marker='x'`

    :usage: plot_kkr(nodes, **kwargs)

    where nodes is a node identifier (the node itself, it's pk or uuid) or a list of node identifiers.

    :note:
        If nodes is a list of nodes then the plots are grouped together if possible.
    """

    def __init__(self, nodes=None, **kwargs):

        # load database if not done already
        from aiida import load_profile
        load_profile()

        if type(nodes)==list:
            from matplotlib.pyplot import show

            node_groups = self.group_nodes(nodes)
            for groupname in list(node_groups.keys()):
                print('\n==================================================================')
                print('Group of nodes: {}\n'.format(groupname))
                # some settings for groups
                if 'noshow' in list(kwargs.keys()): noshow = kwargs.pop('noshow') # this is now removed from kwargs
                noshow = True # always overwrite noshow settings
                if 'only' in list(kwargs.keys()): only = kwargs.pop('only') # this is now removed from kwargs

                # now plot groups one after the other
                self.plot_group(groupname, node_groups, noshow=noshow, nofig=True, **kwargs)

            # finally show all plots
            show()
        elif nodes is not None:
            self.plot_kkr_single_node(nodes, **kwargs)

    ### main wrapper functions ###
    def group_nodes(self, nodes):
        """Go through list of nodes and group them together."""
        groups_dict = {}

        for node in nodes:
            node = self.get_node(node)
            nodeclass = self.classify_and_plot_node(node, return_name_only=True)
            if nodeclass not in list(groups_dict.keys()):
                groups_dict[nodeclass] = []
            groups_dict[nodeclass].append(node)

        return groups_dict

    def plot_group(self, groupname, nodesgroups, **kwargs):
        """Visualize all nodes of one group."""
        from matplotlib.pyplot import figure, subplot, title, xlabel, legend, title
        nodeslist = nodesgroups[groupname]
        # take out label from kwargs since it is overwritten
        if 'label' in list(kwargs.keys()): label = kwargs.pop('label')
        nolegend = False
        if 'nolegend' in list(kwargs.keys()): nolegend = kwargs.pop('nolegend')
        # open a single new figure for each plot here
        if groupname in ['kkr', 'scf']: figure()
        for node in nodeslist:
            print(groupname)
            # open new figure for each plot in these groups
            if groupname in ['eos', 'dos', 'startpot', 'voro']: figure()
            if groupname in ['kkr', 'scf']:
                subplot(2,1,1)
                self.plot_kkr_single_node(node, only='rms', label='pk= {}'.format(node.pk), **kwargs)
                xlabel('') # remove overlapping x label in upper plot
                if not nolegend: legend(fontsize='x-small')
                title('')
                subplot(2,1,2)
                self.plot_kkr_single_node(node, only='neutr', label='pk= {}'.format(node.pk), **kwargs)
                title('')# remove duplicated plot title of lower plot
                if not nolegend: legend(fontsize='x-small')
            else:
                self.plot_kkr_single_node(node, **kwargs)
            print('\n------------------------------------------------------------------\n')

    def plot_kkr_single_node(self, node, **kwargs):
        """ TODO docstring"""

        # determine if some output is printed to stdout
        silent = False
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent') # this is now removed from kwargs
        noshow = False
        if 'noshow' in list(kwargs.keys()):
            noshow = kwargs.pop('noshow') # this is now removed from kwargs

        node = self.get_node(node)

        # print input and output nodes
        if not silent: self.print_clean_inouts(node)

        # classify node and call plotting function
        self.classify_and_plot_node(node, silent=silent, **kwargs)

        if not noshow:
            from matplotlib.pyplot import show
            show()

    def classify_and_plot_node(self, node, return_name_only=False, **kwargs):
        """Find class of the node and call plotting function."""
        # import things
        from pprint import pprint
        from aiida.plugins import DataFactory
        from aiida.orm import CalculationNode

        # basic aiida nodes
        if isinstance(node, DataFactory('structure')):
            if return_name_only: return 'struc'
            self.plot_struc(node, **kwargs)
        elif isinstance(node, DataFactory('dict')):
            if return_name_only: return 'para'
            print('node dict:')
            pprint(node.get_dict())
        elif isinstance(node, DataFactory('remote')):
            if return_name_only: return 'remote'
            print('computer name:', node.get_computer_name())
            print('remote path:', node.get_remote_path())
        elif isinstance(node, DataFactory('folder')):
            if return_name_only: return 'folder'
            print('abs path:')
            pprint(node.get_abs_path())
            print('folder content:')
            pprint(node.get_folder_list())
        # workflows
        elif node.process_label == u'kkr_dos_wc':
            if return_name_only: return 'dos'
            self.plot_kkr_dos(node, **kwargs)
        elif node.process_label == u'kkr_startpot_wc':
            if return_name_only: return 'startpot'
            self.plot_kkr_startpot(node, **kwargs)
        elif node.process_label == u'kkr_scf_wc':
            if return_name_only: return 'scf'
            self.plot_kkr_scf(node, **kwargs)
        elif node.process_label == u'kkr_eos_wc':
            if return_name_only: return 'eos'
            self.plot_kkr_eos(node, **kwargs)
        elif node.process_label == u'kkr_imp_dos_wc':
            if return_name_only: return 'impdos'
            self.plot_kkrimp_dos_wc(node, **kwargs)
        elif node.process_label == u'kkr_imp_wc':
            if return_name_only: return 'imp'
            self.plot_kkrimp_wc(node, **kwargs)
        elif node.process_label == u'kkr_imp_sub_wc':
            if return_name_only: return 'impsub'
            self.plot_kkrimp_sub_wc(node, **kwargs)
        # calculations
        elif node.process_type == u'aiida.calculations:kkr.kkr':
            if return_name_only: return 'kkr'
            self.plot_kkr_calc(node, **kwargs)
        elif node.process_type == u'aiida.calculations:kkr.voro':
            if return_name_only: return 'voro'
            self.plot_voro_calc(node, **kwargs)
        elif node.process_type == u'aiida.calculations:kkr.kkrimp':
            if return_name_only: return 'kkrimp'
            self.plot_kkrimp_calc(node, **kwargs)
        else:
            raise TypeError("input node neither a `Calculation` nor a `WorkChainNode` (i.e. workflow): {} {}".format(type(node), node))

    ### helper functions (structure plot, rms plot, dos plot, data extraction ...) ###

    def get_node(self, node):
        """Get node from pk or uuid"""
        from aiida.orm import load_node, Node
        # load node if pk or uuid is given
        if type(node)==int:
            node = load_node(node)
        elif type(node)==type(''):
            node = load_node(node)
        elif isinstance(node, Node):
            pass
        else:
            raise TypeError("input node should either be the nodes pk (int), it's uuid (str) or the node itself (aiida.orm.Node). Got type(node)={}".format(type(node)))
        return node

    def print_clean_inouts(self, node):
        """print inputs and outputs of nodes without showing 'CALL' and 'CREATE' links in workflows."""
        from pprint import pprint
        # extract inputs and outputs
        inputs = node.get_incoming().all_nodes()
        outputs = node.get_outgoing()
        """
        all_out_labels = outputs.all_link_labels()
        for label in list(all_out_labels):
            try:
                int(label.split('_')[-1])
                has_id = True
            except:
                has_id = False
            # remove 'CALL' and 'CREATE' links
            if 'CALL' in label or 'CREATE' in label or has_id:
                outputs.pop(label)
        """
        outputs = outputs.all_nodes()

        # now print information about node
        print('pk, uuid: {} {}'.format(node.pk, node.uuid))
        print('type:', type(node))
        print('label:', node.label)
        print('description:', node.description)
        try:
            print('process type:', node.process_type)
            print('state:', node.process_state)
        except:
            print('nodes does not have the `process_state` attribute')

        print('\ninputs:')
        pprint(inputs)
        print('\noutputs:')
        pprint(outputs)
        try:
            print('\nexit status: {} ({})'.format(node.exit_status, node.exit_message))
        except:
            pass
        print() # empty line at the end

    def plot_struc(self, node, **kwargs):
        """visualize structure using ase's `view` function"""
        from ase.visualize import view
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        if not isinstance(node, StructureData):
            structure, voro_parent = VoronoiCalculation.find_parent_structure(node)
        else:
            structure = node
        if 'show_empty_atoms' in kwargs:
            show_empty_atoms = kwargs.pop('show_empty_atoms')
        else:
            show_empty_atoms = False
        # check if empty sphere need to be removed for plotting (ase structgure cannot be constructed for alloys or vacancies)
        if structure.has_vacancies:
            print("structure has vacancies, need to remove empty sites for plotting")
            stmp = StructureData(cell=structure.cell)
            for site in structure.sites:
                k = structure.get_kind(site.kind_name)
                pos = site.position
                if not k.has_vacancies:
                    stmp.append_atom(position=pos, symbols=k.symbol)
                elif show_empty_atoms:
                    stmp.append_atom(position=pos, symbols='X')
                else:
                    print("removing atom", site)
            stmp.set_pbc(structure.pbc)
            structure = stmp
        # now construct ase object and use ase's viewer
        ase_atoms = structure.get_ase()
        if 'silent' in kwargs:
            silent = kwargs.pop('silent')
        print("plotting structure using ase's `view` with kwargs={}".format(kwargs))
        view(ase_atoms, **kwargs)

    def dosplot(self, d, natoms, nofig, all_atoms, l_channels, sum_spins, switch_xy, switch_sign_spin2, **kwargs):
        """plot dos from xydata node"""
        from numpy import array, sum, arange
        from matplotlib.pyplot import plot, xlabel, ylabel, gca, figure, legend, fill_between
        import matplotlib as mpl
        from cycler import cycler

        # remove things that will not work for plotting
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')
        if 'filled' in list(kwargs.keys()): filled = kwargs.pop('filled')
        else: filled = False

        # plot only some atoms if 'iatom' is found in input
        show_atoms = []
        if 'iatom' in kwargs:
            show_atoms = kwargs.pop('iatom')
            if type(show_atoms)!=list:
                show_atoms = [show_atoms]
            all_atoms = True # need to set all_atoms to true

        # open new figure
        if not nofig: figure()

        x_all = d.get_x()
        y_all = d.get_y()

        # scale factor for x and/or y
        if 'xscale' in kwargs:
            xscale = kwargs.pop('xscale')
        else:
            xscale = 1.
        if 'yscale' in kwargs:
            yscale = kwargs.pop('yscale')
        else:
            yscale = 1.

        nspin = len(y_all[0][1]) // natoms
        nspin2 = len(y_all[0][1]) // natoms # copy of nspin becaus nspin is reset to 1 if sum_spins is set to True
        if sum_spins: nspin = 1

        xlbl = x_all[0]+' ('+x_all[2]+')'

        #tot only:
        if not l_channels:
            lmax = 1
        else:
            lmax = len(y_all)

        pcycle_default = mpl.rcParams['axes.prop_cycle']
        # change color cycler to match spin up/down colors
        if nspin==2 and (not all_atoms or show_atoms!=[]):
            pcycle_values = pcycle_default.by_key()['color']
            if switch_sign_spin2: # needed for impurity DOS
                n0 = len(show_atoms)
                if n0==0: n0=natoms
                #pcycle_values = array([j for i in range(n0) for j in pcycle_values]).reshape(-1)
                #pcycle_values = list(pcycle_values[:n0]) + list(pcycle_values[:n0])
                pcycle_values = array([[i,i] for i in pcycle_values]).reshape(-1)
            else:
                pcycle_values = array([[i,i] for i in pcycle_values]).reshape(-1)
            pcycle_default = cycler('color', pcycle_values)
        gca().set_prop_cycle(pcycle_default)

        if 'label' in kwargs:
            labels_all = kwargs.pop('label')
            if type(labels_all)!=list:
                labels_all = [labels_all for i in range(natoms)]
        else:
            labels_all = None

        for il in range(lmax):
            y2 = y_all[il]
            # extract label
            ylbl = 'DOS ('+y2[2]+')'
            # take data
            y2 = y2[1].copy()
            y2 = y2.reshape(natoms, nspin2, -1)
            x = x_all[1].copy()

            for ispin in range(nspin):
                if not all_atoms:
                    y = [sum(y2[:,ispin,:], axis=0)] # artificial list so that y[iatom] works later on
                    if sum_spins and nspin2==2:
                        if switch_sign_spin2:
                            y[0] = -y[0] - sum(y2[:,1,:], axis=0)
                        else:
                            y[0] = -y[0] + sum(y2[:,1,:], axis=0)
                    natoms2 = 1
                    yladd = ''
                else:
                    natoms2 = natoms
                    y = y2[:,ispin,:]
                    if sum_spins and nspin2==2:
                        if switch_sign_spin2:
                            y =  y + y2[:,1,:]
                        else:
                            y = -y + y2[:,1,:]

                for iatom in range(natoms2):
                    if iatom in show_atoms or show_atoms==[]:
                        yladd = y_all[il][0].replace('dos ', '')
                        if all_atoms:
                            yladd+=', atom='+str(iatom+1)
                        elif ispin>0:
                            yladd=''
                        if labels_all is not None and ispin==0:
                            yladd = labels_all[iatom]
                        xplt = x[iatom*nspin+ispin] * xscale
                        yplt = y[iatom] * yscale
                        if ispin>0 and switch_sign_spin2:
                            yplt = -yplt
                        yplt = yplt * yscale
                        if not switch_xy:
                            if not filled:
                                plot(xplt, yplt, label=yladd, **kwargs)
                            else:
                                fill_between(xplt, yplt, label=yladd, **kwargs)
                            xlabel(xlbl)
                            ylabel(ylbl)
                        else:
                            if not filled:
                                plot(yplt, xplt, label=yladd, **kwargs)
                            else:
                                fill_between(yplt, xplt, label=yladd, **kwargs)
                            xlabel(ylbl)
                            ylabel(xlbl)

        legend(fontsize='x-small')

    def rmsplot(self, rms, neutr, nofig, ptitle, logscale, only=None, rename_second=None, **kwargs):
        """plot rms and charge neutrality"""
        from numpy import array
        from matplotlib.pylab import figure, plot, twinx, xlabel, ylabel, legend, subplots_adjust, title, gca

        if not nofig: figure()

        # allow to overwrite name for second quantity, plotted on second y axis
        name_second_y = 'charge_neutrality'
        if rename_second is not None:
            name_second_y = rename_second

        if only is None:
            if 'label' not in list(kwargs.keys()):
                label='rms'
            else:
                label=kwargs.pop('label')
            plot(rms, '-xb', label=label)
            ax1 = gca()
            ylabel('rms', color='b')
            xlabel('iteration')
            twinx()
            if logscale: neutr = abs(array(neutr))
            if 'label' not in list(kwargs.keys()):
                label=name_second_y
            else:
                label=kwargs.pop('label')
            plot(neutr,'-or', label=label)
            ax2 = gca()
            ylabel(name_second_y, color='r')
            if logscale:
                ax1.set_yscale('log')
                ax2.set_yscale('log')
        else: # individual plots of rms or neutrality (needed for eos plot)
            if only=='rms':
                plot(rms, **kwargs)
                ylabel('rms')
                xlabel('iteration')
                ax1 = gca()
                if logscale:
                    ax1.set_yscale('log')
            elif only=='neutr':
                if logscale: neutr = abs(array(neutr))
                plot(neutr, **kwargs)
                ylabel(name_second_y)
                xlabel('iteration')
                ax1 = gca()
                if logscale:
                    ax1.set_yscale('log')
            else:
                raise ValueError('`only` can only be `rms or `neutr` but got {}'.format(only))
        title(ptitle)

    def get_rms_kkrcalc(self, node, title=None):
        """extract rms etc from kkr Calculation. Works for both finished and still running Calculations."""
        from aiida.engine import ProcessState
        from aiida.common.folders import SandboxFolder
        from masci_tools.io.common_functions import search_string

        rms, neutr, etot, efermi = [], [], [], []
        ptitle = ''

        if node.process_state == ProcessState.FINISHED:
            if node.is_finished_ok:
                o = node.outputs.output_parameters.get_dict()
                neutr = o['convergence_group'][u'charge_neutrality_all_iterations']
                efermi = o['convergence_group'][u'fermi_energy_all_iterations']
                etot = o['convergence_group'][u'total_energy_Ry_all_iterations']
                rms = o['convergence_group'][u'rms_all_iterations']
                ptitle = 'Time per iteration: ' + str(o['timings_group'].get('Time in Iteration')) + ' s'
        elif node.process_state in [ProcessState.WAITING, ProcessState.FINISHED, ProcessState.RUNNING]:
            # extract info needed to open transport
            c = node.inputs.code
            comp = c.computer
            authinfo = comp.get_authinfo(c.user)
            transport = authinfo.get_transport()

            out_kkr = ''

            # now get contents of out_kkr using remote call of 'cat'
            with SandboxFolder() as tempfolder:
                with tempfolder.open('tempfile', 'w') as f:
                    try:
                        node.outputs.remote_folder.getfile('out_kkr', f.name)
                        has_outfile = True
                    except IOError:
                        has_outfile = False
                if has_outfile:
                    with tempfolder.open('tempfile', 'r') as f:
                        out_kkr = f.readlines()

            # now extract rms, charge neutrality, total energy and value of Fermi energy
            if has_outfile:
                itmp = 0
                while itmp>=0:
                    itmp = search_string('rms', out_kkr)
                    if itmp>=0:
                        tmpline = out_kkr.pop(itmp)
                        tmpval = float(tmpline.split('=')[1].split()[0].replace('D', 'e'))
                        rms.append(tmpval)
                itmp = 0
                while itmp>=0:
                    itmp = search_string('charge neutrality', out_kkr)
                    if itmp>=0:
                        tmpline = out_kkr.pop(itmp)
                        tmpval = float(tmpline.split('=')[1].split()[0].replace('D', 'e'))
                        neutr.append(tmpval)
                itmp = 0
                while itmp>=0:
                    itmp = search_string('TOTAL ENERGY in ryd', out_kkr)
                    if itmp>=0:
                        tmpline = out_kkr.pop(itmp)
                        tmpval = float(tmpline.split(':')[1].split()[0].replace('D', 'e'))
                        etot.append(tmpval)
                itmp = 0
                while itmp>=0:
                    itmp = search_string('E FERMI', out_kkr)
                    if itmp>=0:
                        tmpline = out_kkr.pop(itmp)
                        tmpval = float(tmpline.split('FERMI')[1].split()[0].replace('D', 'e'))
                        efermi.append(tmpval)
        else:
            print('no rms extracted', node.process_state)

        return rms, neutr, etot, efermi, ptitle



    ### Calculations ###

    def plot_kkr_calc(self, node, **kwargs):
        """plot things for a kkr Calculation node"""

        # extract options from kwargs
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        logscale = True
        if 'logscale' in list(kwargs.keys()): logscale = kwargs.pop('logscale')
        only = None
        if 'only' in list(kwargs.keys()): only = kwargs.pop('only')
        silent = False
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')

        #print output
        if not silent:
            from pprint import pprint
            print('results dict (entries with `...` have been removed for this writeout for the sake of shortness):')
            if 'output_parameters' in node.get_outgoing().all_link_labels():
                results_dict = node.get_outgoing().get_node_by_label('output_parameters').get_dict()
                # remove symmetry descriptions from resuts dict before writting output
                if 'symmetries_group' in list(results_dict.keys()): results_dict['symmetries_group']['symmetry_description'] = '...'
                if 'convergence_group' in list(results_dict.keys()):
                    results_dict['convergence_group']['charge_neutrality_all_iterations'] = '...'
                    results_dict['convergence_group']['dos_at_fermi_energy_all_iterations'] = '...'
                    results_dict['convergence_group']['fermi_energy_all_iterations'] = '...'
                    results_dict['convergence_group']['rms_all_iterations'] = '...'
                    results_dict['convergence_group']['total_energy_Ry_all_iterations'] = '...'
                    results_dict['convergence_group']['spin_moment_per_atom_all_iterations'] = '...'
                    results_dict['convergence_group']['orbital_moment_per_atom_all_iterations'] = '...'
                    results_dict['convergence_group']['total_spin_moment_all_iterations'] = '...'
                pprint(results_dict)

        # plot structure
        if strucplot:
            self.plot_struc(node, **kwargs)

        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        else:
            label = None
        try:
            rms, neutr, etot, efermi, ptitle = self.get_rms_kkrcalc(node)
        except KeyError:
            # this happens in case of qdos run
            rms = []

        if len(rms)>1:
            self.rmsplot(rms, neutr, nofig, ptitle, logscale, only, label=label)

        # try to plot dos and qdos data if Calculation was bandstructure or DOS run
        from subprocess import check_output
        from os import listdir
        from numpy import loadtxt, array, where
        from masci_tools.vis.kkr_plot_bandstruc_qdos import dispersionplot
        from masci_tools.vis.kkr_plot_FS_qdos import FSqdos2D
        from masci_tools.vis.kkr_plot_dos import dosplot
        from matplotlib.pyplot import show, figure, title, xticks, xlabel, axvline

        if node.is_finished_ok:
            retlist = node.outputs.retrieved.list_object_names()
            has_dos = 'dos.atom1' in retlist
            has_qvec = 'qvec.dat' in retlist
            has_qdos = False
            
            # remove already automatically set things from kwargs
            if 'ptitle' in list(kwargs.keys()):
                ptitle = kwargs.pop('ptitle')
            else:
                ptitle = 'pk= {}'.format(node.pk)
            if 'newfig' in list(kwargs.keys()): kwargs.pop('newfig')
            
            # qdos
            if has_qvec:
                has_qdos = 'qdos.01.1.dat' in retlist
                if has_qdos:
                    with node.outputs.retrieved.open('qdos.01.1.dat', mode='r') as f:
                        ne = len(set(loadtxt(f)[:,0]))
                        if ne>1 or 'as_e_dimension' in kwargs.keys():
                            try:
                                ef = check_output('grep "Fermi energy" {}'.format(f.name.replace('qdos.01.1.dat', 'output.0.txt')), shell=True) 
                                ef = float(ef.split('=')[2].split()[0])
                            except:
                                # extract Fermi level from parent calculation
                                parent_calc = node.inputs.parent_folder.get_incoming().first().node
                                ef = parent_calc.outputs.output_parameters.get_dict()['fermi_energy']
                            dispersionplot(f, newfig=(not nofig), ptitle=ptitle, logscale=logscale, ef=ef, **kwargs)
                            # add plot labels
                            try:
                                ilbl = node.inputs.kpoints.get_attr('label_numbers')
                                slbl = node.inputs.kpoints.get_attr('labels')
                                ilbl = array(ilbl)
                                slbl = array(slbl)
                                m_overlap = where(abs(ilbl[1:]-ilbl[:-1])== 1)
                                if len(m_overlap[0])>0:
                                    for i in m_overlap[0]:
                                        slbl[i+1] = '\n'+slbl[i+1]
                                xticks(ilbl, slbl)
                                xlabel('')
                                [axvline(i, color='grey', ls=':') for i in ilbl]
                            except:
                                xlabel('id_kpt')
                        else:
                            ef = check_output('grep "Fermi energy" {}'.format(f.name.replace('qdos.01.1.dat', 'output.0.txt')), shell=True) 
                            ef = float(ef.split('=')[2].split()[0])
                            FSqdos2D(f, logscale=logscale, ef=ef, **kwargs)
            
            # dos only if qdos was not plotted already
            if has_dos and not has_qdos:
                with node.outputs.retrieved.open('dos.atom1', mode='r') as f:
                    if not nofig: figure()
                    dosplot(f, **kwargs)
                    title(ptitle)


    def plot_voro_calc(self, node, **kwargs):
        """plot things for a voro Calculation node"""

        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')

        # plot structure
        if strucplot:
            self.plot_struc(node, **kwargs)

        outdict = node.outputs.output_parameters.get_dict()
        # TODO maybe plot some output of voronoi


    def plot_kkrimp_calc(self, node, return_rms=False, return_stot=False, **kwargs):
        """plot things from a kkrimp Calculation node"""
        from numpy import array, ndarray
        from numpy import sqrt, sum

        # read data from output node
        rms_goal, rms = None, []
        if node.is_finished_ok:
            out_para = node.outputs.output_parameters
            out_para_dict = out_para.get_dict()
            out_para_dict['convergence_group']['rms_all_iterations']
            rms = out_para_dict['convergence_group']['rms_all_iterations']
            rms_goal = out_para_dict['convergence_group']['qbound']
            # extract total magnetic moment
            nat = out_para_dict['number_of_atoms_in_unit_cell']
            s = array(out_para_dict['convergence_group']['spin_moment_per_atom_all_iterations'], dtype=float)
            ss = sqrt(sum(s**2, axis=1)).reshape(-1,nat)
            stot = sum(ss, axis=1)

        # now return values
        return_any, return_list = False, []
        if return_rms:
            return_list += [rms, rms_goal]
            return_any = True
        if return_stot:
            return_list += [stot]
            return_any = True
        if return_any:
            return return_list


    def plot_kkrimp_wc(self, node, **kwargs):
        """plot things from a kkrimp_wc workflow"""

        # call imp_sub plotting from here
        from aiida_kkr.workflows import kkr_imp_sub_wc
        sub_wf = [i.node for i in node.get_outgoing(node_class=kkr_imp_sub_wc).all()][0]
        self.plot_kkrimp_sub_wc(sub_wf)


    def plot_kkrimp_sub_wc(self, node, **kwargs):
        """plot things from a kkrimp_sub_wc workflow"""
        from aiida_kkr.calculations import KkrimpCalculation
        from numpy import array
        from matplotlib.pyplot import figure, subplot, axhline, axvline, gca, ylim

        # extract rms from calculations
        impcalcs = [i.node for i in node.get_outgoing(node_class=KkrimpCalculation).all()]
        rms_all, pks_all, stot_all = [], [], []
        rms_goal = None
        for impcalc in impcalcs:
            pks_all.append(impcalc.pk)
            rms_tmp, rms_goal_tmp, stot_tmp = self.plot_kkrimp_calc(impcalc, return_rms=True, return_stot=True)
            rms_all.append(rms_tmp)
            if rms_goal_tmp is not None:
                if rms_goal is not None:
                    rms_goal = min(rms_goal, rms_goal_tmp)
                else:
                    rms_goal = rms_goal_tmp
            stot_all.append(stot_tmp)

        # extract options from kwargs
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        logscale = True
        if 'logscale' in list(kwargs.keys()): logscale = kwargs.pop('logscale')
        if 'subplot' in list(kwargs.keys()):
            subplots = kwargs.pop('subplot')
        else:
            subplots = None
        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        else:
            label = None
        if 'ptitle' in list(kwargs.keys()):
            ptitle = kwargs.pop('ptitle')
        else:
            ptitle = 'pk= {}'.format(node.pk)
        if 'only' in list(kwargs.keys()):
            only = kwargs.pop('only')
        else:
            only = None

        # plotting of convergence properties (rms etc.)
        if len(rms_all)>0:
            # sort rms values and flatten array
            reorder_rms = array(pks_all).argsort()
            rms, niter_calcs, stot = [], [0], []
            for i in array(rms_all)[reorder_rms]:
                rms += list(i)
                niter_calcs.append(len(i)-0.5)
            for i in array(stot_all)[reorder_rms]:
                stot += list(i)
            # now plot
            if len(rms)>0:
                if not nofig:
                    figure()
                if subplots is not None:
                    subplot(subplots[0], subplots[1], subplots[2])
                if rms_goal is not None: axhline(rms_goal, color='grey', ls='--')
                self.rmsplot(rms, stot, nofig=True, ptitle=ptitle, logscale=logscale, only=only, rename_second='sum(spinmom)', label=label)
                # adapt y-limits to take care of showing spin-moment on sensible scale
                if only is None:
                    yl = gca().get_ylim()
                    ylim(yl[0], max(yl[1], 0.1))
                # add lines that indicate different calculations
                tmpsum = 1
                if not nofig and len(niter_calcs)>1:
                    for i in niter_calcs:
                        tmpsum+=i
                        axvline(tmpsum-1, color='k', ls=':')


    def plot_kkrimp_dos_wc(self, node, **kwargs):
        """plot things from a kkrimp_dos workflow node"""

        # try to plot dos and qdos data if Calculation was bandstructure or DOS run
        from os import listdir
        from numpy import loadtxt, array, where
        from masci_tools.vis.kkr_plot_FS_qdos import FSqdos2D
        from masci_tools.vis.kkr_plot_dos import dosplot
        from matplotlib.pyplot import show, figure, title, xticks, xlabel, axvline

        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()): interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()): l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()): sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()): switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')
        if 'switch_sign_spin2' in list(kwargs.keys()): switch_sign_spin2 = kwargs.pop('switch_sign_spin2')
        else: switch_sign_spin2 = True
        if 'yscale' in list(kwargs.keys()): yscale = kwargs.pop('yscale')
        else: yscale = -1

        has_dos = False
        if interpol and 'dos_data_interpol' in node.outputs:
            d = node.outputs.dos_data_interpol
            has_dos = True
        elif 'dos_data' in node.outputs:
            d = node.outputs.dos_data
            has_dos = True

        if has_dos:
            calcnode = [i for i in node.called_descendants if i.process_label=='KkrimpCalculation'][0]
            if calcnode.is_finished_ok:
                natoms = len(calcnode.outputs.output_parameters.get_dict().get('charge_core_states_per_atom'))
                self.dosplot(d, natoms, nofig, all_atoms, l_channels, sum_spins, switch_xy, switch_sign_spin2, yscale=yscale, **kwargs)
                title('pk= {}'.format(node.pk))


    ### workflows ###


    def plot_kkr_dos(self, node, **kwargs):
        """plot outputs of a kkr_dos_wc workflow"""
        from aiida_kkr.calculations.voro import VoronoiCalculation

        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()): interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()): l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()): sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()): switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')

        if node.is_finished_ok:
            if interpol:
                d = node.outputs.dos_data_interpol
            else:
                d = node.outputs.dos_data

            # extract structure (neede in dosplot to extract number of atoms and number of spins)
            struc, voro_parent = VoronoiCalculation.find_parent_structure(node.inputs.remote_data)
            
            # do dos plot after data was extracted
            self.dosplot(d, len(struc.sites), nofig, all_atoms, l_channels, sum_spins, switch_xy, False, **kwargs)


    def plot_kkr_startpot(self, node, **kwargs):
        """plot output of kkr_startpot_wc workflow"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from aiida.common import exceptions
        from matplotlib.pyplot import axvline, legend, title
        from masci_tools.io.common_functions import get_Ry2eV

        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')

        silent = False
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')

        # plot structure
        if strucplot:
            self.plot_struc(node, **kwargs)

        # extract structure (neede in dosplot to extract number of atoms and number of spins)
        struc, voro_parent = VoronoiCalculation.find_parent_structure(node)

        if not silent:
            # print results
            print('results:')
            try:
                res_node = node.outputs.results_vorostart_wc
            except exceptions.NotExistent:
                res_node = None
            if res_node is not None:
                self.plot_kkr_single_node(res_node, noshow=True, silent=True)

        # plot starting DOS

        # follow links until DOS data has been found
        d = None
        d_int = None
        for link_triple in node.get_outgoing().all():
            if link_triple.link_label=='last_doscal_dosdata':
                d = link_triple.node
            elif link_triple.link_label=='last_doscal_dosdata_interpol':
                d_int = link_triple.node
            elif 'CALL_WORK' in link_triple.link_label:
                if link_triple.link_label=='kkr_dos_wc':
                    for link_triple2 in node.get_outgoing().all():
                        if link_triple2.link_label=='dos_data':
                            d = link_triple2.node
                        elif link_triple2.link_label=='dos_data_interpol':
                            d_int = link_triple2.node

        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()): interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()): l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()): sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()): switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')

        if interpol:
            d = d_int

        if d is not None:
            # do dos plot after data was extracted
            self.dosplot(d, len(struc.sites), nofig, all_atoms, l_channels, sum_spins, switch_xy, False, **kwargs)

        # now add lines for emin, core states, EF

        # extract data for dos and energy contour plotting
        if 'last_voronoi_results' in node.get_outgoing().all_link_labels():
            params_dict = node.outputs.last_voronoi_results.get_dict()
        else:
            params_dict = {}
        emin = params_dict.get('emin_minus_efermi', None)
        emin_Ry = params_dict.get('emin', None)
        if emin is not None and params_dict!={}:
            ef_Ry = emin_Ry - params_dict.get('emin_minus_efermi_Ry')
        else:
            ef_Ry = None
        if params_dict!={}: ecore_max = params_dict.get('core_states_group').get('energy_highest_lying_core_state_per_atom', [])

        if d is not None:
            axvline(0, color='k', ls='--', label='EF')
            tit_add = ''
            if emin is not None: axvline(emin, color='r', ls='--', label='emin')
            if ef_Ry is not None and len(ecore_max)>0:
                if abs((ecore_max[0]-ef_Ry)*get_Ry2eV()-emin)<20:
                    axvline((ecore_max[0]-ef_Ry)*get_Ry2eV(), color='b', ls='--', label='ecore_max')
                else:
                    tit_add = '; E_core<=%.2feV'%((ecore_max[0]-ef_Ry)*get_Ry2eV())
                if len(ecore_max)>1:
                    [axvline((i-ef_Ry)*get_Ry2eV(), color='b', ls='--') for i in ecore_max[1:] if abs((i-ef_Ry)*get_Ry2eV()-emin)<20]
            if emin is not None: legend(loc=3, fontsize='x-small')

            title(struc.get_formula()+', starting potential'+tit_add)


    def plot_kkr_scf(self, node, **kwargs):
        """plot outputs of a kkr_scf_wc workflow"""
        from aiida.orm import CalcJobNode, load_node
        from aiida_kkr.calculations.kkr import KkrCalculation
        from numpy import sort
        from matplotlib.pyplot import axvline, axhline, subplot, figure

        # structure plot only if structure is in inputs
        try:
            struc = node.inputs.structure
            strucplot = True
            ptitle = struc.get_formula()
        except:
            strucplot = False
            ptitle = ''

        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        # plot structure
        if strucplot:
            self.plot_struc(struc, **kwargs)

        # next extract information from outputs
        niter_calcs = [0]
        try:
            out_dict = node.outputs.output_kkr_scf_wc_ParameterResults.get_dict()
            neutr = out_dict['charge_neutrality_all_steps']
            rms = out_dict['convergence_values_all_steps']
            rms_goal = node.inputs.wf_parameters.get_dict().get('convergence_criterion')
        except:
            rms_goal = None
            # deal with unfinished workflow
            rms, neutr, etot, efermi = [], [], [], []
            outdict = node.get_outgoing(node_class=CalcJobNode)
            pks_calcs = sort([i.node.pk for i in outdict.all()])
            for pk in pks_calcs:
                node = load_node(pk)
                if node.process_label==u'KkrCalculation':
                    kkrcalc = node
                    rms_tmp, neutr_tmp, etot_tmp, efermi_tmp, ptitle_tmp = self.get_rms_kkrcalc(kkrcalc)
                    if len(rms_tmp)>0:
                        niter_calcs.append(len(rms_tmp))
                        rms += rms_tmp
                        neutr += neutr_tmp

        # extract options from kwargs
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        logscale = True
        if 'logscale' in list(kwargs.keys()): logscale = kwargs.pop('logscale')
        only = None
        if 'only' in list(kwargs.keys()): only = kwargs.pop('only')
        if 'subplot' in list(kwargs.keys()):
            subplots = kwargs.pop('subplot')
        else:
            subplots = None
        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        else:
            label = None
        if 'dos_only' in list(kwargs.keys()):
            dos_only = kwargs.pop('dos_only')
        else:
            dos_only = False

        # extract rms from calculations and plot
        if len(rms)>0 and not dos_only:
            if not nofig:
                figure()
            if subplots is not None:
                subplot(subplots[0], subplots[1], subplots[2])
            self.rmsplot(rms, neutr, True, ptitle, logscale, only, label=label)
            if only == 'rms' and rms_goal is not None: axhline(rms_goal, color='grey', ls='--')
            tmpsum = 1
            if not nofig and len(niter_calcs)>1:
                for i in niter_calcs:
                    tmpsum+=i
                    axvline(tmpsum-1, color='k', ls=':')
            did_plot = True
        else:
            did_plot = False

        # plot DOS

        # follow links until DOS data has been found
        d = None
        d_int = None
        from aiida_kkr.workflows import kkr_dos_wc
        links_dos = node.get_outgoing(node_class=kkr_dos_wc).all()
        if len(links_dos)>0:
            dosnode = links_dos[0].node
            if 'dos_data' in dosnode.outputs:
                d = dosnode.outputs.dos_data
            if 'dos_data_interpol' in dosnode.outputs:
                d_int = dosnode.outputs.dos_data_interpol
        
        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()): interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()): l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()): sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()): switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')

        if interpol:
            d = d_int

        if d is not None:
            # do dos plot after data was extracted
            self.dosplot(d, len(struc.sites), nofig, all_atoms, l_channels, sum_spins, switch_xy, False, **kwargs)

        return did_plot


    def plot_kkr_eos(self, node, **kwargs):
        """plot outputs of a kkr_eos workflow"""
        from numpy import sort, array, where
        from matplotlib.pyplot import figure, title, xlabel, legend, axvline, plot, annotate
        from aiida.orm import load_node
        from aiida_kkr.workflows.voro_start import kkr_startpot_wc
        from ase.eos import EquationOfState

        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')

        # plot structure
        if strucplot:
            self.plot_struc(node, **kwargs)

        # remove unused things from kwargs
        if 'label' in list(kwargs.keys()): label=kwargs.pop('label')
        if 'noshow' in list(kwargs.keys()): kwargs.pop('noshow')
        if 'only' in list(kwargs.keys()): kwargs.pop('only')
        if 'nofig' in list(kwargs.keys()): kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()): kwargs.pop('strucplot')
        silent = False
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')
        nolegend = False
        if 'nolegend' in list(kwargs.keys()): nolegend = kwargs.pop('nolegend')

        # plot convergence behavior
        try:
            results = node.outputs.eos_results
        except:
            silent = True

        if not silent:
            print('results:')
            self.plot_kkr_single_node(results)

        fig_open = False
        plotted_kkr_scf = False
        plotted_kkr_start = False
        outdict = node.get_outgoing()
        for key in outdict.all():
            if key.link_label!='CALL_CALC':
                tmpnode = key.node
                try:
                    tmplabel = tmpnode.process_label
                except:
                    tmplabel = None
                if tmplabel == u'kkr_startpot_wc':
                    self.plot_kkr_startpot(tmpnode, strucplot=False, silent=True, **kwargs)
                    plotted_kkr_start = True
                elif tmplabel == u'kkr_scf_wc':
                    # plot rms
                    did_plot = self.plot_kkr_scf(tmpnode, silent=True, strucplot=False, nofig=fig_open, only='rms', noshow=True, label='pk={}'.format(tmpnode.pk), subplot=(2,1,1), **kwargs) # scf workflow, rms only
                    if did_plot and not fig_open:
                        fig_open = True
                    if did_plot: xlabel('') # remove overlapping x label in upper plot
                    if did_plot and not nolegend: legend(loc=3, fontsize='x-small', ncol=2)
                    # plot charge neutrality
                    self.plot_kkr_scf(tmpnode, silent=True, strucplot=False, nofig=True, only='neutr', noshow=True, label='pk={}'.format(tmpnode.pk), subplot=(2,1,2), **kwargs) # scf workflow, rms only
                    if did_plot: title('') # remove overlapping title
                    if did_plot and  not nolegend: legend(loc=3, fontsize='x-small', ncol=2)
                    plotted_kkr_scf = True

        if not (plotted_kkr_scf or plotted_kkr_start):
            print('found no startpot or kkrstart data to plot')

        # plot eos results
        try:
            # exctract data
            e = array(node.outputs.eos_results.get_dict().get('energies', []))
            v = array(node.outputs.eos_results.get_dict().get('volumes', []))
            fitfunc_gs = node.outputs.eos_results.get_dict().get('gs_fitfunction')
            # now fit and plot
            figure()
            try:
                eos = EquationOfState(v, e, eos=fitfunc_gs)
                v0, e0, B = eos.fit()
                eos.plot()
                axvline(v0, color='k', ls='--')
            except:
                plot(v, e, 'o')

            # add calculation pks to data points
            scalings_all = array(node.outputs.eos_results.get_dict().get('scale_factors_all'))
            scalings = node.outputs.eos_results.get_dict().get('scalings')
            names = sort([name for name in list(node.outputs.eos_results.get_dict().get('sub_workflow_uuids').keys()) if 'kkr_scf' in name])
            pks = array([load_node(node.outputs.eos_results.get_dict().get('sub_workflow_uuids')[name]).pk for name in names])
            mask = []
            for i in range(len(pks)):
                s = scalings[i]
                m = where(scalings_all==s)
                pk = pks[m][0]
                ie = e[m][0]
                iv = v[m][0]
                if not nolegend: annotate(s='pk={}'.format(pk), xy=(iv,ie))

            # investigate fit quality by fitting without first/last datapoint
            if len(e)>4:

                eos = EquationOfState(v[1:-1], e[1:-1], eos=fitfunc_gs)
                v01, e01, B1 = eos.fit()

                # take out smallest data point
                eos = EquationOfState(v[1:], e[1:], eos=fitfunc_gs)
                v01, e01, B1 = eos.fit()

                print('# relative differences to full fit: V0, E0, B (without smallest volume)')
                print('{} {} {}'.format(abs(1-v01/v0), abs(1-e01/e0), abs(1-B1/B)))

                if len(e)>5:
                    # also take out largest volume
                    eos = EquationOfState(v[1:-1], e[1:-1], eos=fitfunc_gs)
                    v02, e02, B2 = eos.fit()

                    print('\n# V0, E0, B (without smallest and largest volume)')
                    print('{} {} {}'.format(abs(1-v02/v0), abs(1-e02/e0), abs(1-B2/B)))
        except:
            pass # do nothing if no eos data there
