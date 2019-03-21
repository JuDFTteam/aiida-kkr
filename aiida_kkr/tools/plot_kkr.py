# -*- coding: utf-8 -*-
"""
contains plot_kkr class for node visualization
"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.4"
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
    :param l_channels: plot l-channels in addition to total DOS (default: True)
    :type l_channels: bool
    :param logscale: plot rms and charge neutrality curves on a log-scale (default: True)
    :type locscale: bool

    additional keyword arguments are passed onto the plotting function which allows, for example,
    to change the markers used in a DOS plot to crosses via `marker='x'`

    :usage: plot_kkr(nodes, **kwargs)

    where nodes is a node identifier (the node itself, it's pk or uuid) or a list of node identifiers.

    :note:
        If nodes is a list of nodes then the plots are grouped together if possible.
    """

    def __init__(self, nodes=None, **kwargs):
        
        # load database if not done already
        from aiida import load_dbenv, is_dbenv_loaded
        if not is_dbenv_loaded():
            load_dbenv()

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
                legend(fontsize='x-small')
                title('')
                subplot(2,1,2)
                self.plot_kkr_single_node(node, only='neutr', label='pk= {}'.format(node.pk), **kwargs)
                title('')# remove duplicated plot title of lower plot
                legend(fontsize='x-small')
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
        from aiida.orm import DataFactory, WorkCalculation, Calculation
        
        # basic aiida nodes
        if isinstance(node, DataFactory('structure')):
            if return_name_only: return 'struc'
            self.plot_struc(node)
        elif isinstance(node, DataFactory('parameter')):
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
        # calculations
        elif node.process_label == u'KkrCalculation':
            if return_name_only: return 'kkr'
            self.plot_kkr_calc(node, **kwargs)
        elif node.process_label == u'VoronoiCalculation':
            if return_name_only: return 'voro'
            self.plot_voro_calc(node, **kwargs)
        elif node.process_label == u'KkrimpCalculation':
            if return_name_only: return 'kkrimp'
            self.plot_kkrimp_calc(node, **kwargs)
        else:
            raise TypeError("input node neither a `Calculation` nor a `WorkCalculation` (i.e. workflow): {}".format(type(node)))
    
    ### helper functions (structure plot, rms plot, dos plot, data extraction ...) ###

    def get_node(self, node):
        """Get node from pk or uuid"""
        from aiida.orm import load_node
        from aiida.orm.node import Node
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
        inputs = node.get_inputs_dict()
        outputs = node.get_outputs_dict()
        for key in list(outputs.keys()):
            try:
                int(key.split('_')[-1])
                has_id = True
            except:
                has_id = False
            # remove 'CALL' and 'CREATE' links
            if 'CALL' in key or 'CREATE' in key or has_id:
                outputs.pop(key)
                
        # now print information about node
        print('type:', type(node))
        print('label:', node.label)
        print('description:', node.description)
        try:
            print('process state:', node.process_label)
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
        from aiida.orm import DataFactory
        StructureData = DataFactory('structure')
        if not isinstance(node, StructureData):
            structure, voro_parent = VoronoiCalculation.find_parent_structure(node)
        else:
            structure = node
        # check if empty sphere need to be removed for plotting (ase structgure cannot be constructed for alloys or vacancies)
        if structure.has_vacancies():
            print("structure has vacancies, need to remove empty sites for plotting")
            stmp = StructureData(cell=structure.cell)
            for site in structure.sites:
                k = structure.get_kind(site.kind_name)
                pos = site.position
                if not k.has_vacancies():
                    stmp.append_atom(position=pos, symbols=k.symbol)
                else:
                    print("removing atom", site)
            stmp.set_pbc(structure.pbc)
            structure = stmp
        # now construct ase object and use ase's viewer
        ase_atoms = structure.get_ase()
        print("plotting structure using ase's `view`")
        view(ase_atoms, **kwargs)
    
    def dosplot(self, d, struc, nofig, all_atoms, l_channels, **kwargs):
        """plot dos from xydata node"""
        from numpy import array, sum
        from matplotlib.pyplot import plot, xlabel, ylabel, gca, figure, legend
        import matplotlib as mpl
        from cycler import cycler
        
        # open new figure
        if not nofig: figure()
        
        x_all = d.get_x()
        y_all = d.get_y()
        
        natoms = len(struc.sites)
        nspin = old_div(len(y_all[0][1]), natoms)
        
        xlbl = x_all[0]+' ('+x_all[2]+')'
    
        #tot only:
        if not l_channels:
            lmax = 1
        else:
            lmax = len(y_all)
            
        pcycle_default = mpl.rcParams['axes.prop_cycle']
        # change color cycler to match spin up/down colors
        if nspin==2 and not all_atoms:
            pcycle_values = pcycle_default.by_key()['color']
            pcycle_values = array([[i,i] for i in pcycle_values]).reshape(-1)
            pcycle_default = cycler('color', pcycle_values)
        gca().set_prop_cycle(pcycle_default)
            
        for il in range(lmax):
            y2 = y_all[il]
            # extract label
            ylbl = 'DOS ('+y2[2]+')'
            # take data
            y2 = y2[1].copy()
            y2 = y2.reshape(natoms, nspin, -1)
            x = x_all[1].copy() 
        
            for ispin in range(nspin):
                if not all_atoms:
                    y = [sum(y2[:,ispin,:], axis=0)] # artificial list so that y[iatom] works later on
                    natoms2 = 1
                    yladd = ''
                else:
                    natoms2 = natoms
                    y = y2[:,ispin,:]
    
                for iatom in range(natoms2):
                    yladd = y_all[il][0].replace('dos ', '')
                    if all_atoms:
                        yladd+=', atom='+str(iatom+1)
                    elif ispin>0:
                        yladd=''
                    xplt = x[iatom*nspin+ispin]
                    yplt = y[iatom]
                    
                    plot(xplt, yplt, label=yladd, **kwargs)
                    xlabel(xlbl)
                    ylabel(ylbl)
        legend(fontsize='x-small')
    
    def rmsplot(self, rms, neutr, nofig, ptitle, logscale, only=None, **kwargs):
        """plot rms and charge neutrality"""
        from numpy import array
        from matplotlib.pylab import figure, plot, twinx, xlabel, ylabel, legend, subplots_adjust, title, gca
        
        if not nofig: figure()

            
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
                label='charge neutrality'
            else:
                label=kwargs.pop('label')
            plot(neutr,'-or', label=label)
            ax2 = gca()
            ylabel('charge neutrality', color='r')
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
                ylabel('neutr')
                xlabel('iteration')
                ax1 = gca()
                if logscale:
                    ax1.set_yscale('log')
            else:
                raise ValueError('`only` can only be `rms or `neutr` but got {}'.format(only))
        title(ptitle)
        subplots_adjust(right=0.85)
    
    def get_rms_kkrcalc(self, node, title=None):
        """extract rms etc from kkr calculation. Works for both finished and still running calculations."""
        from plumpy import ProcessState
        from masci_tools.io.common_functions import search_string
        
        rms, neutr, etot, efermi = [], [], [], []
        ptitle = ''
            
        if node.process_state == ProcessState.FINISHED:
            o = node.out.output_parameters.get_dict()
            neutr = o['convergence_group'][u'charge_neutrality_all_iterations']
            efermi = o['convergence_group'][u'fermi_energy_all_iterations']
            etot = o['convergence_group'][u'total_energy_Ry_all_iterations']
            rms = o['convergence_group'][u'rms_all_iterations']
            ptitle = 'Time per iteration: ' + str(o['timings_group'].get('Time in Iteration')) + ' s'
            """
            pk = o.get('last_calc_nodeinfo').get('pk')
            c = load_node(pk)
            m = c.out.output_parameters.get_dict().get('convergence_group').get('total_spin_moment_all_iterations')
    
            figure(); subplot(1,2,1); plot(rms_all); twinx(); plot(neutr_all,'r'); subplot(1,2,2); plot(m)
            """
        elif node.get_attr('state') in [u'WITHSCHEDULER', u'RETRIEVING']:
            # extract info needed to open transport
            c = node.get_code()
            comp = c.get_computer()
            authinfo = comp.get_authinfo(c.get_user())
            transport = authinfo.get_transport()

            out_kkr = ''
    
            # now get contents of out_kkr using remote call of 'cat'
            with transport as open_transport:
                path = node.get_attr('remote_workdir')
                if 'out_kkr' in transport.listdir(path):
                    out_kkr = transport.exec_command_wait('cat '+path+'/out_kkr')
                    out_kkr = out_kkr[1].split('\n')
                    
            # now extract rms, charge neutrality, total energy and value of Fermi energy
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
            print('no rms extracted')

        return rms, neutr, etot, efermi, ptitle
    
    
    
    ### calculations ###
    
    def plot_kkr_calc(self, node, **kwargs):
        """plot things for a kkr calculation node"""
        
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
            if 'output_parameters' in node.get_outputs_dict():
                results_dict = node.get_outputs_dict().get('output_parameters').get_dict()
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
            self.plot_struc(node)
            
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

        # try to plot dos and qdos data if calculation was bandstructure or DOS run
        from os import listdir
        from numpy import loadtxt, array, where
        from masci_tools.vis.kkr_plot_bandstruc_qdos import dispersionplot
        from masci_tools.vis.kkr_plot_dos import dosplot
        from matplotlib.pyplot import show, figure, title, xticks, xlabel, axvline

        retpath = node.out.retrieved.get_abs_path('')
        has_dos = 'dos.atom1' in listdir(retpath)
        has_qvec = 'qvec.dat' in listdir(retpath)

        # remove already automatically set things from kwargs
        if 'ptitle' in list(kwargs.keys()): 
            ptitle = kwargs.pop('ptitle')
        else:
            ptitle = 'pk= {}'.format(node.pk)
        if 'newfig' in list(kwargs.keys()): kwargs.pop('newfig')

        # qdos
        if has_qvec:
            has_qdos = 'qdos.01.1.dat' in listdir(retpath)
            if has_qdos:
                ne = len(set(loadtxt(retpath+'/qdos.01.1.dat')[:,0]))
                if ne>1:
                    dispersionplot(retpath, newfig=True, ptitle=ptitle, **kwargs)
                    # add plot labels
                    ilbl = node.inp.kpoints.get_attr('label_numbers')
                    slbl = node.inp.kpoints.get_attr('labels')
                    ilbl = array(ilbl)
                    slbl = array(slbl)
                    m_overlap = where(abs(ilbl[1:]-ilbl[:-1])== 1)
                    if len(m_overlap[0])>0:
                        for i in m_overlap[0]:
                            slbl[i+1] = '\n'+slbl[i+1]
                    xticks(ilbl, slbl)
                    xlabel('')
                    [axvline(i, color='grey', ls=':') for i in ilbl]
    
        # dos
        if has_dos:
            figure()
            dosplot(retpath, **kwargs)
            title(ptitle)


    def plot_voro_calc(self, node, **kwargs):
        """plot things for a voro calculation node"""
        
        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
            
        # plot structure
        if strucplot:
            self.plot_struc(node)
            
        outdict = node.out.output_parameters.get_dict()
        # TODO maybe plot some output of voronoi
    
    
    def plot_kkrimp_calc(self, node, **kwargs):
        """plot things from a kkrimp calculation node"""
        print('not implemented yet')
        pass
    
    
    ### workflows ###
    
    
    def plot_kkr_dos(self, node, **kwargs):
        """plot outputs of a kkr_dos_wc workflow"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        
        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels = True, False, True
        if 'interpol' in list(kwargs.keys()): interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()): l_channels = kwargs.pop('l_channels') 
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')
        
        if interpol:
            d = node.out.dos_data_interpol
        else:
            d = node.out.dos_data
            
        # extract structure (neede in dosplot to extract number of atoms and number of spins)
        struc, voro_parent = VoronoiCalculation.find_parent_structure(node.inp.remote_data)
            
        # do dos plot after data was extracted
        self.dosplot(d, struc, nofig, all_atoms, l_channels, **kwargs)
    
    
    def plot_kkr_startpot(self, node, **kwargs):
        """plot output of kkr_startpot_wc workflow"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from matplotlib.pyplot import axvline, legend, title
        from masci_tools.io.common_functions import get_Ry2eV
        
        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        
        silent = False
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')
            
        # plot structure
        if strucplot:
            self.plot_struc(node)
            
        # extract structure (neede in dosplot to extract number of atoms and number of spins)
        struc, voro_parent = VoronoiCalculation.find_parent_structure(node)
            
        if not silent:
            # print results
            print('results:')
            self.plot_kkr_single_node(node.out.results_vorostart_wc, noshow=True, silent=True)
        
        # plot starting DOS

        # follow links until DOS data has been found
        d = None
        d_int = None
        for key, val in node.get_outputs_dict().items():
            if key=='last_doscal_dosdata':
                d = val
            elif key=='last_doscal_dosdata_interpol':
                d_int = val
            elif 'CALL' in key:
                if val.process_label=='kkr_dos_wc':
                    for k2, v2 in val.get_outputs_dict().items():
                        if k2=='dos_data':
                            d = v2
                        elif k2=='dos_data_interpol':
                            d_int = v2
            
        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels = True, False, True
        if 'interpol' in list(kwargs.keys()): interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()): l_channels = kwargs.pop('l_channels') 
        nofig = False
        if 'nofig' in list(kwargs.keys()): nofig = kwargs.pop('nofig')
        
        if interpol:
            d = d_int
            
        if d is not None:
            # do dos plot after data was extracted
            self.dosplot(d, struc, nofig, all_atoms, l_channels, **kwargs)
        
        # now add lines for emin, core states, EF
    
        # extract data for dos and energy contour plotting
        if 'last_voronoi_results' in list(node.get_outputs_dict().keys()):
            params_dict = node.out.last_voronoi_results.get_dict()
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
            if emin is not None: axvline(emin, color='r', ls='--', label='emin')
            if ef_Ry is not None and len(ecore_max)>0:
                axvline((ecore_max[0]-ef_Ry)*get_Ry2eV(), color='b', ls='--', label='ecore_max')
                if len(ecore_max)>1:
                    [axvline((i-ef_Ry)*get_Ry2eV(), color='b', ls='--') for i in ecore_max[1:]]
            if emin is not None: legend(loc=3, fontsize='x-small')
            
            title(struc.get_formula())
    
    
    def plot_kkr_scf(self, node, **kwargs):
        """plot outputs of a kkr_scf_wc workflow"""
        from aiida_kkr.calculations.kkr import KkrCalculation
        from numpy import sort
        from matplotlib.pyplot import axvline, axhline
        
        # structure plot only if structure is in inputs
        try:
            struc = node.inp.structure
            strucplot = True
            ptitle = struc.get_formula()
        except:
            strucplot = False
            ptitle = ''
            
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
        # plot structure
        if strucplot:
            self.plot_struc(struc)
            
        # next extract information from outputs 
        niter_calcs = [0]
        try:
            out_dict = node.out.output_kkr_scf_wc_ParameterResults.get_dict()
            neutr = out_dict['charge_neutrality_all_steps']
            rms = out_dict['convergence_values_all_steps']
            rms_goal = node.inp.wf_parameters.get_dict().get('convergence_criterion')
        except: 
            rms_goal = None
            # deal with unfinished workflow
            rms, neutr, etot, efermi = [], [], [], []
            outdict = node.get_outputs_dict()
            for key in sort(list(outdict.keys())):
                if isinstance(outdict[key], KkrCalculation):
                    kkrcalc = outdict[key]
                    rms_tmp, neutr_tmp, etot_tmp, efermi_tmp, ptitle_tmp = self.get_rms_kkrcalc(kkrcalc)
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
        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        else:
            label = None
    
        # extract rms from calculations and plot
        if len(rms)>0:
            self.rmsplot(rms, neutr, nofig, ptitle, logscale, only, label=label)
            if only == 'rms' and rms_goal is not None: axhline(rms_goal, color='grey', ls='--')
            tmpsum = 0
            if not nofig and len(niter_calcs)>1:
                for i in niter_calcs:
                    tmpsum+=i
                    axvline(tmpsum, color='k', ls=':')
    
    
    def plot_kkr_eos(self, node, **kwargs):
        """plot outputs of a kkr_eos workflow"""
        from numpy import sort, array, where
        from matplotlib.pyplot import figure, subplot, title, xlabel, legend, axvline, plot, annotate
        from aiida.orm import load_node
        from aiida_kkr.workflows.voro_start import kkr_startpot_wc
        from ase.eos import EquationOfState
        
        strucplot = True
        if 'strucplot' in list(kwargs.keys()): strucplot = kwargs.pop('strucplot')
            
        # plot structure
        if strucplot:
            self.plot_struc(node)

        # remove unused things from kwargs
        if 'label' in list(kwargs.keys()): label=kwargs.pop('label')
        if 'noshow' in list(kwargs.keys()): kwargs.pop('noshow')
        if 'only' in list(kwargs.keys()): kwargs.pop('only')
        if 'nofig' in list(kwargs.keys()): kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()): kwargs.pop('strucplot')
        silent = False
        if 'silent' in list(kwargs.keys()): silent = kwargs.pop('silent')
    
        # plot convergence behavior
        outdict = node.get_outputs_dict()
        try:
            results = node.out.eos_results
        except:
            silent = True

        if not silent:
            print('results:')
            self.plot_kkr_single_node(results)

        fig_open = False
        plotted_kkr_scf = False
        plotted_kkr_start = False
        for key in sort(list(outdict.keys())):
            tmpnode = outdict[key]
            try:
                tmplabel = tmpnode.process_label
            except:
                tmplabel = None
            if tmplabel == u'kkr_startpot_wc':
                self.plot_kkr_startpot(tmpnode, strucplot=False, silent=True, **kwargs)
                plotted_kkr_start = True
            elif tmplabel == u'kkr_scf_wc':
                if not fig_open:
                    figure()
                    fig_open = True
                # plot rms
                subplot(2,1,1)
                self.plot_kkr_scf(tmpnode, silent=True, strucplot=False, nofig=True, only='rms', noshow=True, label='pk={}'.format(tmpnode.pk), **kwargs) # scf workflow, rms only
                xlabel('') # remove overlapping x label in upper plot
                legend(loc=3, fontsize='x-small', ncol=2)
                # plot charge neutrality
                subplot(2,1,2)
                self.plot_kkr_scf(tmpnode, silent=True, strucplot=False, nofig=True, only='neutr', noshow=True, label='pk={}'.format(tmpnode.pk), **kwargs) # scf workflow, rms only
                title('')# remove duplicated plot title of lower plot
                legend(loc=3, fontsize='x-small', ncol=2)
                plotted_kkr_scf = True

        if not (plotted_kkr_scf or plotted_kkr_start):
            print('found no startpot or kkrstart data to plot')

        # plot eos results
        try:
            # exctract data
            e = array(node.out.eos_results.get_dict().get('energies', []))
            v = array(node.out.eos_results.get_dict().get('volumes', []))
            fitfunc_gs = node.out.eos_results.get_dict().get('gs_fitfunction')
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
            scalings_all = array(node.out.eos_results.get_dict().get('scale_factors_all'))
            scalings = node.out.eos_results.get_dict().get('scalings')
            names = sort([name for name in list(node.out.eos_results.get_dict().get('sub_workflow_uuids').keys()) if 'kkr_scf' in name])
            pks = array([load_node(node.out.eos_results.get_dict().get('sub_workflow_uuids')[name]).pk for name in names])
            mask = []
            for i in range(len(pks)):
                s = scalings[i]
                m = where(scalings_all==s)
                pk = pks[m][0]
                ie = e[m][0]
                iv = v[m][0]
                annotate(s='pk={}'.format(pk), xy=(iv,ie))

            # investigate fit quality by fitting without first/last datapoint 
            if len(e)>4:
               
                eos = EquationOfState(v[1:-1], e[1:-1], eos=fitfunc_gs)
                v01, e01, B1 = eos.fit()
               
                # take out smallest data point
                eos = EquationOfState(v[1:], e[1:], eos=fitfunc_gs)
                v01, e01, B1 = eos.fit()
           
                print('# relative differences to full fit: V0, E0, B (without smallest volume)')
                print('{} {} {}'.format(abs(1-old_div(v01,v0)), abs(1-old_div(e01,e0)), abs(1-old_div(B1,B))))
               
                if len(e)>5:
                    # also take out largest volume
                    eos = EquationOfState(v[1:-1], e[1:-1], eos=fitfunc_gs)
                    v02, e02, B2 = eos.fit()
           
                    print('\n# V0, E0, B (without smallest and largest volume)')
                    print('{} {} {}'.format(abs(1-old_div(v02,v0)), abs(1-old_div(e02,e0)), abs(1-old_div(B2,B))))
        except:
            pass # do nothing if no eos data there
       
