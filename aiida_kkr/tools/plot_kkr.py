

class plot_kkr():
    """ TODO docstring"""

    def __init__(self, nodes, **kwargs):
        if type(nodes)==list:
            node_groups = self.group_nodes(nodes)
            for groupname in node_groups.keys():
                self.plot_group(node_group[groupname], **kwargs)
        else:
            self.plot_kkr_single_node(nodes, **kwargs)

    ### main wrapper functions ###

    def plot_kkr_single_node(self, node, **kwargs):
        """ TODO docstring"""
        
        # load database if not done already
        from aiida import load_dbenv, is_dbenv_loaded
        if not is_dbenv_loaded():
            load_dbenv()
            
        # import things
        from pprint import pprint
        from aiida.orm.node import Node
        from aiida.orm import load_node, DataFactory, WorkCalculation, Calculation
        
        # determine if some output is printed to stdout
        silent = False
        if 'silent' in kwargs.keys():
            silent = kwargs.pop('silent') # this is now removed from kwargs
        noshow = False
        if 'noshow' in kwargs.keys():
            noshow = kwargs.pop('noshow') # this is now removed from kwargs
            
        # load node if pk or uuid is given
        if type(node)==int:
            node = load_node(node)
        elif type(node)==str:
            node = load_node(node)
        elif isinstance(node, Node):
            pass
        else:
            raise TypeError("input node should either be the nodes pk (int), it's uuid (str) or the node itself (aiida.orm.Node)")
    
        # print input and output nodes
        if not silent: self.print_clean_inouts(node)
        
        # classify node and call plotting function
        
        # basic aiida nodes
        if isinstance(node, DataFactory('structure')):
            self.plot_struc(node)
        elif isinstance(node, DataFactory('parameter')):
            print 'node dict:'
            pprint(node.get_dict())
        elif isinstance(node, DataFactory('remote')):
            print 'computer name:', node.get_computer_name()
            print 'remote path:', node.get_remote_path()
        elif isinstance(node, DataFactory('folder')):
            print 'abs path:'
            pprint(node.get_abs_path())
            print 'folder content:'
            pprint(node.get_folder_list())
        # workflows
        elif node.process_label == u'kkr_dos_wc':
            self.plot_kkr_dos(node, **kwargs)
        elif node.process_label == u'kkr_startpot_wc':
            self.plot_kkr_startpot(node, **kwargs)
        elif node.process_label == u'kkr_scf_wc':
            self.plot_kkr_scf(node, **kwargs)
        elif node.process_label == u'kkr_eos_wc':
            self.plot_kkr_eos(node, **kwargs)
        # calculations
        elif node.process_label == u'KkrCalculation':
            self.plot_kkr_calc(node, **kwargs)
        elif node.process_label == u'VoronoiCalculation':
            self.plot_voro_calc(node, **kwargs)
        elif node.process_label == u'KkrimpCalculation':
            self.plot_kkrimp_calc(node, **kwargs)
        else:
            raise TypeError("input node neither a `Calculation` nor a `WorkCalculation` (i.e. workflow): {}".format(type(node)))
    
        if not noshow:
            from matplotlib.pyplot import show
            show()
    
        ### helper functions (structure plot, rms plot, dos plot, data extraction ...) ###
    
    def print_clean_inouts(self, node):
        from pprint import pprint
        # extract inputs and outputs
        inputs = node.get_inputs_dict()
        outputs = node.get_outputs_dict()
        for key in outputs.keys():
            try:
                int(key.split('_')[-1])
                has_id = True
            except:
                has_id = False
            # remove 'CALL' and 'CREATE' links
            if 'CALL' in key or 'CREATE' in key or has_id:
                outputs.pop(key)
                
        # now print information about node
        print 'type:', type(node)
        print 'label:', node.label
        print 'description:', node.description
        try:
            print 'process state:', node.process_label
            print 'state:', node.process_state
        except:
            print 'nodes does not have the `process_state` attribute'
            
        print '\ninputs:'
        pprint(inputs)
        print '\noutputs:'
        pprint(outputs)
        print # empty line at the end
    
    def plot_struc(self, node, **kwargs):
        """visualize structure using ase's `view` function"""
        from ase.visualize import view
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from aiida.orm import DataFactory
        if not isinstance(node, DataFactory('structure')):
            structure, voro_parent = VoronoiCalculation.find_parent_structure(node)
        else:
            structure = node
        ase_atoms = structure.get_ase()
        print "plotting structure using ase's `view`"
        view(ase_atoms, **kwargs)
    
    def dosplot(self, d, struc, nofig, all_atoms, l_channels, **kwargs):
        """TOTO docstring"""
        from numpy import array, sum
        from matplotlib.pyplot import plot, xlabel, ylabel, gca, figure, legend
        import matplotlib as mpl
        from cycler import cycler
        
        # open new figure
        if not nofig: figure()
        
        x_all = d.get_x()
        y_all = d.get_y()
        
        natoms = len(struc.sites)
        nspin = len(y_all[0][1]) / natoms
        
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
        legend()
    
    def rmsplot(self, rms, neutr, nofig, ptitle, logscale, only=None, **kwargs):
        """TODO docstring"""
        from numpy import array
        from matplotlib.pylab import figure, plot, twinx, xlabel, ylabel, legend, subplots_adjust, title, gca
        
        if not nofig: figure()
            
        if only is None:
            plot(rms, '-xb', label='rms')
            ax1 = gca()
            ylabel('rms', color='b')
            xlabel('iteration')
            twinx()
            if logscale: neutr = abs(array(neutr))
            plot(neutr,'-or', label='charge neutrality')
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
    
    def get_rms_kkrcalc(self, node):
        """TODO docstring"""
        from plumpy import ProcessState
        from masci_tools.io.common_functions import search_string
        
        rms, neutr, etot, efermi = [], [], [], []
            
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
            
            ptitle = ''
    
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
            print 'no rms extracted'
        
        return rms, neutr, etot, efermi, ptitle
    
    
    
        ### calculations ###
    
    def plot_kkr_calc(self, node, **kwargs):
        """TODO docstring"""
        
        # extract options from kwargs
        nofig = False
        if 'nofig' in kwargs.keys(): nofig = kwargs.pop('nofig')
        strucplot = True
        if 'strucplot' in kwargs.keys(): strucplot = kwargs.pop('strucplot')
        logscale = True
        if 'logscale' in kwargs.keys(): logscale = kwargs.pop('logscale')
            
        # plot structure
        if strucplot:
            plot_struc(node)
            
        rms, neutr, etot, efermi, ptitle = self.get_rms_kkrcalc(node)
    
        self.rmsplot(rms, neutr, nofig, ptitle, logscale)
    
    
    def plot_voro_calc(self, node, **kwargs):
        """TODO docstring"""
        
        strucplot = True
        if 'strucplot' in kwargs.keys(): strucplot = kwargs.pop('strucplot')
            
        # plot structure
        if strucplot:
            self.plot_struc(node)
            
        outdict = node.out.output_parameters.get_dict()
        # TODO maybe plot some output of voronoi
    
    
    def plot_kkrimp_calc(self, node, **kwargs):
        print 'not implemented yet'
        pass
    
    
        ### workflows ###
    
    
    def plot_kkr_dos(self, node, **kwargs):
        """TODO: docstring"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        
        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels = False, False, True
        if 'interpol' in kwargs.keys(): interpol = kwargs.pop('interpol')
        if 'all_atoms' in kwargs.keys(): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in kwargs.keys(): l_channels = kwargs.pop('l_channels') 
        nofig = False
        if 'nofig' in kwargs.keys(): nofig = kwargs.pop('nofig')
        if 'strucplot' in kwargs.keys(): strucplot = kwargs.pop('strucplot')
        
        if interpol:
            d = node.out.dos_data_interpol
        else:
            d = node.out.dos_data
            
        # extract structure (neede in dosplot to extract number of atoms and number of spins)
        struc, voro_parent = VoronoiCalculation.find_parent_structure(node.inp.remote_data)
            
        print d, struc
        
        # do dos plot after data was extracted
        self.dosplot(d, struc, nofig, all_atoms, l_channels, **kwargs)
    
    
    def plot_kkr_startpot(self, node, **kwargs):
        """TODO docstring"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from matplotlib.pyplot import axvline, legend, title
        from masci_tools.io.common_functions import get_Ry2eV
        
        strucplot = True
        if 'strucplot' in kwargs.keys(): strucplot = kwargs.pop('strucplot')
            
        # plot structure
        if strucplot:
            self.plot_struc(node)
            
        # extract structure (neede in dosplot to extract number of atoms and number of spins)
        struc, voro_parent = VoronoiCalculation.find_parent_structure(node)
            
        # print results
        print 'results:'
        self.plot_kkr_single_node(node.out.results_vorostart_wc, noshow=True, silent=True)
        
        # plot starting DOS
    
        try:
            d = node.out.last_doscal_dosdata
        except:
            d = None
        try:
            d_int = node.out.last_doscal_dosdata_interpol
        except:
            d_int = None
            
        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels = False, False, True
        if 'interpol' in kwargs.keys(): interpol = kwargs.pop('interpol')
        if 'all_atoms' in kwargs.keys(): all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in kwargs.keys(): l_channels = kwargs.pop('l_channels') 
        nofig = False
        if 'nofig' in kwargs.keys(): nofig = kwargs.pop('nofig')
        
        if interpol:
            d = d_int
            
        if d is not None:
            # do dos plot after data was extracted
            self.dosplot(d, struc, nofig, all_atoms, l_channels, **kwargs)
        
        # now add lines for emin, core states, EF
    
        # extract data for dos and energy contour plotting
        params_dict = node.out.last_voronoi_results.get_dict()
        emin = params_dict.get('emin_minus_efermi')
        emin_Ry = params_dict.get('emin')
        ef_Ry = emin_Ry - params_dict.get('emin_minus_efermi_Ry')
        ecore_max = params_dict.get('core_states_group').get('energy_highest_lying_core_state_per_atom', [])
        
        axvline(0, color='k', ls='--', label='EF')
        axvline(emin, color='r', ls='--', label='emin')
        axvline((ecore_max[0]-ef_Ry)*get_Ry2eV(), color='b', ls='--', label='ecore_max')
        if len(ecore_max)>1:
            [axvline((i-ef_Ry)*get_Ry2eV(), color='b', ls='--') for i in ecore_max[1:]]
        legend(loc=3)
        
        title(struc.get_formula())
    
    
    def plot_kkr_scf(self, node, **kwargs):
        """TODO docstring"""
        from aiida_kkr.calculations.kkr import KkrCalculation
        from numpy import sort
        from matplotlib.pyplot import axvline
        
        # structure plot only if structure is in inputs
        try:
            struc = node.inp.structure
            strucplot = True
            ptitle = struc.get_formula()
        except:
            strucplot = False
            ptitle = ''
            
        if 'strucplot' in kwargs.keys(): strucplot = kwargs.pop('strucplot')
        # plot structure
        if strucplot:
            self.plot_struc(struc)
            
        # next extract information from outputs 
        niter_calcs = [0]
        try:
            out_dict = node.out.output_kkr_scf_wc_ParameterResults.get_dict()
            neutr = out_dict['charge_neutrality_all_steps']
            rms = out_dict['convergence_values_all_steps']
        except: 
            # deal with unfinished workflow
            rms, neutr, etot, efermi = [], [], [], []
            outdict = node.get_outputs_dict()
            for key in sort(outdict.keys()):
                if isinstance(outdict[key], KkrCalculation):
                    kkrcalc = outdict[key]
                    rms_tmp, neutr_tmp, etot_tmp, efermi_tmp, ptitle_tmp = self.get_rms_kkrcalc(kkrcalc)
                    niter_calcs.append(len(rms_tmp))
                    rms += rms_tmp
                    neutr += neutr_tmp
        
        # extract options from kwargs
        nofig = False
        if 'nofig' in kwargs.keys(): nofig = kwargs.pop('nofig')
        logscale = True
        if 'logscale' in kwargs.keys(): logscale = kwargs.pop('logscale')
        only = None
        if 'only' in kwargs.keys(): only = kwargs.pop('only')
    
        # extract rms from calculations and plot
    
        if len(rms)>0:
            self.rmsplot(rms, neutr, nofig, ptitle, logscale, only)
            tmpsum = 0
            if len(niter_calcs)>1:
                for i in niter_calcs:
                    tmpsum+=i
                    axvline(tmpsum, color='k', ls=':')
    
    
    def plot_kkr_eos(self, node, **kwargs):
        """TODO docstring"""
        from numpy import sort
        from matplotlib.pyplot import figure, subplot, title, xlabel
        from aiida_kkr.workflows.voro_start import kkr_startpot_wc
        
        strucplot = True
        if 'strucplot' in kwargs.keys(): strucplot = kwargs.pop('strucplot')
            
        # plot structure
        if strucplot:
            self.plot_struc(node)
    
        outdict = node.get_outputs_dict()
        fig_open = False
        for key in sort(outdict.keys()):
            tmp = outdict[key]
            try:
                if tmp.process_label == u'kkr_startpot_wc':
                    self.plot_kkr_single_node(tmp, silent=True, noshow=True, strucplot=False) # startpot workflow
            except:
                pass
            try:
                if tmp.process_label == u'kkr_scf_wc':
                    if not fig_open:
                        figure()
                        fig_open = True
                    subplot(2,1,1)
                    self.plot_kkr_single_node(tmp, silent=True, strucplot=False, nofig=True, only='rms', noshow=True, **kwargs) # startpot workflow
                    xlabel('')
                    subplot(2,1,2)
                    self.plot_kkr_single_node(tmp, silent=True, strucplot=False, nofig=True, only='neutr', noshow=True, **kwargs) # startpot workflow   
                    title('') 
            except:
                pass


# some tests

#plot_kkr(29116) # structure
#plot_kkr(31521, silent=True) # remote
#plot_kkr(31815, silent=True) # folder
#plot_kkr(31817, silent=True) # parameter
#plot_kkr(31520, silent=False) # voro
#plot_kkr(32958, silent=False, strucplot=False) # kkr, scf run
#plot_kkr(31517, silent=True, strucplot=False) # startpot workflow
#plot_kkr(31819, silent=False, marker='o', interpol=False, strucplot=False) # dos workflow
#plot_kkr(34157, strucplot=False) # kkr_scf
#plot_kkr(31494, silent=False) # eos  workflow
