# -*- coding: utf-8 -*-
"""
contains plot_kkr class for node visualization
"""

import numpy as np
import matplotlib.pyplot as plt
from builtins import object, str
from ..tools.common_workfunctions import get_natyp
from masci_tools.io.common_functions import search_string
import numpy as np

__copyright__ = (u'Copyright (c), 2018, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.7.2'
__contributors__ = ('Philipp Rüßmann')


def get_datetime_from_str(calc, verbose=False):
    """
    Return a datetime object from the last time a calculation was checked by the scheduler.

    Every calculation should have the  'scheduler_lastchecktime' attribute which has the
    following format: '2023-11-08T22:44:13.543215+00:00'.
    This is converted to a datetime object that can be sorted.
    """
    from datetime import datetime
    # get last time stamp of scheduler from calculation attribute
    try:
        last_time_on_computer = calc.attributes['scheduler_lastchecktime']
    except:
        raise ValueError('Failed to get "scheduler_lastchecktime" from calculation.')
    # parse date and time from string
    date = last_time_on_computer.split('T')[0]
    time = last_time_on_computer.split('T')[1].split('.')[0]
    # change format
    datetime_str = date[2:].replace('-', '/') + ' ' + time
    # convert to datetime object
    datetime_object = datetime.strptime(datetime_str, '%y/%m/%d %H:%M:%S')

    if verbose:
        print(datetime_object)  # printed in default format

    #return datetime object of the last time the calculation was checked
    return datetime_object


def get_sorting_indices(calcs):
    """
    Get the sorting index for a list of calculations.

    For each calculation the datetime object of the last time the scheduler checked the
    calculation is extracted. This is then sorted and the sorting index array is returned.
    """
    datetimes = [get_datetime_from_str(calc) for calc in calcs]
    isort = np.array(datetimes).argsort()
    return isort


def remove_empty_atoms(show_empty_atoms, structure, silent=False):
    # check if empty sphere need to be removed for plotting (ase structgure cannot be constructed for alloys or vacancies)
    #print('in remove empty atoms:', structure.has_vacancies, ('X' in [i.kind_name for i in structure.sites]) )
    if structure.has_vacancies or ('X' in [i.kind_name for i in structure.sites]):
        #print("structure has vacancies, need to remove empty sites for plotting")
        from aiida.orm import StructureData
        stmp = StructureData(cell=structure.cell)
        for site in structure.sites:
            k = structure.get_kind(site.kind_name)
            pos = site.position
            if not (k.has_vacancies or 'X' in site.kind_name):
                stmp.append_atom(position=pos, symbols=k.symbol)
            elif show_empty_atoms:
                stmp.append_atom(position=pos, symbols='X')
            elif not silent:
                print('removing atom', site)
        stmp.set_pbc(structure.pbc)
        structure = stmp
    return structure


def _in_notebook():
    """
    Helper function to check if code is executed from within a jupyter notebook
    this is used to change to a different default visualization
    """
    try:
        from IPython import get_ipython
        if get_ipython() is None:
            return False
        if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    return True


def _has_ase_notebook():
    """
    Helper function to check if ase_notebook is installed
    """
    try:
        import ase_notebook
    except ImportError:
        return False
    return True


def _check_tk_gui(static):
    """
    check if tk gui can be openen, otherwise reset static to False
    this is only needed if we are not inside a notebook
    """
    if not _in_notebook() and not static:
        try:
            import tkinter as tk
            window = tk.Tk()
            window.quit()
        except:
            print('cannot open tk gui, fall back to static image')
            static = True

    return static


def save_fig_to_file(kwargs, filename0='plot_kkr_out.png'):
    """
    save the figure as a png file
    look for filename and static in kwargs
    save only if static is True after _check_tk_gui check to make it work in the command line script
    """
    import matplotlib.pyplot as plt
    if not _in_notebook():
        # check if static needs to be enforced
        static = kwargs.get('static', False)
        static = _check_tk_gui(static)
        if static:
            filename = kwargs.get('filename', filename0)
            print('saved static plot to ', filename)
            plt.savefig(filename)


def strucplot_ase_notebook(struc, **kwargs):
    """
    plotting function for aiida structure using ase_notebook visulaization
    """
    from ase_notebook import ViewConfig, AseView  # pylint: disable=import-error

    # extract some setting if given as kwargs
    repeat_uc = kwargs.get('repeat_uc', (1, 1, 1))
    canvas_size = kwargs.get('canvas_size', (300, 300))
    zoom = kwargs.get('zoom', 1.0)
    atom_opacity = kwargs.get('atom_opacity', 0.95)
    static = kwargs.get('static', False)
    rotations = kwargs.get('rotations', '-80x,-20y,-5z')

    if 'show_empty_atoms' in kwargs:
        show_empty_atoms = kwargs.pop('show_empty_atoms')
    else:
        show_empty_atoms = False
    struc = remove_empty_atoms(show_empty_atoms, struc, kwargs.get('silent', False))

    # check if static needs to be enforced
    static = _check_tk_gui(static)

    # set up structure viewer from ase_notebook
    config_dict = {
        'atom_show_label': True,
        'rotations': rotations,
        'show_uc_repeats': True,
        'show_bonds': False,
        'show_unit_cell': True,
        'canvas_size': canvas_size,
        'zoom': zoom,
        'show_axes': True,
        'canvas_background_opacity': 1.00,
        'canvas_color_background': 'white',
        'axes_length': 30,
        'atom_opacity': atom_opacity
    }

    config = ViewConfig(**config_dict)
    ase_view = AseView(config)

    # create ase atoms object from structure
    ase_atoms = struc.get_ase()

    # now create plot
    strucview = None
    if not static and _in_notebook():
        # render view in notebook
        strucview = ase_view.make_render(
            ase_atoms,
            center_in_uc=True,
            create_gui=True,
            repeat_uc=repeat_uc,
            use_atom_arrays=True,
        )
    elif _in_notebook():
        # static plot in notebook (svg)
        strucview = ase_view.make_svg(
            ase_atoms,
            center_in_uc=True,
            repeat_uc=repeat_uc,
        )
    elif static:
        strucview = ase_view.make_svg(
            ase_atoms,
            center_in_uc=True,
            repeat_uc=repeat_uc,
        )
        filename = kwargs.get('filename', 'plot_kkr_out_struc.svg')
        print('saved static plot in svg format to ', filename)
        strucview.saveas(filename)
    else:
        # open gui window
        ase_view.make_gui(
            ase_atoms,
            center_in_uc=True,
            repeat_uc=repeat_uc,
        )

    return strucview


def plot_imp_cluster(kkrimp_calc_node, **kwargs):
    """
    Plot impurity cluster from KkrimpCalculation node

    These kwargs can be used to control the behavior of the plotting tool:

    kwargs = {
        static = False,               # make gui or static (svg) images
        canvas_size = (300, 300),     # size of the canvas
        zoom = 1.0,                   # zoom, set to >1 (<1) to zoom in (out)
        atom_opacity = 0.95,          # set opacity level of the atoms, useful for overlapping atoms
        rotations = "-80x,-20y,-5z",  # rotation in degrees around x,y,z axes
        show_unit_cell = True,        # show the unit cell of the host
        filename = 'plot_kkr_out_impstruc.svg' # filename used for the export of a static svg image
    }

    """
    from aiida.orm import StructureData
    from aiida.common.constants import elements
    from ase_notebook import ViewConfig, AseView  # pylint: disable=import-error
    from aiida_kkr.calculations import VoronoiCalculation
    from aiida_kkr.tools.tools_kkrimp import create_scoef_array
    from masci_tools.io.common_functions import get_alat_from_bravais

    imp_info = kkrimp_calc_node.inputs.impurity_info.get_dict()
    structure0, _ = VoronoiCalculation.find_parent_structure(kkrimp_calc_node.inputs.host_Greenfunction_folder)

    # needed to transform from internal to Angstroem units
    alat = get_alat_from_bravais(np.array(structure0.cell), structure0.pbc[2])

    # extract infos from impurity info node
    rimp_rel = np.array(imp_info.get('Rimp_rel', [[0, 0, 0]]))
    zimp = imp_info.get('Zimp')
    if not (isinstance(zimp, list) or isinstance(zimp, np.ndarray)):
        zimp = [zimp]
    if 'Rimp_rel' in imp_info and 'imp_cls' in imp_info:
        imp_cls = np.array(imp_info['imp_cls'])
    else:
        ilayer = imp_info.get('ilayer_center', 0)
        radius = imp_info['Rcut']
        h = imp_info.get('hcut', -1.)
        vector = imp_info.get('cylinder_orient', [0., 0., 1.])
        i = imp_info.get('ilayer_center', 0)
        imp_cls = np.array(create_scoef_array(structure0, radius, h, vector, i))

    # adapt imp_cls from zimp+rimp_rel and find ghost atoms
    ghost_map = []
    for isite, site in enumerate(imp_cls):
        pos = site[:3]
        dmin = 1e9
        i0 = -1
        for iimp, z in enumerate(zimp):
            dist = np.sqrt(np.sum((rimp_rel[iimp] - pos)**2))
            if dist < dmin:
                dmin = dist
                i0 = iimp
        if dmin < 1e-5:
            ghost_map.append(False)
            imp_cls[isite, 4] = zimp[i0]
        else:
            ghost_map.append(True)
        # position to convert to Ang. units
        imp_cls[isite, :3] *= alat
    ghost_map = np.array(ghost_map)
    # create auxiliary structure
    struc_aux = StructureData(cell=structure0.cell)
    for site in imp_cls:
        struc_aux.append_atom(position=site[:3], symbols=elements[int(site[4])]['symbol'])

    # extract settings from kwargs
    canvas_size = kwargs.get('canvas_size', (300, 300))
    zoom = kwargs.get('zoom', 1.0)
    atom_opacity = kwargs.get('atom_opacity', 0.95)
    static = kwargs.get('static', False)
    rotations = kwargs.get('rotations', '-80x,-20y,-5z')
    show_unit_cell = kwargs.get('show_unit_cell', True)

    # check if static needs to be enforced
    static = _check_tk_gui(static)

    # set up structure viewer from ase_notebook
    config_dict = {
        'atom_show_label': True,
        'rotations': rotations,
        'show_uc_repeats': True,
        'show_bonds': False,
        'show_unit_cell': show_unit_cell,
        'canvas_size': canvas_size,
        'zoom': zoom,
        'show_axes': True,
        'canvas_background_opacity': 0.00,
        'canvas_color_background': 'white',
        'axes_length': 30,
        'atom_opacity': atom_opacity
    }
    config = ViewConfig(**config_dict)
    ase_view_imp = AseView(config)

    # create ase atoms object from auxiliary structure with ghost atoms mapping
    ase_atoms_impcls = struc_aux.get_ase()
    ase_atoms_impcls.set_array('ghost', ghost_map)

    # create plot
    strucview_imp = None
    if not static and _in_notebook():
        strucview_imp = ase_view_imp.make_render(
            ase_atoms_impcls,
            center_in_uc=True,
            create_gui=True,
            use_atom_arrays=True,
        )
    elif _in_notebook():
        strucview_imp = ase_view_imp.make_svg(
            ase_atoms_impcls,
            center_in_uc=True,
        )
    elif static:
        strucview_imp = ase_view_imp.make_svg(
            ase_atoms_impcls,
            center_in_uc=True,
        )
        filename = kwargs.get('filename', 'plot_kkr_out_impstruc.svg')
        print('saved static plot in svg format to ', filename)
        strucview_imp.saveas(filename)
    else:
        ase_view_imp.make_gui(  # pylint: disable=unexpected-keyword-arg
            ase_atoms_impcls,
            center_in_uc=True,
        )

    return strucview_imp


class plot_kkr(object):
    """
    Class grouping all functionality to plot typical nodes (calculations, workflows, ...) of the aiida-kkr plugin.

    :param nodes: node identifier which is to be visualized

    optional arguments:

    :param silent: print information about input node including inputs and outputs (default: False)
    :type silent: bool
    :param strucplot: plot structure using ase’s view function (default: False)
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
    :param debug: activate debug output
    :type debug: bool

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

        # used to keep track of structure plotting
        self.sview = None

        # debug mode
        self.debug = False
        if 'debug' in kwargs:
            self.debug = kwargs.pop('debug')
            print('start plot_kkr')
            print('kwargs:', kwargs)

        # grouping of node if a list of nodes is the input instead of a single node
        groupmode = False
        if type(nodes) == list:
            if len(nodes) > 1:
                groupmode = True
            else:
                nodes = nodes[0]

        if groupmode:
            from matplotlib.pyplot import show

            node_groups = self.group_nodes(nodes)
            for groupname in list(node_groups.keys()):
                print('\n==================================================================')
                print(f'Group of nodes: {groupname}\n')
                # some settings for groups
                if 'noshow' in list(kwargs.keys()):
                    _ = kwargs.pop('noshow')  # this is now removed from kwargs
                if 'only' in list(kwargs.keys()):
                    _ = kwargs.pop('only')  # this is now removed from kwargs
                if 'nofig' in list(kwargs.keys()):
                    _ = kwargs.pop('nofig')  # To revoke the 'nofig' kwarg from the list
                # now plot groups one after the other
                self.plot_group(groupname, node_groups, noshow=True, nofig=True, **kwargs)

            # finally show all plots
            show()
        elif nodes is not None:
            self.plot_kkr_single_node(nodes, **kwargs)

            display = kwargs.get('display', True)
            if display and self.sview is not None:  #self.classify_and_plot_node(nodes, return_name_only=True)=='struc':
                from IPython.display import display
                display(self.sview)

    ### main wrapper functions ###
    def group_nodes(self, nodes):
        """Go through list of nodes and group them together."""
        groups_dict = {}

        for node in nodes:
            node = get_node(node)
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
        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        nolegend = False
        if 'nolegend' in list(kwargs.keys()):
            nolegend = kwargs.pop('nolegend')
        # open a single new figure for each plot here
        if groupname in ['kkr', 'scf']:
            figure()
        for node in nodeslist:
            node = get_node(node)
            print(groupname)
            # open new figure for each plot in these groups
            if groupname in ['eos', 'dos', 'startpot']:
                figure()
            if groupname in ['kkr', 'scf']:
                subplot(2, 1, 1)
                self.plot_kkr_single_node(node, only='rms', label=f'pk= {node.pk}', **kwargs)
                xlabel('')  # remove overlapping x label in upper plot
                if not nolegend:
                    legend(fontsize='x-small')
                title('')
                subplot(2, 1, 2)
                self.plot_kkr_single_node(node, only='neutr', label=f'pk= {node.pk}', **kwargs)
                title('')  # remove duplicated plot title of lower plot
                if not nolegend:
                    legend(fontsize='x-small')
            else:
                self.plot_kkr_single_node(node, **kwargs)
            print('\n------------------------------------------------------------------\n')

    def plot_kkr_single_node(self, node, **kwargs):
        """ TODO docstring"""

        # determine if some output is printed to stdout
        silent = False
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')  # this is now removed from kwargs
        noshow = False
        if 'noshow' in list(kwargs.keys()):
            noshow = kwargs.pop('noshow')  # this is now removed from kwargs

        node = get_node(node)

        # print input and output nodes
        if not silent:
            self.print_clean_inouts(node)

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

        node = get_node(node)

        # basic aiida nodes
        if isinstance(node, DataFactory('structure')):
            if return_name_only:
                return 'struc'
            self.plot_struc(node, **kwargs)
        elif isinstance(node, DataFactory('dict')):
            if return_name_only:
                return 'para'
            print('node dict:')
            pprint(node.get_dict())
        elif isinstance(node, DataFactory('remote')):
            if return_name_only:
                return 'remote'
            print('computer name:', node.get_computer_name())
            print('remote path:', node.get_remote_path())
        elif isinstance(node, DataFactory('folder')):
            if return_name_only:
                return 'folder'
            print('abs path:')
            pprint(node.get_abs_path())
            print('folder content:')
            pprint(node.get_folder_list())
        # workflows
        elif node.process_label == u'kkr_dos_wc':
            if return_name_only:
                return 'dos'
            self.plot_kkr_dos(node, **kwargs)
        elif node.process_label == u'kkr_bs_wc':
            if return_name_only:
                return 'bs'
            self.plot_kkr_bs(node, **kwargs)
        elif node.process_label == u'kkr_startpot_wc':
            if return_name_only:
                return 'startpot'
            self.plot_kkr_startpot(node, **kwargs)
        elif node.process_label == u'kkr_scf_wc':
            if return_name_only:
                return 'scf'
            self.plot_kkr_scf(node, **kwargs)
        elif node.process_label == u'kkr_eos_wc':
            if return_name_only:
                return 'eos'
            self.plot_kkr_eos(node, **kwargs)
        elif node.process_label == u'kkr_imp_dos_wc':
            if return_name_only:
                return 'impdos'
            self.plot_kkrimp_dos_wc(node, **kwargs)
        elif node.process_label == u'kkr_imp_wc':
            if return_name_only:
                return 'imp'
            self.plot_kkrimp_wc(node, **kwargs)
        elif node.process_label == u'kkr_imp_sub_wc':
            if return_name_only:
                return 'impsub'
            self.plot_kkrimp_sub_wc(node, **kwargs)
        # calculations
        elif node.process_type == u'aiida.calculations:kkr.kkr':
            if return_name_only:
                return 'kkr'
            self.plot_kkr_calc(node, **kwargs)
        elif node.process_type == u'aiida.calculations:kkr.voro':
            if return_name_only:
                return 'voro'
            self.plot_voro_calc(node, **kwargs)
        elif node.process_type == u'aiida.calculations:kkr.kkrimp':
            if return_name_only:
                return 'kkrimp'
            self.plot_kkrimp_calc(node, **kwargs)
        elif node.process_type == 'aiida_kkr.workflows._combine_imps.combine_imps_wc':
            if return_name_only:
                return 'combine_imps'
            # extract kkr_imp_sub and plot it
            self.plot_kkrimp_wc(node, **kwargs)
        elif node.process_label == 'KKRnanoCalculation':
            raise TypeError(f'KKRnano input not implemented, yet: {type(node)} {node}')
        else:
            raise TypeError(
                f'input node neither a `Calculation` nor a `WorkChainNode` (i.e. workflow): {type(node)} {node}'
            )

    ### helper functions (structure plot, rms plot, dos plot, data extraction ...) ###

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
        print(f'pk, uuid: {node.pk} {node.uuid}')
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
            print(f'\nexit status: {node.exit_status} ({node.exit_message})')
        except:
            pass
        print()  # empty line at the end

    def plot_struc(self, node, **kwargs):
        """visualize structure using ase's `view` function"""
        from ase.visualize import view
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from aiida.orm import StructureData
        if not isinstance(node, StructureData):
            structure, voro_parent = VoronoiCalculation.find_parent_structure(node)
        else:
            structure = node
        if 'show_empty_atoms' in kwargs:
            show_empty_atoms = kwargs.pop('show_empty_atoms')
        else:
            show_empty_atoms = False
        structure = remove_empty_atoms(show_empty_atoms, structure, kwargs.get('silent', False))

        if _has_ase_notebook() and 'viewer' not in kwargs:
            # by default use ase_notebook if it is available
            self.sview = strucplot_ase_notebook(structure, show_empty_atoms=show_empty_atoms, **kwargs)
        else:
            # use ase's view function instead

            # now construct ase object and use ase's viewer
            ase_atoms = structure.get_ase()
            if 'silent' in kwargs:
                silent = kwargs.pop('silent')
            # remove things that are not understood bt ase's view
            for key in [
                'nofig', 'interpol', 'all_atoms', 'l_channels', 'sum_spins', 'logscale', 'switch_xy', 'iatom', 'label'
            ]:
                if key in kwargs:
                    _ = kwargs.pop(key)
            print(f"plotting structure using ase's `view` with kwargs={kwargs}")

            self.sview = view(ase_atoms, **kwargs)

    def dosplot(self, d, natoms, nofig, all_atoms, l_channels, sum_spins, switch_xy, switch_sign_spin2, **kwargs):
        """plot dos from xydata node"""
        from numpy import array, sum, arange
        from matplotlib.pyplot import plot, xlabel, ylabel, gca, figure, legend, fill_between
        import matplotlib as mpl
        from cycler import cycler

        # remove things that will not work for plotting
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')
        if 'filled' in list(kwargs.keys()):
            filled = kwargs.pop('filled')
        else:
            filled = False

        # plot only some atoms if 'iatom' is found in input
        show_atoms = []
        if 'iatom' in kwargs:
            show_atoms = kwargs.pop('iatom')
            if type(show_atoms) != list:
                show_atoms = [show_atoms]
            all_atoms = True  # need to set all_atoms to true

        # open new figure
        if not nofig:
            figure()

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
        if 'xshift' in kwargs:
            xshift = kwargs.pop('xshift')
        else:
            xshift = 0.

        nspin = len(y_all[0][1]) // natoms
        nspin2 = len(y_all[0][1]) // natoms  # copy of nspin becaus nspin is reset to 1 if sum_spins is set to True
        if sum_spins:
            nspin = 1

        xlbl = x_all[0] + ' (' + x_all[2] + ')'

        #tot only:
        if not l_channels:
            lmax = 1
        else:
            lmax = len(y_all)

        pcycle_default = mpl.rcParams['axes.prop_cycle']
        # change color cycler to match spin up/down colors
        if nspin == 2 and (not all_atoms or show_atoms != []):
            pcycle_values = pcycle_default.by_key()['color']
            if switch_sign_spin2:  # needed for impurity DOS
                n0 = len(show_atoms)
                if n0 == 0:
                    n0 = natoms
                #pcycle_values = array([j for i in range(n0) for j in pcycle_values]).reshape(-1)
                #pcycle_values = list(pcycle_values[:n0]) + list(pcycle_values[:n0])
                pcycle_values = array([[i, i] for i in pcycle_values]).reshape(-1)
            else:
                pcycle_values = array([[i, i] for i in pcycle_values]).reshape(-1)
            pcycle_default = cycler('color', pcycle_values)
        gca().set_prop_cycle(pcycle_default)

        if 'label' in kwargs:
            labels_all = kwargs.pop('label')
            if type(labels_all) != list:
                labels_all = [labels_all for i in range(natoms)]
        else:
            labels_all = None

        for il in range(lmax):
            y2 = y_all[il]
            # extract label
            ylbl = 'DOS (' + y2[2] + ')'
            # take data
            y2 = y2[1].copy()
            y2 = y2.reshape(natoms, nspin2, -1)
            x = x_all[1].copy() + xshift

            for ispin in range(nspin):
                if not all_atoms:
                    y = [sum(y2[:, ispin, :], axis=0)]  # artificial list so that y[iatom] works later on
                    if sum_spins and nspin2 == 2:
                        if switch_sign_spin2:
                            y[0] = -y[0] - sum(y2[:, 1, :], axis=0)
                        else:
                            y[0] = -y[0] + sum(y2[:, 1, :], axis=0)
                    natoms2 = 1
                    yladd = ''
                else:
                    natoms2 = natoms
                    y = y2[:, ispin, :]
                    if sum_spins and nspin2 == 2:
                        if switch_sign_spin2:
                            y = y + y2[:, 1, :]
                        else:
                            y = -y + y2[:, 1, :]

                for iatom in range(natoms2):
                    if iatom in show_atoms or show_atoms == []:
                        yladd = y_all[il][0].replace('dos ', '')
                        if all_atoms:
                            yladd += ', atom=' + str(iatom + 1)
                        if ispin > 0:
                            yladd = ''
                        if labels_all is not None and ispin == 0:
                            yladd = labels_all[iatom]
                        xplt = x[iatom * nspin + ispin] * xscale
                        yplt = y[iatom] * yscale
                        if ispin > 0 and switch_sign_spin2:
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

        if not nofig:
            figure()

        # allow to overwrite name for second quantity, plotted on second y axis
        name_second_y = 'charge neutrality'
        if rename_second is not None:
            name_second_y = rename_second

        if only is None:
            if 'label' not in list(kwargs.keys()):
                label = 'rms'
            else:
                label = kwargs.pop('label')
            plot(rms, '-xb', label=label)
            ax1 = gca()
            ylabel('rms', color='b')
            xlabel('iteration')
            twinx()
            if logscale:
                neutr = abs(array(neutr))
            if 'label' not in list(kwargs.keys()):
                label = name_second_y
            else:
                label = kwargs.pop('label')
            plot(neutr, '-or', label=label)
            ax2 = gca()
            ylabel(name_second_y, color='r')
            if logscale:
                ax1.set_yscale('log')
                ax2.set_yscale('log')
        else:  # individual plots of rms or neutrality (needed for eos plot)
            if only == 'rms':
                plot(rms, **kwargs)
                ylabel('rms')
                xlabel('iteration')
                ax1 = gca()
                if logscale:
                    ax1.set_yscale('log')
            elif only == 'neutr':
                if logscale:
                    neutr = abs(array(neutr))
                plot(neutr, **kwargs)
                ylabel(name_second_y)
                xlabel('iteration')
                ax1 = gca()
                if logscale:
                    ax1.set_yscale('log')
            else:
                raise ValueError(f'`only` can only be `rms or `neutr` but got {only}')
        title(ptitle)

    def get_rms_kkrcalc(self, node, title=None):
        """extract rms etc from kkr Calculation. Works for both finished and still running Calculations."""
        from aiida.engine import ProcessState

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
            rms, neutr, etot, efermi = get_rms_kkrcalc_from_remote(node)
        else:
            print('no rms extracted', node.process_state)

        return rms, neutr, etot, efermi, ptitle

    ### Calculations ###

    def plot_kkr_calc(self, node, **kwargs):
        """plot things for a kkr Calculation node"""

        # extract options from kwargs
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')
        strucplot = False
        if 'strucplot' in list(kwargs.keys()):
            strucplot = kwargs.pop('strucplot')
        logscale = True
        if 'logscale' in list(kwargs.keys()):
            logscale = kwargs.pop('logscale')
        only = None
        if 'only' in list(kwargs.keys()):
            only = kwargs.pop('only')
        silent = False
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')

        #print output
        if not silent:
            from pprint import pprint
            print('results dict (entries with `...` have been removed for this writeout for the sake of shortness):')
            if 'output_parameters' in node.get_outgoing().all_link_labels():
                results_dict = node.get_outgoing().get_node_by_label('output_parameters').get_dict()
                # remove symmetry descriptions from resuts dict before writting output
                if 'symmetries_group' in list(results_dict.keys()):
                    results_dict['symmetries_group']['symmetry_description'] = '...'
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

        if len(rms) > 1:
            self.rmsplot(rms, neutr, nofig, ptitle, logscale, only, label=label)
            # maybe save as file
            save_fig_to_file(kwargs, 'plot_kkr_out_rms.png')

        # try to plot dos and qdos data if Calculation was bandstructure or DOS run
        from subprocess import check_output
        from os import listdir
        from numpy import loadtxt, array, where
        from masci_tools.io.common_functions import open_general
        from masci_tools.vis.kkr_plot_bandstruc_qdos import dispersionplot
        from masci_tools.vis.kkr_plot_FS_qdos import FSqdos2D
        from masci_tools.vis.kkr_plot_dos import dosplot
        from matplotlib.pyplot import show, figure, title, xticks, xlabel, axvline
        from aiida import __version__ as aiida_core_version

        if node.is_finished_ok:
            retlist = node.outputs.retrieved.list_object_names()
            has_dos = 'dos.atom1' in retlist
            has_qvec = 'qvec.dat' in retlist
            has_qdos = False

            # remove already automatically set things from kwargs
            if 'ptitle' in list(kwargs.keys()):
                ptitle = kwargs.pop('ptitle')
            else:
                ptitle = f'pk= {node.pk}'
            if 'newfig' in list(kwargs.keys()):
                kwargs.pop('newfig')

            # qdos
            if has_qvec:
                has_qdos = 'qdos.01.1.dat' in retlist or 'qdos.001.1.dat' in retlist
                if has_qdos:
                    qdos_filenames, ne = get_qdos_filenames(node)

                    if ne > 1 or 'as_e_dimension' in list(kwargs.keys()):
                        a0 = get_a0_from_node(node)
                        ef = get_ef_from_parent(node)
                        data = get_qdos_data_from_node(node, qdos_filenames)
                        data_all = (a0, data, None, None, ef)
                        dispersionplot(
                            data_all=data_all, newfig=(not nofig), ptitle=ptitle, logscale=logscale, **kwargs
                        )
                        # add plot labels
                        try:
                            labels = node.inputs.kpoints.labels
                            ilbl = array([int(i[0]) for i in labels])
                            slbl = array([i[1] for i in labels])
                            m_overlap = where(abs(ilbl[1:] - ilbl[:-1]) == 1)
                            if len(m_overlap[0]) > 0:
                                for i in m_overlap[0]:
                                    slbl[i + 1] = '\n' + slbl[i + 1]
                            xticks(ilbl, slbl)
                            xlabel('')
                            [axvline(i, color='grey', ls=':') for i in ilbl]
                        except:
                            xlabel('id_kpt')

                        # maybe save as file
                        save_fig_to_file(kwargs, 'plot_kkr_out_bs.png')
                    else:
                        if int(aiida_core_version.split('.')[0]) < 2:
                            ef = get_ef_from_parent(node)
                            with node.outputs.retrieved.open('qvec.dat') as f:
                                FSqdos2D(f.name.replace('qvec.dat', ''), logscale=logscale, ef=ef, **kwargs)
                            # maybe save as file
                            save_fig_to_file(kwargs, 'plot_kkr_out_FS.png')
                        else:
                            raise ValueError('FS plot does not work with aiida-core>=2.0')

            # dos only if qdos was not plotted already
            if has_dos and not has_qdos:
                with node.outputs.retrieved.open('dos.atom1', mode='r') as f:
                    if not nofig:
                        figure()
                    dosplot(f, **kwargs)
                    title(ptitle)
                    # maybe save as file
                    save_fig_to_file(kwargs, 'plot_kkr_out_dos.png')

    def plot_voro_calc(self, node, **kwargs):
        """plot things for a voro Calculation node"""

        strucplot = False
        if 'strucplot' in list(kwargs.keys()):
            strucplot = kwargs.pop('strucplot')

        # plot structure
        if strucplot:
            self.plot_struc(node, **kwargs)

        # TODO maybe plot some output of voronoi
        #outdict = node.outputs.output_parameters.get_dict()

    def plot_kkrimp_calc(self, node, return_rms=False, return_stot=False, plot_rms=True, **kwargs):
        """plot things from a kkrimp Calculation node"""

        if self.debug:
            print('in plot_kkrimp_calc')
            print('kwargs:', kwargs)

        # plot impurity cluster
        if kwargs.get('strucplot', False):
            if _has_ase_notebook():
                self.sview = plot_imp_cluster(node, **kwargs)
            else:
                print('Cannot plot impurity structure because ase_notebook is not installed')
        # remove plotting-exclusive keys from kwargs
        for k in ['static', 'canvas_size', 'zoom', 'atom_opacity', 'rotations', 'show_unit_cell', 'strucplot']:
            if k in kwargs:
                kwargs.pop(k)

        # read data from output node
        rms_goal, rms = None, []
        if node.is_finished_ok:
            out_para = node.outputs.output_parameters
            out_para_dict = out_para.get_dict()
            out_para_dict['convergence_group']['rms_all_iterations']
            rms = out_para_dict['convergence_group']['rms_all_iterations']
            rms_goal = out_para_dict['convergence_group']['qbound']

            # extract total magnetic moment
            nspin = out_para_dict['nspin']
            if nspin > 1:
                try:
                    nat = out_para_dict['number_of_atoms_in_unit_cell']
                    s = np.array(out_para_dict['convergence_group']['total_spin_moment_all_iterations'][1], dtype=float)
                    ss = np.sqrt(np.sum(s**2, axis=1)).reshape(-1, nat)
                    stot = np.sum(ss, axis=1)
                except:
                    stot = None
            else:
                stot = None
        else:
            stot = None

        # make rms plot
        if plot_rms:
            if 'ptitle' in list(kwargs.keys()):
                ptitle = kwargs.pop('ptitle')
            else:
                ptitle = f'pk= {node.pk}'

            self.make_kkrimp_rmsplot([rms], [stot], [node], rms_goal, ptitle, **kwargs)

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
        sub_wf = [i.node for i in node.get_outgoing(node_class=kkr_imp_sub_wc).all()]
        if len(sub_wf) > 0:
            self.plot_kkrimp_sub_wc(sub_wf[0], **kwargs)

    def plot_kkrimp_sub_wc(self, node, **kwargs):
        """plot things from a kkrimp_sub_wc workflow"""
        from aiida_kkr.calculations import KkrimpCalculation

        if self.debug:
            print('in plot_kkrimp_sub_wc')
            print('kwargs:', kwargs)

        impcalcs = [i.node for i in node.get_outgoing(node_class=KkrimpCalculation).all()]

        # plot impurity cluster
        if len(impcalcs) > 0 and kwargs.get('strucplot', False):
            if _has_ase_notebook():
                self.sview = plot_imp_cluster(impcalcs[0], **kwargs)
            else:
                print('Cannot plot impurity structure because ase_notebook is not installed')
        # remove plotting-exclusive keys from kwargs
        for k in ['static', 'canvas_size', 'zoom', 'atom_opacity', 'rotations', 'show_unit_cell', 'strucplot']:
            if k in kwargs:
                kwargs.pop(k)

        # extract rms from calculations
        rms_all, stot_all = [], []
        rms_goal = None
        for impcalc in impcalcs:
            rms_tmp, rms_goal_tmp, stot_tmp = self.plot_kkrimp_calc(
                impcalc, return_rms=True, return_stot=True, plot_rms=False
            )
            rms_all.append(rms_tmp)
            if rms_goal_tmp is not None:
                if rms_goal is not None:
                    rms_goal = min(rms_goal, rms_goal_tmp)
                else:
                    rms_goal = rms_goal_tmp
            stot_all.append(stot_tmp)

        if 'ptitle' in list(kwargs.keys()):
            ptitle = kwargs.pop('ptitle')
        else:
            ptitle = f'pk= {node.pk}'

        self.make_kkrimp_rmsplot(rms_all, stot_all, impcalcs, rms_goal, ptitle, **kwargs)

    def make_kkrimp_rmsplot(self, rms_all, stot_all, list_of_impcalcs, rms_goal, ptitle, **kwargs):
        """
        plot rms and total spin moment of kkrimp calculation or series of kkrimp calculations
        """

        from numpy import array
        from matplotlib.pyplot import figure, subplot, axhline, axvline, gca, ylim

        # extract options from kwargs
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')
        logscale = True
        if 'logscale' in list(kwargs.keys()):
            logscale = kwargs.pop('logscale')
        if 'subplot' in list(kwargs.keys()):
            subplots = kwargs.pop('subplot')
        else:
            subplots = None
        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        else:
            label = None
        if 'only' in list(kwargs.keys()):
            only = kwargs.pop('only')
        else:
            only = None

        # plotting of convergence properties (rms etc.)
        if len(rms_all) > 0:
            # sort rms values and flatten array
            reorder_rms = get_sorting_indices(list_of_impcalcs)
            rms, niter_calcs, stot = [], [0], []
            rms_all_sorted = [rms_all[i] for i in reorder_rms]
            for i in rms_all_sorted:
                rms += list(i)
                niter_calcs.append(len(i) - 0.5)
            stot_sorted = [stot_all[i] for i in reorder_rms]
            for i in stot_sorted:
                if i is not None:
                    stot += list(i)
            # now plot
            if len(rms) > 0:
                if not nofig:
                    figure()
                if subplots is not None:
                    subplot(subplots[0], subplots[1], subplots[2])
                if rms_goal is not None:
                    axhline(rms_goal, color='grey', ls='--')
                self.rmsplot(
                    rms,
                    stot,
                    nofig=True,
                    ptitle=ptitle,
                    logscale=logscale,
                    only=only,
                    rename_second='sum(spinmom)',
                    label=label
                )
                # adapt y-limits to take care of showing spin-moment on sensible scale
                if only is None:
                    yl = gca().get_ylim()
                    ylim(yl[0], max(yl[1], 0.1))
                # add lines that indicate different calculations
                tmpsum = 1
                if not nofig and len(niter_calcs) > 1:
                    for i in niter_calcs:
                        tmpsum += i
                        axvline(tmpsum - 1, color='k', ls=':')
                # maybe save as file
                save_fig_to_file(kwargs, 'plot_kkr_out_rms.png')

    def plot_kkrimp_dos_wc(self, node, **kwargs):
        """plot things from a kkrimp_dos workflow node"""

        # try to plot dos and qdos data if Calculation was bandstructure or DOS run
        from os import listdir
        from numpy import loadtxt, array, where
        from masci_tools.vis.kkr_plot_FS_qdos import FSqdos2D
        from masci_tools.vis.kkr_plot_dos import dosplot
        from matplotlib.pyplot import show, figure, title, xticks, xlabel, axvline

        if self.debug:
            print('in plot_kkrimp_dos_wc')
            print('kwargs:', kwargs)

        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        ptitle = None
        if 'ptitle' in list(kwargs.keys()):
            ptitle = kwargs.pop('ptitle')
        if 'interpol' in list(kwargs.keys()):
            interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()):
            all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()):
            l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()):
            sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()):
            switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')
        if 'switch_sign_spin2' in list(kwargs.keys()):
            switch_sign_spin2 = kwargs.pop('switch_sign_spin2')
        else:
            switch_sign_spin2 = True
        if 'yscale' in list(kwargs.keys()):
            yscale = kwargs.pop('yscale')
        else:
            yscale = -1

        has_dos = False
        if interpol and 'dos_data_interpol' in node.outputs:
            d = node.outputs.dos_data_interpol
            has_dos = True
        elif 'dos_data' in node.outputs:
            d = node.outputs.dos_data
            has_dos = True

        if has_dos:
            calcnode = [i for i in node.called_descendants if i.process_label == 'KkrimpCalculation'][0]

            # plot impurity cluster
            if kwargs.get('strucplot', True):
                self.sview = plot_imp_cluster(calcnode, **kwargs)
            # remove plotting-exclusive keys from kwargs
            for k in ['static', 'canvas_size', 'zoom', 'atom_opacity', 'rotations', 'show_unit_cell', 'strucplot']:
                if k in kwargs:
                    kwargs.pop(k)

            if calcnode.is_finished_ok:
                natoms = len(calcnode.outputs.output_parameters.get_dict().get('charge_core_states_per_atom'))
                from aiida_kkr.tools import find_parent_structure
                natoms = get_natyp(find_parent_structure(calcnode))
                self.dosplot(
                    d,  # pylint: disable=possibly-used-before-assignment
                    natoms,
                    nofig,
                    all_atoms,
                    l_channels,
                    sum_spins,
                    switch_xy,
                    switch_sign_spin2,
                    yscale=yscale,
                    **kwargs
                )
                if ptitle is None:
                    title(f'pk= {node.pk}')
                else:
                    title(ptitle)
                # maybe save as file
                save_fig_to_file(kwargs, 'plot_kkr_out_dos.png')

    ### workflows ###

    def plot_kkr_dos(self, node, **kwargs):
        """plot outputs of a kkr_dos_wc workflow"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from matplotlib.pylab import title

        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()):
            interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()):
            all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()):
            l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()):
            sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()):
            switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()):
            strucplot = kwargs.pop('strucplot')
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')

        if node.is_finished_ok:
            if interpol:
                d = node.outputs.dos_data_interpol
            else:
                d = node.outputs.dos_data

            # extract structure (neede in dosplot to extract number of atoms and number of spins)
            struc, voro_parent = VoronoiCalculation.find_parent_structure(node.inputs.remote_data)

            # do dos plot after data was extracted
            if 'ptitle' in list(kwargs.keys()):
                ptitle = kwargs.pop('ptitle')
            else:
                ptitle = f'pk= {node.pk}'
            self.dosplot(d, get_natyp(struc), nofig, all_atoms, l_channels, sum_spins, switch_xy, False, **kwargs)
            title(ptitle)
            # maybe save as file
            save_fig_to_file(kwargs, 'plot_kkr_out_dos.png')

    def plot_kkr_bs(self, node, **kwargs):
        import matplotlib.pyplot as plt
        if node.is_finished_ok:
            BSF = node.outputs.BS_Data.get_array('BlochSpectralFunction')
            eng = node.outputs.BS_Data.get_array('energy_points')
            Kpts = node.outputs.BS_Data.get_array('Kpts')
            k_label = node.outputs.BS_Data.extras['k-labels']

            ixlbl = [int(i) for i in k_label.keys()]
            sxlbl = [i for i in k_label.values()]
            j = 0
            for i in ixlbl[:-1]:
                if (ixlbl[j + 1] - i) < 2:
                    sxlbl[j + 1] = str(sxlbl[j]) + '|' + str(sxlbl[j + 1])
                    sxlbl[j] = ''
                j += 1

            y, x = np.mgrid[slice(0, len(eng) + 1, 1), slice(0, len(Kpts[:, 0]) + 1, 1)]

            eng_extend = np.ones(len(eng[:]) + 1)

            eng = eng[::-1]
            eng_extend[:-1] = np.sort(eng)
            eng_extend[-1] = eng_extend[-2]

            y = np.array(y, float)
            for i in range(len(x[0, :])):
                y[:, i] = eng_extend

            nofig = kwargs.get('nofig', False)
            if not nofig:
                fig = plt.figure(figsize=(5, 5))

            # maybe change the colormap
            cmap = kwargs.get('cmap', plt.cm.viridis)  # pylint: disable=no-member

            # maybe scale and shift the x values
            xscale = kwargs.get('xscale', 1.0)
            xshift = kwargs.get('xshift', 0.0)
            x = (x + xshift) * xscale

            # maybe scale and shift the y values
            yscale = kwargs.get('yscale', 1.0)
            yshift = kwargs.get('yshift', 0.0)
            y = (y + yshift) * yscale

            # now create the plot
            plt.pcolormesh(x, y, np.log(abs(BSF.T)), cmap=cmap, edgecolor='face', rasterized=True)
            plt.ylabel('E-E_F (eV)')
            plt.xlabel('')

            # contol limits of the color scale
            clim = kwargs.get('clim', None)
            if clim is not None:
                plt.clim(clim[0], clim[1])
            else:
                # fix lower bound
                plt.clim(-6)

            show_cbar = kwargs.get('show_cbar', True)
            if show_cbar:
                plt.colorbar()

            plt.title(f'band structure from kkr_bs_wc (pk= {node.pk})')

            plt.xticks(ixlbl, sxlbl)
            plt.axhline(0, color='red', ls=':', lw=2)

            # maybe save as file
            save_fig_to_file(kwargs, 'plot_kkr_out_bs.png')

    def plot_kkr_startpot(self, node, **kwargs):
        """plot output of kkr_startpot_wc workflow"""
        from aiida_kkr.calculations.voro import VoronoiCalculation
        from aiida.common import exceptions
        from matplotlib.pyplot import axvline, legend, title
        from masci_tools.io.common_functions import get_Ry2eV

        strucplot = False
        if 'strucplot' in list(kwargs.keys()):
            strucplot = kwargs.pop('strucplot')

        silent = False
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')

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
            if link_triple.link_label == 'last_doscal_dosdata':
                d = link_triple.node
            elif link_triple.link_label == 'last_doscal_dosdata_interpol':
                d_int = link_triple.node
            elif 'CALL_WORK' in link_triple.link_label:
                if link_triple.link_label == 'kkr_dos_wc':
                    for link_triple2 in node.get_outgoing().all():
                        if link_triple2.link_label == 'dos_data':
                            d = link_triple2.node
                        elif link_triple2.link_label == 'dos_data_interpol':
                            d_int = link_triple2.node

        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()):
            interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()):
            all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()):
            l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()):
            sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()):
            switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')

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
        if emin is not None and params_dict != {}:
            ef_Ry = emin_Ry - params_dict.get('emin_minus_efermi_Ry')
        else:
            ef_Ry = None
        if params_dict != {}:
            ecore_max = params_dict.get('core_states_group').get('energy_highest_lying_core_state_per_atom', [])

        if d is not None:
            axvline(0, color='k', ls='--', label='EF')
            tit_add = ''
            if emin is not None:
                axvline(emin, color='r', ls='--', label='emin')
            if ef_Ry is not None and len(ecore_max) > 0:  # pylint: disable=possibly-used-before-assignment
                if abs((ecore_max[0] - ef_Ry) * get_Ry2eV() - emin) < 20:
                    axvline((ecore_max[0] - ef_Ry) * get_Ry2eV(), color='b', ls='--', label='ecore_max')
                else:
                    tit_add = f'; E_core<={(ecore_max[0] - ef_Ry) * get_Ry2eV():.2f}eV'
                if len(ecore_max) > 1:
                    [
                        axvline((i - ef_Ry) * get_Ry2eV(), color='b', ls='--')
                        for i in ecore_max[1:]
                        if abs((i - ef_Ry) * get_Ry2eV() - emin) < 20
                    ]
            if emin is not None:
                legend(loc=3, fontsize='x-small')

            title(struc.get_formula() + ', starting potential' + tit_add)

    def plot_kkr_scf(self, node, **kwargs):
        """plot outputs of a kkr_scf_wc workflow"""
        from aiida.orm import CalcJobNode, load_node
        from aiida_kkr.calculations.kkr import KkrCalculation
        from numpy import sort
        from matplotlib.pyplot import axvline, axhline, subplot, figure, title

        # structure plot only if structure is in inputs
        if 'structure' in node.inputs:
            struc = node.inputs.structure
        else:
            from aiida_kkr.tools import find_parent_structure
            struc = find_parent_structure(node.inputs.remote_data)

        try:
            ptitle = struc.get_formula()
        except:
            ptitle = ''

        strucplot = False
        if 'strucplot' in list(kwargs.keys()):
            strucplot = kwargs.pop('strucplot')

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
                if node.process_label == u'KkrCalculation':
                    kkrcalc = node
                    rms_tmp, neutr_tmp, etot_tmp, efermi_tmp, ptitle_tmp = self.get_rms_kkrcalc(kkrcalc)
                    if len(rms_tmp) > 0:
                        niter_calcs.append(len(rms_tmp))
                        rms += rms_tmp
                        neutr += neutr_tmp

        # extract options from kwargs
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')
        logscale = True
        if 'logscale' in list(kwargs.keys()):
            logscale = kwargs.pop('logscale')
        only = None
        if 'only' in list(kwargs.keys()):
            only = kwargs.pop('only')
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
        if len(rms) > 0 and not dos_only:
            if not nofig:
                figure()
            if subplots is not None:
                subplot(subplots[0], subplots[1], subplots[2])
            self.rmsplot(rms, neutr, True, ptitle, logscale, only, label=label)
            if only == 'rms' and rms_goal is not None:
                axhline(rms_goal, color='grey', ls='--')
            tmpsum = 1
            if not nofig and len(niter_calcs) > 1:
                for i in niter_calcs:
                    tmpsum += i
                    axvline(tmpsum - 1, color='k', ls=':')
            # maybe save as file
            save_fig_to_file(kwargs, 'plot_kkr_out_rms.png')
            did_plot = True
        else:
            did_plot = False

        # plot DOS

        # follow links until DOS data has been found
        d = None
        d_int = None
        from aiida_kkr.workflows import kkr_dos_wc
        links_dos = node.get_outgoing(node_class=kkr_dos_wc).all()
        if len(links_dos) > 0:
            dosnode = links_dos[0].node
            if 'dos_data' in dosnode.outputs:
                d = dosnode.outputs.dos_data
            if 'dos_data_interpol' in dosnode.outputs:
                d_int = dosnode.outputs.dos_data_interpol

        # extract all options that should not be passed on to plot function
        interpol, all_atoms, l_channels, sum_spins, switch_xy = True, False, True, False, False
        if 'interpol' in list(kwargs.keys()):
            interpol = kwargs.pop('interpol')
        if 'all_atoms' in list(kwargs.keys()):
            all_atoms = kwargs.pop('all_atoms')
        if 'l_channels' in list(kwargs.keys()):
            l_channels = kwargs.pop('l_channels')
        if 'sum_spins' in list(kwargs.keys()):
            sum_spins = kwargs.pop('sum_spins')
        if 'switch_xy' in list(kwargs.keys()):
            switch_xy = kwargs.pop('switch_xy')
        nofig = False
        if 'nofig' in list(kwargs.keys()):
            nofig = kwargs.pop('nofig')

        if interpol:
            d = d_int

        if d is not None:
            # do dos plot after data was extracted
            if 'ptitle' in list(kwargs.keys()):
                ptitle = kwargs.pop('ptitle')
            else:
                ptitle = f'pk= {node.pk}'
            self.dosplot(d, len(struc.sites), nofig, all_atoms, l_channels, sum_spins, switch_xy, False, **kwargs)
            plt.title(ptitle)
            # maybe save as file
            save_fig_to_file(kwargs, 'plot_kkr_out_dos.png')

        return did_plot

    def plot_kkr_eos(self, node, **kwargs):
        """plot outputs of a kkr_eos workflow"""
        from numpy import sort, array, where
        from matplotlib.pyplot import figure, title, xlabel, legend, axvline, plot, annotate
        from aiida.orm import load_node
        from aiida_kkr.workflows.voro_start import kkr_startpot_wc
        from ase.eos import EquationOfState

        strucplot = False
        if 'strucplot' in list(kwargs.keys()):
            strucplot = kwargs.pop('strucplot')

        # plot structure
        if strucplot:
            self.plot_struc(node, **kwargs)

        # remove unused things from kwargs
        if 'label' in list(kwargs.keys()):
            label = kwargs.pop('label')
        if 'noshow' in list(kwargs.keys()):
            kwargs.pop('noshow')
        if 'only' in list(kwargs.keys()):
            kwargs.pop('only')
        if 'nofig' in list(kwargs.keys()):
            kwargs.pop('nofig')
        if 'strucplot' in list(kwargs.keys()):
            kwargs.pop('strucplot')
        silent = False
        if 'silent' in list(kwargs.keys()):
            silent = kwargs.pop('silent')
        nolegend = False
        if 'nolegend' in list(kwargs.keys()):
            nolegend = kwargs.pop('nolegend')

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
            if key.link_label != 'CALL_CALC':
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
                    did_plot = self.plot_kkr_scf(
                        tmpnode,
                        silent=True,
                        strucplot=False,
                        nofig=fig_open,
                        only='rms',
                        noshow=True,
                        label=f'pk={tmpnode.pk}',
                        subplot=(2, 1, 1),
                        **kwargs
                    )  # scf workflow, rms only
                    if did_plot and not fig_open:
                        fig_open = True
                    if did_plot:
                        xlabel('')  # remove overlapping x label in upper plot
                    if did_plot and not nolegend:
                        legend(loc=3, fontsize='x-small', ncol=2)
                    # plot charge neutrality
                    self.plot_kkr_scf(
                        tmpnode,
                        silent=True,
                        strucplot=False,
                        nofig=True,
                        only='neutr',
                        noshow=True,
                        label=f'pk={tmpnode.pk}',
                        subplot=(2, 1, 2),
                        **kwargs
                    )  # scf workflow, rms only
                    if did_plot:
                        title('')  # remove overlapping title
                    if did_plot and not nolegend:
                        legend(loc=3, fontsize='x-small', ncol=2)
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
            names = sort([
                name for name in list(node.outputs.eos_results.get_dict().get('sub_workflow_uuids').keys())
                if 'kkr_scf' in name
            ])
            pks = array([
                load_node(node.outputs.eos_results.get_dict().get('sub_workflow_uuids')[name]).pk for name in names
            ])
            mask = []
            for i in range(len(pks)):
                s = scalings[i]
                m = where(scalings_all == s)
                pk = pks[m][0]
                ie = e[m][0]
                iv = v[m][0]
                if not nolegend:
                    annotate(text=f'pk={pk}', xy=(iv, ie))

            # investigate fit quality by fitting without first/last datapoint
            if len(e) > 4:

                eos = EquationOfState(v[1:-1], e[1:-1], eos=fitfunc_gs)
                v01, e01, B1 = eos.fit()

                # take out smallest data point
                eos = EquationOfState(v[1:], e[1:], eos=fitfunc_gs)
                v01, e01, B1 = eos.fit()

                print('# relative differences to full fit: V0, E0, B (without smallest volume)')
                print(f'{abs(1 - v01 / v0)} {abs(1 - e01 / e0)} {abs(1 - B1 / B)}')

                if len(e) > 5:
                    # also take out largest volume
                    eos = EquationOfState(v[1:-1], e[1:-1], eos=fitfunc_gs)
                    v02, e02, B2 = eos.fit()

                    print('\n# V0, E0, B (without smallest and largest volume)')
                    print(f'{abs(1 - v02 / v0)} {abs(1 - e02 / e0)} {abs(1 - B2 / B)}')
        except:
            pass  # do nothing if no eos data there


def get_node(node):
    """Get node from pk or uuid"""
    from aiida.orm import load_node, Node
    # load node if pk or uuid is given
    if type(node) == int:
        node = load_node(node)
    elif type(node) == type(''):
        node = load_node(node)
    elif isinstance(node, Node):
        pass
    else:
        raise TypeError(
            "input node should either be the nodes pk (int), it's uuid (str) or the node itself (aiida.orm.Node). Got type(node)={}"
            .format(type(node))
        )
    return node


def get_rms_kkrcalc_from_remote(node, **kwargs):
    """
    connect to remote folder and extract the rms etc from the out_kkr file
    if kwargs are given, a dict is retuned for each of the search keys
    """
    from aiida.common.folders import SandboxFolder
    from masci_tools.io.common_functions import search_string
    # extract info needed to open transport
    #c = node.inputs.code
    #comp = c.computer
    #authinfo = comp.get_authinfo(c.user)
    #transport = authinfo.get_transport()
    transport = node.get_transport()

    out_kkr = ''
    rms, neutr, etot, efermi = [], [], [], []
    return_dict = {}

    # now get contents of out_kkr using remote call of 'cat'
    with SandboxFolder() as tempfolder:
        with tempfolder.open('tempfile', 'w') as f:
            try:
                node.outputs.remote_folder.getfile('out_kkr', f.name)
                has_outfile = True
            except:
                has_outfile = False
        if has_outfile:
            with tempfolder.open('tempfile', 'r') as f:
                out_kkr = f.readlines()

    # now extract rms, charge neutrality, total energy and value of Fermi energy
    if has_outfile:
        itmp = 0
        while itmp >= 0:
            itmp = search_string('rms', out_kkr)
            if itmp >= 0:
                tmpline = out_kkr.pop(itmp)
                tmpval = float(tmpline.split('=')[1].split()[0].replace('D', 'e'))
                rms.append(tmpval)
        itmp = 0
        while itmp >= 0:
            itmp = search_string('charge neutrality', out_kkr)
            if itmp >= 0:
                tmpline = out_kkr.pop(itmp)
                tmpval = float(tmpline.split('=')[1].split()[0].replace('D', 'e'))
                neutr.append(tmpval)
        itmp = 0
        while itmp >= 0:
            itmp = search_string('TOTAL ENERGY in ryd', out_kkr)
            if itmp >= 0:
                tmpline = out_kkr.pop(itmp)
                tmpval = float(tmpline.split(':')[1].split()[0].replace('D', 'e'))
                etot.append(tmpval)
        itmp = 0
        while itmp >= 0:
            itmp = search_string('E FERMI', out_kkr)
            if itmp >= 0:
                tmpline = out_kkr.pop(itmp)
                tmpval = float(tmpline.split('FERMI')[1].split()[0].replace('D', 'e'))
                efermi.append(tmpval)
        # search for all additinoal keys
        for key in kwargs:
            return_dict[key] = []
            itmp = 0
            while itmp >= 0:
                itmp = search_string(key, out_kkr)
                if itmp >= 0:
                    tmpline = out_kkr.pop(itmp)
                    tmpval = float(tmpline.split(key)[1])
                return_dict[key].append(tmpval)

    if return_dict == {}:
        return rms, neutr, etot, efermi
    else:
        return rms, neutr, etot, efermi, return_dict


def get_ef_from_parent(node):
    """Extract Fermi level from parent calculation"""

    try:
        parent_calc = node.inputs.parent_folder.get_incoming().first().node
        try:
            # parent is KkrCalc
            ef = parent_calc.outputs.output_parameters.get_dict()['fermi_energy']
        except:
            # parent is scf workflow
            ef = parent_calc.outputs.last_calc_out['fermi_energy']
    except:
        with node.outputs.retrieved.open('output.0.txt', mode='r') as file_handle:
            txt = file_handle.readlines()
            iline = search_string('Fermi energy =', txt)
            if iline >= 0:
                ef = txt[iline].split('=')[1]
                ef = float(ef.split()[0])
            else:
                ef = None
        if ef is None:
            raise ValueError('error loading Fermi energy from outfile, retry extracting from parent')

    return ef


def get_a0_from_node(node):
    """Extract factor 2pi/alat from inputcard"""
    with node.outputs.retrieved.open('inputcard') as _f:
        txt = _f.readlines()
        txt = [line for line in txt if 'ALATBASIS' in line]
    if len(txt) == 0:
        raise ValueError('ALATBASIS not found in inputcard')
    alat = float(txt[0].split('=')[1].split()[0])
    a0 = 2 * np.pi / alat / 0.52918
    return a0


def get_qdos_filenames(node):
    """get the sorted list of qdos filenames and find number of energy points"""
    # get qdos filenames
    retrieved = node.outputs.retrieved
    qdos_filenames = np.sort([i for i in retrieved.list_object_names() if 'qdos.' in i])

    # read number of energy points
    with retrieved.open(qdos_filenames[0]) as _f:
        ne = len(set(np.loadtxt(_f)[:, 0]))

    return qdos_filenames, ne


def get_qdos_data_from_node(node, qdos_filenames):
    """Extract the qdos data (summed over all atoms)"""
    from aiida import orm

    if 'saved_dispersion_data_uuid' in node.extras:
        # only return the existing numpy array if it exists
        data = orm.load_node(node.extras['saved_dispersion_data_uuid']).get_array('data')
    else:
        # read qdos data and store as extra
        for i, fname in enumerate(qdos_filenames):
            # read data from qdos file
            with node.outputs.retrieved.open(fname) as _f:
                tmp = np.loadtxt(_f)
            # sum up all contributions
            if i == 0:
                data = tmp
            else:
                if len(fname.split('.')[1]) == 2:
                    # old data format that also kept the imaginary part
                    # of the energy point
                    data[:, 5:] += tmp[:, 5:]
                else:
                    # new data format
                    data[:, 4:] += tmp[:, 4:]

        # now store as extra
        # TODO: make this nicer and connect to provenance graph with calcfunction
        data_node = orm.ArrayData()
        data_node.set_array('data', data)
        data_node.store()
        node.set_extra('saved_dispersion_data_uuid', data_node.uuid)

    return data
