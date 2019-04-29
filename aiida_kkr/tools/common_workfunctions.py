# -*- coding: utf-8 -*-
"""
Here workfunctions and normal functions using aiida-stuff (typically used
within workfunctions) are collected.
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from aiida.common.exceptions import InputValidationError
from aiida.engine import calcfunction
from aiida.plugins import DataFactory
from masci_tools.io.kkr_params import kkrparams
from masci_tools.io.common_functions import open_general
from six.moves import range

#define aiida structures from DataFactory of aiida
Dict = DataFactory('dict')

# keys that are used by aiida-kkr some something else than KKR parameters
_ignored_keys = ['ef_set', 'use_input_alat']

@calcfunction
def update_params_wf(parameternode, updatenode):
    """
    Work function to update a KKR input parameter node.
    Stores new node in database and creates a link from old parameter node to new node
    Returns updated parameter node using update_params function

    :note: Input nodes need to be valid aiida Dict objects.

    :param parameternode: Input aiida Dict node cotaining KKR specific parameters
    :param updatenode: Input aiida Dict node containing a dictionary with the parameters that are supposed to be changed.

    :note: If 'nodename' is contained in dict of updatenode the string corresponding to this key will be used as nodename for the new node. Otherwise a default name is used
    :note: Similar for 'nodedesc' which gives new node a description

    :example: updated_params = Dict(dict={'nodename': 'my_changed_name', 'nodedesc': 'My description text', 'EMIN': -1, 'RMAX': 10.})
              new_params_node = update_params_wf(input_node, updated_params)
    """
    updatenode_dict = updatenode.get_dict()
    if 'nodename' in list(updatenode_dict.keys()):
        # take nodename out of dict (should only contain valid KKR parameter)
        nodename = updatenode_dict.pop('nodename')
    else:
        nodename = None
    if 'nodedesc' in list(updatenode_dict.keys()):
        # take nodename out of dict (should only contain valid KKR parameter later on)
        nodedesc = updatenode_dict.pop('nodedesc')
    else:
        nodedesc = None

    # do nothing if updatenode is empty
    if len(list(updatenode_dict.keys()))==0:
        print('Input node is empty, do nothing!')
        raise InputValidationError('Nothing to store in input')
    #
    new_parameternode = update_params(parameternode, nodename=nodename,
                                      nodedesc=nodedesc, **updatenode_dict)
    return new_parameternode


def update_params(node, nodename=None, nodedesc=None, **kwargs):
    """
    Update parameter node given with the values given as kwargs.
    Returns new node.

    :param node: Input parameter node (needs to be valid KKR input parameter node).
    :param **kwargs: Input keys with values as in kkrparams.
    :param linkname: Input linkname string. Give link from old to new node a name .
                     If no linkname is given linkname defaults to 'updated parameters'

    :return: parameter node

    :example usage: OutputNode = KkrCalculation.update_params(InputNode, EMIN=-1, NSTEPS=30)

    :note: Keys are set as in kkrparams class. Check documentation of kkrparams for further information.
    :note: By default nodename is 'updated KKR parameters' and description contains list of changed
    """
    # check if node is a valid KKR parameters node
    if not isinstance(node, Dict):
        print('Input node is not a valid Dict node')
        raise InputValidationError('update_params needs valid parameter node as input')

    #initialize temporary kkrparams instance containing all possible KKR parameters
    params = kkrparams()

    # extract input dict from node
    inp_params = node.get_dict()

    # check if input dict contains only values for KKR parameters
    for key in inp_params:
        if key not in list(params.values.keys()) and key not in _ignored_keys:
            print('Input node contains unvalid key "{}"'.format(key))
            raise InputValidationError('unvalid key "{}" in input parameter node'.format(key))

    # copy values from input node
    for key in inp_params:
        value = inp_params[key]
        params.set_value(key, value, silent=True)

    # to keep track of changed values:
    changed_params = {}

    # check if values are given as **kwargs (otherwise return input node)
    if len(kwargs)==0:
        print('No additional input keys given, return input node')
        return node.clone()
    else:
        for key in kwargs:
            # check if value of 'key' should be set (either because it differs from old para node or because it was not set at all)
            update_value = False
            if key in list(inp_params.keys()):
                if kwargs[key] != inp_params[key]:
                    update_value = True
            else:
                update_value = True
            if update_value:
                params.set_value(key, kwargs[key], silent=True)
                changed_params[key] = kwargs[key]

    if len(list(changed_params.keys()))==0:
        print('No keys have been changed, return input node')
        return node.clone()

    # set linkname with input or default value
    if nodename is None or type(nodename) is not str:
        nodename = 'updated KKR parameters'
    if nodedesc is None or type(nodedesc) is not str:
        nodedesc = 'changed parameters: {}'.format(changed_params)


    # create new node
    ParaNode = Dict(dict=params.values)
    ParaNode.label = nodename
    ParaNode.description = nodedesc

    return ParaNode

#TODO implment VCA functionality
# maybe one starts from a calculation closest to the VCA case and slowly
# increase ZATOM which violates the _do_never_modify rule in KKR calculation
# this should then create a new structure and modify the old potential accordingly
# general rule: Nover destroy the data provenance!!!
@calcfunction
def prepare_VCA_structure_wf():
    pass

def prepare_VCA_structure():
    pass


#TODO implement 2D input helper
# a helper workfunction would be nice to create the vacuum region etc. for 2D calculation
@calcfunction
def prepare_2Dcalc_wf():
    pass

def prepare_2Dcalc():
    pass


def test_and_get_codenode(codenode, expected_code_type, use_exceptions=False):
    """
    Pass a code node and an expected code (plugin) type. Check that the
    code exists, is unique, and return the Code object.

    :param codenode: the name of the code to load (in the form label@machine)
    :param expected_code_type: a string with the plugin that is expected to
      be loaded. In case no plugins exist with the given name, show all existing
      plugins of that type
    :param use_exceptions: if True, raise a ValueError exception instead of
      calling sys.exit(1)
    :return: a Code object

    :example usage: from kkr_scf workflow::

        if 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                error = ("The code you provided for voronoi  does not "
                         "use the plugin kkr.voro")
                self.control_end_wc(error)
    """
    import sys
    from aiida.common.exceptions import NotExistent
    from aiida.orm import Code


    try:
        if codenode is None:
            raise ValueError
        code = codenode
        if code.get_input_plugin_name() != expected_code_type:
            raise ValueError
    except (NotExistent, ValueError):
        from aiida.orm.querybuilder import QueryBuilder
        qb = QueryBuilder()
        qb.append(Code,
                  filters={'attributes.input_plugin':
                               {'==': expected_code_type}},
                  project='*')

        valid_code_labels = ["{}@{}".format(c.label, c.get_computer().name)
                             for [c] in qb.all()]

        if valid_code_labels:
            msg = ("Pass as further parameter a valid code label.\n"
                   "Valid labels with a {} executable are:\n".format(
                expected_code_type))
            msg += "\n".join("* {}".format(l) for l in valid_code_labels)

            if use_exceptions:
                raise ValueError(msg)
            else:
                print(msg, file=sys.stderr)
                sys.exit(1)
        else:
            msg = ("Code not valid, and no valid codes for {}.\n"
                   "Configure at least one first using\n"
                   "    verdi code setup".format(
                expected_code_type))
            if use_exceptions:
                raise ValueError(msg)
            else:
                print(msg, file=sys.stderr)
                sys.exit(1)

    return code


def get_inputs_kkr(code, remote, options, label='', description='', parameters=None, serial=False, imp_info=None):
    """
    Get the input for a voronoi calc.
    Wrapper for KkrProcess setting structure, code, options, label, description etc.
    :param code: a valid KKRcode installation (e.g. input from Code.get_from_string('codename@computername'))
    :param remote: remote directory of parent calculation (Voronoi or previous KKR calculation)

    """
    from aiida_kkr.calculations.kkr import KkrCalculation

    # then reuse common inputs setter
    builder = get_inputs_common(KkrCalculation, code, remote, None, options, label,
                               description, parameters, serial, imp_info)

    return builder



def get_inputs_kkrimporter(code, remote, options, label='', description='', parameters=None, serial=False):
    """
    Get the input for a voronoi calc.
    Wrapper for KkrProcess setting structure, code, options, label, description etc.
    """
    from aiida_kkr.calculations.kkr import KkrCalculation
    KkrProcess = KkrCalculation.process()

    # then reuse common inputs setter
    inputs = get_inputs_common(KkrProcess, code, remote, None, options, label,
                               description, parameters, serial)

    return inputs


def get_inputs_voronoi(code, structure, options, label='', description='', params=None, serial=True):
    """
    Get the input for a voronoi calc.
    Wrapper for VoronoiProcess setting structure, code, options, label, description etc.
    """
    # get process for VoronoiCalculation
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # then reuse common inputs setter all options
    builder = get_inputs_common(VoronoiCalculation, code, None, structure, options, label,
                               description, params, serial)

    return builder


def get_inputs_kkrimp(code, options, label='', description='', parameters=None, serial=False, imp_info=None, host_GF=None, imp_pot=None, kkrimp_remote=None):
    """
    Get the input for a kkrimp calc.
    Wrapper for KkrimpProcess setting structure, code, options, label, description etc.
    :param code: a valid KKRimpcode installation (e.g. input from Code.get_from_string('codename@computername'))
    TBD
    """

    from aiida_kkr.calculations.kkrimp import KkrimpCalculation

    # then reuse common inputs setter
    builder = get_inputs_common(KkrimpCalculation, code, None, None, options, label,
                               description, parameters, serial, imp_info, host_GF, imp_pot, kkrimp_remote)

    return builder


def get_inputs_common(calculation, code, remote, structure, options, label, description, params, serial, imp_info=None, host_GF=None, imp_pot=None, kkrimp_remote=None):
    """
    Base function common in get_inputs_* functions for different codes
    """
    inputs = calculation.get_builder()

    if structure:
        inputs.structure = structure

    if remote:
        inputs.parent_folder = remote

    if code:
        inputs.code = code

    if params:
        inputs.parameters = params

    if not options:
        options = {}

    #for key, val in options.iteritems():
    #    if val==None:
    #        #leave them out, otherwise the dict schema won't validate
    #        continue
    #    else:
    #        inputs.options[key] = val

    if description:
        inputs.metadata.description = description
    else:
        inputs.metadata.description = ''

    if label:
        inputs.metadata.label = label
    else:
        inputs.metadata.label = ''

    if serial:
        options['withmpi'] = False # for now
        options['resources'] = {"num_machines": 1}

    if options:
        inputs.metadata.options = options
    '''
    options = {
    "max_wallclock_seconds": int,
    "resources": dict,
    "custom_scheduler_commands": unicode,
    "queue_name": basestring,
    "computer": Computer,
    "withmpi": bool,
    "mpirun_extra_params": Any(list, tuple),
    "import_sys_environment": bool,
    "environment_variables": dict,
    "priority": unicode,
    "max_memory_kb": int,
    "prepend_text": unicode,
    "append_text": unicode}
    '''

    # for kkrimp calculations
    if imp_info is not None:
        inputs.impurity_info = imp_info

    if host_GF is not None:
        inputs.host_Greenfunction_folder = host_GF

    if imp_pot is not None:
        inputs.impurity_potential = imp_pot

    if kkrimp_remote is not None:
        inputs.parent_calc_folder = kkrimp_remote

    return inputs


def get_parent_paranode(remote_data):
    """
    Return the input parameter of the parent calulation giving the remote_data node
    """
    inp_calc = remote_data.get_incoming(link_label_filter='remote_folder').first().node
    inp_para = inp_calc.get_incoming(link_label_filter='parameters').first().node
    return inp_para


def generate_inputcard_from_structure(parameters, structure, input_filename, parent_calc=None, shapes=None, isvoronoi=False, use_input_alat=False, vca_structure=False):
    """
    Takes information from parameter and structure data and writes input file 'input_filename'

    :param parameters: input parameters node containing KKR-related input parameter
    :param structure: input structure node containing lattice information
    :param input_filename: input filename, typically called 'inputcard'

    optional arguments
    :param parent_calc: input parent calculation node used to determine if EMIN
                        parameter is automatically overwritten (from voronoi output)
                        or not
    :param shapes: input shapes array (set automatically by
                   aiida_kkr.calculations.Kkrcaluation and shall not be overwritten)
    :param isvoronoi: tell whether or not the parameter set is for a voronoi calculation or kkr calculation (have different lists of mandatory keys)
    :param use_input_alat: True/False, determines whether the input alat value is taken or the new alat is computed from the Bravais vectors


    :note: assumes valid structure and parameters, i.e. for 2D case all necessary
           information has to be given. This is checked with function
           'check_2D_input' called in aiida_kkr.calculations.Kkrcaluation
    """

    from aiida.common.constants import elements as PeriodicTableElements
    from numpy import array
    from masci_tools.io.kkr_params import kkrparams
    from masci_tools.io.common_functions import get_Ang2aBohr, get_alat_from_bravais
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # initialize list of warnings
    warnings = []

    #list of globally used constants
    a_to_bohr = get_Ang2aBohr()


    # Get the connection between coordination number and element symbol
    # maybe do in a differnt way

    _atomic_numbers = {data['symbol']: num for num,
                    data in PeriodicTableElements.items()}

    # KKR wants units in bohr
    bravais = array(structure.cell)*a_to_bohr
    alat_input = parameters.get_dict().get('ALATBASIS')
    if use_input_alat and alat_input is not None:
        alat = alat_input
    else:
        alat = get_alat_from_bravais(bravais, is3D=structure.pbc[2])
        wmess = 'found alat in input parameters, this will trigger scaling of RMAX, GMAX and RCLUSTZ!'
        print('WARNING: '+wmess)
        warnings.append(wmess)
    bravais = bravais/alat

    sites = structure.sites
    naez = len(sites)
    positions = []
    charges = []
    weights = [] # for CPA
    isitelist = [] # counter sites array for CPA
    isite = 0
    for site in sites:
        pos = site.position
        #TODO maybe convert to rel pos and make sure that type is right for script (array or tuple)
        abspos = array(pos)*a_to_bohr/alat # also in units of alat
        positions.append(abspos)
        isite += 1
        sitekind = structure.get_kind(site.kind_name)
        for ikind in range(len(sitekind.symbols)):
            site_symbol = sitekind.symbols[ikind]
            if sitekind.is_alloy:
                wght = sitekind.weights[ikind]
            else:
                wght = 1.
            if not sitekind.has_vacancies:
                zatom_tmp = _atomic_numbers[site_symbol]
            else:
                zatom_tmp = 0.0
            if vca_structure and ikind>0 and not isvoronoi:
                # for VCA case take weighted average (only for KKR code, voronoi code uses zatom of first site for dummy calculation)
                zatom  = zatom*wght_last + zatom_tmp*wght
                # also reset weight to 1
                wght = 1.
            else:
                zatom = zatom_tmp
                if vca_structure and isvoronoi:
                    wght = 1.

            wght_last = wght # for VCA mode

            # make sure that for VCA only averaged position is written (or first for voronoi code)
            if ( (vca_structure and ((len(sitekind.symbols)==1) or
                                     (not isvoronoi and ikind==1) or
                                     (isvoronoi and ikind==0)))
                 or (not vca_structure) ):
                charges.append(zatom)
                weights.append(wght)
                isitelist.append(isite)

    weights = array(weights)
    isitelist = array(isitelist)
    charges = array(charges)
    positions = array(positions)

    # workaround for voronoi calculation with Zatom=83 (Bi potential not there!)
    if isvoronoi:
        from numpy import where
        mask_replace_Bi_Pb = where(charges==83)
        if len(mask_replace_Bi_Pb[0])>0:
            charges[mask_replace_Bi_Pb] = 82
            wmess = 'Bi potential not available, using Pb instead!!!'
            print('WARNING: '+wmess)
            warnings.append(wmess)


    ######################################
    # Prepare keywords for kkr from input structure

    # get parameter dictionary
    input_dict = parameters.get_dict()

    # remove special keys that are used for special cases but are not part of the KKR parameter set
    for key in _ignored_keys:
        if input_dict.get(key) is not None:
            wmess = 'automatically removing value of key {}'.format(key)
            print('WARNING: '+wmess)
            warnings.append(wmess)
            input_dict.pop(key)

    # get rid of structure related inputs that are overwritten from structure input
    for key in ['BRAVAIS', 'ALATBASIS', 'NAEZ', '<ZATOM>', '<RBASIS>', 'CARTESIAN']:
        if input_dict.get(key) is not None:
            wmess = 'automatically removing value of key {}'.format(key)
            print('WARNING: '+wmess)
            warnings.append(wmess)
            input_dict.pop(key)

    # automatically rescale RMAX, GMAX, RCLUSTZ, RCLUSTXY which are scaled with the lattice constant
    if alat_input is not None:
        if input_dict.get('RMAX') is not None:
            wmess = 'rescale RMAX: {}'.format(alat_input/alat)
            print('WARNING: '+wmess)
            warnings.append(wmess)
            input_dict['RMAX'] = input_dict['RMAX']*alat_input/alat
        if input_dict.get('GMAX') is not None:
            wmess='rescale GMAX: {}'.format(1/(alat_input/alat))
            print('WARNING: '+wmess)
            warnings.append(wmess)
            input_dict['GMAX'] = input_dict['GMAX']*1/(alat_input/alat)
        if input_dict.get('RCLUSTZ') is not None:
            wmess='rescale RCLUSTZ: {}'.format(alat_input/alat)
            print('WARNING: '+wmess)
            warnings.append(wmess)
            input_dict['RCLUSTZ'] = input_dict['RCLUSTZ']*alat_input/alat
        if input_dict.get('RCLUSTXY') is not None:
            wmess='rescale RCLUSTXY: {}'.format(alat_input/alat)
            print('WARNING: '+wmess)
            warnings.append(wmess)
            input_dict['RCLUSTXY'] = input_dict['RCLUSTXY']*alat_input/alat

    # empty kkrparams instance (contains formatting info etc.)
    if not isvoronoi:
        params = kkrparams()
    else:
        params = kkrparams(params_type='voronoi')

    # for KKR calculation set EMIN automatically from parent_calc (always in res.emin of voronoi and kkr) if not provided in input node
    if ('EMIN' not in list(input_dict.keys()) or input_dict['EMIN'] is None) and parent_calc is not None:
        wmess='Overwriting EMIN with value from parent calculation {}'.format(parent_calc)
        print('WARNING: '+wmess)
        warnings.append(wmess)
        if parent_calc.process_class == VoronoiCalculation:
            emin = parent_calc.outputs.output_parameters.get_dict().get('emin')
        else:
            emin = parent_calc.outputs.output_parameters.get_dict().get('energy_contour_group').get('emin')
        print('Setting emin:',emin, 'is emin None?',emin is None)
        params.set_value('EMIN', emin)

    # overwrite keywords with input parameter
    for key in list(input_dict.keys()):
        params.set_value(key, input_dict[key], silent=True)

    # Write input to file (the parameters that are set here are not allowed to be modfied externally)
    params.set_multiple_values(BRAVAIS=bravais, ALATBASIS=alat, NAEZ=naez,
                               ZATOM=charges, RBASIS=positions, CARTESIAN=True)
    # for CPA case:
    if len(weights)>naez:
        natyp = len(weights)
        params.set_value('NATYP', natyp)
        params.set_value('<CPA-CONC>', weights)
        params.set_value('<SITE>', isitelist)
    else:
        natyp = naez

    # write shapes (extracted from voronoi parent automatically in kkr calculation plugin)
    if shapes is not None:
        params.set_value('<SHAPE>', shapes)

    # change input values of 2D input to new alat:
    rbl = params.get_value('<RBLEFT>')
    rbr = params.get_value('<RBRIGHT>')
    zper_l = params.get_value('ZPERIODL')
    zper_r = params.get_value('ZPERIODR')
    if rbl is not None: params.set_value('<RBLEFT>', array(rbl)*a_to_bohr/alat)
    if rbr is not None: params.set_value('<RBRIGHT>', array(rbr)*a_to_bohr/alat)
    if zper_l is not None: params.set_value('ZPERIODL', array(zper_l)*a_to_bohr/alat)
    if zper_r is not None: params.set_value('ZPERIODR', array(zper_r)*a_to_bohr/alat)

    # write inputfile
    params.fill_keywords_to_inputfile(output=input_filename)

    nspin = params.get_value('NSPIN')

    newsosol = False
    if 'NEWSOSOL' in params.get_value('RUNOPT'):
        newsosol = True

    return natyp, nspin, newsosol, warnings


def check_2Dinput_consistency(structure, parameters):
    """
    Check if structure and parameter data are complete and matching.

    :param input: structure, needs to be a valid aiida StructureData node
    :param input: parameters, needs to be valid aiida Dict node

    returns (False, errormessage) if an inconsistency has been found, otherwise return (True, '2D consistency check complete')
    """
    # default is bulk, get 2D info from structure.pbc info (periodic boundary contitions)
    is2D = False
    if not all(structure.pbc):
        # check periodicity, assumes finite size in z-direction
        if structure.pbc != (True, True, False):
            return (False, "Structure.pbc is neither (True, True, True) for bulk nor (True, True, False) for surface calculation!")
        is2D = True

    # check for necessary info in 2D case
    inp_dict = parameters.get_dict()
    set_keys = [i for i in list(inp_dict.keys()) if inp_dict[i] is not None]
    has2Dinfo = True
    for icheck in ['INTERFACE', '<NRBASIS>', '<RBLEFT>', '<RBRIGHT>', 'ZPERIODL', 'ZPERIODR', '<NLBASIS>']:
        if icheck not in set_keys:
            has2Dinfo = False
    if has2Dinfo and not inp_dict['INTERFACE'] and is2D:
        return (False, "'INTERFACE' parameter set to False but structure is 2D")

    if has2Dinfo!=is2D:
        if is2D:
            return (False, "2D info given in parameters but structure is 3D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}".format(is2D, has2Dinfo, set_keys))
        else:
            return (False, "3D info given in parameters but structure is 2D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}".format(is2D, has2Dinfo, set_keys))

    # if everything is ok:
    return (True, "2D consistency check complete")


def structure_from_params(parameters):
    """
    Construct aiida structure out of kkr parameter set (if ALATBASIS, RBASIS, ZATOM etc. are given)

    :param input: parameters, kkrparams object with structure information set (e.g. extracted from read_inputcard function)

    :returns: success, boolean to determine if structure creatoin was successful
    :returns: structure, an aiida StructureData object
    """
    from masci_tools.io.common_functions import get_aBohr2Ang
    from aiida.common.constants import elements as PeriodicTableElements
    from masci_tools.io.kkr_params import kkrparams
    from numpy import array

    #check input
    if not isinstance(parameters, kkrparams):
        raise InputValidationError('input parameters needs to be a "kkrparams" instance!')

    # initialize some stuff
    StructureData = DataFactory('structure')
    is_complete = True
    for icheck in ['<ZATOM>', '<RBASIS>', 'BRAVAIS', 'ALATBASIS']:
        if parameters.get_value(icheck) is None:
            is_complete = False

    # set natyp
    natyp = parameters.get_value('NATYP')
    naez = parameters.get_value('NAEZ')
    if natyp is None:
        if naez is None:
            is_complete = False
        else:
            natyp = naez

    # check if all necessary info for 2D calculation is there
    if parameters.get_value('INTERFACE'):
        for icheck in ['<NRBASIS>', '<RBLEFT>', '<RBRIGHT>', 'ZPERIODL', 'ZPERIODR', '<NLBASIS>']:
            if parameters.get_value(icheck) is None:
                is_complete = False

    # check CPA case
    if natyp != naez:
        for icheck in ['<SITE>', '<CPA-CONC>']:
            if parameters.get_value(icheck) is None:
                is_complete = False

    if not is_complete:
        return is_complete, StructureData()

    # extract cell using BRAVAIS and ALATBASIS and create empty structure
    alat = parameters.get_value('ALATBASIS')
    cell = array(parameters.get_value('BRAVAIS')) * alat * get_aBohr2Ang()
    struc = StructureData(cell=cell)

    # extract sites with positions, charges/Atom labels, weights
    pos_all = array(parameters.get_value('<RBASIS>')) # positions in units of alat
    if not parameters.get_value('CARTESIAN'):
        # convert from internal to cartesian coordinates
        for isite in range(len(pos_all)):
            tmp_pos = pos_all[isite]
            pos_all[isite] = tmp_pos[0]*cell[0]+tmp_pos[1]*cell[1]+tmp_pos[2]*cell[2] # cell already contains alat factor to convert to Ang. units
    else:
        pos_all = pos_all * alat * get_aBohr2Ang() # now positions are in Ang. units

    # extract atom numbers
    zatom_all = parameters.get_value('<ZATOM>')
    # convert to list if input contains a single entry only
    if type(zatom_all) != list:
        zatom_all = [zatom_all]
        pos_all = [pos_all]

    #extract weights and sites for CPA calculations
    if natyp==naez:
        weights = [1. for i in range(natyp)]
        sites = list(range(1,natyp+1))
    else:
        weights = parameters.get_value('<CPA-CONC>')
        sites = parameters.get_value('<SITE>')

    # fill structure from zatom, weights and sites information
    for isite in sites:
        pos = pos_all[sites.index(isite)]
        weight = weights[sites.index(isite)]
        if abs(zatom_all[isite-1]-int(zatom_all[isite-1]))>10**-4:
            # TODO deal with VCA (non-integer zatom)
            print('VCA not implemented yet, stopping here!')
            raise NotImplementedError('VCA functionality not implemented')

        if zatom_all[isite-1]<1:
            symbol = 'H'
            weight = 0.0
            struc.append_atom(position=pos, symbols='H', weights=0.0, mass=1.0)
        else:
            symbol = PeriodicTableElements.get(zatom_all[isite-1]).get('symbol')
            struc.append_atom(position=pos, symbols=symbol, weights=weight)

    # set correct pbc for 2D case
    if parameters.get_value('INTERFACE'):
        struc.set_pbc((True, True, False))

    # finally return structure
    return is_complete, struc

@calcfunction
def neworder_potential_wf(settings_node, parent_calc_folder, **kwargs) : #, parent_calc_folder2=None):
    """
    Workfunction to create database structure for aiida_kkr.tools.modify_potential.neworder_potential function
    A temporary file is written in a Sandbox folder on the computer specified via
    the input computer node before the output potential is stored as SingleFileData
    in the Database.

    :param settings_node: settings for the neworder_potentail function (Dict)
    :param parent_calc_folder: parent calculation remote folder node where the input
        potential is retreived from (RemoteData)
    :param parent_calc_folder2: *optional*, parent calculation remote folder node where
        the second input potential is retreived from in case 'pot2' and 'replace_newpos'
        are also set in settings_node (RemoteData)

    :returns: output_potential node (SingleFileData)

    .. note::

        The settings_node dictionary needs to be of the following form::

            settings_dict = {'pot1': '<filename_input_potential>',  'out_pot': '<filename_output_potential>', 'neworder': [list of intended order in output potential]}

        Optional entries are::

            'pot2': '<filename_second_input_file>'
            'replace_newpos': [[position in neworder list which is replace with potential from pot2, position in pot2 that is chosen for replacement]]
            'label': 'label_for_output_node'
            'description': 'longer_description_for_output_node'
    """
    import os
    from aiida_kkr.tools.tools_kkrimp import modify_potential
    from aiida.common.folders import SandboxFolder
    from aiida.common.exceptions import UniquenessError
    from aiida.orm import CalcJobNode
    from aiida.plugins import DataFactory

    if 'parent_calc_folder2' in list(kwargs.keys()):
        parent_calc_folder2=kwargs.get('parent_calc_folder2', None)
    else:
        parent_calc_folder2=None

    # get aiida data types used here
    Dict = DataFactory('dict')
    RemoteData = DataFactory('remote')
    SingleFileData = DataFactory('singlefile')

    # check input consistency
    if not isinstance(settings_node, Dict):
        raise InputValidationError('settings_node needs to be a valid aiida Dict node')
    if not isinstance(parent_calc_folder, RemoteData):
        raise InputValidationError('parent_calc_folder needs to be a valid aiida RemoteData node')
    if parent_calc_folder2 is not None and not isinstance(parent_calc_folder2, RemoteData):
        raise InputValidationError('parent_calc_folder2 needs to be a valid aiida RemoteData node')

    settings_dict = settings_node.get_dict()
    pot1 = settings_dict.get('pot1', None)
    if pot1 is None:
        raise InputValidationError('settings_node_dict needs to have key "pot1" containing the filename of the input potential')
    out_pot = settings_dict.get('out_pot', None)
    if out_pot is None:
        raise InputValidationError('settings_node_dict needs to have key "out_pot" containing the filename of the input potential')
    neworder = settings_dict.get('neworder', None)
    if neworder is None:
        raise InputValidationError('settings_node_dict needs to have key "neworder" containing the list of new positions')
    pot2 = settings_dict.get('pot2', None)
    replace_newpos = settings_dict.get('replace_newpos', None)

    # Create Sandbox folder for generation of output potential file
    # and construct output potential
    with SandboxFolder() as tempfolder:
        # Get abolute paths of input files from parent calc and filename
        parent_calcs = parent_calc_folder.get_incoming(node_class=CalcJobNode).all()
        n_parents = len(parent_calcs)
        if n_parents != 1:
            raise UniquenessError(
                    "Input RemoteData is child of {} "
                    "calculation{}, while it should have a single parent"
                    "".format(n_parents, "" if n_parents == 0 else "s"))
        else:
            parent_calc = parent_calcs[0].node
        pot1_fhandle = parent_calc.outputs.retrieved.open(pot1)

        # extract nspin from parent calc's input parameter node
        nspin = parent_calc.inputs.parameters.get_dict().get('NSPIN')
        neworder_spin = []
        for iatom in neworder:
            for ispin in range(nspin):
                neworder_spin.append(iatom*nspin+ispin)
        neworder = neworder_spin

        # Copy optional files?
        if pot2 is not None and parent_calc_folder2 is not None:
            parent_calcs = parent_calc_folder2.get_incoming(node_class=CalcJobNode).all()
            n_parents = len(parent_calcs)
            if n_parents != 1:
                raise UniquenessError(
                        "Input RemoteData of parent_calc_folder2 is child of {} "
                        "calculation{}, while it should have a single parent"
                        "".format(n_parents, "" if n_parents == 0 else "s"))
            else:
                parent_calc = parent_calcs[0].node
            pot2_fhandle = parent_calc.outputs.retrieved.open(pot2)
        else:
            pot2_fhandle = None

        # change file path to Sandbox folder accordingly
        out_pot_fhandle = tempfolder.open(out_pot, u'w')

        # run neworder_potential function
        modify_potential().neworder_potential(pot1_fhandle, out_pot_fhandle, neworder, potfile_2=pot2_fhandle,
                                              replace_from_pot2=replace_newpos)

        # store output potential to SingleFileData
        output_potential_sfd_node = SingleFileData(file=tempfolder.open(out_pot))

        lbl = settings_dict.get('label', None)
        if lbl is not None:
            output_potential_sfd_node.label = lbl
        desc = settings_dict.get('description', None)
        if desc is not None:
            output_potential_sfd_node.description = desc

        #TODO create shapefun sfd node accordingly
        """
        out_shape_path =

        output_shapefun_sfd_node = SingleFileData(file=out_shape_path)

        lbl2 = settings_dict.get('label_shape', None)
        if lbl2 is None and lbl is not None:
            lbl2 = lbl
        if lbl2 is not None:
            output_shapefun_sfd_node.label = lbl2
        desc2 = settings_dict.get('description_shape', None)
        if desc2 is None and desc is not None:
            desc2 = desc
        if desc2 is not None:
            output_shapefun_sfd_node.description = desc2

        return output_potential_sfd_node, output_shapefun_sfd_node
        """
        return output_potential_sfd_node



def vca_check(structure, parameters):
    """

    """
    nsites = 0
    for site in structure.sites:
        sitekind = structure.get_kind(site.kind_name)
        nsites += len(sitekind.symbols)
    # VCA mode if CPAINFO = [-1,-1] first
    try:
        if parameters.get_dict().get('CPAINFO')[0]<0:
            params_vca_mode= True
        else:
            params_vca_mode = False
    except:
        params_vca_mode = False
    # check if structure supports VCA mode
    vca_structure = False
    if params_vca_mode:
        if nsites>len(structure.sites):
            vca_structure = True

    return vca_structure


def kick_out_corestates(potfile, potfile_out, emin):
    """
    Read potential file and kick out all core states that lie higher than emin.
    If no core state lies higher than emin then the output potential will be the same as the input potential
    :param potfile: input potential
    :param potfile_out: output potential where some core states are kicked out
    :param emin: minimal energy above which all core states are kicked out from potential
    :returns: number of lines that have been deleted
    """
    from masci_tools.io.common_functions import get_corestates_from_potential
    from numpy import where, array

    # read core states
    nstates, energies, lmoments = get_corestates_from_potential(potfile)

    with open_general(potfile) as f:
        # read potential file
        txt = f.readlines()

        #get start of each potential part
        istarts = [iline for iline in range(len(txt)) if 'POTENTIAL' in txt[iline]]
        all_lines = list(range(len(txt))) # index array

        # change list of core states
        for ipot in range(len(nstates)):
            if nstates[ipot]>0:
                m = where(energies[ipot]>emin)
                if len(m[0])>0:
                    istart = istarts[ipot]
                    # change number of core states in potential
                    print(txt[istart+6])
                    txt[istart+6] = '%i 1\n'%(nstates[ipot]-len(m[0]))
                    # now remove energy line accordingly
                    for ie_out in m[0][::-1]:
                        m_out = where(array(all_lines)==istart+6+ie_out+1)[0][0]
                        e_out = all_lines.pop(m_out)

        # find number of deleted lines
        num_deleted = len(txt)-len(all_lines)

        if num_deleted>0:
            # write output potential
            with open_general(potfile_out, u'w') as f2:
                txt_new = list(array(txt)[all_lines])
                f2.writelines(txt_new)

        # return number of lines that were deleted
        return num_deleted


@calcfunction
def kick_out_corestates_wf(potential_sfd, emin):
    """
    Workfunction that kicks out all core states from single file data potential that are higher than emin.
    :param potential_sfd: SingleFileData type of potential
    :param emin: Energy threshold above which all core states are removed from potential (Float)
    :returns: potential without core states higher than emin (SingleFileData)
    """
    from aiida.common.folders import SandboxFolder
    from aiida.plugins import DataFactory

    SingleFileData = DataFactory('singlefile')

    with SandboxFolder() as tmpdir:
        with tmpdir.open('potential_deleted_core_states', 'w') as potfile_out: 
            with potential_sfd.open(potential_sfd.filename) as potfile_in:
                num_deleted = kick_out_corestates(potfile_in, potfile_out, emin)
                if num_deleted>0:
                    potential_nocore_sfd = SingleFileData(file=potfile_out)
                    potential_nocore_sfd.store()

    # return potential
    if num_deleted>0:
        return potential_nocore_sfd
    else:
        return potential_sfd.clone()
