# -*- coding: utf-8 -*-
"""
Here workfunctions and normal functions using aiida-stuff (typically used
within workfunctions) are collected.
"""

from aiida.common.exceptions import InputValidationError
from aiida.engine import calcfunction
from aiida.orm import Dict
from masci_tools.io.kkr_params import kkrparams
from builtins import str

# keys that are used by aiida-kkr some something else than KKR parameters
_ignored_keys = ['ef_set', 'use_input_alat', '<NEWVERSION_BDG>', '<DECOUPLE_SPINS_CHEBY>']
_ignored_keys += [i.upper() for i in _ignored_keys]


@calcfunction
def update_params_wf(parameternode, updatenode, **link_inputs):
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
    if len(list(updatenode_dict.keys())) == 0:
        print('Input node is empty, do nothing!')
        raise InputValidationError('Nothing to store in input')
    #
    new_parameternode = update_params(parameternode, nodename=nodename, nodedesc=nodedesc, **updatenode_dict)
    return new_parameternode


def update_params(node, nodename=None, nodedesc=None, strict=False, **kwargs):
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
    :note: If kwargs contain the key `add_direct`, then no kkrparams instance is used and no checks are performed but the dictionary is filled directly!
    :note: By default nodename is 'updated KKR parameters' and description contains list of changed
    """
    # check if node is a valid KKR parameters node
    if not isinstance(node, Dict):
        print('Input node is not a valid Dict node')
        raise InputValidationError('update_params needs valid parameter node as input')

    # check if add_direct is in kwargs (shortcuts checks of kkrparams by not using the kkrparams class to set the dict)
    add_direct = False
    if 'add_direct' in list(kwargs.keys()):
        add_direct = kwargs.pop('add_direct')

    # initialize temporary kkrparams instance containing all possible KKR parameters
    if not add_direct:
        params = kkrparams()
    else:
        params = {}

    # extract input dict from node
    inp_params = node.get_dict()

    # check if input dict contains only values for KKR parameters
    if not add_direct:
        remove_keys = []
        for key in inp_params:
            if key not in list(params.values.keys()) and key not in _ignored_keys:
                print(f'WARNING: Input node contains invalid key "{key}"')
                if strict:
                    raise InputValidationError(f'invalid key "{key}" in input parameter node')
                else:
                    # print a warning and remove the key
                    print(f'ignore this key/value pair: {key}: {inp_params.get(key)}')
                    remove_keys.append(key)
        for key in remove_keys:
            inp_params.pop(key)

    # copy values from input node
    for key in inp_params:
        value = inp_params[key]
        if not add_direct:
            params.set_value(key, value, silent=True)
        else:
            params[key] = value

    # to keep track of changed values:
    changed_params = {}

    # check if values are given as **kwargs (otherwise return input node)
    if len(kwargs) == 0:
        print('No additional input keys given, return input node')
        return node.clone()
    for key in kwargs:
        # check if value of 'key' should be set (either because it differs from old para node or because it was not set at all)
        update_value = False
        if key in list(inp_params.keys()):
            if kwargs[key] != inp_params[key]:
                update_value = True
        else:
            update_value = True
        if update_value:
            if not add_direct:
                params.set_value(key, kwargs[key], silent=True)
            else:
                params[key] = kwargs[key]
            changed_params[key] = kwargs[key]

    if len(list(changed_params.keys())) == 0:
        print('No keys have been changed, return input node')
        return node.clone()

    # set linkname with input or default value
    if nodename is None or not isinstance(nodename, str):
        nodename = 'updated KKR parameters'
    if nodedesc is None or not isinstance(nodedesc, str):
        nodedesc = f'changed parameters: {changed_params}'

    # create new node
    if not add_direct:
        ParaNode = Dict(params.values)
    else:
        ParaNode = Dict(params)
    ParaNode.label = nodename
    ParaNode.description = nodedesc

    return ParaNode


# TODO implment VCA functionality
# maybe one starts from a calculation closest to the VCA case and slowly
# increase ZATOM which violates the _do_never_modify rule in KKR calculation
# this should then create a new structure and modify the old potential accordingly
# general rule: Nover destroy the data provenance!!!


@calcfunction
def prepare_VCA_structure_wf():
    pass


def prepare_VCA_structure():
    pass


# TODO implement 2D input helper
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
        qb.append(Code, filters={'attributes.input_plugin': {'==': expected_code_type}}, project='*')

        valid_code_labels = [f'{c.label}@{c.computer.label}' for [c] in qb.all()]

        if valid_code_labels:
            msg = (
                'Pass as further parameter a valid code label.\n'
                'Valid labels with a {} executable are:\n'.format(expected_code_type)
            )
            msg += '\n'.join(f'* {label}' for label in valid_code_labels)

            if use_exceptions:
                raise ValueError(msg)
            else:
                print(msg, file=sys.stderr)
                sys.exit(1)
        else:
            msg = (
                'Code not valid, and no valid codes for {}.\n'
                'Configure at least one first using\n'
                '    verdi code setup'.format(expected_code_type)
            )
            if use_exceptions:
                raise ValueError(msg)
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
    builder = get_inputs_common(
        KkrCalculation, code, remote, None, options, label, description, parameters, serial, imp_info
    )

    return builder


def get_inputs_kkrimporter(code, remote, options, label='', description='', parameters=None, serial=False):
    """
    Get the input for a voronoi calc.
    Wrapper for KkrProcess setting structure, code, options, label, description etc.
    """
    from aiida_kkr.calculations.kkr import KkrCalculation
    KkrProcess = KkrCalculation.process()

    # then reuse common inputs setter
    inputs = get_inputs_common(KkrProcess, code, remote, None, options, label, description, parameters, serial)

    return inputs


def get_inputs_voronoi(code, structure, options, label='', description='', params=None, serial=True, parent_KKR=None):
    """
    Get the input for a voronoi calc.
    Wrapper for VoronoiProcess setting structure, code, options, label, description etc.
    """
    # get process for VoronoiCalculation
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # then reuse common inputs setter all options
    if structure is not None:
        # for 'normal' case starting from structure
        builder = get_inputs_common(
            VoronoiCalculation, code, None, structure, options, label, description, params, serial
        )
    else:
        # for parent_KKR feature used to increase lmax which cannot have 'structure' in inputs
        builder = get_inputs_common(
            VoronoiCalculation, code, None, None, options, label, description, params, serial, parent_KKR=parent_KKR
        )

    return builder


def get_inputs_kkrimp(
    code,
    options,
    label='',
    description='',
    parameters=None,
    serial=False,
    imp_info=None,
    host_GF=None,
    imp_pot=None,
    kkrimp_remote=None,
    host_GF_Efshift=None
):
    """
    Get the input for a kkrimp calc.
    Wrapper for KkrimpProcess setting structure, code, options, label, description etc.
    :param code: a valid KKRimpcode installation (e.g. input from Code.get_from_string('codename@computername'))
    TBD
    """

    from aiida_kkr.calculations.kkrimp import KkrimpCalculation

    # then reuse common inputs setter
    builder = get_inputs_common(
        KkrimpCalculation, code, None, None, options, label, description, parameters, serial, imp_info, host_GF,
        imp_pot, kkrimp_remote, host_GF_Efshift
    )

    return builder


def get_inputs_common(
    calculation,
    code,
    remote,
    structure,
    options,
    label,
    description,
    params,
    serial,
    imp_info=None,
    host_GF=None,
    imp_pot=None,
    kkrimp_remote=None,
    host_GF_Efshift=None,
    **kwargs
):
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
        _sched = code.computer.scheduler_type
    else:
        _sched = None
    if params:
        inputs.parameters = params

    if not options:
        options = {}

    if description:
        inputs.metadata.description = description
    else:
        inputs.metadata.description = ''

    if label:
        inputs.metadata.label = label
    else:
        inputs.metadata.label = ''

    if serial:
        # check for old aiida name (e.g. "slurm") and new aiida name ("core.slurm") of the scheduler
        if _sched in ['core.slurm', 'slurm', 'core.pbspro', 'pbspro']:
            # overwrite settings for serial run
            options['withmpi'] = False
            options['resources'] = {'num_machines': 1, 'tot_num_mpiprocs': 1}
        if _sched in ['core.sge', 'sge']:
            options['withmpi'] = False
            options['resources'] = {'parallel_env': 'smpslots', 'tot_num_mpiprocs': 1}
    else:
        # otherwise assume MPI parallelism if not given in input options
        if 'withmpi' not in list(options.keys()):
            options['withmpi'] = True

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

    if host_GF_Efshift is not None:
        inputs.host_Greenfunction_folder_Efshift = host_GF_Efshift

    if imp_pot is not None:
        inputs.impurity_potential = imp_pot

    if kkrimp_remote is not None:
        inputs.parent_calc_folder = kkrimp_remote

    # add additional inputs
    for link_label, node in kwargs.items():
        inputs[link_label] = node

    return inputs


def get_parent_paranode(remote_data):
    """
    Return the input parameter of the parent calculation giving the remote_data node
    """
    inp_calc = remote_data.get_incoming(link_label_filter='remote_folder').first().node
    inp_para = inp_calc.get_incoming(link_label_filter='parameters').first().node
    return inp_para


def generate_inputcard_from_structure(
    parameters,
    structure,
    input_filename,
    parent_calc=None,
    shapes=None,
    isvoronoi=False,
    use_input_alat=False,
    vca_structure=False
):
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
                   aiida_kkr.calculations.Kkrcalculation and shall not be overwritten)
    :param isvoronoi: tell whether or not the parameter set is for a voronoi calculation or kkr calculation (have different lists of mandatory keys)
    :param use_input_alat: True/False, determines whether the input alat value is taken or the new alat is computed from the Bravais vectors


    :note: assumes valid structure and parameters, i.e. for 2D case all necessary
           information has to be given. This is checked with function
           'check_2D_input' called in aiida_kkr.calculations.Kkrcalculation
    """

    from aiida.common.constants import elements as PeriodicTableElements
    from numpy import array
    from masci_tools.io.kkr_params import kkrparams
    from masci_tools.io.common_functions import get_Ang2aBohr, get_alat_from_bravais
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # initialize list of warnings
    warnings = []

    # list of globally used constants
    a_to_bohr = get_Ang2aBohr()

    # Get the connection between coordination number and element symbol
    # maybe do in a different way

    _atomic_numbers = {data['symbol']: num for num, data in PeriodicTableElements.items()}

    # KKR wants units in bohr
    bravais = array(structure.cell) * a_to_bohr
    alat_input = parameters.get_dict().get('ALATBASIS')
    if use_input_alat and alat_input is not None:
        alat = alat_input
        wmess = 'found alat in input parameters, this will trigger scaling of RMAX, GMAX and RCLUSTZ!'
        print(f'WARNING: {wmess}')
        warnings.append(wmess)
    else:
        alat = get_alat_from_bravais(bravais, is3D=structure.pbc[2])
    bravais = bravais / alat

    sites = structure.sites
    naez = len(sites)
    positions = []
    charges = []
    weights = []  # for CPA
    isitelist = []  # counter sites array for CPA
    isite = 0
    for site in sites:
        pos = site.position
        # TODO maybe convert to rel pos and make sure that type is right for script (array or tuple)
        abspos = array(pos) * a_to_bohr / alat  # also in units of alat
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
            if vca_structure and ikind > 0 and not isvoronoi:
                # for VCA case take weighted average (only for KKR code, voronoi code uses zatom of first site for dummy calculation)
                zatom = zatom * wght_last + zatom_tmp * wght
                # also reset weight to 1
                wght = 1.
            else:
                zatom = zatom_tmp
                if vca_structure and isvoronoi:
                    wght = 1.

            wght_last = wght  # for VCA mode

            # make sure that for VCA only averaged position is written (or first for voronoi code)
            if ((
                vca_structure and ((len(sitekind.symbols) == 1) or (not isvoronoi and ikind == 1) or
                                   (isvoronoi and ikind == 0))
            ) or (not vca_structure)):
                charges.append(zatom)
                weights.append(wght)
                isitelist.append(isite)

    weights = array(weights)
    isitelist = array(isitelist)
    charges = array(charges)
    positions = array(positions)

    ######################################
    # Prepare keywords for kkr from input structure

    # get parameter dictionary
    input_dict = parameters.get_dict()

    # remove special keys that are used for special cases but are not part of the KKR parameter set
    for key in _ignored_keys:
        if input_dict.get(key) is not None:
            wmess = f'automatically removing value of key {key}'
            print('WARNING: ' + wmess)
            warnings.append(wmess)
            input_dict.pop(key)

    # get rid of structure related inputs that are overwritten from structure input
    for key in ['BRAVAIS', 'ALATBASIS', 'NAEZ', '<ZATOM>', '<RBASIS>', 'CARTESIAN']:
        if input_dict.get(key) is not None:
            wmess = f'automatically removing value of key {key}'
            print('WARNING: ' + wmess)
            warnings.append(wmess)
            input_dict.pop(key)

    # automatically rescale RMAX, GMAX, RCLUSTZ, RCLUSTXY which are scaled with the lattice constant
    if alat_input is not None:
        if input_dict.get('RMAX') is not None:
            wmess = f'rescale RMAX: {alat_input / alat}'
            print('WARNING: ' + wmess)
            warnings.append(wmess)
            input_dict['RMAX'] = input_dict['RMAX'] * alat_input / alat
        if input_dict.get('GMAX') is not None:
            wmess = f'rescale GMAX: {1 / (alat_input / alat)}'
            print('WARNING: ' + wmess)
            warnings.append(wmess)
            input_dict['GMAX'] = input_dict['GMAX'] * 1 / (alat_input / alat)
        if input_dict.get('RCLUSTZ') is not None:
            wmess = f'rescale RCLUSTZ: {alat_input / alat}'
            print('WARNING: ' + wmess)
            warnings.append(wmess)
            input_dict['RCLUSTZ'] = input_dict['RCLUSTZ'] * alat_input / alat
        if input_dict.get('RCLUSTXY') is not None:
            wmess = f'rescale RCLUSTXY: {alat_input / alat}'
            print('WARNING: ' + wmess)
            warnings.append(wmess)
            input_dict['RCLUSTXY'] = input_dict['RCLUSTXY'] * alat_input / alat

    # empty kkrparams instance (contains formatting info etc.)
    if not isvoronoi:
        params = kkrparams()
    else:
        params = kkrparams(params_type='voronoi')

    # for KKR calculation set EMIN automatically from parent_calc (always in res.emin of voronoi and kkr) if not provided in input node
    if (('EMIN' not in list(input_dict.keys()) or input_dict['EMIN'] is None) and parent_calc is not None):
        wmess = f'Overwriting EMIN with value from parent calculation {parent_calc}'
        print('WARNING: ' + wmess)
        warnings.append(wmess)
        if parent_calc.process_class == VoronoiCalculation:
            emin = parent_calc.outputs.output_parameters.get_dict().get('emin')
        else:
            emin = parent_calc.outputs.output_parameters.get_dict().get('energy_contour_group').get('emin')
        print('Setting emin:', emin, 'is emin None?', emin is None)
        params.set_value('EMIN', emin)

    # overwrite keywords with input parameter
    for key in list(input_dict.keys()):
        params.set_value(key, input_dict[key], silent=True)

    # Write input to file (the parameters that are set here are not allowed to be modfied externally)
    params.set_multiple_values(
        BRAVAIS=bravais, ALATBASIS=alat, NAEZ=naez, ZATOM=charges, RBASIS=positions, CARTESIAN=True
    )
    # for CPA case:
    if len(weights) > naez:
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
    if rbl is not None:
        params.set_value('<RBLEFT>', array(rbl) * a_to_bohr / alat)
    if rbr is not None:
        params.set_value('<RBRIGHT>', array(rbr) * a_to_bohr / alat)
    if zper_l is not None:
        params.set_value('ZPERIODL', array(zper_l) * a_to_bohr / alat)
    if zper_r is not None:
        params.set_value('ZPERIODR', array(zper_r) * a_to_bohr / alat)

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
            return (
                False,
                'Structure.pbc is neither (True, True, True) for bulk nor (True, True, False) for surface calculation!'
            )
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

    if has2Dinfo != is2D:
        if is2D:
            return (
                False,
                '2D info given in parameters but structure is 3D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}'
                .format(is2D, has2Dinfo, set_keys)
            )
        return (
            False,
            '3D info given in parameters but structure is 2D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}'
            .format(is2D, has2Dinfo, set_keys)
        )

    # if everything is ok:
    return (True, '2D consistency check complete')


def structure_from_params(parameters):
    """
    Construct aiida structure out of kkr parameter set (if ALATBASIS, RBASIS, ZATOM etc. are given)

    :param input: parameters, kkrparams object with structure information set (e.g. extracted from read_inputcard function)

    :returns: success, boolean to determine if structure creatoin was successful
    :returns: structure, an aiida StructureData object
    """
    from masci_tools.io.common_functions import get_aBohr2Ang
    from aiida.common.constants import elements as PeriodicTableElements
    from aiida.orm import StructureData
    from masci_tools.io.kkr_params import kkrparams
    from numpy import array

    # check input
    if not isinstance(parameters, kkrparams):
        raise InputValidationError('input parameters needs to be a "kkrparams" instance!')

    # initialize some stuff
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

    # extract atom numbers
    zatom_all = parameters.get_value('<ZATOM>')

    # extract sites with positions, charges/Atom labels, weights
    # positions in units of alat
    pos_all = array(parameters.get_value('<RBASIS>'))
    if len(pos_all.shape) == 1:
        pos_all = array([pos_all])
        zatom_all = [zatom_all]
    if not parameters.get_value('CARTESIAN'):
        # convert from internal to cartesian coordinates
        for isite, tmp_pos in enumerate(pos_all):
            # cell already contains alat factor to convert to Ang. units
            pos_all[isite] = tmp_pos[0]*cell[0] + \
                tmp_pos[1]*cell[1] + tmp_pos[2]*cell[2]
    else:
        pos_all = pos_all * alat * get_aBohr2Ang()  # now positions are in Ang. units

    # convert to list if input contains a single entry only
    if not isinstance(zatom_all, list):
        zatom_all = [zatom_all]
        pos_all = [pos_all]

    # extract weights and sites for CPA calculations
    if natyp == naez:
        weights = [1. for i in range(natyp)]
        sites = list(range(1, natyp + 1))
    else:
        weights = parameters.get_value('<CPA-CONC>')
        sites = parameters.get_value('<SITE>')

    # fill structure from zatom, weights and sites information
    for isite in sites:
        pos = pos_all[sites.index(isite)]
        weight = weights[sites.index(isite)]
        if abs(zatom_all[isite - 1] - int(zatom_all[isite - 1])) > 10**-4:
            # TODO deal with VCA (non-integer zatom)
            print('VCA not implemented yet, stopping here!')
            raise NotImplementedError('VCA functionality not implemented')

        if zatom_all[isite - 1] < 1:
            symbol = 'X'
            struc.append_atom(position=pos, symbols='X', weights=weight)
        else:
            symbol = PeriodicTableElements.get(zatom_all[isite - 1]).get('symbol')
            struc.append_atom(position=pos, symbols=symbol, weights=weight)

    # set correct pbc for 2D case
    if parameters.get_value('INTERFACE'):
        struc.set_pbc((True, True, False))

    # finally return structure
    return is_complete, struc


def vca_check(structure, parameters):
    """

    """
    nsites = 0
    for site in structure.sites:
        sitekind = structure.get_kind(site.kind_name)
        nsites += len(sitekind.symbols)
    # VCA mode if CPAINFO = [-1,-1] first
    try:
        if parameters.get_dict().get('CPAINFO')[0] < 0:
            params_vca_mode = True
        else:
            params_vca_mode = False
    except:
        params_vca_mode = False
    # check if structure supports VCA mode
    vca_structure = False
    if params_vca_mode:
        if nsites > len(structure.sites):
            vca_structure = True

    return vca_structure


def get_username(computer):
    """
    set upload dir (get the remote username and try 5 times if there was a connection error
    """
    import time
    try_trans = 0
    while try_trans < 5:
        try_trans += 1
        try:
            with computer.get_transport() as transport:
                remote_user = transport.whoami()
        except:
            # this means we have some ssh connection error, thus we wait 5 seconds before we try again
            remote_user = None
            time.sleep(5)
        if remote_user is not None:
            break
    # check if username was extracted correctly and raise an error otherwise
    if remote_user is None:
        raise ValueError('Error getting the username from the computer!')

    return remote_user


def get_natyp(structure):
    """Count number of atom types (>NAEZ for CPA) for the structure"""
    counter = 0  # for CPA
    for site in structure.sites:
        sitekind = structure.get_kind(site.kind_name)
        for ikind in range(len(sitekind.symbols)):
            counter += 1
    return counter
