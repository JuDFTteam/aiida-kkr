# -*- coding: utf-8 -*-
"""
Input plug-in for a KKR calculation.
"""
from __future__ import print_function, absolute_import
from __future__ import unicode_literals
import os
import numpy as np
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode, load_node, RemoteData, Dict, StructureData, KpointsData, Bool
from .voro import VoronoiCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import InputValidationError, ValidationError
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.exceptions import UniquenessError
from aiida_kkr.tools.common_workfunctions import (
    generate_inputcard_from_structure, check_2Dinput_consistency, update_params_wf, vca_check, kick_out_corestates
)
from masci_tools.io.common_functions import get_alat_from_bravais, get_Ang2aBohr
from aiida_kkr.tools.tools_kkrimp import make_scoef, write_scoef_full_imp_cls
from masci_tools.io.kkr_params import __kkr_default_params__, kkrparams
import six
from six.moves import range

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.11.8'
__contributors__ = ('Jens Broeder', 'Philipp Rüßmann')


class KkrCalculation(CalcJob):
    """
    AiiDA calculation plugin for a KKR calculation.
    """

    # calculation plugin version
    _CALCULATION_PLUGIN_VERSION = __version__

    # Default input and output files
    _DEFAULT_INPUT_FILE = 'inputcard'  # will be shown with inputcat
    _DEFAULT_OUTPUT_FILE = 'out_kkr'  # verdi shell output will be shown with outputcat

    # same as _DEFAULT_OUTPUT_FILE: piped output of kkr execution to this file
    _OUTPUT_FILE_NAME = _DEFAULT_OUTPUT_FILE

    # List of mandatory input files
    _INPUT_FILE_NAME = _DEFAULT_INPUT_FILE
    _POTENTIAL = 'potential'

    # List of optional input files (may be mandatory for some settings in inputcard)
    _SHAPEFUN = 'shapefun'  # mandatory if nonspherical calculation
    _SCOEF = 'scoef'  # mandatory for KKRFLEX calculation and some functionalities
    _NONCO_ANGLES = 'nonco_angle.dat'  # mandatory if noncollinear directions are used that are not (theta, phi)= (0,0) for all atoms
    _NONCO_ANGLES_IMP = 'nonco_angle_imp.dat'  # mandatory for GREENIMP option (scattering code)
    _SHAPEFUN_IMP = 'shapefun_imp'  # mandatory for GREENIMP option (scattering code)
    _POTENTIAL_IMP = 'potential_imp'  # mandatory for GREENIMP option (scattering code)

    # List of output files that should always be present
    _OUT_POTENTIAL = 'out_potential'
    _OUTPUT_0_INIT = 'output.0.txt'
    _OUTPUT_000 = 'output.000.txt'
    _OUTPUT_2 = 'output.2.txt'
    _OUT_TIMING_000 = 'out_timing.000.txt'
    _NONCO_ANGLES_OUT = 'nonco_angle_out.dat'
    _NONCO_ANGLES_ALL_ITER = 'nonco_angle_out_all_iter.dat'

    # special files (some runs)
    # DOS files
    _COMPLEXDOS = 'complex.dos'
    _DOS_ATOM = 'dos.atom%i'
    _LMDOS = 'lmdos.%2i.%i.dat'
    # qdos files
    _QVEC = 'qvec.dat'
    _QDOS_ATOM = 'qdos.%3i.%i.dat'
    _QDOS_ATOM_OLD = 'qdos.%2i.%i.dat'
    _QDOS_SX = 'qdos_sx.%2i.dat'
    _QDOS_SY = 'qdos_sy.%2i.dat'
    _QDOS_SZ = 'qdos_sz.%2i.dat'
    # kkrflex files for impurity calculation
    _KKRFLEX_GREEN = 'kkrflex_green'
    _KKRFLEX_TMAT = 'kkrflex_tmat'
    _KKRFLEX_ATOMINFO = 'kkrflex_atominfo'
    _KKRFLEX_INTERCELL_REF = 'kkrflex_intercell_ref'
    _KKRFLEX_INTERCELL_CMOMS = 'kkrflex_intercell_cmoms'
    _ALL_KKRFLEX_FILES = [
        _KKRFLEX_GREEN,
        _KKRFLEX_TMAT,
        _KKRFLEX_ATOMINFO,
        _KKRFLEX_INTERCELL_REF,
        _KKRFLEX_INTERCELL_CMOMS,
    ]
    # Jij files
    _Jij_ATOM = 'Jij.atom%0.5i'
    _SHELLS_DAT = 'shells.dat'
    # deci-out and decimation
    _DECIFILE = 'decifile'

    # template.product entry point defined in setup.json
    _default_parser = 'kkr.kkrparser'

    # list of keywords that are not allowed to be modified (new calculation
    # starting from structure and voronoi run is needed instead):
    _do_never_modify = [
        'ALATBASIS',
        'BRAVAIS',
        'NAEZ',
        '<RBASIS>',
        'CARTESIAN',
        'INTERFACE',
        '<NLBASIS>',
        '<RBLEFT>',
        'ZPERIODL',
        '<NRBASIS>',
        '<RBRIGHT>',
        'ZPERIODR',
        'KSHAPE',
        '<SHAPE>',
        '<ZATOM>',
        'NATYP',
        '<SITE>',
        '<CPA-CONC>',
        '<KAOEZL>',
        '<KAOEZR>',
    ]
    #TODO implement workfunction to modify structure (e.g. to use VCA)

    # small number used to check for equivalence
    _eps = 10**-12

    @classmethod
    def define(cls, spec):
        """
        Init internal parameters at class load time
        """
        # reuse base class function
        super(KkrCalculation, cls).define(spec)

        # now define input files and parser
        spec.input(
            'metadata.options.parser_name',
            valid_type=six.string_types,
            default=cls._default_parser,
            non_db=True,
        )
        spec.input(
            'metadata.options.input_filename',
            valid_type=six.string_types,
            default=cls._DEFAULT_INPUT_FILE,
            non_db=True
        )
        spec.input(
            'metadata.options.output_filename',
            valid_type=six.string_types,
            default=cls._DEFAULT_OUTPUT_FILE,
            non_db=True
        )

        # define input nodes (optional ones have required=False)
        spec.input(
            'parameters',
            valid_type=Dict,
            required=True,
            help='Use a node that specifies the input parameters',
        )
        spec.input(
            'parent_folder',
            valid_type=RemoteData,
            required=True,
            help="""Use a remote or local repository folder as parent folder
            (also for restarts and similar). It should contain all the  needed
            files for a KKR calc, only edited files should be uploaded from the
            repository."""
        )
        spec.input(
            'impurity_info',
            valid_type=Dict,
            required=False,
            help="""Use a Parameter node that specifies properties for a following
            impurity calculation (e.g. setting of impurity cluster in scoef
            file that is automatically created)."""
        )
        spec.input(
            'kpoints',
            valid_type=KpointsData,
            required=False,
            help="""Use a KpointsData node that specifies the kpoints for which a
            bandstructure (i.e. 'qdos') calculation should be performed."""
        )
        spec.input(
            'initial_noco_angles',
            valid_type=Dict,
            required=False,
            help="""Initial non-collinear angles for the magnetic moments of
            the impurities. These values will be written into the
            `kkrflex_angle` input file of KKRimp.
            The Dict node should be of the form
            initial_noco_angles = Dict(dict={
                'theta': [theta_at1, theta_at2, ..., theta_atN],
                # list theta values in degrees (0..180)
                'phi': [phi_at1, phi_at2, ..., phi_atN],
                # list phi values in degrees (0..360)
                'fix_dir': [True/False at_1, ..., True/False at_N]
                # list of booleans indicating if the direction of the magnetic
                # moment should be fixed or is allowed relax (True means keep the
                # direction of the magnetic moment fixed)
            })
            Note: The length of the theta, phi and fix_dir lists have to be
            equal to the number of atoms.
            """
        )
        spec.input(
            'deciout_parent',
            valid_type=RemoteData,
            required=False,
            help="""KkrCalculation RemoteData folder from deci-out calculation"""
        )

        # define outputs
        spec.output(
            'output_parameters',
            valid_type=Dict,
            required=True,
            help='results of the KKR calculation',
        )
        spec.default_output_node = 'output_parameters'

        # define exit codes, also used in parser
        spec.exit_code(
            301,
            'ERROR_NO_OUTPUT_FILE',
            message='KKR output file not found',
        )
        spec.exit_code(
            302,
            'ERROR_KKR_PARSING_FAILED',
            message='KKR parser retuned an error',
        )
        spec.exit_code(
            303,
            'ERROR_NO_SHAPEFUN_FOUND',
            message='Could not find shapefun from voronoi parent',
        )

    def prepare_for_submission(self, tempfolder):
        """
        Create input files.

            :param tempfolder: aiida.common.folders.Folder subclass where
                the plugin should put all its files.
            :param inputdict: dictionary of the input nodes as they would
                be returned by get_inputs_dict
        """

        has_parent = False
        local_copy_list = []

        # get mandatory input nodes
        parameters = self.inputs.parameters
        code = self.inputs.code
        parent_calc_folder = self.inputs.parent_folder

        # now check for optional nodes

        # for GF writeout
        if 'impurity_info' in self.inputs:
            imp_info = self.inputs.impurity_info
            found_imp_info = True
        else:
            imp_info = None
            found_imp_info = False

        # for qdos funcitonality
        if 'kpoints' in self.inputs:
            kpath = self.inputs.kpoints
            found_kpath = True
        else:
            found_kpath = False

        # extract parent calculation
        parent_calcs = parent_calc_folder.get_incoming(node_class=CalcJobNode)
        n_parents = len(parent_calcs.all_link_labels())
        if n_parents != 1:
            raise UniquenessError(
                'Input RemoteData is child of {} '
                'calculation{}, while it should have a single parent'
                ''.format(n_parents, '' if n_parents == 0 else 's')
            )
            # TODO change to exit code
        if n_parents == 1:
            parent_calc = parent_calcs.first().node
            has_parent = True

        # check if parent is either Voronoi or previous KKR calculation
        #self._check_valid_parent(parent_calc)

        # extract parent input parameter dict for following check
        try:
            parent_inp_dict = parent_calc.inputs.parameters.get_dict()
        except:
            self.logger.error('Failed trying to find input parameter of parent {}'.format(parent_calc))
            raise InputValidationError('No parameter node found of parent calculation.')

        # check if no keys are illegally overwritten (i.e. compare with keys in self._do_never_modify)
        for key in list(parameters.get_dict().keys()):
            value = parameters.get_dict()[key]
            #self.logger.info("Checking {} {}".format(key, value))
            if value is not None:
                if key in self._do_never_modify:
                    oldvalue = parent_inp_dict[key]
                    try:
                        if oldvalue is None and key in __kkr_default_params__:
                            oldvalue = __kkr_default_params__.get(key)
                        if value == oldvalue:
                            values_eqivalent = True
                        else:
                            values_eqivalent = False
                            # check if values match up to certain numerical accuracy
                            if type(value) == float:
                                if abs(value - oldvalue) < self._eps:
                                    values_eqivalent = True
                            elif type(value) == list or type(value) == np.ndarray:
                                tmp_value, tmp_oldvalue = np.array(value).reshape(-1), np.array(oldvalue).reshape(-1)
                                values_eqivalent_tmp = []
                                for ival in range(len(tmp_value)):
                                    if abs(tmp_value[ival] - tmp_oldvalue[ival]) < self._eps:
                                        values_eqivalent_tmp.append(True)
                                    else:
                                        values_eqivalent_tmp.append(False)
                                if all(values_eqivalent_tmp) and len(value) == len(oldvalue):
                                    values_eqivalent = True
                    except:
                        raise InputValidationError(
                            f'Error while trying to compare old and new values with key={key} in do_never_modify list, oldval={oldvalue}; newval={value}'
                        )
                    if not values_eqivalent:
                        self.logger.error(
                            f'You are trying to set keyword {key} = {value} but this is not allowed since the structure would be modified. Please use a suitable workfunction instead.'
                        )
                        raise InputValidationError(
                            f'You are trying to modify a keyword that is not allowed to be changed! (key={key}, oldvalue={oldvalue}, newvalue={value})'
                        )

        #TODO check for remote folder (starting from folder data not implemented yet)
        # if voronoi calc check if folder from db given, or get folder from rep.
        # Parent calc does not has to be on the same computer.
        # so far we copy every thing from local computer ggf if kkr we want to copy remotely

        # get StructureData node from Parent if Voronoi
        structure = None
        self.logger.info('KkrCalculation: Get structure node from voronoi parent')
        try:
            structure, voro_parent = VoronoiCalculation.find_parent_structure(parent_calc)
        except:
            self.logger.error(f'KkrCalculation: Could not get structure from Voronoi parent ({parent_calc}).')
            raise ValidationError(f'Cound not find structure node from parent {parent_calc}')

        # for VCA: check if input structure and parameter node define VCA structure
        vca_structure = vca_check(structure, parameters)

        # check whether or not the alat from the input parameters are used
        # (this enters as a scaling factor for some parameters)
        use_alat_input = parameters.get_dict().get('use_input_alat', False)
        use_alat_input = parameters.get_dict().get('USE_INPUT_ALAT', use_alat_input)

        # prepare scoef file if impurity_info was given
        write_scoef = False
        runopt = parameters.get_dict().get('RUNOPT', None)
        kkrflex_opt = False

        if runopt is not None and 'KKRFLEX' in runopt:
            kkrflex_opt = True
        if kkrflex_opt:
            write_scoef = True
        elif found_imp_info:
            self.logger.info('Found impurity_info in inputs of the calculation, automatically add runopt KKRFLEX')
            write_scoef = True
            runopt = parameters.get_dict().get('RUNOPT', [])
            runopt.append('KKRFLEX')
            parameters = update_params_wf(
                parameters,
                Dict(
                    dict={
                        'RUNOPT': runopt,
                        'nodename': 'update_KKRFLEX',
                        'nodedesc': 'Update Parameter node with KKRFLEX runopt'
                    }
                )
            )

        if found_imp_info and write_scoef:

            imp_info_dict = imp_info.get_dict()

            # find alat input if needed
            if use_alat_input:
                alat_input = parameters.get_dict().get('ALATBASIS', None) / get_Ang2aBohr()
                self.logger.info('alat_input is ' + str(alat_input))
            else:
                self.logger.info('alat_input is None')
                alat_input = None

            # create scoef file
            if 'imp_cls' not in imp_info_dict:
                # this means cluster is found from parameters in imp_info

                # extract cluster settings
                Rcut = imp_info_dict.get('Rcut', None)
                hcut = imp_info_dict.get('hcut', -1.)
                cylinder_orient = imp_info_dict.get('cylinder_orient', [0., 0., 1.])
                ilayer_center = imp_info_dict.get('ilayer_center', 0)
                for i in range(len(cylinder_orient)):
                    try:
                        len(cylinder_orient[i])
                        vec_shape = False
                    except TypeError:
                        vec_shape = True

                # some consistency checks
                if ilayer_center > len(structure.sites) - 1:
                    raise IndexError(
                        f'Index of the reference site is out of range! Possible values: 0 to {len(structure.sites) - 1}.'
                    )
                elif Rcut < 0:
                    raise ValueError('Cutoff radius has to be positive!')
                if not vec_shape or len(cylinder_orient) != 3:
                    raise TypeError(
                        f'Input orientation vector ({cylinder_orient}) has the wrong shape! It needs to be a 3D-vector!'
                    )

                # now write scoef file
                print('Input parameters for make_scoef read in correctly!')
                with tempfolder.open(self._SCOEF, 'w') as scoef_file:
                    make_scoef(
                        structure,
                        Rcut,
                        scoef_file,
                        hcut,
                        cylinder_orient,
                        ilayer_center,
                        alat_input,
                    )

            else:

                # this means the full imp cluster is given in the input
                # TODO add some consistency checks with structure etc.
                print('Write scoef from imp_cls input!', len(imp_info.get_dict().get('imp_cls')))
                with tempfolder.open(self._SCOEF, 'w') as scoef_file:
                    if alat_input is not None:
                        alat = get_alat_from_bravais(np.array(structure.cell), structure.pbc[2])
                        rescale_alat = alat / alat_input
                        self.report(f'INFO: rescaling imp cls due to alat_input: {rescale_alat}')
                    else:
                        rescale_alat = None
                    write_scoef_full_imp_cls(imp_info, scoef_file, rescale_alat)

        elif write_scoef:
            # if we end up here there is a problem with the input
            self.logger.info('Need to write scoef file but no impurity_info given!')
            raise ValidationError('Found RUNOPT KKRFLEX but no impurity_info in inputs')

        # Check for 2D case
        twoDimcheck, msg = check_2Dinput_consistency(structure, parameters)
        if not twoDimcheck:
            raise InputValidationError(msg)

        # set shapes array either from parent voronoi run or read from inputcard in kkrimporter calculation
        if parent_calc.process_label == 'VoronoiCalculation' or parent_calc.process_label == 'KkrCalculation':
            # get shapes array from voronoi parent
            shapes = voro_parent.outputs.output_parameters.get_dict().get('shapes')
        else:
            # extract shapes from input parameters node constructed by kkrimporter calculation
            shapes = voro_parent.inputs.parameters.get_dict().get('<SHAPE>')
        self.logger.info('Extracted shapes: {}'.format(shapes))

        # qdos option, ensure low T, E-contour, qdos run option and write qvec.dat file
        if found_kpath:
            parameters = self._prepare_qdos_calc(parameters, kpath, structure, tempfolder, use_alat_input)

        # write nonco_angle.dat file and adapt RUNOPTS if needed (i.e. add FIXMOM if directions are not relaxed)
        if 'initial_noco_angles' in self.inputs:
            parameters = self._use_initial_noco_angles(parameters, structure, tempfolder)

        # activate decimation mode and copy decifile from deciout parent
        if 'deciout_parent' in self.inputs:
            parameters = self._use_decimation(parameters, tempfolder)

        # Prepare inputcard from Structure and input parameter data
        with tempfolder.open(self._INPUT_FILE_NAME, u'w') as input_file:
            natom, nspin, newsosol, warnings_write_inputcard = generate_inputcard_from_structure(
                parameters,
                structure,
                input_file,
                parent_calc,
                shapes=shapes,
                vca_structure=vca_structure,
                use_input_alat=use_alat_input
            )

        #################
        # Decide what files to copy based on settings to the code (e.g. KKRFLEX option needs scoef)
        if has_parent:
            # copy the right files #TODO check first if file, exists and throw
            # warning, now this will throw an error
            outfolder = parent_calc.outputs.retrieved

            copylist = []
            if parent_calc.process_class == KkrCalculation:
                copylist = [self._OUT_POTENTIAL]
                # TODO ggf copy remotely from remote node if present ...

            elif parent_calc.process_class == VoronoiCalculation:
                copylist = [parent_calc.process_class._SHAPEFUN]
                # copy either overwrite potential or voronoi output potential
                # (voronoi caclualtion retreives only one of the two)
                if parent_calc.process_class._POTENTIAL_IN_OVERWRITE in outfolder.list_object_names():
                    copylist.append(parent_calc.process_class._POTENTIAL_IN_OVERWRITE)
                else:
                    copylist.append(parent_calc.process_class._OUT_POTENTIAL_voronoi)

            #change copylist in case the calculation starts from an imported calculation
            else:  #if parent_calc.process_class == KkrImporterCalculation:
                if self._OUT_POTENTIAL in outfolder.list_object_names():
                    copylist.append(self._OUT_POTENTIAL)
                else:
                    copylist.append(self._POTENTIAL)
                if self._SHAPEFUN in outfolder.list_object_names():
                    copylist.append(self._SHAPEFUN)

            # create local_copy_list from copylist and change some names automatically
            for file1 in copylist:
                # deal with special case that file is written to another name
                if (
                    file1 == 'output.pot' or file1 == self._OUT_POTENTIAL or (
                        parent_calc.process_class == VoronoiCalculation and
                        file1 == parent_calc.process_class._POTENTIAL_IN_OVERWRITE
                    )
                ):
                    filename = self._POTENTIAL
                else:
                    filename = file1
                # now add to copy list
                local_copy_list.append((outfolder.uuid, file1, filename))

                # add shapefun file from voronoi parent if needed
                if self._SHAPEFUN not in copylist:
                    try:
                        struc, voro_parent = VoronoiCalculation.find_parent_structure(parent_calc)
                    except ValueError:
                        return self.exit_codes.ERROR_NO_SHAPEFUN_FOUND
                    # copy shapefun from retrieved of voro calc
                    voro_retrieved = voro_parent.outputs.retrieved
                    local_copy_list.append((voro_retrieved.uuid, VoronoiCalculation._SHAPEFUN, self._SHAPEFUN))

            # check if core state lie within energy contour and take them out if needed
            local_copy_list = self._kick_out_corestates_kkrhost(local_copy_list, tempfolder)

            # for set-ef option (needs to be done AFTER kicking out core states):
            ef_set = parameters.get_dict().get('ef_set', None)
            ef_set = parameters.get_dict().get('EF_SET', ef_set)
            self.report(
                f"efset: {ef_set}  efset_1: {parameters.get_dict().get('ef_set')} efset_2: {parameters.get_dict().get('EF_SET')} params: {parameters.get_dict()}"
            )
            if ef_set is not None:
                local_copy_list = self._set_ef_value_potential(ef_set, local_copy_list, tempfolder)

            # TODO different copy lists, depending on the keywors input
            print('local copy list: {}'.format(local_copy_list))
            self.report('local copy list: {}'.format(local_copy_list))

        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = []

        # TODO retrieve list needs some logic, retrieve certain files,
        # only if certain input keys are specified....
        calcinfo.retrieve_list = [
            self._DEFAULT_OUTPUT_FILE,
            self._INPUT_FILE_NAME,
            self._SCOEF,
            self._NONCO_ANGLES_OUT,
            self._NONCO_ANGLES_ALL_ITER,
            self._OUT_POTENTIAL,
            self._OUTPUT_0_INIT,
            self._OUTPUT_000,
            self._OUTPUT_2,
            self._OUT_TIMING_000,
        ]

        # for special cases add files to retireve list:

        # 1. dos calculation, add *dos* files if NPOL==0
        retrieve_dos_files = False
        if 'NPOL' in list(parameters.get_dict().keys()):
            if parameters.get_dict()['NPOL'] == 0:
                retrieve_dos_files = True
        if 'TESTOPT' in list(parameters.get_dict().keys()):
            testopts = parameters.get_dict()['TESTOPT']
            if testopts is not None:
                stripped_test_opts = [i.strip() for i in testopts]
                if 'DOS' in stripped_test_opts:
                    retrieve_dos_files = True
        if retrieve_dos_files:
            print('adding files for dos output', self._COMPLEXDOS, self._DOS_ATOM, self._LMDOS)
            add_files = [self._COMPLEXDOS]
            for iatom in range(natom):
                add_files.append(self._DOS_ATOM % (iatom + 1))
                for ispin in range(nspin):
                    add_files.append((self._LMDOS % (iatom + 1, ispin + 1)).replace(' ', '0'))
            calcinfo.retrieve_list += add_files

        # 2. KKRFLEX calculation
        retrieve_kkrflex_files = False
        if 'RUNOPT' in list(parameters.get_dict().keys()):
            runopts = parameters.get_dict()['RUNOPT']
            if runopts is not None:
                stripped_run_opts = [i.strip() for i in runopts]
                if 'KKRFLEX' in stripped_run_opts:
                    retrieve_kkrflex_files = True
        if retrieve_kkrflex_files:
            if self.inputs.retrieve_kkrflex.value:
                # retrieve all kkrflex files
                add_files = self._ALL_KKRFLEX_FILES
            else:
                # do not retrieve kkrflex_tmat and kkrflex_green, they are kept on the remote and used from there
                add_files = [self._KKRFLEX_ATOMINFO, self._KKRFLEX_INTERCELL_REF, self._KKRFLEX_INTERCELL_CMOMS]
            print('adding files for KKRFLEX output', add_files)
            calcinfo.retrieve_list += add_files

        # 3. qdos claculation
        retrieve_qdos_files = False
        if 'RUNOPT' in list(parameters.get_dict().keys()):
            runopts = parameters.get_dict()['RUNOPT']
            if runopts is not None:
                stripped_run_opts = [i.strip() for i in runopts]
                if 'qdos' in stripped_run_opts:
                    retrieve_qdos_files = True
        if retrieve_qdos_files:
            print('adding files for qdos output', self._QDOS_ATOM, self._QVEC)
            add_files = [self._QVEC]
            for iatom in range(natom):
                for ispin in range(nspin):
                    add_files.append((self._QDOS_ATOM % (iatom + 1, ispin + 1)).replace(' ', '0'))
                    add_files.append((self._QDOS_ATOM_OLD % (iatom + 1, ispin + 1)).replace(' ', '0')
                                     )  # try to retrieve both old and new version of the files
                # retrieve also qdos_sx,y,z files if written out
                add_files.append((self._QDOS_SX % (iatom + 1)).replace(' ', '0'))
                add_files.append((self._QDOS_SY % (iatom + 1)).replace(' ', '0'))
                add_files.append((self._QDOS_SZ % (iatom + 1)).replace(' ', '0'))
            calcinfo.retrieve_list += add_files

        # 4. Jij calculation
        retrieve_Jij_files = False
        if 'RUNOPT' in list(parameters.get_dict().keys()):
            runopts = parameters.get_dict()['RUNOPT']
            if runopts is not None:
                stripped_run_opts = [i.strip() for i in runopts]
                if 'XCPL' in stripped_run_opts:
                    retrieve_Jij_files = True
        if retrieve_Jij_files:
            add_files = [self._SHELLS_DAT] + [self._Jij_ATOM % iatom for iatom in range(1, natom + 1)]
            print('adding files for Jij output', add_files)
            calcinfo.retrieve_list += add_files

        # 5. deci-out
        retrieve_decifile = False
        if 'RUNOPT' in list(parameters.get_dict().keys()):
            runopts = parameters.get_dict()['RUNOPT']
            if runopts is not None:
                stripped_run_opts = [i.strip() for i in runopts]
                if 'deci-out' in stripped_run_opts:
                    retrieve_decifile = True
        if retrieve_decifile:
            add_files = [self._DECIFILE]
            print('adding files for deci-out', add_files)
            calcinfo.retrieve_list += add_files

        # now set calcinfo and return
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.code_uuid = code.uuid
        codeinfo.stdout_name = self._DEFAULT_OUTPUT_FILE
        calcinfo.codes_info = [codeinfo]

        return calcinfo

    def _set_parent_remotedata(self, remotedata):
        """
        Used to set a parent remotefolder in the restart of fleur.
        """
        if not isinstance(remotedata, RemoteData):
            raise ValueError('remotedata must be a RemoteData')

        # complain if another remotedata is already found
        input_remote = self.get_inputs(node_type=RemoteData)
        if input_remote:
            raise ValidationError('Cannot set several parent calculation to a KKR calculation')

        self.use_parent_folder(remotedata)

    def _set_ef_value_potential(self, ef_set, local_copy_list, tempfolder):
        """
        Set EF value ef_set in the potential file.
        """
        self.report('local copy list before change: {}'.format(local_copy_list))
        self.report("found 'ef_set' in parameters: change EF of potential to this value")

        # first read old potential

        if self._POTENTIAL in tempfolder.get_content_list():
            has_potfile = True
        else:
            has_potfile = False
        self.report(f'has_potfile? {has_potfile}')

        txt = []
        if has_potfile:
            # this is the case if we kicked out core states before
            with tempfolder.open(self._POTENTIAL, 'r') as potfile:
                # read potential
                txt = potfile.readlines()

        if not has_potfile or len(txt) == 0:
            # this is the case when we take the potential from an existing folder
            potcopy_info = [i for i in local_copy_list if i[2] == self._POTENTIAL][0]
            with load_node(potcopy_info[0]).open(potcopy_info[1]) as potfile:
                # remove previous output potential from copy list
                local_copy_list.remove(potcopy_info)
                # read potential
                txt = potfile.readlines()

        self.report(f'len(potfile)? {len(txt)}')

        # now change value of Fermi level in potential text
        potstart = []
        for iline in range(len(txt)):
            line = txt[iline]
            if 'exc:' in line:
                potstart.append(iline)
        for ipotstart in potstart:
            self.report(f'set ef {ef_set} in potential starting in line {ipotstart}')
            tmpline = txt[ipotstart + 3]
            tmpline = tmpline.split()
            newline = '%10.5f%20.14f%20.14f\n' % (float(tmpline[0]), ef_set, float(tmpline[-1]))

            txt[ipotstart + 3] = newline

        # now (over)writing potential file in tempfolder with changed Fermi energy
        with tempfolder.open(self._POTENTIAL, 'w') as pot_new_ef:
            # write new file
            pot_new_ef.writelines(txt)
            # now this directory (tempfolder) contains the updated potential file
            # thus it is not needed to put it in the local copy list anymore

        # return updated local_copy_list
        return local_copy_list

    def _kick_out_corestates_kkrhost(self, local_copy_list, tempfolder):
        """
        Compare value of core states from potential file in local_copy_list with EMIN
        and kick corestate out of potential if they lie inside the energy contour.
        """

        # read EMIN value from inputcard
        params = kkrparams()
        with tempfolder.open(self._INPUT_FILE_NAME) as input_file:
            params.read_keywords_from_inputcard(input_file)
        emin = params.get_value('EMIN')

        # run kick_out_corestates routine to remove core states that lie above emin
        potcopy_info = [i for i in local_copy_list if i[2] == self._POTENTIAL][0]
        with tempfolder.open(self._POTENTIAL, 'w') as potfile_out:
            with load_node(potcopy_info[0]).open(potcopy_info[1]) as potfile_in:
                num_deleted = kick_out_corestates(potfile_in, potfile_out, emin)

        # remove changed potential from local copy list (already in tempfolder without overlapping core states)
        if num_deleted > 0:
            local_copy_list.remove(potcopy_info)
        else:
            # remove temporarily created file
            tempfolder.remove_path(self._POTENTIAL)

        # return updated local_copy_list
        return local_copy_list

    def _prepare_qdos_calc(self, parameters, kpath, structure, tempfolder, use_alat_input):
        """
        prepare a qdos (i.e. bandstructure) calculation, can only be done if k-points are given in input
        Note: this changes some settings in the parameters to ensure a DOS contour and low smearing temperature
        Also the qvec.dat file is written here.
        """
        # check qdos settings
        change_values = []
        runopt = parameters.get_dict().get('RUNOPT')
        if runopt is None:
            runopt = []
        runopt = [i.strip() for i in runopt]
        if 'qdos' not in runopt:
            runopt.append('qdos')
            change_values.append(['RUNOPT', runopt])
        tempr = parameters.get_dict().get('TEMPR')
        if tempr is None or tempr > 100.:
            change_values.append(['TEMPR', 50.])
        N1 = parameters.get_dict().get('NPT1')
        if N1 is None or N1 > 0:
            change_values.append(['NPT1', 0])
        N2 = parameters.get_dict().get('NPT2')
        if N2 is None:
            change_values.append(['NPT2', 100])
        N3 = parameters.get_dict().get('NPT3')
        if N3 is None or N3 > 0.:
            change_values.append(['NPT3', 0])
        NPOL = parameters.get_dict().get('NPOL')
        if NPOL is None or NPOL > 0.:
            change_values.append(['NPOL', 0])
        parameters = _update_params(parameters, change_values)
        # write qvec.dat file
        kpath_array = kpath.get_kpoints(cartesian=True)
        # convert automatically to internal units
        alat = get_alat_from_bravais(np.array(structure.cell), is3D=structure.pbc[2]) * get_Ang2aBohr()
        if use_alat_input:
            alat_input = parameters.get_dict().get('ALATBASIS')
        else:
            alat_input = alat
        kpath_array = kpath_array * (alat_input / alat) / get_Ang2aBohr() / (2 * np.pi / alat)
        # now write file
        qvec = ['%i\n' % len(kpath_array)]
        qvec += ['%e %e %e\n' % (kpt[0], kpt[1], kpt[2]) for kpt in kpath_array]
        with tempfolder.open(self._QVEC, 'w') as qvecfile:
            qvecfile.writelines(qvec)

        return parameters

    def _use_initial_noco_angles(self, parameters, structure, tempfolder):
        """
        Set starting values for non-collinear calculation (writes nonco_angle.dat to tempfolder).
        Adapt FIXMOM runopt according to fix_dir input in initial_noco_angle input node
        """
        self.report('Found `initial_noco_angles` input node, writing nonco_angle.dat file')

        # extract fix_dir flag and set FIXMOM RUNOPT in parameters accordingly
        fix_dir = self.inputs.initial_noco_angles['fix_dir']
        natom = len(structure.sites)
        if len(fix_dir) != natom:
            raise InputValidationError(
                'Error: `fix_dir` list in `initial_noco_angles` input node needs to have the same length as number of atoms!'
            )

        change_values = []
        runopt = parameters.get_dict().get('RUNOPT')
        if runopt is None:
            runopt = []
        runopt = [i.strip() for i in runopt]
        if all(fix_dir) and 'FIXMOM' not in runopt:
            runopt.append('FIXMOM')
            change_values.append(['RUNOPT', runopt])
        elif not all(fix_dir) and 'FIXMOM' in runopt:
            runopt.pop('FIXMOM')
            change_values.append(['RUNOPT', runopt])
        parameters = _update_params(parameters, change_values)

        # extract theta and phi values from input node
        thetas = self.inputs.initial_noco_angles['theta']
        if len(thetas) != natom:
            raise InputValidationError(
                'Error: `theta` list in `initial_noco_angles` input node needs to have the same length as number of atoms!'
            )
        phis = self.inputs.initial_noco_angles['phi']
        if len(phis) != natom:
            raise InputValidationError(
                'Error: `phi` list in `initial_noco_angles` input node needs to have the same length as number of atoms!'
            )

        # now write kkrflex_angle file
        with tempfolder.open(self._NONCO_ANGLES, 'w') as noco_angle_file:
            for iatom in range(natom):
                theta, phi = thetas[iatom], phis[iatom]
                # check consistency
                if theta < 0. or theta > 180.:
                    raise InputValidationError(
                        f'Error: theta value out of range (0..180): iatom={iatom}, theta={theta}'
                    )
                # write line
                noco_angle_file.write(f'   {theta}    {phi}    {fix_dir[iatom]}\n')

        return parameters

    def _use_decimation(self, parameters, tempfolder):
        """
        Activate decimation mode and copy decifile from output of deciout_parent calculation
        """
        self.report('Found `deciout_parent` input node, activae decimation mode')

        # check if deciout parent calculation was proper deci-out calculation
        deciout_parent = self.inputs.deciout_parent
        parent_calcs = deciout_parent.get_incoming(node_class=CalcJobNode)
        n_parents = len(parent_calcs.all_link_labels())
        if n_parents != 1:
            raise UniquenessError(
                'Input RemoteData is child of {} '
                'calculation{}, while it should have a single parent'
                ''.format(n_parents, '' if n_parents == 0 else 's')
            )
            # TODO change to exit code
        parent_calc = parent_calcs.first().node
        deciout_retrieved = parent_calc.outputs.retrieved
        if self._DECIFILE not in deciout_retrieved.list_object_names():
            raise InputValidationError('Error: deciout_parent does not contain decifile!')

        # add 'DECIMATE' flag, decifile and NSTEPS=1
        change_values = []
        runopt = parameters.get_dict().get('RUNOPT')
        if runopt is None:
            runopt = []
        runopt = [i.strip() for i in runopt]
        runopt.append('DECIMATE')
        change_values.append(['RUNOPT', runopt])
        change_values.append(['FILES', [self._POTENTIAL, self._SHAPEFUN]])  # needed to make DECIFILE work
        change_values.append(['DECIFILES', ['vacuum', self._DECIFILE]])  # works only for right continuation for now!
        change_values.append(['NSTEPS', 1])  # decimation works only in one-shot mode
        parameters = _update_params(parameters, change_values)

        # now write kkrflex_angle file
        with deciout_retrieved.open(self._DECIFILE, 'r') as decifile_handle:
            decifile_txt = decifile_handle.readlines()
        with tempfolder.open(self._DECIFILE, 'w') as decifile_handle:
            decifile_handle.writelines(decifile_txt)

        return parameters


def _update_params(parameters, change_values):
    """
    change parameters node from change_values list of key value pairs
    Retrun input parameter node if change_values list is empty
    """
    if change_values != []:
        new_params = {}
        #{'nodename': 'changed_params_qdos', 'nodedesc': 'Changed parameters to mathc qdos mode. Changed values: {}'.format(change_values)}
        for key, val in parameters.get_dict().items():
            new_params[key] = val
        for key, val in change_values:
            new_params[key] = val
        new_params_node = Dict(dict=new_params)
        #parameters = update_params_wf(parameters, new_params_node)
        parameters = new_params_node
    return parameters
