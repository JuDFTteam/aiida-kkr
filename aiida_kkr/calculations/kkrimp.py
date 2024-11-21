# -*- coding: utf-8 -*-
"""
Input plug-in for a KKRimp calculation.
"""

from aiida.engine import CalcJob
from aiida.orm import CalcJobNode, Dict, RemoteData, SinglefileData, Bool
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError, UniquenessError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from masci_tools.io.kkr_params import kkrparams
from .voro import VoronoiCalculation
from .kkr import KkrCalculation
from aiida_kkr.tools.tools_kkrimp import modify_potential, make_scoef, write_scoef_full_imp_cls, get_imp_info_from_parent
from aiida_kkr.tools.common_workfunctions import get_username
from aiida_kkr.tools.ldau import get_ldaupot_text
from aiida_kkr.tools.imp_cluster_tools import get_scoef_single_imp
from masci_tools.io.common_functions import search_string, get_ef_from_potfile
import os
import tarfile
import numpy as np

__copyright__ = (u'Copyright (c), 2018, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.10.2'
__contributors__ = (u'Philipp Rüßmann', u'Fabian Bertoldo')

#TODO: implement 'ilayer_center' consistency check


class KkrimpCalculation(CalcJob):
    """
    AiiDA calculation plugin for a KKRimp calculation.
    """

    # calculation plugin version
    _CALCULATION_PLUGIN_VERSION = __version__

    # List of mandatory input files
    _CONFIG = u'config.cfg'
    _POTENTIAL = u'potential'
    _KKRFLEX_GREEN = KkrCalculation._KKRFLEX_GREEN
    _KKRFLEX_TMAT = KkrCalculation._KKRFLEX_TMAT
    _KKRFLEX_ATOMINFO = KkrCalculation._KKRFLEX_ATOMINFO
    _KKRFLEX_INTERCELL_REF = KkrCalculation._KKRFLEX_INTERCELL_REF
    _KKRFLEX_INTERCELL_CMOMS = KkrCalculation._KKRFLEX_INTERCELL_CMOMS

    # List of optional input files (may be mandatory for some setting in inputputcard)
    _SHAPEFUN = u'shapefun'
    _KKRFLEX_ANGLE = u'kkrflex_angle'
    _KKRFLEX_LLYFAC = u'kkrflex_llyfac'
    _KKRFLEX_SOCFAC = u'kkrflex_spinorbitperatom'
    _KKRFLEX_RIMPSHIFT = u'kkrflex_rimpshift'

    # full list of kkrflex files
    _ALL_KKRFLEX_FILES = KkrCalculation._ALL_KKRFLEX_FILES
    _ALL_KKRFLEX_FILES.append(_KKRFLEX_ANGLE)
    _ALL_KKRFLEX_FILES.append(_KKRFLEX_LLYFAC)
    _ALL_KKRFLEX_FILES.append(_KKRFLEX_SOCFAC)

    # List of output files that should always be present (are always retrieved)
    _OUT_POTENTIAL = u'out_potential'
    _OUTPUT_000 = u'out_log.000.txt'
    _OUT_TIMING_000 = u'out_timing.000.txt'
    _OUT_ENERGYSP = u'out_energysp_eV'
    _OUT_ENERGYSP_PER_ATOM = u'out_energysp_per_atom_eV'
    _OUT_ENERGYTOT = u'out_energytotal_eV'
    _OUT_ENERGYTOT_PER_ATOM = u'out_energytotal_per_atom_eV'

    # List of output files that are retrieved if special conditions are fulfilled
    _OUT_JIJMAT = u'out_Jijmatrix'
    _OUT_JIJ_OF_E_BASE = 'out_Jijmatrix_Eres_IE*'
    _OUT_LDOS = u'out_ldos.atom=*'
    _OUT_LDOS_INTERPOL = u'out_ldos.interpol.atom=*'
    _OUT_LMDOS = u'out_lmdos.atom=*'
    _OUT_LMDOS_INTERPOL = u'out_lmdos.interpol.atom=*'
    _OUT_MAGNETICMOMENTS = u'out_magneticmoments'
    _OUT_ORBITALMOMENTS = u'out_orbitalmoments'
    _LDAUPOT = 'ldaupot'

    # template.product entry point defined in setup.json
    _default_parser = u'kkr.kkrimpparser'

    # default names used within aiida (verdi calculation out(in)putcat)
    _OUTPUT_FILE_NAME = u'out_kkrimp'
    _INPUT_FILE_NAME = _CONFIG
    _DEFAULT_OUTPUT_FILE = _OUTPUT_FILE_NAME  # will be shown with outputcat
    _DEFAULT_INPUT_FILE = _INPUT_FILE_NAME  # will be shown with inputcat

    # name of tarfile which is created by parser after successful parsing (to reduce amount of data stored in repo)
    _FILENAME_TAR = 'output_all.tar.gz'
    _DIRNAME_GF_UPLOAD = 'kkrflex_green_upload'

    @classmethod
    def define(cls, spec):
        """
        Init internal parameters at class load time
        """

        # reuse base class function
        super(KkrimpCalculation, cls).define(spec)
        # now define input files and parser
        spec.input('metadata.options.parser_name', valid_type=str, default=cls._default_parser, non_db=True)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE, non_db=True)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE, non_db=True)
        # define input nodes (optional ones have required=False)
        spec.input(
            'parameters',
            valid_type=Dict,
            required=False,
            help='Use a node that specifies the input parameters (calculation settings).'
        )
        spec.input(
            'host_Greenfunction_folder',
            valid_type=RemoteData,
            required=True,
            help=
            'Use a node that specifies the host KKR calculation contaning the host Green function and tmatrix (KkrCalculation with impurity_info input).'
        )
        spec.input(
            'host_Greenfunction_folder_Efshift',
            valid_type=RemoteData,
            required=False,
            help=
            'Use a node that specifies the host KKR calculation contaning the host Green function and tmatrix with Fermi level shift (used to set Fermi level).'
        )
        spec.input(
            'impurity_potential',
            valid_type=SinglefileData,
            required=False,
            help='Use a node that contains the input potential.'
        )
        spec.input(
            'parent_calc_folder',
            valid_type=RemoteData,
            required=False,
            help='Use a node that specifies a parent KKRimp calculation.'
        )
        spec.input(
            'impurity_info',
            valid_type=Dict,
            required=False,
            help='Use a parameter node that specifies properties for a immpurity calculation.'
        )
        spec.input(
            'settings_LDAU',
            valid_type=Dict,
            required=False,
            help="""
Settings for running a LDA+U calculation. The Dict node should be of the form
    settings_LDAU = Dict({'iatom=0':{
        'L': 3,         # l-block which gets U correction (1: p, 2: d, 3: f-electrons)
        'U': 7.,        # U value in eV
        'J': 0.75,      # J value in eV
        'Eref_EF': 0.,  # reference energy in eV relative to the Fermi energy. This is the energy where the projector wavefunctions are calculated (should be close in energy where the states that are shifted lie (e.g. for Eu use the Fermi energy))
    }})

Note: you can add multiple entries like the one for iatom==0 in this example. The atom index refers to the corresponding atom in the impurity cluster.
"""
        )
        spec.input(
            'initial_noco_angles',
            valid_type=Dict,
            required=False,
            help="""
Initial non-collinear angles for the magnetic moments of the impurities. These values will be written into the `kkrflex_angle` input file of KKRimp.
The Dict node should be of the form
    initial_noco_angles = Dict({
        'theta': [theta_at1, theta_at2, ..., theta_atN], # list theta values in degrees (0..180)
        'phi': [phi_at1, phi_at2, ..., phi_atN],         # list phi values in degrees (0..360)
        'fix_dir': [True, False, ..., True/False],       # list of booleans indicating of the direction of the magentic moment should be fixed or is allowed to be updated (True means keep the direction of the magnetic moment fixed)
    })

Note: The length of the theta, phi and fix_dir lists have to be equal to the number of atoms in the impurity cluster.
"""
        )
        spec.input(
            'rimpshift',
            valid_type=Dict,
            required=False,
            help="""
Shift for atoms in the impurity cluster used in U-transformation.

The Dict node should be of the form
    rimpshift = Dict({'shifts': [[0., 0., 0.], ... ]})
where the shifts are given in atomic units (i.e. the internal KKR units).

Note: The length of the 'shifts' attribute should be an array with three numbers indicating the shift for each atom in the impurity cluster.
"""
        )

        # define outputs
        spec.output('output_parameters', valid_type=Dict, required=True, help='results of the KKRimp calculation')
        spec.default_output_node = 'output_parameters'
        # define exit codes, also used in parser
        spec.exit_code(301, 'ERROR_NO_RETRIEVED_FOLDER', message='Retrieved folder of KKRimp calculation not found.')
        spec.exit_code(302, 'ERROR_PARSING_KKRIMPCALC', message='KKRimp parser returned an error.')
        #TBD

    def prepare_for_submission(self, tempfolder):
        """
        Create input files.

            :param tempfolder: aiida.common.folders.Folder subclass where
                the plugin should put all its files.
            :param inputdict: dictionary of the input nodes as they would
                be returned by get_inputs_dict
        """
        # Check inputdict and extrace nodes
        tmp = self._check_and_extract_input_nodes(tempfolder)
        parameters, code, imp_info = tmp[0], tmp[1], tmp[2]
        kkrflex_file_paths, shapefun_path, shapes = tmp[3], tmp[4], tmp[5]
        host_parent_calc, params_host, impurity_potential = tmp[6], tmp[7], tmp[8]
        parent_calc_folder, structure = tmp[9], tmp[10]

        # Prepare input files for KKRimp calculation
        # 1. fill kkr params for KKRimp, write config file and eventually also kkrflex_llyfac file
        self._extract_and_write_config(
            parent_calc_folder, params_host, parameters, tempfolder,
            kkrflex_file_paths[KkrCalculation._KKRFLEX_ATOMINFO]
        )
        # 2. write shapefun from impurity info and host shapefun and copy imp. potential
        self._get_pot_and_shape(
            imp_info, shapefun_path, shapes, impurity_potential, parent_calc_folder, tempfolder, structure
        )
        # 3. change kkrflex_atominfo to match impurity case
        self._change_atominfo(imp_info, kkrflex_file_paths, tempfolder)

        # prepare copy and retrieve lists
        local_copy_list = []  # potential already in tempfolder which is copied in full
        for filename in list(kkrflex_file_paths.keys()):
            if filename != self._KKRFLEX_ATOMINFO:
                src_path = kkrflex_file_paths[filename]
                local_copy_list.append((src_path.uuid, filename, filename))

        retrieve_list = [
            self._OUTPUT_FILE_NAME, self._CONFIG, self._KKRFLEX_ATOMINFO, self._KKRFLEX_ANGLE, self._KKRFLEX_LLYFAC,
            self._KKRFLEX_SOCFAC, self._OUT_POTENTIAL, self._OUTPUT_000, self._OUT_TIMING_000,
            self._OUT_ENERGYSP_PER_ATOM, self._OUT_ENERGYTOT_PER_ATOM
        ]

        # maybe create kkrflex_rimpshift file
        self._write_kkrflex_rimpshift(tempfolder, parameters)

        # extract run and test options (these change retrieve list in some cases)
        allopts = self.get_run_test_opts(parameters)

        # retrieve l(m)dos files
        retrieve_list = self.add_lmdos_files_to_retrieve(tempfolder, allopts, retrieve_list, kkrflex_file_paths)

        # change retrieve list for Chebychev solver
        retrieve_list = self.adapt_retrieve_tmatnew(tempfolder, allopts, retrieve_list)

        # write ldaupot file and change retrieved list for LDA+U calculation
        # this is triggered with having the settings_LDAU dict node in input.
        retrieve_list = self.init_ldau(tempfolder, retrieve_list, parent_calc_folder)

        # add Jij files to retrieve list if calculation if in Jij mode
        retrieve_list = self.add_jij_files(tempfolder, retrieve_list)

        # change local and remote copy list if GF is found on remote machine
        remote_symlink_list, local_copy_list = self.get_remote_symlink(local_copy_list)

        # Prepare CalcInfo to be returned to aiida (e.g. retreive_list etc.)
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = []
        calcinfo.remote_symlink_list = remote_symlink_list
        calcinfo.retrieve_list = retrieve_list

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo

    #####################################################
    # helper functions

    def _get_and_verify_hostfiles(self, tempfolder):
        """
        Check inputdict for host_Greenfunction_folder and extract impurity_info, paths to kkrflex-files and path of shapefun file

        :param inputdict: input dictionary containing all input nodes to KkrimpCalculation
        :returns:
            * imp_info: Dict node containing impurity information like position, Z_imp, cluster size, etc.
            * kkrflex_file_paths: dict of absolute file paths for the kkrflex files
            * shapefun_path: absolute path of the shapefunction file in the host calculation (needed to construct shapefun_imp)
            * shapes: mapping array of atoms to shapes (<SHAPE> input)
        :note: shapefun_path is None if host_Greenfunction calculation was not full-potential
        :raises:
            * InputValidationError, if inputdict does not contain 'host_Greenfunction'
            * InputValidationError, if host_Greenfunction_folder not of right type
            * UniquenessError, if host_Greenfunction_folder does not have exactly one parent
            * InputValidationError, if host_Greenfunction does not have an input node impurity_info
            * InputValidationError, if host_Greenfunction was not a KKRFLEX calculation
        """

        # get mandatory input nodes (extract host_Greenfunction_folder)
        host_parent = self.inputs.host_Greenfunction_folder

        # extract parent calculation
        parent_calcs = host_parent.get_incoming(node_class=CalcJob)
        n_parents = len(parent_calcs.all_link_labels())
        if n_parents != 1:
            raise UniquenessError(
                'Input RemoteData is child of {} '
                'calculation{}, while it should have a single parent'
                ''.format(n_parents, '' if n_parents == 0 else 's')
            )
        else:
            parent_calc = parent_calcs.first().node
        # extract impurity_info
        found_impurity_inputnode = False
        found_host_parent = False
        if 'impurity_info' in self.inputs:
            imp_info_inputnode = self.inputs.impurity_info
            if not isinstance(imp_info_inputnode, Dict):
                raise InputValidationError('impurity_info not of type Dict')
            found_impurity_inputnode = True
            found_host_parent = True
        imp_info = get_imp_info_from_parent(parent_calc)
        if imp_info is None:
            raise InputValidationError('host_Greenfunction calculation does not have an input node impurity_info')

        # if impurity input is seperate input, check if it is the same as
        # the one from the parent calc (except for 'Zimp'). If that's not the
        # case, raise an error
        if found_impurity_inputnode and found_host_parent:
            check_consistency_imp_info = False
            if 'imp_cls' in imp_info_inputnode.get_dict().keys():
                input_imp_cls_arr = np.array(imp_info_inputnode.get_dict()['imp_cls'])

                try:
                    parent_imp_cls_arr = np.array(imp_info.get_dict()['imp_cls'])
                except:
                    host_structure, _ = VoronoiCalculation.find_parent_structure(parent_calc)
                    parent_imp_cls_arr = get_scoef_single_imp(host_structure, imp_info)

                is_identical = np.array_equal(
                    np.round(input_imp_cls_arr[:, 0:4], 5), np.round(parent_imp_cls_arr[:, 0:4], 5)
                )

                if is_identical:
                    check_consistency_imp_info = True
                else:
                    self.report('impurity_info node from input and from previous GF calculation are NOT compatible!.')
                    self.report(f'impurity_info node from input: {input_imp_cls_arr[:, 0:4]}')
                    self.report(f'impurity_info node from previous GF calculation: {parent_imp_cls_arr[:, 0:4]}')

            #TODO: implement also 'ilayer_center' check
            elif imp_info_inputnode.get_dict().get('Rcut') == imp_info.get_dict().get('Rcut'):

                check_consistency_imp_info = True
                try:
                    if (
                        imp_info_inputnode.get_dict().get('hcut') == imp_info.get_dict().get('hcut') and
                        imp_info_inputnode.get_dict().get('cylinder_orient')
                        == imp_info.get_dict().get('cylinder_orient') and
                        imp_info_inputnode.get_dict().get('Rimp_rel') == imp_info.get_dict().get('Rimp_rel')
                    ):

                        self.report('impurity_info node from input and from previous GF calculation are compatible')
                        check_consistency_imp_info = True
                    else:
                        self.report(
                            'impurity_info node from input and from previous GF calculation are NOT compatible!. '
                            'Please check your impurity_info nodes for consistency.'
                        )
                        check_consistency_imp_info = False
                except AttributeError:
                    self.report(
                        'Non default values of the impurity_info node from input and from previous '
                        "GF calculation are compatible. Default values haven't been checked"
                    )
                    check_consistency_imp_info = True
            else:
                self.report(
                    'impurity_info node from input and from previous GF calculation are NOT compatible!. '
                    'Please check your impurity_info nodes for consistency.'
                )
                check_consistency_imp_info = False
            if check_consistency_imp_info:
                imp_info = imp_info_inputnode
            else:
                raise InputValidationError('impurity_info nodes (input and GF calc) are not compatible')

        # check if host parent was KKRFLEX calculation
        hostfolder = parent_calc.outputs.retrieved
        with hostfolder.open(KkrCalculation._DEFAULT_INPUT_FILE) as fhandle:
            # use read_keywords_from_inputcard of kkrparams class
            params_host_calc = kkrparams(params_type='kkr')
            params_host_calc.read_keywords_from_inputcard(inputcard=fhandle)

        if 'RUNOPT' not in list(params_host_calc.get_dict().keys()):
            host_ok = False
        elif 'KKRFLEX' not in params_host_calc.get_dict().get('RUNOPT', []):
            host_ok = False
        else:
            host_ok = True

        if not host_ok:
            raise InputValidationError('host_Greenfunction calculation was not a KKRFLEX run')

        # extract information from Efshift host GF input node (not mandatory)
        if 'host_Greenfunction_folder_Efshift' in self.inputs:
            host_parent_Efshift = self.inputs.host_Greenfunction_folder_Efshift
            parent_calcs_Efshift = host_parent_Efshift.get_incoming(node_class=CalcJob)
            parent_calc_Efshift = parent_calcs_Efshift.first().node
            hostfolder_Efshift = parent_calc_Efshift.outputs.retrieved
        else:
            hostfolder_Efshift = None

        kkrflex_file_paths = {}
        for filename in self._ALL_KKRFLEX_FILES:
            if filename in hostfolder.list_object_names():
                kkrflex_file_paths[filename] = hostfolder
            # take tmat and green file from Fermi level overwrite directory (second GF_writeout calculation)
            if hostfolder_Efshift is not None and filename in [self._KKRFLEX_TMAT, self._KKRFLEX_GREEN]:
                if filename in hostfolder_Efshift.list_object_names():
                    kkrflex_file_paths[filename] = hostfolder_Efshift

        # extract shapes array from parameters read from inputcard
        shapes = params_host_calc.get_dict().get('<SHAPE>', None)
        if shapes is None:
            # fallback if SHAPES is not explicityl set (assume each atom has it's own shapefun
            shapes = range(1, params_host_calc.get_dict().get('NAEZ') + 1)
        elif type(shapes) == int:
            shapes = [shapes]

        # extract input structure and voro_parent to get shapefun in next step
        try:
            structure, voro_parent = VoronoiCalculation.find_parent_structure(parent_calc)
        except:
            raise InputValidationError('No structure node found from host GF parent')

        # extract shapefun path for read-in
        shapefun_path = {}
        if VoronoiCalculation._SHAPEFUN in voro_parent.outputs.retrieved.list_object_names():
            shapefun_path = voro_parent.outputs.retrieved
        else:
            shapefun_path = None

        return imp_info, kkrflex_file_paths, shapefun_path, shapes, parent_calc, params_host_calc, structure

    def _check_and_extract_input_nodes(self, tempfolder):
        """
        Extract input nodes from inputdict and check consitency of input nodes
        :param inputdict: dict of inputnodes
        :returns:
            * parameters (aiida_kkr.tools.kkr_params.kkrparams), optional: parameters of KKRimp that end up in config.cfg
            * code (KKRimpCodeNode): code of KKRimp on some machine
            * imp_info (DictNode): parameter node of the impurity information, extracted from host_parent_calc
            * kkrflex_file_paths (dict): dictionary of {filenames: absolute_path_to_file} for the kkrflex-files
            * shapfun_path (str): absolute path of the shapefunction of the host parent calculation
            * host_parent_calc (KkrCalculation): node of the parent host calculation where the kkrflex-files were created
            * impurity_potential (SinglefileData): single file data node containing the starting potential for the impurity calculation
            * parent_calc_folder (RemoteData): remote directory of a parent KKRimp calculation
        """

        # get mandatory input nodes (extract code)
        code = self.inputs.code

        # now check for optional nodes
        if 'parameters' in self.inputs:
            parameters = self.inputs.parameters
        else:
            parameters = None
        if parameters is not None:  # convert to kkrparams instance
            parameters = kkrparams(params_type='kkrimp', **parameters.get_dict())

        # get hostfiles
        imp_info, kkrflex_file_paths, shapfun_path, shapes, host_parent_calc, params_host, structure = self._get_and_verify_hostfiles(
            tempfolder
        )

        # check impurity potential or parent calculation input
        # impurity_potential
        if 'impurity_potential' in self.inputs:
            impurity_potential = self.inputs.impurity_potential
            found_imp_pot = True
        else:
            impurity_potential = None
            found_imp_pot = False
        # parent calculation folder
        if 'parent_calc_folder' in self.inputs:
            parent_calc_folder = self.inputs.parent_calc_folder
            found_parent_calc = True
        else:
            parent_calc_folder = None
            found_parent_calc = False
        # consistency checks
        if not found_parent_calc and not found_imp_pot:
            raise InputValidationError(
                'Neither impurity_potential nor parent_calc_folder specified for this calculation.\n'
                'Please provide either impurity_potential or parent_calc_folder.'
            )
        elif found_parent_calc and found_imp_pot:
            raise InputValidationError(
                'Both impurity_potential and parent_calc_folder specified for this calculation.\n'
                'Please provide one one, i.e. either impurity_potential or parent_calc_folder.'
            )

        # Done checking inputs, returning...
        return parameters, code, imp_info, kkrflex_file_paths, shapfun_path, shapes, host_parent_calc, params_host, impurity_potential, parent_calc_folder, structure

    def _extract_and_write_config(self, parent_calc_folder, params_host, parameters, tempfolder, GFhost_folder):
        """
        fill kkr params for KKRimp and write config file
        also writes kkrflex_llyfac file if Lloyd is used in the host system
        """

        runflag, params_kkrimp = self._initialize_kkrimp_params(params_host, parameters, GFhost_folder, tempfolder)

        # overwrite keys if found in parent_calc (previous KKRimp calculation)
        # here `parent_calc_folder` is the `remote` output node of the previous KKRimp calculation
        if parent_calc_folder is not None:
            parent_calc = parent_calc_folder.get_incoming().get_node_by_label('remote_folder')
            params_parent = parent_calc.get_incoming().get_node_by_label('parameters')
        else:
            params_parent = None
        if params_parent is not None:
            params_parent = kkrparams(params_type='kkrimp', **params_parent.get_dict())
            for (key, val) in params_parent.get_set_values():
                self._check_key_setting_consistency(params_kkrimp, key, val)
                params_kkrimp.set_value(key, val)

        # finally overwrite from input parameters
        if parameters is not None:
            for (key, val) in parameters.get_set_values():
                self._check_key_setting_consistency(params_kkrimp, key, val)
                params_kkrimp.set_value(key, val)

        # special run mode: calculation of Jijs
        if parameters.get_value('CALCJIJMAT') is not None and parameters.get_value('CALCJIJMAT') == 1:
            self.report('Found CALCJIJMAT=1: trigger JIJ mode which overwrites IMIX, MIXFAC, SCFSTEPS and RUNFLAGs')
            runflag = self._activate_jij_calc(runflag, params_kkrimp, GFhost_folder, tempfolder)

        # write kkrflex_angle file
        # DOES NOT WORK TOGETHER WITH JIJ mode!!!
        if 'initial_noco_angles' in self.inputs:
            self.report('Found `initial_noco_angles` input node, writing kkrflex_angle file')
            self._write_kkrflex_angle(parameters, GFhost_folder, tempfolder)

        # write config.cfg
        with tempfolder.open(self._CONFIG, u'w') as config_file:
            params_kkrimp.fill_keywords_to_inputfile(output=config_file)

    def _change_atominfo(self, imp_info, kkrflex_file_paths, tempfolder):
        """
        change kkrflex_atominfo to match impurity case
        """
        # read in atominfo file as it is written out
        with kkrflex_file_paths[KkrCalculation._KKRFLEX_ATOMINFO].open(KkrCalculation._KKRFLEX_ATOMINFO) as file:
            atominfo = file.readlines()

        #TODO implement logic to extract this info from imp_info
        replace_zatom_imp = []

        # read scoef for comparison with Rimp_rel
        scoef = []
        with tempfolder.open(KkrCalculation._SCOEF, u'r') as scoeffile:
            n_rows = len(scoeffile.readlines()) - 1
        with tempfolder.open(KkrCalculation._SCOEF, u'r') as scoeffile:
            if n_rows > 1:
                scoef = np.loadtxt(scoeffile, skiprows=1)[:, :3]
            else:
                scoef = np.loadtxt(scoeffile, skiprows=1)[:3]

        # find replaceZimp list from Zimp and Rimp_rel
        imp_info_dict = imp_info.get_dict()
        Zimp_list = imp_info_dict.get(u'Zimp')
        if type(Zimp_list) != list:
            Zimp_list = [Zimp_list]  # fast fix for cases when Zimp is not a list but a single value
        Rimp_rel_list = imp_info_dict.get(u'Rimp_rel', [[0, 0, 0]])
        self.report(f'DEBUG: Rimp_rel_list: {Rimp_rel_list}.')

        for iatom in range(len(Zimp_list)):
            rtmp = np.array(Rimp_rel_list[iatom])[:3]
            self.report(f'INFO: Rimp_rel {iatom}, {rtmp}')
            if n_rows > 1:
                diff = np.sqrt(np.sum((rtmp - scoef)**2, axis=1))
            else:
                diff = np.sqrt(np.sum((rtmp - scoef)**2, axis=0))
            Zimp = Zimp_list[iatom]
            ipos_replace = np.where(diff == diff.min())[0][0]
            replace_zatom_imp.append([ipos_replace, Zimp])

        for (iatom, zimp) in replace_zatom_imp:
            tmp = atominfo[iatom + 4].split()
            x, y, z = float(tmp[0]), float(tmp[1]), float(tmp[2])
            zatom = float(tmp[3])
            virt, remove, lmax = int(tmp[4]), int(tmp[5]), int(tmp[6])
            zatom = zimp
            tmp = u' %24.16f %24.16f %24.16f %5.2f %4i %4i %4i\n' % (x, y, z, zatom, virt, remove, lmax)
            atominfo[iatom + 4] = tmp

        # write atominfo file
        with tempfolder.open(self._KKRFLEX_ATOMINFO, u'w') as atominfofile:
            atominfofile.writelines(atominfo)

    def _get_pot_and_shape(
        self, imp_info, shapefun, shapes, impurity_potential, parent_calc_folder, tempfolder, structure
    ):
        """
        write shapefun from impurity info and host shapefun and copy imp. potential

        returns: file handle to potential file
        """

        imp_info_dict = imp_info.get_dict()

        # create scoef file
        if 'imp_cls' not in imp_info_dict:
            # this means cluster is found from parameters in imp_info

            # extract cluster settings
            Rcut = imp_info_dict.get('Rcut', None)
            hcut = imp_info_dict.get('hcut', -1.)
            cylinder_orient = imp_info_dict.get('cylinder_orient', [0., 0., 1.])
            ilayer_center = imp_info_dict.get('ilayer_center', 0)

            # now write scoef file
            self.report('Input parameters for make_scoef read in correctly!')
            with tempfolder.open(KkrCalculation._SCOEF, 'w') as scoef_file:
                make_scoef(structure, Rcut, scoef_file, hcut, cylinder_orient, ilayer_center)

        else:
            # this means the full imp cluster is given in the input

            self.report(f"Write scoef from imp_cls input! {len(imp_info_dict.get('imp_cls'))}")
            with tempfolder.open(KkrCalculation._SCOEF, 'w') as scoef_file:
                write_scoef_full_imp_cls(imp_info, scoef_file)

        # now create impurity shapefun (reads scoef file)
        with tempfolder.open(KkrCalculation._SCOEF, u'r') as scoef_file:
            with tempfolder.open(KkrimpCalculation._SHAPEFUN, u'w') as shapefun_new:
                with shapefun.open(KkrimpCalculation._SHAPEFUN, u'r') as shapefun_file:
                    shapelen = len(shapefun_file.readlines())
                    if shapelen > 1:
                        modify_potential().shapefun_from_scoef(scoef_file, shapefun_file, shapes, shapefun_new)

        # get path of tempfolder
        with tempfolder.open('.dummy', 'w') as tmpfile:
            tempfolder_path = os.path.dirname(tmpfile.name)

        # find path to input potential
        if impurity_potential is not None:
            potfile_name, potfile_folder = impurity_potential.filename, impurity_potential
        elif parent_calc_folder is not None:
            self.report(
                f'parent_calc_folder {parent_calc_folder} {parent_calc_folder.get_incoming().all_link_labels()}'
            )
            retrieved = parent_calc_folder.get_incoming(node_class=CalcJobNode
                                                        ).first().node.get_outgoing().get_node_by_label('retrieved')
            self.report(f'potfile {retrieved} {self._OUT_POTENTIAL}')

            # extract file from host's tarball (extract to tempfolder and use from there, this way the unnessesary files are deleted once submission is done)
            tar_filenames = []
            if self._FILENAME_TAR in retrieved.list_object_names():
                # get path of tarfile
                with retrieved.open(self._FILENAME_TAR) as tf:
                    tfpath = tf.name
                # extract file from tarfile of retrieved to tempfolder
                with tarfile.open(tfpath) as tf:
                    tar_filenames = [ifile.name for ifile in tf.getmembers()]
                    self.report(f'extract potfile from tarball {tar_filenames}')
                    filename = self._OUT_POTENTIAL
                    if filename in tar_filenames:
                        tf.extract(filename, tempfolder_path)  # extract to tempfolder
            else:  # otherwise copy from retrieved to tempfolder (rest of calculation needs files to be in tempfolder)
                filename = self._OUT_POTENTIAL
                if filename in retrieved.list_object_names():
                    with tempfolder.open(filename, u'w') as newfile:
                        with retrieved.open(filename, u'r') as oldfile:
                            newfile.writelines(oldfile.readlines())
                    self.report('copied potfile from retrieved to tempfolder')

            # now out_potential is in tempfolder (either copied or extracted) and can be copied from there
            potfile_name, potfile_folder = self._OUT_POTENTIAL, tempfolder
        else:
            raise InputValidationError(
                'ERROR in _get_pot_and_shape: neither impurity potential nor parent_calc_folder given!'
            )

        # copy potential to correct input name ('potential')
        with potfile_folder.open(potfile_name, u'r') as oldfile:
            with tempfolder.open(self._POTENTIAL, u'w') as newfile:
                newfile.writelines(oldfile.readlines())

        # remove old potential file if still present (saves unnecessary copying)
        if self._OUT_POTENTIAL in os.listdir(tempfolder_path):
            os.remove(os.path.join(tempfolder_path, self._OUT_POTENTIAL))
        if '.dummy' in os.listdir(tempfolder_path):
            os.remove(os.path.join(tempfolder_path, '.dummy'))

    def _get_natom(self, tempfolder):
        """Get the number of atoms in the impurity cluster from kkrflex_atominfo file"""
        with tempfolder.open(self._KKRFLEX_ATOMINFO) as file:
            atominfo = file.readlines()
        itmp = search_string('NATOM', atominfo)
        if itmp >= 0:
            natom = int(atominfo[itmp + 1].split()[0])
        else:
            raise ValueError('Could not extract NATOM value from kkrflex_atominfo')
        return natom

    def _initialize_kkrimp_params(self, params_host, parameters, GFhost_folder, tempfolder):
        """Initialize KKRimp parameters and set keys that are the same as in the host calculation"""

        # initialize kkrimp parameter set with default values
        params_kkrimp = kkrparams(
            params_type='kkrimp',
            NPAN_LOGPANELFAC=2,
            RADIUS_MIN=-1,
            NCOLL=0,
            SPINORBIT=0,
            SCFSTEPS=1,
            IMIX=0,
            MIXFAC=0.05,
            ITDBRY=20,
            BRYMIX=0.05,
            QBOUND=10**-7,
            RUNFLAG=[],
            TESTFLAG=[],
            HFIELD=[0.0, 0],
            CALCFORCE=0,
            CALCJIJMAT=0,
            CALCORBITALMOMENT=0,
            ICST=2
        )

        # keys that are being overwritten from host calculation settings
        keys_overwrite = [
            'NSPIN', 'KVREL', 'XC', 'INS', 'ICST', 'RADIUS_LOGPANELS', 'NPAN_EQ', 'NPAN_LOG', 'NCHEB', 'QBOUND'
        ]
        for key in keys_overwrite:
            if key == 'XC':
                key0 = 'KEXCOR'
            elif key == 'RADIUS_LOGPANELS':
                key0 = 'R_LOG'
            elif key == 'MIXFAC':
                key0 = 'STRMIX'
            else:
                key0 = key
            val = params_host.get_value(key0)
            if val is not None:
                params_kkrimp.set_value(key, val)

        # settings for SOC solver
        use_cheby = False
        runopts = params_host.get_value('RUNOPT')
        if runopts is not None and 'NEWSOSOL' in runopts:
            use_cheby = True
        if (
            params_host.get_value('<USE_CHEBYCHEV_SOLVER>') is not None and
            params_host.get_value('<USE_CHEBYCHEV_SOLVER>')
        ):
            use_cheby = True
        if use_cheby:
            params_kkrimp.set_multiple_values(NCOLL=1, SPINORBIT=1, CALCORBITALMOMENT=1, TESTFLAG=['tmatnew'])
        else:
            params_kkrimp.set_multiple_values(NCOLL=0, SPINORBIT=0, CALCORBITALMOMENT=0, TESTFLAG=[])

        # SOC solver but with SOC strength scaled to zero
        if use_cheby:
            # maybe deactivate SOC in KKRimp calculation
            self._set_nosoc(params_host, GFhost_folder, tempfolder)

        # extract input RUNFLAGS
        runflag = None
        if parameters is not None:
            for (key, val) in parameters.get_set_values():
                if key == 'RUNFLAG':
                    runflag = list(val)
        if runflag is None:
            runflag = []

        # special settings
        runopts = params_host.get_value('RUNOPT')
        if 'SIMULASA' in runopts or (params_kkrimp.get_value('NCOLL') > 0 and params_kkrimp.get_value('INS') == 0):
            runflag.append('SIMULASA')
        # take care of LLYsimple (i.e. Lloyd in host system)
        if 'LLOYD' in runopts:
            runflag = self._use_lloyd(runflag, GFhost_folder, tempfolder)

        # now set runflags
        params_kkrimp.set_value('RUNFLAG', runflag)

        return runflag, params_kkrimp

    def _set_nosoc(self, params_host, GFhost_folder, tempfolder):
        """Check if host is a noSOC calculation and then set the kkrflex_spinorbitperatom accordingly"""

        # check if noSOC mode was used in the host run
        nosoc = False
        if (params_host.get_value('<SET_CHEBY_NOSOC>') is not None and params_host.get_value('<SET_CHEBY_NOSOC>')):
            nosoc = True
        if params_host.get_value('<SOCSCL>') is not None:
            # TODO make this more flexible, now only SOC on/off is possible
            socscale = params_host.get_value('<SOCSCL>')
            try:
                for scl in socscale:
                    if scl < 1e-14:
                        nosoc = True
            except TypeError:
                if socscale < 1e-14:
                    nosoc = True

        # maybe deactivate SOC by writing the kkrflex_spinorbitperatom file
        if nosoc:
            # extract NATOM from atominfo file
            natom = self._get_natom(GFhost_folder)

            # set soc scale factor (should be integer at the moment!)
            socscale = [0 for i in range(natom)]

            # now write kkrflex_spinorbitperatom file
            with tempfolder.open(self._KKRFLEX_SOCFAC, 'w') as kkrflex_socfac:
                for socscl in socscale:
                    kkrflex_socfac.write(f' {socscl}\n')

    def _use_lloyd(self, runflag, GFhost_folder, tempfolder):
        """Use the LLYsimple version of KKRimp code with the average renormalization factor from the host calculation"""
        # add runflag for imp code
        runflag.append('LLYsimple')
        # also extract renormalization factor and create kkrflex_llyfac file (contains one value only)
        with GFhost_folder.open('output.000.txt') as f:
            txt = f.readlines()
            iline = search_string('RENORM_LLY: Renormalization factor of total charge', txt)
            if iline >= 0:
                llyfac = txt[iline].split()[-1]
                # now write kkrflex_llyfac to tempfolder where later on config file is also written
                with tempfolder.open(self._KKRFLEX_LLYFAC, 'w') as f2:
                    f2.writelines([llyfac])
        return runflag

    def _activate_jij_calc(self, runflag, params_kkrimp, GFhost_folder, tempfolder):
        """Adapt runoptions to use Jij calculation and set inputs for Jij run"""

        # settings in config file
        runflag.append('force_angles')
        # take care of LDA+U
        if 'settings_LDAU' in self.inputs:
            # this prevents mixing LDAU potential in between iterations
            runflag.append('freezeldau')

        # now add runflags
        params_kkrimp.set_multiple_values(IMIX=0, MIXFAC=0., SCFSTEPS=3, RUNFLAG=runflag)

        # for DOS mode add flag to writeout Jij info energy resolved (this will be retrieved if )
        testflag = params_kkrimp.get_value('TESTFLAG')
        testflag.append('Jij(E)')
        params_kkrimp.set_value('TESTFLAG', testflag)

        # extract NATOM from atominfo file
        natom = self._get_natom(GFhost_folder)

        # now write kkrflex_angle file
        with tempfolder.open(self._KKRFLEX_ANGLE, 'w') as kkrflex_angle_file:
            for istep in range(3):
                for iatom in range(natom):
                    if istep == 0:
                        kkrflex_angle_file.write(f'   0.0    0.0    1\n')
                    elif istep == 1:
                        kkrflex_angle_file.write(f'  90.0    0.0    1\n')
                    else:
                        kkrflex_angle_file.write(f'  90.0   90.0    1\n')

        return runflag

    def _write_kkrflex_angle(self, parameters, GFhost_folder, tempfolder):
        """Create the kkrflex_angle file in tempfolder"""

        # check if calculation is no Jij run
        if parameters.get_value('CALCJIJMAT') is not None and parameters.get_value('CALCJIJMAT') == 1:
            raise InputValidationError('ERROR: angles cannot be set if Jij mode is chosen!')

        # extract NATOM from atominfo file
        natom = self._get_natom(GFhost_folder)

        # extract values from input node
        thetas = self.inputs.initial_noco_angles['theta']
        if len(thetas) != natom:
            raise InputValidationError(
                'Error: `theta` list in `initial_noco_angles` input node needs to have the same length as number of atoms in the impurity cluster!'
            )
        phis = self.inputs.initial_noco_angles['phi']
        if len(phis) != natom:
            raise InputValidationError(
                'Error: `phi` list in `initial_noco_angles` input node needs to have the same length as number of atoms in the impurity cluster!'
            )
        fix_dirs = self.inputs.initial_noco_angles['fix_dir']
        if len(fix_dirs) != natom:
            raise InputValidationError(
                'Error: `fix_dir` list in `initial_noco_angles` input node needs to have the same length as number of atoms in the impurity cluster!'
            )

        # now write kkrflex_angle file
        with tempfolder.open(self._KKRFLEX_ANGLE, 'w') as kkrflex_angle_file:
            for iatom in range(natom):
                theta, phi, fix_dir = thetas[iatom], phis[iatom], fix_dirs[iatom]
                # check consistency
                if theta < 0 or theta > 180:
                    raise InputValidationError(
                        f'Error: theta value out of range (0..180): iatom={iatom}, theta={theta}'
                    )
                if phi < 0 or phi > 360:
                    raise InputValidationError(f'Error: phi value out of range (0..360): iatom={iatom}, phi={phi}')
                # convert fix_dir to integer if given as boolean
                if isinstance(fix_dir, bool):
                    fix_dir = (1 if fix_dir else 0)
                # write line
                kkrflex_angle_file.write(f'   {theta}    {phi}    {fix_dir}\n')

    def _write_kkrflex_rimpshift(self, tempfolder, parameters):
        """Create the kkrflex_rimpshift file in tempfolder for U-transformation"""

        if 'rimpshift' in self.inputs:
            # check if calculation is no Jij run
            if parameters.get_value('LATTICE_RELAX') is not None and parameters.get_value('LATTICE_RELAX') != 1:
                raise InputValidationError('ERROR: "LATTICE_RELAX" in the input parameters needs to be set to 1.')

            # extract NATOM from atominfo file
            natom = self._get_natom(tempfolder)

            # extract values from input node
            rimpshift = self.inputs.rimpshift['shifts']
            if len(rimpshift) != natom:
                raise InputValidationError(
                    'Error: `shifts` list in `rimpshift` input node needs to have the same length as number of atoms in the impurity cluster!'
                )

            # now write kkrflex_rimpshift file
            with tempfolder.open(self._KKRFLEX_RIMPSHIFT, 'w') as kkrflex_rimpshift_file:
                for iatom in range(natom):
                    shift = rimpshift[iatom]
                    # write line
                    kkrflex_rimpshift_file.write(f'   {shift[0]}    {shift[1]}    {shift[2]}\n')

    def _check_key_setting_consistency(self, params_kkrimp, key, val):
        """
        Check if key/value pair that is supposed to be set is not in conflict with previous settings of parameters in params_kkrimp
        """
        param_ok = True

        #TODO implement checks

        if not param_ok:
            raise ValueError(f'Trying to set key "{key}" with value "{val}" which is in conflict to previous settings!')

    def get_run_test_opts(self, parameters):
        """Extract run and test options from input parameters"""
        runopts = parameters.get_dict().get('RUNFLAG')
        if runopts is None:
            runopts = []
        testopts = parameters.get_dict().get('TESTFLAG')
        if testopts is None:
            testopts = []
        allopts = runopts + testopts
        return allopts

    def add_lmdos_files_to_retrieve(self, tempfolder, allopts, retrieve_list, kkrflex_file_paths):
        """Add DOS files to retrieve list"""

        if 'lmdos' in allopts or 'ldos' in allopts:
            # check if mode is Jij
            with tempfolder.open(self._CONFIG) as file:
                config = file.readlines()
            itmp = search_string('CALCJIJMAT', config)
            if itmp >= 0:
                calcjijmat = int(config[itmp].split()[-1])
            else:
                raise ValueError('Could not extract CALCJIJMAT value from config.cfg')

            # add DOS output files accordingly (names have '*' ending to catch all files for atoms
            retrieve_list.append(self._OUT_LDOS)
            retrieve_list.append(self._OUT_LDOS_INTERPOL)
            retrieve_list.append(self._OUT_LMDOS)
            retrieve_list.append(self._OUT_LMDOS_INTERPOL)
            # add Jij of E file if Jij mode
            if calcjijmat > 0:
                with kkrflex_file_paths[self._KKRFLEX_TMAT].open(self._KKRFLEX_TMAT, 'r') as f:
                    txt = []
                    for iline in range(3):
                        txt.append(f.readline())
                    nepts = int(txt[1].split()[3])
                for ie in range(1, nepts + 1):
                    retrieve_list.append(self._OUT_JIJ_OF_E)  # energy resolved values

        return retrieve_list

    def adapt_retrieve_tmatnew(self, tempfolder, allopts, retrieve_list):
        """Add out_magneticmoments and orbitalmoments files to retrieve list"""

        # extract NSPIN value
        with tempfolder.open(self._CONFIG) as file:
            config = file.readlines()
        itmp = search_string('NSPIN', config)
        if itmp >= 0:
            nspin = int(config[itmp].split()[-1])
        else:
            raise ValueError('Could not extract NSPIN value from config.cfg')

        # change retrieve list
        if 'tmatnew' in allopts and nspin > 1:
            retrieve_list.append(self._OUT_MAGNETICMOMENTS)
            with tempfolder.open(self._CONFIG) as file:
                outorb = file.readlines()
            itmp = search_string('CALCORBITALMOMENT', outorb)
            if itmp >= 0:
                calcorb = int(outorb[itmp].split()[-1])
            else:
                calcorb = 0
            if calcorb == 1:
                retrieve_list.append(self._OUT_ORBITALMOMENTS)

        return retrieve_list

    def init_ldau(self, tempfolder, retrieve_list, parent_calc_folder):
        """
        Check if settings_LDAU is in input and set up LDA+U calculation. Reuse old ldaupot of parent_folder contains a file ldaupot.
        """

        # first check if settings_LDAU is in inputs
        if 'settings_LDAU' not in self.inputs:
            # do nothing
            return retrieve_list
        else:
            # this means we need to set up LDA+U

            #TODO: check consistency of settings_LDAU

            # add ldaupot to retrieve and local copy lists
            retrieve_list.append(self._LDAUPOT)

            # add runopt LDA+U
            with tempfolder.open(self._CONFIG) as file:
                config = file.readlines()
            itmp = search_string('RUNFLAG', config)
            if itmp >= 0:
                config[itmp] = config[itmp].replace('\n', ' LDA+U \n')
                # overwrite config.cfg
                with tempfolder.open(self._CONFIG, 'w') as file:
                    file.writelines(config)
            else:
                raise ValueError('Could not find RUNFLAG in config.cfg to add LDA+U option')

            # now create ldaupot file
            self.create_or_update_ldaupot(parent_calc_folder, tempfolder)

        return retrieve_list

    def create_or_update_ldaupot(self, parent_calc_folder, tempfolder):
        """
        Writes ldaupot to tempfolder.

        If parent_calc_folder is found and it contains an onld ldaupot, we reuse the values for wldau, uldau and phi from there.
        """

        # extract Fermi energy from potential file
        with tempfolder.open(self._POTENTIAL) as potfile:
            ef_Ry = get_ef_from_potfile(potfile)

        # extract NATOM from atominfo file
        natom = self._get_natom(tempfolder)

        # get old ldaupot file
        reuse_old_ldaupot = self.get_old_ldaupot(parent_calc_folder, tempfolder)

        # settings dict (defines U, J etc.)
        ldau_settings = self.inputs.settings_LDAU.get_dict()

        if reuse_old_ldaupot:
            # reuse wldau, ildau and phi from old ldaupot file
            # Attention the first number needs to be non-zero
            txt = get_ldaupot_text(ldau_settings, ef_Ry, natom, initialize=False)
            # now we read the old file
            with tempfolder.open(self._LDAUPOT + '_old', 'r') as ldaupot_file:
                txt0 = ldaupot_file.readlines()
                # find start of wldau etc.
                ii = 0
                for line in txt0:
                    if 'wldau' in line:
                        break
                    ii += 1
                txt0 = txt0[ii:]

            #remove last line (is replace from txt0)
            txt.pop(-1)
            # put new header and old bottom together
            newtxt = txt + ['\n'] + txt0
        else:
            # initialize ldaupot file
            # Attention: here the first number needs to be 0 which triggers generating initial values in KKRimp
            newtxt = get_ldaupot_text(ldau_settings, ef_Ry, natom, initialize=True)

        # now write to file
        with tempfolder.open(self._LDAUPOT, 'w') as out_filehandle:
            out_filehandle.writelines(newtxt)

    def get_old_ldaupot(self, parent_calc_folder, tempfolder):
        """
        Copy old ldaupot from retrieved of parent or extract from tarball.
        If no parent_calc_folder is present this step is skipped.
        """

        has_ldaupot = False

        if parent_calc_folder is not None:
            retrieved = parent_calc_folder.get_incoming(node_class=CalcJobNode
                                                        ).first().node.get_outgoing().get_node_by_label('retrieved')
            # extract file from host's tarball (extract to tempfolder and use from there, this way the unnessesary files are deleted once submission is done)
            has_ldaupot = self.get_ldaupot_from_retrieved(retrieved, tempfolder)

            return has_ldaupot

    def get_remote_symlink(self, local_copy_list):
        """Check if host GF is found on remote machine and reuse from there"""
        remote_symlink_list = []

        # extract remote computer information
        code = self.inputs.code
        comp = code.computer
        # set upload dir (get the remote username and try 5 times if there was a connection error
        remote_user = get_username(comp)
        workdir = comp.get_workdir().format(username=remote_user)
        GFpath_remote = os.path.join(workdir, self._DIRNAME_GF_UPLOAD)

        self.report(f'local_copy_list: {local_copy_list}')

        # extract GF information from retrieved folder of host GF calc
        uuid_GF_calc = local_copy_list[0][0]
        GF_local_copy_info = [i for i in local_copy_list if i[1] == self._KKRFLEX_GREEN]
        TM_local_copy_info = [i for i in local_copy_list if i[1] == self._KKRFLEX_TMAT]
        if len(TM_local_copy_info) > 0:
            TM_local_copy_info = TM_local_copy_info[0]
        else:
            TM_local_copy_info = None
        if len(GF_local_copy_info) > 0:
            GF_local_copy_info = GF_local_copy_info[0]
        else:
            GF_local_copy_info = None

        # open transport to remote computer
        with comp.get_transport() as connection:
            # do this for GMAT if it was found in retreived folder, otherwise we must assume to find it on the remote
            filename = self._KKRFLEX_GREEN
            GF_remote_path = os.path.join(GFpath_remote, uuid_GF_calc, filename)
            # check if file exists on remote
            if connection.isfile(GF_remote_path):
                # remove GF from local copy list and add to remote symlink list
                if GF_local_copy_info is not None:
                    local_copy_list.remove(GF_local_copy_info)
                remote_symlink_list.append((comp.uuid, GF_remote_path, filename))

            # do the same for TMAT
            filename = self._KKRFLEX_TMAT
            TM_remote_path = os.path.join(GFpath_remote, uuid_GF_calc, filename)
            # check if file exists on remote
            if connection.isfile(TM_remote_path):
                # remove TMAT from local copy list and add to remote symlink list
                if TM_local_copy_info is not None:
                    local_copy_list.remove(TM_local_copy_info)
                remote_symlink_list.append((comp.uuid, TM_remote_path, filename))

        # print symlink and local copy list (for debugging purposes)
        self.report(f'local_copy_list: {local_copy_list}')
        self.report(f'symlink_list: {remote_symlink_list}')

        # now return updated remote_symlink and local_copy lists
        return remote_symlink_list, local_copy_list

    @classmethod
    def get_ldaupot_from_retrieved(self, retrieved, tempfolder):
        """
        Extract ldaupot from output of KKRimp retreived to tempfolder.
        The extracted file in tempfolder will be named ldaupot_old.

        returns True of ldaupot was found, otherwise returns False
        """
        tar_filenames = []
        has_ldaupot = False
        if self._FILENAME_TAR in retrieved.list_object_names():
            # get path of tempfolder
            with tempfolder.open('.dummy', 'w') as tmpfile:
                tempfolder_path = os.path.dirname(tmpfile.name)

            # get path of tarfile
            with retrieved.open(self._FILENAME_TAR) as tf:
                tfpath = tf.name

            # extract file from tarfile of retrieved to tempfolder
            with tarfile.open(tfpath) as tf:
                tar_filenames = [ifile.name for ifile in tf.getmembers()]
                if self._LDAUPOT in tar_filenames:
                    has_ldaupot = True
                    tf.extract(self._LDAUPOT, tempfolder_path)  # extract to tempfolder
            # copy to new filename
            if has_ldaupot:
                with tempfolder.open(self._LDAUPOT + '_old', u'w') as newfile:
                    with tempfolder.open(self._LDAUPOT, u'r') as oldfile:
                        newfile.writelines(oldfile.readlines())

            # remove dummy file
            if '.dummy' in os.listdir(tempfolder_path):
                os.remove(os.path.join(tempfolder_path, '.dummy'))

        else:  # otherwise copy from retrieved to tempfolder (rest of calculation needs files to be in tempfolder)
            if self._LDAUPOT in retrieved.list_object_names():
                has_ldaupot = True
                with tempfolder.open(self._LDAUPOT + '_old', u'w') as newfile:
                    with retrieved.open(self._LDAUPOT, u'r') as oldfile:
                        newfile.writelines(oldfile.readlines())

        return has_ldaupot

    def add_jij_files(self, tempfolder, retrieve_list):
        """
        check if KkrimpCalculation is in Jij mode and add OUT_JIJMAT to retrieve list if needed
        """

        # check if Jij mode is set in config file and add OUT_JIJMAT to retrieved list if needed
        with tempfolder.open(self._CONFIG) as file:
            config = file.readlines()
            itmp = search_string('CALCJIJMAT', config)
            if itmp >= 0:
                calcjijmat = int(config[itmp].split()[-1])
                if calcjijmat > 0:
                    retrieve_list.append(self._OUT_JIJMAT)
            else:
                raise ValueError('Could not extract CALCJIJMAT value from config.cfg')

        return retrieve_list
