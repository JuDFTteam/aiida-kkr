# -*- coding: utf-8 -*-
"""
Input plug-in for a voronoi calculation.
"""
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode, Dict, StructureData, RemoteData, SinglefileData
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida_kkr.tools.common_workfunctions import generate_inputcard_from_structure, check_2Dinput_consistency, vca_check
from aiida.common.exceptions import UniquenessError
import os

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.5.3'
__contributors__ = ('Jens Broeder', 'Philipp Rüßmann')


class VoronoiCalculation(CalcJob):
    """
    AiiDA calculation plugin for a voronoi calculation (creation of starting potential and shapefun).
    """

    ####################
    # File names etc.
    ####################
    # calculation plugin version
    _CALCULATION_PLUGIN_VERSION = __version__
    # Default input and output files
    _DEFAULT_INPUT_FILE = 'inputcard'  # will be shown with inputcat
    _DEFAULT_OUTPUT_FILE = 'out_voronoi'  #'shell output will be shown with outputca
    # List of mandatory input files
    _INPUT_FILE_NAME = 'inputcard'
    # List of output files that should always be present
    _OUTPUT_FILE_NAME = 'out_voronoi'
    # template.product entry point defined in setup.json
    _default_parser = 'kkr.voroparser'
    # File names
    _ATOMINFO = 'atominfo.txt'
    _RADII = 'radii.dat'
    _SHAPEFUN = 'shapefun'
    _VERTICES = 'vertices.dat'
    _OUT_POTENTIAL_voronoi = 'output.pot'
    _POTENTIAL_IN_OVERWRITE = 'overwrite_potential'

    @classmethod
    def define(cls, spec):
        """
        define internals and inputs / outputs of calculation
        """
        # reuse base class (i.e. CalcJob) functions
        super(VoronoiCalculation, cls).define(spec)
        # now define input files and parser
        spec.input('metadata.options.parser_name', valid_type=str, default=cls._default_parser, non_db=True)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE, non_db=True)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE, non_db=True)
        # define input nodes (optional ones have required=False)
        spec.input('parameters', valid_type=Dict, help='Use a node that specifies the input parameters')
        spec.input(
            'structure',
            valid_type=StructureData,
            required=False,
            help='Use a node that specifies the input crystal structure'
        )
        spec.input(
            'parent_KKR',
            valid_type=RemoteData,
            required=False,
            help='Use a node that specifies a parent KKR calculation'
        )
        spec.input(
            'potential_overwrite',
            valid_type=SinglefileData,
            required=False,
            help='Use a node that specifies the potential which is used instead of the voronoi output potential'
        )
        spec.input(
            'shapefun_overwrite',
            valid_type=SinglefileData,
            required=False,
            help='Use a node that specifies the shapefun which is used instead of the voronoi output'
        )
        # define outputs
        spec.output('output_parameters', valid_type=Dict, required=True, help='results of the calculation')
        spec.default_output_node = 'output_parameters'
        # define exit codes, also used in parser
        spec.exit_code(301, 'ERROR_NO_OUTPUT_FILE', message='Voronoi output file not found')
        spec.exit_code(302, 'ERROR_VORONOI_PARSING_FAILED', message='Voronoi parser retuned an error')

    def prepare_for_submission(self, tempfolder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param tempfolder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        # Check inputdict
        parameters = self.inputs.parameters

        if 'structure' in self.inputs:
            structure = self.inputs.structure
            found_structure = True
        else:
            found_structure = False

        vca_structure = False
        if found_structure:
            # for VCA: check if input structure and parameter node define VCA structure
            vca_structure = vca_check(structure, parameters)  # pylint: disable=used-before-assignment

        code = self.inputs.code

        # check if a parent folder containing a potential file (out_potential) is given
        if 'parent_KKR' in self.inputs:
            parent_calc_folder = self.inputs.parent_KKR
            found_parent = True
        else:
            found_parent = False

        if found_parent:
            # check if parent is either Voronoi or previous KKR calculation
            overwrite_potential, parent_calc = self._check_valid_parent(parent_calc_folder)  # pylint: disable=possibly-used-before-assignment

            #cross check if no structure was given and extract structure from parent
            if found_structure and not vca_structure:
                raise InputValidationError(
                    'parent_KKR and structure found in input. '
                    'Can only use either parent_KKR or structure in input.'
                )
            else:
                structure_remote_KKR, voro_parent = self.find_parent_structure(parent_calc)
                if not vca_structure:
                    structure = structure_remote_KKR
                else:
                    # check consistency of input vca structure and structure  from remote KKR folder
                    # TODO check consistency
                    pass
        else:
            overwrite_potential = False
            if not found_structure:
                raise InputValidationError('Neither structure nor parent_KKR specified for this '
                                           'calculation')

        # check if overwrite potential is given explicitly
        if 'potential_overwrite' in self.inputs:
            potfile_overwrite = self.inputs.potential_overwrite
            has_potfile_overwrite = True
        else:
            has_potfile_overwrite = False

        if has_potfile_overwrite:
            overwrite_potential = True
            if not found_structure:
                raise InputValidationError(
                    'Input structure needed for this calculation '
                    "(using 'potential_overwrite' input node)"
                )

        ###################################
        # Check for 2D case
        twoDimcheck, msg = check_2Dinput_consistency(structure, parameters)
        if not twoDimcheck:
            raise InputValidationError(msg)

        # Prepare inputcard from Structure and input parameter data
        with tempfolder.open(self._INPUT_FILE_NAME, u'w') as input_file:
            try:
                use_alat_input = parameters.get_dict().get('use_input_alat', False)
                use_alat_input = parameters.get_dict().get('USE_INPUT_ALAT', use_alat_input)
                natom, nspin, newsosol, warnings_write_inputcard = generate_inputcard_from_structure(
                    parameters,
                    structure,
                    input_file,
                    isvoronoi=True,
                    vca_structure=vca_structure,
                    use_input_alat=use_alat_input,
                )
            except ValueError as e:
                raise InputValidationError(f'Input Dict not consistent: {e}')

        # Decide what files to copy
        local_copy_list = []
        if overwrite_potential:
            # copy the right files #TODO check first if file, exists and throw
            # warning, now this will throw an error
            if found_parent and self._is_KkrCalc(parent_calc):
                outfolder = parent_calc.outputs.retrieved  # copy from remote folder
                copylist = [parent_calc.process_class._OUT_POTENTIAL]
            elif has_potfile_overwrite:
                outfolder = potfile_overwrite  # copy from potential sfd  # pylint: disable=possibly-used-before-assignment
                copylist = [potfile_overwrite.filename]
            else:
                copylist = []

            for file1 in copylist:
                filename = file1
                if (found_parent or has_potfile_overwrite) and file1 == copylist[0]:
                    filename = self._POTENTIAL_IN_OVERWRITE
                local_copy_list.append((outfolder.uuid, file1, filename))  # pylint: disable=possibly-used-before-assignment

            # add shapefun to overwrite
            if 'shapefun_overwrite' in self.inputs:
                shapefun_overwrite = self.inputs.shapefun_overwrite
                filename = shapefun_overwrite.filename
                local_copy_list.append((shapefun_overwrite.uuid, filename, 'shapefun_overwrite'))

        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = [
            self._OUTPUT_FILE_NAME,
            self._ATOMINFO,
            self._RADII,
            self._SHAPEFUN,
            self._VERTICES,
            self._INPUT_FILE_NAME,
        ]

        # pass on overwrite potential if this was given in input
        # (KkrCalculation checks if this file is there and takes this file instead of _OUT_POTENTIAL_voronoi
        #  if given)
        if overwrite_potential:
            calcinfo.retrieve_list += [self._POTENTIAL_IN_OVERWRITE]
        else:
            calcinfo.retrieve_list += [self._OUT_POTENTIAL_voronoi]

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo

    def _check_valid_parent(self, parent_calc_folder):
        """
        Check that calc is a valid parent for a FleurCalculation.
        It can be a VoronoiCalculation, KKRCalculation
        """
        overwrite_pot = False

        # extract parent calculation
        parent_calcs = parent_calc_folder.get_incoming(node_class=CalcJobNode)
        n_parents = len(parent_calcs.all_link_labels())
        if n_parents != 1:
            raise UniquenessError(
                'Input RemoteData is child of {} '
                'calculation{}, while it should have a single parent'
                ''.format(n_parents, '' if n_parents == 0 else 's')
            )
        else:
            parent_calc = parent_calcs.first().node
            overwrite_pot = True

        if ((not self._is_KkrCalc(parent_calc))):
            raise ValueError('Parent calculation must be a KkrCalculation')

        return overwrite_pot, parent_calc

    def _is_KkrCalc(self, calc):
        """
        check if calc contains the file out_potential
        """
        is_KKR = False
        if calc.process_type == 'aiida.calculations:kkr.kkr':
            retrieved_node = calc.get_retrieved_node()
            if 'out_potential' in retrieved_node.list_object_names():
                is_KKR = True

        return is_KKR

    @classmethod
    def find_parent_structure(self, parent_folder):
        """
        Find the Structure node recuresively in chain of parent calculations (structure node is input to voronoi calculation)
        This is a copy of the find_parent_structure that moved to tools.find_parent to keep backwards compatibility.
        """
        from aiida_kkr.tools.find_parent import find_parent_structure
        return find_parent_structure(parent_folder)
