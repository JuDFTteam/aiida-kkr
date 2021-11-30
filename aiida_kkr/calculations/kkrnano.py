# -*- coding: utf-8 -*-
"""
Input plug-in for a KKRnano calculation.
"""
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode, Dict, RemoteData
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.common.exceptions import UniquenessError
from aiida_kkr.tools import find_parent_structure

__copyright__ = (u'Copyright (c), 2021, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.1'
__contributors__ = ('Markus Struckmann', 'Philipp Rüßmann')


class KKRnanoCalculation(CalcJob):
    """
    AiiDA calculation plugin for a KKRnano calculation
    """

    ####################
    # File names etc.
    ####################
    # calculation plugin version
    _CALCULATION_PLUGIN_VERSION = __version__
    # Default input and output files
    _DEFAULT_INPUT_FILE = 'input.conf'  # will be shown with inputcat
    _DEFAULT_OUTPUT_FILE = 'out'  #'shell output will be shown with outputcat
    # template.product entry point defined in setup.json
    _DEFAULT_PARSER = 'kkr.kkrnanoparser'
    # File names
    _SHAPEFUN = 'shapefun'
    _POTENTIAL = 'potential'
    _OUT_POTENTIAL = 'out_potential'
    _RBASIS = 'rbasis.xyz'

    _DEFAULT_KKRNANO_PARA = {}

    @classmethod
    def define(cls, spec):
        """
        define internals and inputs / outputs of calculation
        """
        # reuse base class (i.e. CalcJob) functions
        super(KKRnanoCalculation, cls).define(spec)
        # now define input files and parser
        spec.inputs['metadata']['options']['parser_name'].default = cls._DEFAULT_PARSER
        spec.inputs['metadata']['options']['input_filename'].default = cls._DEFAULT_INPUT_FILE
        spec.inputs['metadata']['options']['output_filename'].default = cls._DEFAULT_OUTPUT_FILE

        # define input nodes (optional ones have required=False)
        spec.input('parameters', valid_type=Dict, required=False,
                   default= lambda: Dict(dict=cls._DEFAULT_KKRNANO_PARA),
                   help='Dict node that specifies the input parameters for KKRnano (k-point density etc.)')
        spec.input(
            'parent_calc',
            valid_type=RemoteData,
            required=True,
            help='Use a node that specifies a parent KKRnano or voronoi calculation'
        )

        # define outputs
        spec.output('output_parameters', valid_type=Dict, required=True, help='results of the calculation')
        spec.default_output_node = 'output_parameters'

        # define exit codes, also used in parser
        spec.exit_code(301, 'ERROR_NO_OUTPUT_FILE', message='KKRnano output file not found')
        spec.exit_code(302, 'ERROR_PARSING_FAILED', message='KKRnano parser retuned an error')


    def prepare_for_submission(self, tempfolder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param tempfolder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        # Check inputdict
        parameters = self.inputs.parameters
        structure = find_parent_structure(self.inputs.parent_calc)
        code = self.inputs.code

        # Prepare inputcard from Structure and input parameter data
        with tempfolder.open(self._DEFAULT_INPUT_FILE, u'w') as input_file_handle:
            self._write_input_file(input_file_handle, parameters, structure)
        with tempfolder.open(self._RBASIS, u'w') as rbasis_handle:
            self._write_rbasis(rbasis_handle, structure)
        
        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = self._get_local_copy_list(self.inputs.parent_calc)
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = [
            # TODO fill retrieve list with binary output file
        ]

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._DEFAULT_OUTPUT_FILE
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo


    def _write_input_file(self, input_file_handle, parameters, structure):
        """write the input.conf file for KKRnano"""
        input_file_handle.writelines(['TO BE FILLED\n'])

    def _write_rbasis(self, rbasis_handle, structure):
        """write the rbasis.xyz file for KKRnano"""
        rbasis_handle.writelines(['TO BE FILLED\n'])

    def _get_local_copy_list(self, parent_calc_remote):
        """find files to copy from the parent_folder's retrieved to the iniput of a KKRnano calculation"""

        parent_calc = parent_calc_remote.get_incoming(node_class=CalcJobNode).first().node
        retrieved = parent_calc.outputs.retrieved

        local_copy_list = []
        if not self._is_KkrnanoCalc(parent_calc):
            # copy input potential from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._OUT_POTENTIAL_voronoi, self._POTENTIAL)]

            # copy shapefun from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._SHAPEFUN, self._SHAPEFUN)]

        else:  #if parent_calc.process_class == KkrImporterCalculation:
            #TODO add functionality to start from KKRnano parent
            raise NotImplementedError('starting from KKRnano parent not implemented yet.') 

        return local_copy_list


    def _check_valid_parent(self, parent_calc_folder):
        """
        Check that calc is a valid parent for a KKRnano.
        It can be a VoronoiCalculation, KKRnanoCalculation
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

        if ((not self._is_KkrnanoCalc(parent_calc))):
            raise ValueError('Parent calculation must be a KkrCalculation')

        return overwrite_pot, parent_calc

    def _is_KkrnanoCalc(self, calc):
        """
        check if calc contains the file out_potential
        """
        is_KKR = False
        if calc.process_type == 'aiida.calculations:kkr.kkrnano':
            retrieved_node = calc.get_retrieved_node()
            if 'out_potential' in retrieved_node.list_object_names():
                is_KKR = True

        return is_KKR
