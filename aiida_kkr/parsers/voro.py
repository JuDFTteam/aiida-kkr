# -*- coding: utf-8 -*-

from aiida.parsers.parser import Parser
from aiida.orm import Dict
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.exceptions import InputValidationError, NotExistent
from masci_tools.io.parsers.voroparser_functions import parse_voronoi_output
from aiida_kkr.tools.context import open_context_to_stack, open_files_in_context
import os
from contextlib import ExitStack

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.4.0'
__contributors__ = ('Jens Broeder', 'Philipp Rüßmann')


class VoronoiParser(Parser):
    """
    Parser class for parsing output of voronoi code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of Voronoi_Parser
        """
        # these files should be present after success of voronoi
        self._default_files = {
            'outfile': VoronoiCalculation._OUTPUT_FILE_NAME,
            'atominfo': VoronoiCalculation._ATOMINFO,
            'radii': VoronoiCalculation._RADII
        }

        self._ParserVersion = __version__

        # reuse init of base class
        super(VoronoiParser, self).__init__(calc)

    # pylint: disable=protected-access
    def parse(self, debug=False, **kwargs):
        """
        Parse output data folder, store results in database.

        :param retrieved: a dictionary of retrieved nodes, where
          the key is the link name
        :returns: nothing if everything is fine or an exit code defined in the voronoi calculation class
        """

        success = False
        node_list = ()

        # Get retrieved folders
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = out_folder.list_object_names()

        # we need at least the output file name as defined in calcs.py
        if VoronoiCalculation._OUTPUT_FILE_NAME not in list_of_files:
            self.logger.error(f"Output file '{VoronoiCalculation._OUTPUT_FILE_NAME}' not found")
            return self.exit_codes.ERROR_NO_OUTPUT_FILE

        # Parse voronoi output, results that are stored in database are in out_dict

        # get path to files and catch errors if files are not present
        file_errors = []
        potfile_name = VoronoiCalculation._OUT_POTENTIAL_voronoi
        if potfile_name not in out_folder.list_object_names():
            # cover case where potfile is overwritten from input to voronoi calculation
            if VoronoiCalculation._POTENTIAL_IN_OVERWRITE in out_folder.list_object_names():
                potfile_name = VoronoiCalculation._POTENTIAL_IN_OVERWRITE
            else:
                file_errors.append(
                    'Critical error! Neither potfile {}  not {} '
                    'was found'.format(potfile_name, VoronoiCalculation._POTENTIAL_IN_OVERWRITE)
                )
                potfile_name = None

        outfile_name = VoronoiCalculation._OUTPUT_FILE_NAME
        if outfile_name not in out_folder.list_object_names():
            file_errors.append(f'Critical error! outfile not found {outfile_name}')
            outfile_name = None

        atominfo_name = VoronoiCalculation._ATOMINFO
        if atominfo_name not in out_folder.list_object_names():
            file_errors.append(f'Critical error! atominfo not found {atominfo_name}')
            atominfo_name = None

        radii_name = VoronoiCalculation._RADII
        if radii_name not in out_folder.list_object_names():
            file_errors.append(f'Critical error! radii not found {radii_name}')
            radii_name = None

        inputfile_name = VoronoiCalculation._INPUT_FILE_NAME
        if inputfile_name not in out_folder.list_object_names():
            file_errors.append(f'Critical error! inputfile not found {inputfile_name}')
            inputfile_name = None

        # initialize out_dict and parse output files
        out_dict = {'parser_version': self._ParserVersion}
        out_dict['calculation_plugin_version'] = VoronoiCalculation._CALCULATION_PLUGIN_VERSION
        # TODO add job description, compound name, calculation title

        with ExitStack() as stack:
            # open files only if they exist
            # openend files are added to the context manager stack
            # so that at the end of the parsing all files are closed
            # when the context manager is exited
            (potfile, outfile, atominfo, radii, inputfile) = open_files_in_context(
                stack, out_folder, potfile_name, outfile_name, atominfo_name, radii_name, inputfile_name
            )
            success, msg_list, out_dict = parse_voronoi_output(
                out_dict, outfile, potfile, atominfo, radii, inputfile, debug=debug
            )

        # add file open errors to parser output of error messages
        for f_err in file_errors:
            msg_list.append(f_err)
        out_dict['parser_errors'] = msg_list

        # create output node and link
        self.out('output_parameters', Dict(dict=out_dict))

        if not success:
            return self.exit_codes.ERROR_VORONOI_PARSING_FAILED
