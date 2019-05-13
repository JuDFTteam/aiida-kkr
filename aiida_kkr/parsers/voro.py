#-*- coding: utf-8 -*-

from __future__ import absolute_import
from aiida.parsers.parser import Parser
from aiida.orm import Dict
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.exceptions import InputValidationError
from masci_tools.io.parsers.voroparser_functions import parse_voronoi_output
import os


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.3"
__contributors__ = ("Jens Broeder", "Philipp Rüßmann")


class VoronoiParser(Parser):
    """
    Parser class for parsing output of voronoi code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of Voronoi_Parser
        """
        # these files should be present after success of voronoi
        self._default_files = {'outfile': VoronoiCalculation._OUTPUT_FILE_NAME,
                               'atominfo': VoronoiCalculation._ATOMINFO,
                               'radii': VoronoiCalculation._RADII}

        self._ParserVersion = __version__

        #reuse init of base class
        super(VoronoiParser, self).__init__(calc)

    # pylint: disable=protected-access
    def parse(self, **kwargs):
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
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = out_folder._repository.list_object_names()

        # we need at least the output file name as defined in calcs.py
        if VoronoiCalculation._OUTPUT_FILE_NAME not in list_of_files:
            self.logger.error("Output file '{}' not found".format(VoronoiCalculation._OUTPUT_FILE_NAME))
            return self.exit_codes.ERROR_NO_OUTPUT_FILE

        #Parse voronoi output, results that are stored in database are in out_dict

        # get path to files and catch errors if files are not present
        file_errors = []
        if VoronoiCalculation._OUT_POTENTIAL_voronoi in out_folder.list_object_names():
            potfile = out_folder.open(VoronoiCalculation._OUT_POTENTIAL_voronoi)
        else:
            # cover case where potfile is overwritten from input to voronoi calculation
            if VoronoiCalculation._POTENTIAL_IN_OVERWRITE in out_folder.list_object_names():
                potfile = out_folder.open(VoronoiCalculation._POTENTIAL_IN_OVERWRITE)
            else:
                file_errors.append("Critical error! Neither potfile {}  not {} "
                                   "was found".format(VoronoiCalculation._OUT_POTENTIAL_voronoi,
                                                      VoronoiCalculation._POTENTIAL_IN_OVERWRITE))
                potfile = 'file_not_found'
        if VoronoiCalculation._OUTPUT_FILE_NAME in out_folder.list_object_names():
            outfile = out_folder.open(VoronoiCalculation._OUTPUT_FILE_NAME)
        else:
            file_errors.append("Critical error! outfile not found {}".format(VoronoiCalculation._OUTPUT_FILE_NAME))
            outfile = 'file_not_found'
        if VoronoiCalculation._ATOMINFO in out_folder.list_object_names():
            atominfo = out_folder.open(VoronoiCalculation._ATOMINFO)
        else:
            file_errors.append("Critical error! atominfo not found {}".format(VoronoiCalculation._ATOMINFO))
            atominfo = 'file_not_found'
        if VoronoiCalculation._RADII in out_folder.list_object_names():
            radii = out_folder.open(VoronoiCalculation._RADII)
        else:
            file_errors.append("Critical error! radii not found {}".format(VoronoiCalculation._RADII))
            radii = 'file_not_found'
        if VoronoiCalculation._INPUT_FILE_NAME in out_folder.list_object_names():
            inputfile = out_folder.open(VoronoiCalculation._INPUT_FILE_NAME)
        else:
            file_errors.append("Critical error! inputfile not found {}".format(VoronoiCalculation._INPUT_FILE_NAME))
            inputfile = 'file_not_found'
        # initialize out_dict and parse output files
        out_dict = {'parser_version': self._ParserVersion}
        out_dict['calculation_plugin_version'] = VoronoiCalculation._CALCULATION_PLUGIN_VERSION
        #TODO add job description, compound name, calculation title
        success, msg_list, out_dict = parse_voronoi_output(out_dict, outfile,
                                                           potfile, atominfo,
                                                           radii, inputfile)
        # add file open errors to parser output of error messages
        for f_err in file_errors:
            msg_list.append(f_err)
        out_dict['parser_errors'] = msg_list

        #create output node and link
        self.out('output_parameters', Dict(dict=out_dict))

        if not success:
            return self.exit_codes.ERROR_VORONOI_PARSING_FAILED
