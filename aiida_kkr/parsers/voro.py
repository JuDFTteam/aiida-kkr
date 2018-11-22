#-*- coding: utf-8 -*-

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.exceptions import InputValidationError
from masci_tools.io.parsers.voroparser_functions import parse_voronoi_output



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
        # check for valid input
        if not isinstance(calc, VoronoiCalculation):
            raise InputValidationError("Input calc must be a Voronoi Calculation")

        # these files should be present after success of voronoi
        self._default_files = {'outfile': calc._OUTPUT_FILE_NAME, 
                               'atominfo': calc._ATOMINFO, 
                               'radii': calc._RADII}
        
        self._ParserVersion = __version__

        #reuse init of base class
        super(VoronoiParser, self).__init__(calc)

    # pylint: disable=protected-access
    def parse_with_retrieved(self, retrieved):
        """
        Parse output data folder, store results in database.

        :param retrieved: a dictionary of retrieved nodes, where
          the key is the link name
        :returns: a tuple with two values ``(bool, node_list)``, 
          where:

          * ``bool``: variable to tell if the parsing succeeded
          * ``node_list``: list of new nodes to be stored in the db
            (as a list of tuples ``(link_name, node)``)
        """
        
        success = False
        node_list = ()

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error("No retrieved folder found")
            return success, node_list

        # check what is inside the folder
        list_of_files = out_folder.get_folder_list()

        # we need at least the output file name as defined in calcs.py
        if self._calc._OUTPUT_FILE_NAME not in list_of_files:
            self.logger.error("Output file '{}' not found".format(self._calc._OUTPUT_FILE_NAME))
            return success, node_list
        
        #Parse voronoi output, results that are stored in database are in out_dict
        
        # get path to files and catch errors if files are not present
        file_errors = []
        try:
            potfile = out_folder.get_abs_path(self._calc._OUT_POTENTIAL_voronoi)
        except OSError:
            # cover case where potfile is overwritten from input to voronoi calculation
            try:
                potfile = out_folder.get_abs_path(self._calc._POTENTIAL_IN_OVERWRITE)
            except OSError:
                file_errors.append("Critical error! Neither potfile {}  not {} "
                                   "was found".format(self._calc._OUT_POTENTIAL_voronoi, 
                                                      self._calc._POTENTIAL_IN_OVERWRITE))
                potfile = 'file_not_found'
        try:
            outfile = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
        except OSError:
            file_errors.append("Critical error! outfile not found {}".format(self._calc._OUTPUT_FILE_NAME))
            outfile = 'file_not_found'
        try:
            atominfo = out_folder.get_abs_path(self._calc._ATOMINFO)
        except OSError:  
            file_errors.append("Critical error! atominfo not found {}".format(self._calc._ATOMINFO))
            atominfo = 'file_not_found'
        try:
            radii = out_folder.get_abs_path(self._calc._RADII)
        except OSError:  
            file_errors.append("Critical error! radii not found {}".format(self._calc._RADII))
            radii = 'file_not_found'
        try:
            inputfile = out_folder.get_abs_path(self._calc._INPUT_FILE_NAME)
        except OSError:
            file_errors.append("Critical error! inputfile not found {}".format(self._calc._INPUT_FILE_NAME))
            inputfile = 'file_not_found'
        # initialize out_dict and parse output files
        out_dict = {'parser_version': self._ParserVersion}
        out_dict['calculation_plugin_version'] = self._calc._CALCULATION_PLUGIN_VERSION
        #TODO add job description, compound name, calculation title
        success, msg_list, out_dict = parse_voronoi_output(out_dict, outfile, 
                                                           potfile, atominfo, 
                                                           radii, inputfile)
        # add file open errors to parser output of error messages
        for f_err in file_errors: 
            msg_list.append(f_err)
        out_dict['parser_errors'] = msg_list
        
        #create output node and link
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        
        return success, node_list
    
     
