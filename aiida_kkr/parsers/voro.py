#-*- coding: utf-8 -*-

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.exceptions import InputValidationError
from aiida_kkr.parsers.voroparser_functions import parse_voronoi_output



__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
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
            self.logger.error("Output file not found")
            return success, node_list
        
        #Parse voronoi output, results that are stored in database are in out_dict
        out_dict = {'parser_version': self._ParserVersion}
        outfile = 
        potfile = 
        atominfo = 
        radii = 
        success, msg_list, out_dict = parse_voronoi_output(out_dict, outfile, potfile, atominfo, radii)
        
        out_dict['parser_warnings'] = msg_list
        
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        
        return success, node_list
    
     