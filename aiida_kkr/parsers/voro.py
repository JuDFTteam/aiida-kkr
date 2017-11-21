#-*- coding: utf-8 -*-

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_kkr.tools.kkrcontrol import check_voronoi_output

class VoronoiParser(Parser):
    """
    Parser class for parsing output of voronoi code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of Fleur_inputgenParser
        """
        # check for valid input
        #if not isinstance(calc, FleurinputgenCalculation):
        #    raise FleurOutputParsingError(
        #        "Input calc must be a FleurInpgenCalculation")

        # these files should be at least present after success of inpgen
        #self._default_files = {calc._OUTPUT_FILE_NAME, calc._INPXML_FILE_NAME}
        #self._other_files = {calc._SHELLOUT_FILE_NAME}

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
        try:
            outfile = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
            potfile = out_folder.get_abs_path(self._calc._OUT_POTENTIAL_voronoi)
            emin = check_voronoi_output(potfile, outfile)
            out_dict = {'EMIN' : emin}
        except:
            self.logger.error("Error parsing output of voronoi: 'EMIN'")
            return success, node_list
            
        
        #TODO: parse version info (out_voronoi)
        #TODO: RMT0/Rout (radii.dat)
        #TODO: Rout/distNN (radii.dat)
        #TODO: Ncls (out_voronoi)
        #TODO: icls-list (out_voronoi)
        #TODO: NclsSites-list (out_voronoi)
        #TODO: atom volumes vs total volume, (out_voronoi)
        #TODO: Jellstart/Genpotstart+potserial (out_voronoi)
        #TODO: <SHAPE> array (atominfo): should always be set to have correct clusters
        


        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        success = True

        return success, node_list
