# -*- coding: utf-8 -*-
"""
Parser for the KKR Code.
The parser should never fail, but it should catch 
all errors and warnings and show them to the user.
"""

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData


class KkrParser(Parser):
    """
    Parser class for parsing output of KKR code..
    """

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

        try: # TODO: open here or in subroutine?
            with open(out_folder.get_abs_path(
                    self._calc._OUTPUT_FILE_NAME)) as f:
                outtxt = f.readlines()
                result_dict = parse_kkr_outputfile(outtxt)
        except ValueError:
            self.logger.error("Error parsing the output json")
            return success, node_list
        
        out_dict = result_dict
        
        #TODO other keys to add...
        
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        success = True

        return success, node_list

def parse_kkr_outputfile(outfile):
    """
    Parser method for the kkr outfile. It returns a dictionary with results
    """
    # TODO This still needs to be overworked.
    # If we want to store the arrays we have to store them as array data
    # If only certain results are imported they have to be distilled and stored in the output node
    
    from scipy import array
    
    #parse out-file
    #get rms error and other relevant output-data
    results = {}
    #outfi = open(outfile, 'r') # TODO try block
    outtxt = outfile#outfile.readlines()
    try:
        results['rms'] = array([float(outtxt[iline].split('=')[1].replace('D', 'E')) for iline in range(len(outtxt)) if 'average rms-error' in outtxt[iline]])[-1]
    except:
        results['rms'] = None
    try:
        results['charge_neutrality'] = array([float(outtxt[iline].split('=')[1]) for iline in range(len(outtxt)) if 'charge neutrality in unit cell' in outtxt[iline]])[-1]
    except:
        results['charge_neutrality'] = None
    try:
        results['total_magnetic_moment'] = array([float(outtxt[iline].split('=')[1]) for iline in range(len(outtxt)) if 'TOTAL mag. moment in unit cell' in outtxt[iline]])[-1]
    except:
        results['total_magnetic_moment'] = None
    try:
        results['EF'] = array([float(outtxt[iline].split('FERMI')[1].split()[0]) for iline in range(len(outtxt)) if 'E FERMI' in outtxt[iline]])[-1]
    except:
        results['EF'] = None
    try:
        results['DOS_EF'] = array([float(outtxt[iline].split('=')[1]) for iline in range(len(outtxt)) if 'DOS(E_F)' in outtxt[iline]])[-1]
    except:
        results['DOS_EF'] = None
    try:
        results['total_energy'] = array([float(outtxt[iline].split(':')[1]) for iline in range(len(outtxt)) if 'TOTAL ENERGY in ryd.' in outtxt[iline]])[-1]
    except:
        results['total_energy'] = None
        
    return results