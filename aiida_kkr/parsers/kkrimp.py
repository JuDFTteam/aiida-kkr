# -*- coding: utf-8 -*-
"""
Parser for the KKR-impurity Code.
The parser should never fail, but it should catch 
all errors and warnings and show them to the user.
"""

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_kkr.calculations.kkrimp import KkrimpCalculation
from aiida.common.exceptions import InputValidationError
from masci_tools.io.parsers.kkrparser_functions import check_error_category
from aiida_kkr.tools.tools_kkrimp import kkrimp_parser_functions


__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = ("Philipp Rüßmann")


class KkrimpParser(Parser):
    """
    Parser class for parsing output of the KKRimp code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of KkrimpParser
        """
    
        # check for valid input
        if not isinstance(calc, KkrimpCalculation):
            raise InputValidationError("Input calc must be a KkrimpCalculation")
        
        self._ParserVersion = __version__

        #reuse init of base class
        super(KkrimpParser, self).__init__(calc)
        

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

        file_errors = []
        files = {}

        # Parse output files of KKRimp calculation
        
        # first get path to files and catch errors if files are not present
        # append tupels (error_category, error_message) where error_category is 
        # 1: critical error, always leads to failing of calculation
        # 2: warning, is inspected and checked for consistency with read-in 
        #    out_dict values (e.g. nspin, newsosol, ...)
        
        # we need at least the output file name as defined in calcs.py
        if self._calc._DEFAULT_OUTPUT_FILE not in list_of_files:
            msg = "Output file '{}' not found in list of files: {}".format(self._calc._DEFAULT_OUTPUT_FILE, list_of_files)
            self.logger.error(msg)
        try:
            filepath = out_folder.get_abs_path(self._calc._DEFAULT_OUTPUT_FILE)
            files['outfile'] = filepath
        except OSError:
            file_errors.append((1, msg))
            files['outfile'] = None
        # collect list of files
        try:
            fname = self._calc._OUTPUT_000
            filepath = out_folder.get_abs_path(fname)
            files['out_log'] = filepath
        except OSError:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_log'] = None
        try:
            fname = self._calc._OUT_POTENTIAL
            filepath = out_folder.get_abs_path(fname)
            files['out_pot'] = filepath
        except OSError:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_pot'] = None
        try:
            fname = self._calc._OUT_TIMING_000
            filepath = out_folder.get_abs_path(fname)
            files['out_timing'] = filepath
        except OSError:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_timing'] = None
        try:
            fname = self._calc._OUT_ENERGYSP_PER_ATOM
            filepath = out_folder.get_abs_path(fname)
            files['out_enersp_at'] = filepath
        except OSError:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_enersp_at'] = None
        try:
            fname = self._calc._OUT_ENERGYTOT_PER_ATOM
            filepath = out_folder.get_abs_path(fname)
            files['out_enertot_at'] = filepath
        except OSError:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_enertot_at'] = None
        try:
            fname = self._calc._KKRFLEX_LLYFAC
            filepath = out_folder.get_abs_path(fname)
            files['kkrflex_llyfac'] = filepath
        except OSError:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['kkrflex_llyfac'] = None
        try:
            fname = self._calc._KKRFLEX_ANGLE
            filepath = out_folder.get_abs_path(fname)
            files['kkrflex_angles'] = filepath
        except OSError:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['kkrflex_angles'] = None
        try:
            fname = self._calc._OUT_MAGNETICMOMENTS
            filepath = out_folder.get_abs_path(fname)
            files['out_spinmoms'] = filepath
        except OSError:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['out_spinmoms'] = None
        try:
            fname = self._calc._OUT_ORBITALMOMENTS
            filepath = out_folder.get_abs_path(fname)
            files['out_orbmoms'] = filepath
        except OSError:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['out_orbmoms'] = None
        
        # now parse file output
        out_dict = {'parser_version': self._ParserVersion, 
                    'calculation_plugin_version': self._calc._CALCULATION_PLUGIN_VERSION}
        
        success, msg_list, out_dict = kkrimp_parser_functions().parse_kkrimp_outputfile(out_dict, files)
        
        out_dict['parser_errors'] = msg_list
         # add file open errors to parser output of error messages
        for (err_cat, f_err) in file_errors: 
            if err_cat == 1:
                msg_list.append(f_err)
            elif check_error_category(err_cat, f_err, out_dict):
                msg_list.append(f_err)
            else:
                if 'parser_warnings' not in list(out_dict.keys()):
                    out_dict['parser_warnings'] = []
                out_dict['parser_warnings'].append(f_err.replace('Error', 'Warning'))
        out_dict['parser_errors'] = msg_list
        
        #create output node and link
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
            
        return success, node_list
    
