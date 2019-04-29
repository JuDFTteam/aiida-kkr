# -*- coding: utf-8 -*-
"""
Parser for the KKR-impurity Code.
The parser should never fail, but it should catch
all errors and warnings and show them to the user.
"""

from __future__ import absolute_import
from aiida.parsers.parser import Parser
from aiida.orm import Dict
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

        self._ParserVersion = __version__

        #reuse init of base class
        super(KkrimpParser, self).__init__(calc)


    # pylint: disable=protected-access
    def parse(self, **kwargs):
        """
        Parse output data folder, store results in database.

        :param retrieved: a dictionary of retrieved nodes, where
          the key is the link name
        """
        success = False
        node_list = ()

        # Check that the retrieved folder is there
        try:
            out_folder = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = out_folder._repository.list_object_names()

        file_errors = []
        files = {}

        # Parse output files of KKRimp calculation

        # first get path to files and catch errors if files are not present
        # append tupels (error_category, error_message) where error_category is
        # 1: critical error, always leads to failing of calculation
        # 2: warning, is inspected and checked for consistency with read-in
        #    out_dict values (e.g. nspin, newsosol, ...)

        if KkrimpCalculation._DEFAULT_OUTPUT_FILE not in list_of_files:
            msg = "Output file '{}' not found in list of files: {}".format(KkrimpCalculation._DEFAULT_OUTPUT_FILE, list_of_files)
        
    
        if KkrimpCalculation._DEFAULT_OUTPUT_FILE in out_folder.list_object_names():
            outfile = out_folder.open(KkrimpCalculation._DEFAULT_OUTPUT_FILE)
            files['outfile'] = outfile
        else:
            file_errors.append((1,msg))
            outfile = None
            
        fname = KkrimpCalculation._OUTPUT_000
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_log'] = filepath
        else:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_log'] = None
        fname = KkrimpCalculation._OUT_POTENTIAL
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_pot'] = filepath
        else:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_pot'] = None
        fname = KkrimpCalculation._OUT_TIMING_000
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_timing'] = filepath
        else:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_timing'] = None 
        fname = KkrimpCalculation._OUT_ENERGYSP_PER_ATOM
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_enersp_at'] = filepath
        else:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_enersp_at'] = None
        fname = KkrimpCalculation._OUT_ENERGYTOT_PER_ATOM
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_enertot_at'] = filepath
        else:
            file_errors.append((1, "Critical error! file '{}' not found ".format(fname)))
            files['out_enertot_at'] = None
        fname = KkrimpCalculation._KKRFLEX_LLYFAC
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['kkrflex_llyfac'] = filepath
        else:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['kkrflex_llyfac'] = None
        fname = KkrimpCalculation._KKRFLEX_ANGLE
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['kkrflex_angles'] = filepath
        else:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['kkrflex_angles'] = None
        fname = KkrimpCalculation._OUT_MAGNETICMOMENTS
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_spinmoms'] = filepath
        else:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['out_spinmoms'] = None
        fname = KkrimpCalculation._OUT_ORBITALMOMENTS
        if fname in out_folder.list_object_names():
            filepath = out_folder.open(fname)
            files['out_orbmoms'] = filepath
        else:
            file_errors.append((2, "Warning! file '{}' not found ".format(fname)))
            files['out_orbmoms'] = None


        # now parse file output
        out_dict = {'parser_version': self._ParserVersion,
                    'calculation_plugin_version': KkrimpCalculation._CALCULATION_PLUGIN_VERSION}

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
        self.out('output_parameters', Dict(dict=out_dict))

