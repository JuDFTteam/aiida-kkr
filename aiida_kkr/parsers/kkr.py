# -*- coding: utf-8 -*-
"""
Parser for the KKR Code.
The parser should never fail, but it should catch
all errors and warnings and show them to the user.
"""

from aiida import __version__ as aiida_core_version
from aiida.parsers.parser import Parser
from aiida.orm import Dict
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida.common.exceptions import InputValidationError, NotExistent
from masci_tools.io.parsers.kkrparser_functions import parse_kkr_outputfile, check_error_category
from masci_tools.io.common_functions import search_string
from aiida_kkr.tools.context import open_files_in_context
from contextlib import ExitStack

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.8.0'
__contributors__ = ('Jens Broeder', u'Philipp Rüßmann')


class KkrParser(Parser):
    """
    Parser class for parsing output of KKR code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of KkrParser
        """

        # needed for KKRimporter parser
        self.icrit = 0

        self._ParserVersion = __version__

        # reuse init of base class
        super(KkrParser, self).__init__(calc)

    # pylint: disable=protected-access

    def parse(self, debug=False, **kwargs):
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

        # Get retrieved folders
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = out_folder.list_object_names()

        # we need at least the output file name as defined in calcs.py
        if KkrCalculation._DEFAULT_OUTPUT_FILE not in list_of_files:
            msg = f"Output file '{KkrCalculation._DEFAULT_OUTPUT_FILE}' not found in list of files: {list_of_files}"
            if self.icrit == 0:  # this check turns this off for the KKRimporter calculation
                self.logger.error(msg)
                return self.exit_codes.ERROR_NO_OUTPUT_FILE

        # determine whether or not everything should be parsed or not (e.g. qdos option)
        skip_mode = False
        only_000_present = False
        with out_folder.open(KkrCalculation._INPUT_FILE_NAME) as file:
            txt = file.readlines()
            itmp = search_string('RUNOPT', txt)
            if itmp >= 0:
                runopts = txt[itmp + 1].lower()
                if 'qdos' in runopts:
                    skip_mode = True
                if 'kkrflex' in runopts:
                    only_000_present = True

        # now collect the rest of the files
        file_errors = []

        # Parse output files of KKR calculation
        outfile_name = KkrCalculation._DEFAULT_OUTPUT_FILE
        if outfile_name not in out_folder.list_object_names():
            file_errors.append((1 + self.icrit, msg))
            outfile_name = None

        # get path to files and catch errors if files are not present
        # append tupels (error_category, error_message) where error_category is
        # 1: critical error, always leads to failing of calculation
        # 2: warning, is inspected and checked for consistency with read-in
        #    out_dict values (e.g. nspin, newsosol, ...)
        outfile_0init_name = KkrCalculation._OUTPUT_0_INIT
        if outfile_0init_name not in out_folder.list_object_names():
            file_errors.append((1 + self.icrit, f'Critical error! OUTPUT_0_INIT not found {outfile_0init_name}'))
            outfile_0init_name = None
        outfile_000_name = KkrCalculation._OUTPUT_000
        if outfile_000_name not in out_folder.list_object_names():
            file_errors.append((1 + self.icrit, f'Critical error! OUTPUT_000 not found {outfile_000_name}'))
            outfile_000_name = None
        outfile_2_name = KkrCalculation._OUTPUT_2
        if outfile_2_name not in out_folder.list_object_names():
            if not only_000_present:
                file_errors.append((1 + self.icrit, f'Critical error! OUTPUT_2 not found {outfile_2_name}'))
                outfile_2_name = None
            else:
                outfile_2_name = outfile_000_name
        potfile_out_name = KkrCalculation._OUT_POTENTIAL
        if potfile_out_name not in out_folder.list_object_names():
            file_errors.append((1 + self.icrit, f'Critical error! OUT_POTENTIAL not found {potfile_out_name}'))
            potfile_out_name = None
        timing_file_name = KkrCalculation._OUT_TIMING_000
        if timing_file_name not in out_folder.list_object_names():
            file_errors.append((1 + self.icrit, f'Critical error! OUT_TIMING_000  not found {timing_file_name}'))
            timing_file_name = None
        nonco_out_file_name = KkrCalculation._NONCO_ANGLES_OUT
        if nonco_out_file_name not in out_folder.list_object_names():
            file_errors.append((2, f'Error! NONCO_ANGLES_OUT not found {nonco_out_file_name}'))
            nonco_out_file_name = None

        out_dict = {
            'parser_version': self._ParserVersion,
            'calculation_plugin_version': KkrCalculation._CALCULATION_PLUGIN_VERSION
        }

        # TODO job title, compound description

        # open all files in context managers
        with ExitStack() as stack:
            # open files only if they exist
            # openend files are added to the context manager stack
            # so that at the end of the parsing all files are closed
            # when the context manager is exited
            (outfile, outfile_0init, outfile_000, outfile_2, potfile_out, timing_file,
             nonco_out_file) = open_files_in_context(
                 stack, out_folder, outfile_name, outfile_0init_name, outfile_000_name, outfile_2_name,
                 potfile_out_name, timing_file_name, nonco_out_file_name
             )

            # then parse the output
            out_dict = {}
            success, msg_list, out_dict = parse_kkr_outputfile(
                out_dict,
                outfile,
                outfile_0init,
                outfile_000,
                timing_file,
                potfile_out,
                nonco_out_file,
                outfile_2,
                skip_readin=skip_mode,
                debug=debug
            )

        # try to parse with other combinations of files to minimize parser errors
        if self.icrit != 0:
            self.logger.info(f'msg_list0: {msg_list}')
            # try second combination of files
            out_dict2 = out_dict.copy()
            success2, msg_list2, out_dict2 = parse_kkr_outputfile(
                out_dict2,
                outfile_2,
                outfile_0init,
                outfile_000,
                timing_file,
                potfile_out,
                nonco_out_file,
                outfile_2,
                skip_readin=skip_mode
            )
            self.logger.info(f'msg_list1: {msg_list2}')
            if len(msg_list2) < len(msg_list):  # overwrite parser outputs if fewer errors
                self.logger.info('take output of parser run 1')
                success, msg_list, out_dict = success2, msg_list2, out_dict2
            # try third combination of files
            out_dict2 = out_dict.copy()
            success2, msg_list2, out_dict2 = parse_kkr_outputfile(
                out_dict2,
                outfile_000,
                outfile_0init,
                outfile_000,
                timing_file,
                potfile_out,
                nonco_out_file,
                outfile_2,
                skip_readin=skip_mode
            )
            self.logger.info(f'msg_list2: {msg_list2}')
            if len(msg_list2) < len(msg_list):  # overwrite parser outputs if fewer errors
                self.logger.info('take output of parser run 2')
                success, msg_list, out_dict = success2, msg_list2, out_dict2

        # TODO discriminate between real errors and warnings
        msg_list = [i for i in msg_list if 'single particle energies' not in i]
        if len(msg_list) == 0:
            success = True

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

        # create output node and link
        self.out('output_parameters', Dict(dict=out_dict))

        if self.icrit != 0 and not success:  # overwrite behavior with KKRimporter
            success = True  # set automatically to True even if only partial output was parsed
            msg = 'Automatically returned success=True for KKR importer although some parsing errors occurred'
            self.logger.warning(msg)

        if not success:
            return self.exit_codes.ERROR_KKR_PARSING_FAILED
