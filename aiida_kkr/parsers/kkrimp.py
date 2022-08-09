# -*- coding: utf-8 -*-
"""
Parser for the KKR-impurity Code.
The parser should never fail, but it should catch
all errors and warnings and show them to the user.
"""

from __future__ import absolute_import
import tarfile
import os
from aiida import __version__ as aiida_core_version
from aiida.orm import Dict
from aiida.common.exceptions import InputValidationError, NotExistent
from aiida.parsers.parser import Parser
from aiida_kkr.calculations.kkrimp import KkrimpCalculation
from aiida_kkr.tools.context import open_files_in_context
from masci_tools.io.parsers.kkrparser_functions import check_error_category
from masci_tools.io.parsers.kkrimp_parser_functions import KkrimpParserFunctions
from six.moves import range
from pprint import pprint
from contextlib import ExitStack

__copyright__ = (u'Copyright (c), 2018, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.5.0'
__contributors__ = ('Philipp Rüßmann')


class KkrimpParser(Parser):
    """
    Parser class for parsing output of the KKRimp code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of KkrimpParser
        """

        self._ParserVersion = __version__

        # reuse init of base class
        super(KkrimpParser, self).__init__(calc)

    # pylint: disable=protected-access

    def parse(self, debug=False, **kwargs):
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
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = out_folder.list_object_names()

        file_errors = []
        files = {}

        # Parse output files of KKRimp calculation

        # first get path to files and catch errors if files are not present
        # append tupels (error_category, error_message) where error_category is
        # 1: critical error, always leads to failing of calculation
        # 2: warning, is inspected and checked for consistency with read-in
        #    out_dict values (e.g. nspin, newsosol, ...)

        # first go through list of critical files that should always be there
        critical_files = [
            ('outfile', KkrimpCalculation._DEFAULT_OUTPUT_FILE),
            ('out_log', KkrimpCalculation._OUTPUT_000),
            ('out_pot', KkrimpCalculation._OUT_POTENTIAL),
            ('out_timing', KkrimpCalculation._OUT_TIMING_000),
            ('out_enersp_at', KkrimpCalculation._OUT_ENERGYSP_PER_ATOM),
            ('out_enertot_at', KkrimpCalculation._OUT_ENERGYTOT_PER_ATOM),
        ]
        for keyname, fname in critical_files:
            self._check_file_existance(files, keyname, fname, 1, file_errors)

        additional_files = [
            ('kkrflex_llyfac', KkrimpCalculation._KKRFLEX_LLYFAC),
            ('kkrflex_angles', KkrimpCalculation._KKRFLEX_ANGLE),
            ('out_spinmoms', KkrimpCalculation._OUT_MAGNETICMOMENTS),
            ('out_orbmoms', KkrimpCalculation._OUT_ORBITALMOMENTS),
        ]
        for keyname, fname in additional_files:
            self._check_file_existance(files, keyname, fname, 2, file_errors)

        if debug:
            pprint(files)

        # now parse file output
        out_dict = {
            'parser_version': self._ParserVersion,
            'calculation_plugin_version': KkrimpCalculation._CALCULATION_PLUGIN_VERSION
        }

        # open all files in context managers
        with ExitStack() as stack:
            # open files only if they exist
            # openend files are added to the context manager stack
            # so that at the end of the parsing all files are closed
            # when the context manager is exited
            fhandles = open_files_in_context(stack, out_folder, *files.values())
            named_file_handles = dict([[key, fhandles[ik]] for ik, key in enumerate(files.keys())])

            # now we can parse the output files
            success, msg_list, out_dict = KkrimpParserFunctions().parse_kkrimp_outputfile(
                out_dict, named_file_handles, debug=debug
            )

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

        # cleanup after parsing (only if parsing was successful), only works below aiida-core v2.0
        if success:
            if int(aiida_core_version.split('.')[0]) < 2:
                # check if we should do the cleanup or not
                cleanup_outfiles = False
                if 'cleanup_outfiles' in self.node.inputs:
                    cleanup_outfiles = self.node.inputs.cleanup_outfiles.value
                if cleanup_outfiles:
                    # reduce size of timing file
                    self.cleanup_outfiles(files['out_timing'], ['Iteration number', 'time until scf starts'])
                    # reduce size of out_log file
                    self.cleanup_outfiles(files['out_log'], ['Iteration Number'])
                    # delete completely parsed output files and create a tar ball to reduce size
                    self.remove_unnecessary_files()
                    self.final_cleanup()
        else:
            return self.exit_codes.ERROR_PARSING_KKRIMPCALC

    def _check_file_existance(self, files, keyname, fname, icrit, file_errors):
        """Check if file called `fname` exists and then
        add it to the `files` dict with a given `keyname`.
        The icrit index determines how critical it is if a file
        is not found (1=critical error, 2=only a warning).
        """
        if fname in self.retrieved.list_object_names():
            files[keyname] = fname
        else:
            # file not
            if icrit == 1:
                crit_level = 'Critical file error!'
            elif icrit == 2:
                crit_level = 'Warning!'
            else:
                raise ValueError('icrit should be either 1 or 2')
            file_errors.append((icrit, crit_level + f" File '{fname}' not found."))
            files[keyname] = None

    def cleanup_outfiles(self, fileidentifier, keyslist):
        """open file and remove unneeded output"""
        if fileidentifier is not None:
            lineids = []
            with self.retrieved.open(fileidentifier) as tfile:
                txt = tfile.readlines()
                for iline in range(len(txt)):
                    for key in keyslist:  # go through all keys
                        if key in txt[iline]:  # add line id to list if key has been found
                            lineids.append(iline)
            # rewrite file deleting the middle part
            if len(lineids) > 1:  # cut only if more than one iteration was found
                txt = txt[:lineids[0]] + \
                    ['# ... [removed output except for last iteration] ...\n'] + \
                    txt[lineids[-1]:]
                with self.retrieved.open(fileidentifier, 'w') as tfilenew:
                    tfilenew.writelines(txt)

    def remove_unnecessary_files(self):
        """
        Remove files that are not needed anymore after parsing
        The information is completely parsed (i.e. in outdict of calculation)
        and keeping the file would just be a duplication.
        """
        # first delete unused files (completely in parsed output)
        files_to_delete = [
            KkrimpCalculation._OUT_ENERGYSP_PER_ATOM, KkrimpCalculation._OUT_ENERGYTOT_PER_ATOM,
            KkrimpCalculation._SHAPEFUN
        ]
        for fileid in files_to_delete:
            if fileid in self.retrieved.list_object_names():
                self.retrieved.delete_object(fileid, force=True)

    def final_cleanup(self):
        """Create a tarball of the rest."""

        # short name for retrieved folder
        ret = self.retrieved

        # Now create tarball of output
        #
        # check if output has been packed to tarfile already
        # only if tarfile is not there we create the output tar file
        if KkrimpCalculation._FILENAME_TAR not in ret.list_object_names():
            # first create dummy file which is used to extract the full path that is given to tarfile.open
            with ret.open(KkrimpCalculation._FILENAME_TAR, 'w') as f:
                filepath_tar = f.name

            # now create tarfile and loop over content of retrieved directory
            to_delete = []
            with tarfile.open(filepath_tar, 'w:gz') as tf:
                for f in ret.list_object_names():
                    with ret.open(f) as ftest:
                        filesize = os.stat(ftest.name).st_size
                        ffull = ftest.name
                    if (
                        f != KkrimpCalculation._FILENAME_TAR  # ignore tar file
                        and filesize > 0  # ignore empty files
                        # ignore files starting with '.' like '.nfs...'
                        and f[0] != '.'
                    ):
                        tf.add(ffull, arcname=os.path.basename(ffull))
                        to_delete.append(f)

            # finally delete files that have been added to tarfile
            for f in to_delete:
                ret.delete_object(f, force=True)
