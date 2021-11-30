# -*- coding: utf-8 -*-

from __future__ import absolute_import
from aiida.parsers.parser import Parser
from aiida.orm import Dict
from aiida_kkr.calculations.kkrnano import KKRnanoCalculation
import os

__copyright__ = (u'Copyright (c), 2021, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.1'
__contributors__ = ('Markus Struckmann', 'Philipp Rüßmann')


class KKRnanoParser(Parser):
    """
    Parser class for parsing output of the KKRnano code
    """

    def __init__(self, calc):
        """
        Initialize 
        """
        # these files should be present after success of voronoi
        self._default_files = {
            'stdout': KKRnanoCalculation._DEFAULT_OUTPUT_FILE
        }

        self._ParserVersion = __version__

        # reuse init of base class
        super(KKRnanoParser, self).__init__(calc)

    # pylint: disable=protected-access
    def parse(self, debug=False, **kwargs):
        """
        Parse output data folder, store results in database.

        :param retrieved: a dictionary of retrieved nodes, where
          the key is the link name
        :returns: nothing if everything is fine or an exit code defined in the voronoi calculation class
        """

        success = True
        node_list = ()

        # initialize out_dict and parse output files
        out_dict = {'parser_version': self._ParserVersion}
        
        #TODO fill parsed output dict

        # create output node and link
        self.out('output_parameters', Dict(dict=out_dict))

        if not success:
            return self.exit_codes.ERROR_PARSING_FAILED
