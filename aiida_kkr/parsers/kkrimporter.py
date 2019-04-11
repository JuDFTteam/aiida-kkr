# -*- coding: utf-8 -*-
"""
Parser for the KKR imprter, slight modification to KKr parser (dealing of missing output files).
The parser should never fail, but it should catch 
all errors and warnings and show them to the user.
"""

from __future__ import absolute_import
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.parsers.kkr import KkrParser
from aiida.common.exceptions import InputValidationError


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = ("Philipp Rüßmann")


class KkrImporterParser(KkrParser):
    """
    Parser class for parsing output of KKR code after import
    """

    def __init__(self, calc):
        """
        Initialize the instance of KkrParser
        """
        # check for valid input
        if not isinstance(calc, KkrCalculation):
            raise InputValidationError("Input calc must be a KkrCalculation")
        
        self._ParserVersion = __version__

        #reuse init of base class but select icrit=1 (determines how missing files are interpreted)
        self.icrit = 1
        super(KkrParser, self).__init__(calc)
