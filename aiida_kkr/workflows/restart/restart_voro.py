#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow that defines the restart
functionality for the VoronoiCalculation.
"""

from aiida import orm
from .base_restart import CalculationBaseWorkChain
from aiida_kkr.calculations import VoronoiCalculation

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = u'Philipp Rüßmann'


class VoronoiCalculationBaseWorkChain(CalculationBaseWorkChain):
    """Workchain to run a Voronoi calculation with automated error handling and restarts"""
    _workflowversion = __version__
    _process_class = VoronoiCalculation

    _max_queue_nodes = 1
    _max_queue_wallclock_sec = 3600

    # list of exit codes we cannot handle automatically
    _exit_codes_nohandler = [
        VoronoiCalculation.exit_codes.ERROR_NO_OUTPUT_FILE,
        VoronoiCalculation.exit_codes.ERROR_VORONOI_PARSING_FAILED,
        VoronoiCalculation.exit_codes.ERROR_OPENING_OUTPUTS,
        VoronoiCalculation.exit_codes.ERROR_CALCULATION_FAILED,
    ]
    # common exit codes for which we know what to do
    _exit_code_memory = VoronoiCalculation.exit_codes.ERROR_NOT_ENOUGH_MEMORY
    _exit_code_timelimit = VoronoiCalculation.exit_codes.ERROR_TIME_LIMIT

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)

    # additional calculation specific process handlers go here
    # @process_handler(...)
    # def _handle_something(self, calculation):
    #     """
    #     If calculation fails due to time some reason ...
    #     """
    #     ...

    #     return ProcessHandlerReport(True)
