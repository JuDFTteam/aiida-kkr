#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow that defines the restart
functionality for the VoronoiCalculation.
"""

from aiida import orm
from aiida_kkr.workflows.restart.base_restart import CalculationBaseWorkChain
from aiida_kkr.calculations import VoronoiCalculation
from aiida_kkr.tools import kkrparams
from aiida.engine.processes.workchains.utils import process_handler, ProcessHandlerReport

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = u'Philipp Rüßmann'


class VoronoiCalculationBaseWorkChain(CalculationBaseWorkChain):
    """Workchain to run a Voronoi calculation with automated error handling and restarts"""

    # set process class (needed by BaseRestartWorkChain)
    _process_class = VoronoiCalculation

    # change default value (voronoi runs in serial and is fast)
    _max_queue_nodes = 1
    _max_queue_wallclock_sec = 3600
    _optimize_resources = False

    # list of exit codes we cannot handle automatically
    _exit_codes_nohandler = [
        VoronoiCalculation.exit_codes.ERROR_NO_OUTPUT_FILE,  # pylint: disable=no-member
        VoronoiCalculation.exit_codes.ERROR_VORONOI_PARSING_FAILED,  # pylint: disable=no-member
        VoronoiCalculation.exit_codes.ERROR_OPENING_OUTPUTS,  # pylint: disable=no-member
        VoronoiCalculation.exit_codes.ERROR_CALCULATION_FAILED,  # pylint: disable=no-member
    ]
    # common exit codes for which we know what to do
    _exit_code_memory = VoronoiCalculation.exit_codes.ERROR_NOT_ENOUGH_MEMORY  # pylint: disable=no-member
    _exit_code_timelimit = VoronoiCalculation.exit_codes.ERROR_TIME_LIMIT  # pylint: disable=no-member

    # some boilerplate code needed to set the proper exit codes to process_handler decorators
    @process_handler(priority=1, exit_codes=_exit_codes_nohandler)
    def _handle_general_error(self, calculation):
        super()._handle_general_error(calculation)

    @process_handler(priority=10, exit_codes=_exit_code_memory)
    def _handle_not_enough_memory(self, calculation):
        super()._handle_not_enough_memory(calculation)

    @process_handler(priority=20, exit_codes=_exit_code_timelimit)
    def _handle_time_limits(self, calculation):
        super()._handle_time_limits(calculation)

    # additional calculation specific process handlers go here

    @process_handler(priority=30, exit_codes=VoronoiCalculation.exit_codes.ERROR_NACLSD_TOO_SMALL)  # pylint: disable=no-member
    def _handle_naclsd(self, calculation):
        """Handle NACLSD too small error"""
        try:
            # update to smaller RCLUSTZ (here we try decreasing it by 20%)
            params = kkrparams(**self.ctx.inputs.parameters)
            params.set_multiple_values(RCLUSTZ=self.ctx.inputs.parameters['RCLUSTZ'] * 0.8)
            self.ctx.inputs.parameters = orm.Dict(params)
            return ProcessHandlerReport(True)
        except:
            return ProcessHandlerReport(True, self.exit_codes.ERROR_SOMETHING_WENT_WRONG)  # pylint: disable=no-member

    # @process_handler(...)
    # def _handle_something(self, calculation):
    #     """
    #     If calculation fails due to time some reason ...
    #     """
    #     ...

    #     return ProcessHandlerReport(True)
