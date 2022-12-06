#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow that defines the restart
functionality for the VoronoiCalculation.
"""

from aiida import orm
from aiida_kkr.workflows.restart.base_restart import CalculationBaseWorkChain
from aiida_kkr.calculations import KkrCalculation
from aiida_kkr.tools import kkrparams
from aiida.engine.processes.workchains.utils import process_handler, ProcessHandlerReport

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = u'Philipp Rüßmann'


class KkrCalculationBaseWorkChain(CalculationBaseWorkChain):
    """Workchain to run a KKRhost calculation with automated error handling and restarts"""

    # set process class (needed by BaseRestartWorkChain)
    _process_class = KkrCalculation

    # list of exit codes we cannot handle automatically
    _exit_codes_nohandler = [
        KkrCalculation.exit_codes.ERROR_NO_OUTPUT_FILE,  # pylint: disable=no-member
        KkrCalculation.exit_codes.ERROR_KKR_PARSING_FAILED,  # pylint: disable=no-member
        KkrCalculation.exit_codes.ERROR_OPENING_OUTPUTS,  # pylint: disable=no-member
        KkrCalculation.exit_codes.ERROR_CALCULATION_FAILED,  # pylint: disable=no-member
        KkrCalculation.exit_codes.ERROR_NO_SHAPEFUN_FOUND,  # pylint: disable=no-member
    ]
    # common exit codes for which we know what to do
    _exit_code_memory = KkrCalculation.exit_codes.ERROR_NOT_ENOUGH_MEMORY  # pylint: disable=no-member
    _exit_code_timelimit = KkrCalculation.exit_codes.ERROR_TIME_LIMIT  # pylint: disable=no-member

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

    @process_handler(priority=30, exit_codes=KkrCalculation.exit_codes.ERROR_RLOG_TOO_SMALL)  # pylint: disable=no-member
    def _handle_rlog_too_small(self, calculation):
        """
        Increase R_LOG parameter
        """
        try:
            with calculation.outputs.retrieved.open('out_kkr') as _f:
                txt = _f.readlines()
            rlog = [float(i.split()[-1]) for i in txt if 'Rmesh(IRMIN' in i][0] + 0.05
            params = kkrparams(**self.ctx.inputs.parameters)
            params.set_multiple_values(R_LOG=rlog)
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
