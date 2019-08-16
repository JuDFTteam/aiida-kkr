#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow that defines the restart
caclulation functionality.
"""

from aiida_kkr.calculations import KkrimpCalculation, VoronoiCalculation, KkrCalculation

__copyright__ = (u"Copyright (c), 2019, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = u"Philipp Rüßmann"


# similar to base restart workflow of aiida-quantumespresso

class base_restart_calc(WorkChain):
    """
    Use this base workflow in workflows that directly sumbit calculations (not workflows)

    This workflow defines the following functions that can be used to encapsulate functions that submit calculations:
      * cls.prepare_new_restart_calc
      * cls.inspect_restart_calc_done

    You can set the behavior of the restart block by setting these values (defaults after the `=`:
      * _max_iter_restart_calc = 2
      * _clean_workdir_restart_calc = False

    Then you should make sure to run prepare_new_restart_calc each time you want to do a new calculation (with restarts).
    
    Make sure self.ctx.last_calc is being set to the calculation you want to be able to resubmit.

    If you want to cleanup all called claculations after running you should have a call to `cleanup_all_calcs_work` in your workflow.

    usage example:
        class my_workflow(base_restart_calc):

        ...

        spec.outline(
            ...
            #here the calculation is submitted:
            cls.prepare_new_restart_calc,
            while_(cls.inspect_restart_calc_done)(
                cls.my_custom_run_function),
            # now the calculation should have finished, eventually after being restarted to overcome cluster problems
            ...
            # if we want to cleanup the work directories (`my_workflow` class needs to 
            # `set _clean_workdir_restart_calc` to True if the work directories are supposed to be removed!)
            cls.cleanup_all_calcs_work
        )
    """

    _workflowversion = __version__

    _max_iter_restart_calc = 2
    _clean_workdir_restart_calc = False

    def __init__(self, *args, **kwargs):
        super(base_restart_calc, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(base_restart_calc, cls).define(spec)

        # define general exit codes for restart of calculations
        spec.exit_code(100, 'ERROR_ITERATION_RETURNED_NO_CALCULATION', 
                       message='the run_calculation step did not successfully add a calculation node to the context')

    def prepare_new_restart_calc(self):
        """initialize counters etc. used within inspect_restart_calc_done"""
        self.ctx.restart_iteration = 0
        self.ctx.last_calc_restart_failed = False

    def inspect_restart_calc_done(self):
        """check if the calculation should be submitted"""
        # initial go
        if self.ctx.restart_iteration == 0:
            self.ctx.restart_iteration += 1
            return True

        #  stop after _max_iter_restart_calc iterations
        if self.ctx.restart_iteration>self._max_iter_restart_calc:
            self.cleanup_calc_work()
            return False

        # extract last calculation
        try:
            last_calculation = self.ctx.last_calc
        except:
            self.report('')
            return self.exit_codes.ERROR_ITERATION_RETURNED_NO_CALCULATION

        # now inspect last calculation
        if last_calculation.is_finished_ok:
            # calculation already succesful
            self.cleanup_calc_work()
            return False
        else:
            # check if this is the second consecutive failure
            if self.ctx.last_calc_restart_failed = True:
                self.cleanup_calc_work()
                return False

            # here we know already that the current calculation has failed
            self.ctx.last_calc_restart_failed = True

            # TODO: handle other known errors (maybe pause calculation until some problem is fixed)
            # TODO: maybe capture out of walltime or out of memory issues here as well (could automatically increase the options)

            # if nothing else was caught, run the calculation
            self.ctx.restart_iteration += 1
            return True

    # this is probably not working yet (check again with aiida-quantumespresso plugin)
    def cleanup_all_calcs_work(self):
        # cleanup work directories of all called calculations
        if self._clean_workdir_restart_calc:
            cleaned_calcs = []

            for called_descendant in self.calc.called_descendants:
                if isinstance(called_descendant, CalcJobNode):
                    try:
                        called_descendant.outputs.remote_folder._clean()
                        cleaned_calcs.append(called_descendant.pk)
                    except (IOError, OSError, KeyError):
                        pass
    
            if cleaned_calcs:
                self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))

