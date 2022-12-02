#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow that defines the restart
functionality for calculations of the AiiDA-KKR package.
"""

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import BaseRestartWorkChain, while_
from aiida.engine.processes.workchains.utils import process_handler, ProcessHandlerReport

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__contributors__ = u'Philipp Rüßmann'


class CalculationBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Calculation with automated error handling and restarts"""
    # calculation process class, needs to be overwritten from calculation restart plugin!
    _process_class = None

    # these exit codes have to be overwritten by the calculation restart plugin!
    # list of exit codes we cannot handle automatically
    _exit_codes_nohandler = None  # can be a list of exit codes
    # common exit codes for which we know what to do
    _exit_code_memory = None
    _exit_code_timelimit = None

    # default values for limits, can be overwritten
    _max_queue_nodes = 20
    _max_queue_wallclock_sec = 86400  # 24h

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)

        # inputs
        spec.expose_inputs(cls._process_class, namespace='calc')
        # additionally some parameters for retries:
        spec.input(
            'add_comp_para',
            valid_type=orm.Dict,
            default=lambda: orm.Dict(
                dict={
                    'only_even_MPI': False,
                    'forbid_single_mpi': False,
                    'max_queue_nodes': cls._max_queue_nodes,
                    'max_queue_wallclock_sec': cls._max_queue_wallclock_sec
                }
            ),
            help='Gives additional control over computational parameters'
            'only_even_MPI: set to true if you want to suppress odd number of MPI processes in parallelisation.'
            'This might speedup a calculation for machines having even number of sockets per node.'
            'max_queue_nodes: maximal number of nodes allowed on the remote machine.'
            'max_queue_wallclock_sec: maximal wallclock time allowed on the remote machine.'
        )

        # outputs
        spec.expose_outputs(cls._process_class)

        # exit codes
        spec.exit_code(397, 'ERROR_MEMORY_HANDLER_FAILED', message='Cannot use a larger timelimit.')
        spec.exit_code(398, 'ERROR_TIMELIMIT_HANDLER_FAILED', message='Cannot use a larger timelimit.')
        spec.exit_code(
            399,
            'ERROR_SOMETHING_WENT_WRONG',
            message='Calculation failed and BaseWorkChain has no strategy to resolve this'
        )

        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

    def validate_inputs(self):
        """
        Validate inputs that might depend on each other and cannot be validated by the spec.
        Also define dictionary `inputs` in the context, that will contain the inputs for the
        calculation that will be launched in the `run_calculation` step.
        """
        self.ctx.inputs = AttributeDict(self.exposed_inputs(self._process_class, namespace='calc'))

        self.ctx.max_queue_nodes = self.inputs.add_comp_para['max_queue_nodes']
        self.ctx.max_queue_wallclock_sec = self.inputs.add_comp_para['max_queue_wallclock_sec']

        input_options = self.ctx.inputs.metadata.options
        self.ctx.optimize_resources = input_options.pop('optimize_resources', True)
        self.ctx.inputs.metadata.options = input_options

        if 'description' in self.inputs:
            self.ctx.inputs.metadata.description = self.inputs.description
        else:
            self.ctx.inputs.metadata.description = ''
        if 'label' in self.inputs:
            self.ctx.inputs.metadata.label = self.inputs.label
        else:
            self.ctx.inputs.metadata.label = ''

        if not self.ctx.optimize_resources:
            self.ctx.can_be_optimised = False  # set this for handlers to not change resources
            return

        resources_input = self.ctx.inputs.metadata.options['resources']
        try:
            self.ctx.num_machines = int(resources_input['num_machines'])
            self.ctx.num_mpiprocs_per_machine = int(resources_input['num_mpiprocs_per_machine'])
        except KeyError:
            self.ctx.can_be_optimised = False
            self.report('WARNING: Computation resources were not optimised.')
        else:
            try:
                self.ctx.num_cores_per_mpiproc = int(resources_input['num_cores_per_mpiproc'])
                self.ctx.use_omp = True
                self.ctx.suggest_mpi_omp_ratio = self.ctx.num_mpiprocs_per_machine / self.ctx.num_cores_per_mpiproc
            except KeyError:
                self.ctx.num_cores_per_mpiproc = 1
                self.ctx.use_omp = False
                self.ctx.suggest_mpi_omp_ratio = 1

    @process_handler(priority=1, exit_codes=_exit_codes_nohandler)
    def _handle_general_error(self, calculation):
        """
        Calculation failed for unknown reason.
        """
        self.ctx.restart_calc = calculation
        self.ctx.is_finished = True
        self.report('Calculation failed for a reason that can not be resolved automatically')
        self.results()
        return ProcessHandlerReport(True, self.exit_codes.ERROR_SOMETHING_WENT_WRONG)

    @process_handler(priority=10, exit_codes=_exit_code_memory)
    def _handle_not_enough_memory(self, calculation):
        """
        Calculation failed due to lack of memory.
        """

        if not self.ctx.can_be_optimised:
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True
            self.report(
                'I am not allowed to optimize your settings. Consider providing at least'
                'num_machines and num_mpiprocs_per_machine'
            )
            self.results()
            return ProcessHandlerReport(True, self.exit_codes.ERROR_MEMORY_HANDLER_FAILED)

        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.report(
            'Calculation failed due to lack of memory, I resubmit it with twice larger'
            ' amount of computational nodes and smaller MPI/OMP ratio'
        )

        # increase number of nodes
        propose_nodes = self.ctx.num_machines * 2
        if propose_nodes > self.ctx.max_queue_nodes:
            propose_nodes = self.ctx.max_queue_nodes
        self.ctx.num_machines = propose_nodes

        return ProcessHandlerReport(True)

    @process_handler(priority=20, exit_codes=_exit_code_timelimit)
    def _handle_time_limits(self, calculation):
        """
        If calculation fails due to time limits, we simply resubmit it.
        """
        from aiida.common.exceptions import NotExistent

        # if previous calculation failed for the same reason, do not restart
        prev_calculation_remote_link = calculation.get_incoming(node_class=orm.CalcJobNode).first()
        if prev_calculation_remote_link is not None:
            prev_calculation_remote = prev_calculation_remote_link.node
            prev_calculation_status = prev_calculation_remote.creator.exit_status
            if prev_calculation_status in self._process_class.get_exit_statuses(['ERROR_TIME_LIMIT']):
                self.ctx.is_finished = True
                self.results()
                return ProcessHandlerReport(True)

        self.report('Calculation failed due to time limits, I restart it from where it ended')

        # increase wallclock time
        if 'max_wallclock_seconds' in self.ctx.inputs.metadata.options:
            propose_wallclock = self.ctx.inputs.metadata.options['max_wallclock_seconds'] * 2
            if propose_wallclock > self.ctx.max_queue_wallclock_sec:
                propose_wallclock = self.ctx.max_queue_wallclock_sec
            self.inputs.calc.metadata.options['max_wallclock_seconds'] = propose_wallclock
        else:
            return ProcessHandlerReport(True, self.exit_codes.ERROR_TIMELIMIT_HANDLER_FAILED)

        # increase number of nodes
        propose_nodes = self.ctx.num_machines * 2
        if propose_nodes > self.ctx.max_queue_nodes:
            propose_nodes = self.ctx.max_queue_nodes
        self.ctx.num_machines = propose_nodes

        remote = calculation.get_outgoing().get_node_by_label('remote_folder')

        # # resubmit providing inp.xml and cdn from the remote folder
        # self.ctx.is_finished = False
        # if _is_remote_reusable(self.ctx.inputs, calculation):
        #     if 'fleurinpdata' in self.ctx.inputs:
        #         del self.ctx.inputs.fleurinpdata
        #     self.ctx.inputs.parent_folder = remote

        return ProcessHandlerReport(True)


# def _is_remote_reusable(inputs, calculation):
#     """
#     Check whether the remote folder of the given calculation
#     can be resubmitted
#     """
#     can_use_remote = False
#     #If no charge density file is available to restart from the calculation will except
#     #with a not nice error message. So we can only reuse the charge density if these files are available
#     retrieved_filenames = calculation.get_outgoing().get_node_by_label('retrieved').list_object_names()
#     if any(file in retrieved_filenames for file in (
#             'cdn_last.hdf',
#             'cdn1',
#     )):
#         can_use_remote = True

#     if 'fleurinpdata' in inputs:
#         modes = inputs.fleurinpdata.get_fleur_modes()
#         if modes['force_theorem'] or modes['dos'] or modes['band']:
#             # in modes listed above it makes no sense copying cdn.hdf
#             can_use_remote = False
#     # without fleurinp it is harder to extract modes in this case
#     # - simply try to reuse cdn.hdf and hope it works

#     return can_use_remote
