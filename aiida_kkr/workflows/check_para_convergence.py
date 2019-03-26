#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""

from __future__ import absolute_import
from aiida.plugins import DataFactory

__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.0"
__contributors__ = u"Philipp Rüßmann"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
Dict = DataFactory('dict')

class kkr_check_para_wc(WorkChain):
    """
    Workchain a set of KKR calculations checking the convergence.
    """

    _workflowversion = "0.0.0"
    _wf_default = {
                   }
    _options_default = {'queue_name' : '',                        # Queue name to submit jobs too
                        'resources': {"num_machines": 1},         # resources to allowcate for the job
                        'max_wallclock_seconds' : 60*60,          # walltime after which the job gets killed (gets parsed to KKR)
                        'use_mpi' : False,                        # execute KKR with mpi or without
                        'custom_scheduler_commands' : ''          # some additional scheduler commands
                        }

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_check_para_wc, cls).define(spec)

        # Here the structure of the workflow is defined
        spec.outline(
            #TODO implement
        )
        #spec.dynamic_output()
