#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""

from __future__ import absolute_import
from aiida.orm import RemoteData, StructureData, Dict

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0'
__contributors__ = u'Philipp Rüßmann'


class kkr_check_mag_wc(WorkChain):
    """
    Workchain used to initialize a magnetic KKR calculation starting from the
    remoteData node of a previous calculation (either Voronoi or KKR).
    """

    _workflowversion = __version__
    _wf_default = {  # DOS params: end of energy contour
    }
    _options_default = {
        'queue_name': '',  # Queue name to submit jobs too
        'resources': {
            'num_machines': 1
        },  # resources to allowcate for the job
        'max_wallclock_seconds': 60 * 60,  # walltime after which the job gets killed (gets parsed to KKR)
        'withmpi': False,  # execute KKR with mpi or without
        'custom_scheduler_commands': ''  # some additional scheduler commands
    }

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_maginit_wc, cls).define(spec)

        # Here the structure of the workflow is defined
        spec.outline(
            #TODO impmement
        )
