#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""

from aiida.orm import Code, DataFactory
from aiida.work.workchain import WorkChain, while_, if_, ToContext
from aiida.work.run import submit, run
from aiida.work import workfunction as wf
from aiida.common.datastructures import calc_states
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from masci_tools.io.kkr_params import kkrparams


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.2"
__contributors__ = u"Philipp Rüßmann"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()

class kkr_check_para_wc(WorkChain):
    """
    Workchain a set of KKR calculations checking the convergence. 
    Starts from the remoteData node of a previous KKR calculation.

    :param wf_parameters: (ParameterData), Workchain specifications
    :param remote_data: (RemoteData), from a KKR calculation
    :param kkr: (Code)

    :return results_kkr_check_para_wc: (ParameterData), Information of workflow results 
        like Success, last result node, ...
    """

    _workflowversion = "0.1.0"
    _wf_default = {'dos_params' : {"nepts": 61,              # DOS params: number of points in contour
                                   "tempr": 200,             # DOS params: temperature
                                   "emin": -1,               # DOS params: start of energy contour
                                   "emax": 1}                 # DOS params: end of energy contour
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
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict=cls._wf_default))
        spec.input("remote_data", valid_type=RemoteData, required=False)
        spec.input("kkr", valid_type=Code, required=True)

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,
            if_(cls.validate_input)(
                cls.run_voronoi,
                cls.get_dos),
            cls.run_kkr,
            cls.inspect_kkr,
            while_(cls.condition)(
                cls.run_kkr,
                cls.inspect_kkr,
                cls.get_res),
            cls.return_results
        )
        #spec.dynamic_output()


    
