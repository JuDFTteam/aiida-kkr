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
from aiida.work.process_registry import ProcessRegistry
from aiida.common.datastructures import calc_states
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.tools.kkr_params import kkrparams


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = u"Philipp Rüßmann"


StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()
VoronoiProcess = VoronoiCalculation.process()

class kkr_startpot_wc(WorkChain):
    """
    Workchain to start a KKR calculation by running voronoi and getting the 
    starting DOS for first checks on the validity of the input setting. 
    Starts from a structure together with a KKR parameter node.

    :param wf_parameters: (ParameterData), Workchain specifications
    :param structure: (StructureData), aiida structure node to begin 
        calculation from (needs to contain vacancies, if KKR needs empty spheres)
    :param kkr: (Code)
    :param voronoi: (Code)

    :return result_kkr_startpot_wc: (ParameterData), Information of workflow results 
        like Success, last result node, dos array data
    """

    _workflowversion = "0.1.0"
    _wf_default = {'queue_name' : '',                        # Queue name to submit jobs too
                   'resources': {"num_machines": 1},         # resources to allowcate for the job
                   'walltime_sec' : 60*60,                   # walltime after which the job gets killed (gets parsed to KKR)
                   'mpirun' : False,                         # execute KKR with mpi or without
                   'custom_scheduler_commands' : '',         # some additional scheduler commands 
                   'dos_params' : {"nepts": 61,              # DOS params: number of points in contour
                                   "tempr": 200,             # DOS params: temperature
                                   "emin": -1,               # DOS params: start of energy contour
                                   "emax": 1}                 # DOS params: end of energy contour
                   }

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow. 
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_startpot_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict=cls._wf_default))
        spec.input("structure", valid_type=StructureData, required=True)
        spec.input("kkr", valid_type=Code, required=True)
        spec.input("voronoi", valid_type=Code, required=True)

        # Here the structure of the workflow is defined
        spec.outline(
            # initialize workflow
            cls.start,
            # check input and run voronoi
            if_(cls.validate_input)(
                cls.run_voronoi),
            # check voronoi output and create starting DOS using dos sub-workflow
            if_(cls.condition)(cls.get_dos),
            # collect results and return
            cls.return_results
        )

    