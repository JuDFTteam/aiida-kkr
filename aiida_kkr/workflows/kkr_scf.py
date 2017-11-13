#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for converging a kkr calculation and
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


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.27"
__contributors__ = "Jens Broeder"


RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KkrProcess = KkrCalculation.process()
VoronoiProcess = VoronoiCalculation.process()

class kkr_scf_wc(WorkChain):
    """
    Workchain for converging a KKR calculation (SCF).

    It converges the charge density, potential?.
    Two paths are possible:

    (1) Start from a structure and run a voronoi calculation first,
    optional with calc_parameters
    (2) Start from an existing Voronoi or KKR calculation, with a remoteData

    :param wf_parameters: (ParameterData), Workchain Spezifications
    :param structure: (StructureData), Crystal structure
    :param calc_parameters: (ParameterData), Vornoi/Kkr Parameters
    :param remote_data: (RemoteData), from a KKR, or Vornoi calculation
    :param voronoi: (Code)
    :param kkr: (Code)

    :return output_kkr_scf_wc_para: (ParameterData), Information of workflow results 
        like Success, last result node, list with convergence behavior

    minimum input example:
    1. Code1, Code2, Structure, (Parameters), (wf_parameters)
    2. Code2, remote_data, (Parameters), (wf_parameters)

    maximum input example:
    1. Code1, Code2, Structure, Parameters
        wf_parameters: {
                        'queue' : String,
                        'resources' : dict({"num_machines": int, "num_mpiprocs_per_machine" : int})
                        'walltime' : int}
    2. Code2, (remote-data), wf_parameters as in 1.

    Hints:
    1. This workflow does not work with local codes!
    """

    _workflowversion = "0.1.0"
    _wf_default = {'kkr_runmax': 4,              # Maximum number of kkr jobs/starts (defauld iterations per start)
                   'density_criterion' : 0.00002,  # Stop if charge denisty is converged below this value
                   'energy_criterion' : 0.002,     # if converge energy run also this total energy convergered below this value
                   'converge_density' : True,      # converge the charge density
                   'converge_energy' : False,      # converge the total energy (usually converged before density)
                   'resue' : True,                 # AiiDA fastforwarding (currently not there yet)
                   'queue_name' : '',              # Queue name to submit jobs too
                   'resources': {"num_machines": 1},# resources to allowcate for the job
                   'walltime_sec' : 60*60,          # walltime after which the job gets killed (gets parsed to fleur)
                   'serial' : False,                # execute fleur with mpi or without
                   'custom_scheduler_commands' : ''}

    @classmethod
    def define(cls, spec):
        super(kkr_scf_wc, cls).define(spec)
        spec.input("wf_parameters", valid_type=ParameterData, required=False,
                   default=ParameterData(dict={'kkr_runmax': 4,
                                               'resources': {"num_machines": 1},
                                               'walltime_sec': 60*60,
                                               'queue_name': '',
                                               'serial' : True,
                                               'custom_scheduler_commands' : ''}))
        spec.input("structure", valid_type=StructureData, required=False)
        spec.input("calc_parameters", valid_type=ParameterData, required=False)
        #spec.input("settings", valid_type=ParameterData, required=False)
        spec.input("remote_data", valid_type=RemoteData, required=False)
        spec.input("voronoi", valid_type=Code, required=False)
        spec.input("kkr", valid_type=Code, required=True)

        spec.outline(
            cls.start,
            if_(cls.validate_input)(
                cls.run_voronoi),
            cls.run_kkr,
            while_(cls.condition)(
                cls.run_kkr,
                cls.inspect_kkr,
                cls.get_res),
            cls.return_results
        )
        #spec.dynamic_output()

    def start(self):
        """
        init context and some parameters
        """
        self.report('INFO: started KKR convergence workflow version {}\n'
                    'INFO: Workchain node identifiers: {}'
                    ''.format(self._workflowversion, ProcessRegistry().current_calc_node))

        ####### init    #######

        # internal para /control para
        self.ctx.last_calc = None
        self.ctx.loop_count = 0
        self.ctx.calcs = []
        self.ctx.abort = False

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()

        if wf_dict == {}:
            wf_dict = self._wf_default

        self.ctx.serial = wf_dict.get('serial', True)

        # set values, or defaults
        self.ctx.max_number_runs = wf_dict.get('fleur_runmax', 4)
        self.ctx.resources = wf_dict.get('resources', {"num_machines": 1})
        self.ctx.walltime_sec = wf_dict.get('walltime_sec', 60*60)
        self.ctx.queue = wf_dict.get('queue_name', '')
        self.ctx.custom_scheduler_commands = wf_dict.get('custom_scheduler_commands', '')
        self.ctx.description_wf = self.inputs.get('_description', '') + '|kkr_scf_wc|'
        self.ctx.label_wf = self.inputs.get('_label', 'kkr_scf_wc')

        # return para/vars
        self.ctx.successful = True
        self.ctx.distance = []
        self.ctx.warnings = []
        self.ctx.errors = []
        self.ctx.formula = ''

    def validate_input(self):
        """
        # validate input and find out which path (1, or 2) to take
        # return True means run inpgen if false run fleur directly
        """
        run_voronoi = True
        inputs = self.inputs

        if 'structure' in inputs:
            if not 'voronoi' in inputs:
                error = 'ERROR: StructureData was provided, but no voronoi code was provided'
                self.ctx.errors.append(error)
                self.control_end_wc(error)
        elif 'remote_data' in inputs:
            run_voronoi = False
        else:
            error = 'ERROR: No StructureData nor remote_data was provided as Input'
            self.ctx.errors.append(error)
            self.control_end_wc(error)

        if 'voronoi' in inputs:
            try:
                test_and_get_codenode(inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                error = ("The code you provided for voronoi  does not "
                         "use the plugin kkr.voro")
                self.control_end_wc(error)

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                error = ("The code you provided for kkr does not "
                         "use the plugin kkr.kkr")
                self.control_end_wc(error)

        # maybe ckeck here is unessesary...
        wf_dict = self.inputs.wf_parameters.get_dict()

        if wf_dict == {}:
            wf_dict = self._wf_default

        return run_voronoi


    def run_voronoi(self):
        """
        run the voronoi
        """
        structure = self.inputs.structure
        self.ctx.formula = structure.get_formula()
        label = 'scf: voronoi'
        description = '{} vornoi on {}'.format(self.ctx.description_wf, self.ctx.formula)

        voronoicode = self.inputs.voronoi
        if 'calc_parameters' in self.inputs:
            params = self.inputs.calc_parameters
        else:
            params = None # TODO use the defaults

        options = {"max_wallclock_seconds": self.ctx.walltime_sec,
                   "resources": self.ctx.resources,
                   "queue_name" : self.ctx.queue}

        inputs = get_inputs_voronoi(structure, voronoicode, options, label, description, params=params)
        self.report('INFO: run voronoi')
        future = submit(VoronoiProcess, **inputs)

        return ToContext(voronoi=future, last_calc=future)

    def get_kkr_param(self):
        """
        Always provides the right parameters for a KKR calc, 
        if needed extracts things from calculation before and changes the input
        """
        # prepare parameter node if needed from last kkr or voronoi output
        if 'calc_parameters' in self.inputs:
            params = self.inputs.calc_parameters
        else:
            params = None # TODO, currently will fail.  
        # TODO for now, but this should (take the defaults?) overwrite user input and overwrite some keys
        # which we gain out ouf the voronoi calc
        
        # check if user provided something
        # check if voronoi was run and last calc is a voronoi calculation
        # overwrite things from voronoi
        # ggf overwrite things from user input
        
        # if calculation before was a KKR run, check if things need to be changed..

        
        return params
        
    def run_kkr(self):
        """
        run a KKR calculation
        """

        if self.ctx['last_calc']:
            remote = self.ctx['last_calc'].out.remote_folder
        elif 'remote_data' in self.inputs:
            remote = self.inputs.remote_data
        else:
            # Something went wrong?
            remote = None
        

        params = self.get_kkr_param()   
        
        label = ' '
        description = ' '
        if self.ctx.formula:
            label = 'scf: kkr run {}'.format(self.ctx.loop_count+1)
            description = '{} kkr run {} on {}'.format(self.ctx.description_wf, self.ctx.loop_count+1, self.ctx.formula)
        else:
            label = 'scf: kkr run {}'.format(self.ctx.loop_count+1)
            description = '{} kkr run {}, voronoi given'.format(self.ctx.description_wf, self.ctx.loop_count+1)

        code = self.inputs.kkr
        options = {"max_wallclock_seconds": self.ctx.walltime_sec,
                   "resources": self.ctx.resources,
                   "queue_name" : self.ctx.queue}#,
        if self.ctx.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.ctx.custom_scheduler_commands
        inputs = get_inputs_kkr(code, remote, options, label, description, parameters=params, serial=self.ctx.serial)


        future = submit(KkrProcess, **inputs)
        self.ctx.loop_count = self.ctx.loop_count + 1
        self.report('INFO: run KKR number: {}'.format(self.ctx.loop_count))
        self.ctx.calcs.append(future)

        return ToContext(last_calc=future)


    def inspect_kkr(self):
        """
        Analyse the results of the previous KKR Calculation,
        checking whether it finished successfully or if not troubleshoot the
        cause and adapt the input parameters accordingly before
        restarting, or abort if unrecoverable error was found
        """
        #expected_states = [calc_states.FINISHED, calc_states.FAILED, 
        #                   calc_states.SUBMISSIONFAILED]
        #print(self.ctx['last_calc'])
        #self.report('I am in inspect_fleur')
        try:
            calculation = self.ctx.last_calc
        except AttributeError:
            self.ctx.successful = False
            error = 'ERROR: Something went wrong I do not have a last calculation'
            self.control_end_wc(error)
            return
        calc_state = calculation.get_state()
        #self.report('the state of the last calculation is: {}'.format(calc_state))

        if calc_state != calc_states.FINISHED:
            #TODO kill workflow in a controled way, call return results, or write a end_routine
            self.ctx.successful = False
            self.ctx.abort = True
            error = ('ERROR: Last KKR calculation failed somehow it is '
                    'in state {}'.format(calc_state))
            self.control_end_wc(error)
            return
        elif calc_state == calc_states.FINISHED:
            pass

    def inspect_voronoi(self):
        """
        Analyse if everything was ok with the voronoi calculation
        """
        pass
        
    def get_res(self):
        """
        Check how the last KKR calculation went
        Parse some results.
        """

        if self.ctx.successful:
            #self.report('last calc successful = {}'.format(self.ctx.successful))
            last_calc = self.ctx.last_calc


    def condition(self):
        """
        check convergence condition
        """

        #density_converged = False
        #energy_converged = False
        # TODO do a test first if last_calculation was successful, otherwise,
        # 'output_parameters' wont exist.
        
        return False #TODO:  So war for testing


    def return_results(self):
        """
        return the results of the calculations
        This shoudl run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """
        try:
            last_calc_uuid = self.ctx.last_calc.uuid
        except AttributeError:
            last_calc_uuid = None
        try: # if something failed, we still might be able to retrieve something
            last_calc_out = self.ctx.last_calc.out['output_parameters']
            last_calc_out_dict = last_calc_out.get_dict()
        except AttributeError:
            last_calc_out = None
            last_calc_out_dict = {}



        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['material'] = self.ctx.formula
        outputnode_dict['loop_count'] = self.ctx.loop_count
        #outputnode_dict['iterations_total'] = last_calc_out_dict.get('number_of_iterations_total', None)
        #outputnode_dict['distance_charge'] = last_calc_out_dict.get('charge_density', None)
        #outputnode_dict['distance_charge_all'] = self.ctx.distance
        #outputnode_dict['total_energy'] = last_calc_out_dict.get('energy_hartree', None)
        #outputnode_dict['total_energy_all'] = self.ctx.total_energy
        #outputnode_dict['distance_charge_units'] = 'me/bohr^3'
        #outputnode_dict['total_energy_units'] = 'Htr'
        outputnode_dict['warnings'] = self.ctx.warnings
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['last_calc_uuid'] = last_calc_uuid
        # maybe also store some information about the formula
        #also lognotes, which then can be parsed from subworkflow too workflow, list of calculations involved (pks, and uuids),
        #This node should contain everything you wish to plot, here iteration versus, total energy and distance.

        if self.ctx.successful:
            self.report('STATUS: Done, the convergence criteria are reached.\n'
                        'INFO: The charge density of the KKR calculation pk= '
                        'converged after {} KKR runs and {} iterations to {} '
                        '"me/bohr^3" \n'.format(self.ctx.loop_count,
                                       last_calc_out_dict.get('number_of_iterations_total', None),
                                       last_calc_out_dict.get('charge_density', None)))

        else: # Termination ok, but not converged yet...
            if self.ctx.abort: # some error occured, donot use the output.
                self.report('STATUS/ERROR: I abort, see logs and '
                            'erros/warning/hints in output_kkr_scf_wc_para')
            else:
                self.report('STATUS/WARNING: Done, the maximum number of runs '
                            'was reached or something failed.\n INFO: The '
                            'charge density of the KKR calculation pk= '
                            'after {} KKR runs and {} iterations is {} "me/bohr^3"\n'
                            ''.format(self.ctx.loop_count,
                            last_calc_out_dict.get('number_of_iterations_total', None),
                            last_calc_out_dict.get('charge_density', None)))

        #also lognotes, which then can be parsed from subworkflow too workflow, list of calculations involved (pks, and uuids),
        #This node should contain everything you wish to plot, here iteration versus, total energy and distance.


        outputnode_t = ParameterData(dict=outputnode_dict)
        
         # this is unsafe so far, because last_calc_out could not exist...
        if last_calc_out:
            outdict = create_scf_result_node(outpara=outputnode_t, last_calc_out=last_calc_out)
        else:
            outdict = create_scf_result_node(outpara=outputnode_t)

        if last_calc_out:
            outdict['last_kkr_calc_output'] = last_calc_out

        # TODO return other nodes? depending what calculation was run

        #outdict['output_scf_wc_para'] = outputnode
        for link_name, node in outdict.iteritems():
            self.out(link_name, node)


    def handle_kkr_failure(self):
        """
        handle a failure of a KKR calculation
        """
        return
    
    def handle_voronoi_failure(self):
        """
        Handle a failure of voronoi
        """
        return

    def control_end_wc(self, errormsg):
        """
        Controled way to shutdown the workchain. will initalize the output nodes
        """
        self.ctx.successful = False
        self.ctx.abort = True
        self.report(errormsg) # because return_results still fails somewhen
        self.return_results()
        #self.abort_nowait(errormsg)
        self.abort(errormsg)

'''
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=('SCF with FLEUR. workflow to'
                 ' converge the chargedensity and optional the total energy.'))
    parser.add_argument('--wf_para', type=ParameterData, dest='wf_parameters',
                        help='The pseudopotential family', required=False)
    parser.add_argument('--structure', type=StructureData, dest='structure',
                        help='The crystal structure node', required=False)
    parser.add_argument('--calc_para', type=ParameterData, dest='calc_parameters',
                        help='Parameters for the FLEUR calculation', required=False)
    parser.add_argument('--fleurinp', type=FleurInpData, dest='fleurinp',
                        help='FleurinpData from which to run the FLEUR calculation', required=False)
    parser.add_argument('--remote', type=RemoteData, dest='remote_data',
                        help=('Remote Data of older FLEUR calculation, '
                        'from which files will be copied (broyd ...)'), required=False)
    parser.add_argument('--inpgen', type=Code, dest='inpgen',
                        help='The inpgen code node to use', required=False)
    parser.add_argument('--fleur', type=Code, dest='fleur',
                        help='The FLEUR code node to use', required=True)

    args = parser.parse_args()
    res = run(fleur_scf_wc, 
              wf_parameters=args.wf_parameters,
              structure=args.structure,
              calc_parameters=args.calc_parameters,
              fleurinp=args.fleurinp,
              remote_data=args.remote_data,
              inpgen = args.inpgen,
              fleur=args.fleur)
'''


@wf
def create_scf_result_node(**kwargs):
    """
    This is a pseudo wf, to create the rigth graph structure of AiiDA.
    This wokfunction will create the output node in the database.
    It also connects the output_node to all nodes the information commes from.
    So far it is just also parsed in as argument, because so far we are to lazy
    to put most of the code overworked from return_results in here.
    """
    for key, val in kwargs.iteritems():
        if key == 'outpara': #  should be alwasys there
            outpara = val
    outdict = {}
    outputnode = outpara.copy()
    outputnode.label = 'output_kkr_scf_wc_para'
    outputnode.description = ('Contains self-consistency results and '
                             'information of an kkr_scf_wc run.')

    outdict['output_kkr_scf_wc_para'] = outputnode
    # copy, because we rather produce the same node twice then have a circle in the database for now...
    #output_para = args[0]
    #return {'output_eos_wc_para'}
    return outdict



def test_and_get_codenode(codenode, expected_code_type, use_exceptions=False):
    """
    Pass a code node and an expected code (plugin) type. Check that the
    code exists, is unique, and return the Code object.

    :param codenode: the name of the code to load (in the form label@machine)
    :param expected_code_type: a string with the plugin that is expected to
      be loaded. In case no plugins exist with the given name, show all existing
      plugins of that type
    :param use_exceptions: if True, raise a ValueError exception instead of
      calling sys.exit(1)
    :return: a Code object
    """
    import sys
    from aiida.common.exceptions import NotExistent
    from aiida.orm import Code


    try:
        if codenode is None:
            raise ValueError
        code = codenode
        if code.get_input_plugin_name() != expected_code_type:
            raise ValueError
    except (NotExistent, ValueError):
        from aiida.orm.querybuilder import QueryBuilder
        qb = QueryBuilder()
        qb.append(Code,
                  filters={'attributes.input_plugin':
                               {'==': expected_code_type}},
                  project='*')

        valid_code_labels = ["{}@{}".format(c.label, c.get_computer().name)
                             for [c] in qb.all()]

        if valid_code_labels:
            msg = ("Pass as further parameter a valid code label.\n"
                   "Valid labels with a {} executable are:\n".format(
                expected_code_type))
            msg += "\n".join("* {}".format(l) for l in valid_code_labels)

            if use_exceptions:
                raise ValueError(msg)
            else:
                print >> sys.stderr, msg
                sys.exit(1)
        else:
            msg = ("Code not valid, and no valid codes for {}.\n"
                   "Configure at least one first using\n"
                   "    verdi code setup".format(
                expected_code_type))
            if use_exceptions:
                raise ValueError(msg)
            else:
                print >> sys.stderr, msg
                sys.exit(1)

    return code

    
def get_inputs_kkr(code, remote, options, label='', description='', parameters=None, serial=False):
    """
    get the input for a KKR calc
    """
    inputs = KkrProcess.get_inputs_template()
    #print('Template fleur {} '.format(inputs))
    if remote:
        inputs.parent_folder = remote
    if code:
        inputs.code = code
    if parameters:
        inputs.parameters = parameters
    for key, val in options.iteritems():
        if val==None:
            continue
        else:
            inputs._options[key] = val

    if description:
        inputs['_description'] = description
    else:
        inputs['_description'] = ''
    if label:
        inputs['_label'] = label
    else:
        inputs['_label'] = ''

    if serial:
        inputs._options.withmpi = False # for now
        inputs._options.resources = {"num_machines": 1}


    '''
    options = {
    "max_wallclock_seconds": int,
    "resources": dict,
    "custom_scheduler_commands": unicode,
    "queue_name": basestring,
    "computer": Computer,
    "withmpi": bool,
    "mpirun_extra_params": Any(list, tuple),
    "import_sys_environment": bool,
    "environment_variables": dict,
    "priority": unicode,
    "max_memory_kb": int,
    "prepend_text": unicode,
    "append_text": unicode}
    '''
    return inputs


def get_inputs_voronoi(structure, voronoicode, options, label='', description='', params=None, serial=True):
    """
    get the input for a voronoi calc
    """
    inputs = VoronoiProcess.get_inputs_template()
    #print('Template voronoi {} '.format(inputs))

    if structure:
        inputs.structure = structure
    if voronoicode:
        inputs.code = voronoicode
    if params:
        inputs.parameters = params
    for key, val in options.iteritems():
        if val==None:
            #leave them out, otherwise the dict schema won't validate
            continue
        else:
            inputs._options[key] = val

    if description:
        inputs['_description'] = description
    else:
        inputs['_description'] = ''

    if label:
        inputs['_label'] = label
    else:
        inputs['_label'] = ''

    if serial:
        inputs._options.withmpi = False # for now
        inputs._options.resources = {"num_machines": 1}

    return inputs