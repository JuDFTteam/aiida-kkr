#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""
from masci_tools.io.kkr_params import kkrparams
from aiida import orm
from aiida.engine import WorkChain, if_, ToContext
from aiida.engine import submit
from aiida_kkr.tools.common_workfunctions import (
    test_and_get_codenode,
    get_parent_paranode,
    update_params_wf,
    get_inputs_kkr,
)
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.engine import CalcJob, calcfunction
from aiida.common.exceptions import InputValidationError
from aiida_kkr.tools.save_output_nodes import create_out_dict_node
from aiida_kkr.workflows.bs import set_energy_params

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.8.0'
__contributors__ = u'Philipp Rüßmann'


class kkr_dos_wc(WorkChain):
    """
    Workchain a DOS calculation with KKR starting from the remoteData node
    of a previous calculation (either Voronoi or KKR).

    :param wf_parameters: (Dict); Workchain specifications
    :param options: (Dict); specifications for the computer
    :param remote_data: (RemoteData), mandatory; from a KKR or Vornoi calculation
    :param kkr: (Code), mandatory; KKR code running the dos calculation

    :return result_kkr_dos_wc: (Dict), Information of workflow results
        like Success, last result node, list with convergence behavior
    """

    _workflowversion = __version__
    _wf_label = 'kkr_dos_wc'
    _wf_description = 'Workflow for a KKR dos calculation starting either from a structure with automatic voronoi calculation or a valid RemoteData node of a previous calculation.'
    _wf_default = {
        'nepts': 96,  # number of points in contour
        'tempr': 200.,  # K, smearing temperature
        'emin': -10.,  # eV, rel to EF, start of energy contour
        'emax': 5.,  # eV       end of energy contour
        'kmesh': [30, 30, 30],  # kmesh for DOS calculation (typically higher than in scf contour)
        'RCLUSTZ': None,  # cluster radiu, only used if a value is set
    }
    _options_default = {
        'queue_name': '',  # Queue name to submit jobs to
        # resources to allowcate for the job
        'resources': {
            'num_machines': 1
        },
        # walltime after which the job gets killed (gets parsed to KKR)
        'max_wallclock_seconds': 60 * 60,
        'withmpi': True,  # execute KKR with mpi or without
        'custom_scheduler_commands': '',  # some additional scheduler commands
    }

    # intended to guide user interactively in setting up a valid wf_params node
    @classmethod
    def get_wf_defaults(cls, silent=False):
        """
        Print and return _wf_defaults dictionary.

        Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """
        if not silent:
            print(f'Version of workflow: {cls._workflowversion}')
        return cls._wf_default

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow.
        """
        # Take input of the workflow or use defaults defined above
        super(kkr_dos_wc, cls).define(spec)
        spec.input(
            'wf_parameters',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
        )
        spec.input(
            'options',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
        )
        spec.input(
            'remote_data',
            valid_type=orm.RemoteData,
            required=True,
        )
        spec.input(
            'kkr',
            valid_type=orm.Code,
            required=True,
        )

        # define outputs
        spec.output(
            'results_wf',
            valid_type=orm.Dict,
            required=True,
        )
        spec.output(
            'dos_data',
            valid_type=orm.XyData,
            required=False,
        )
        spec.output(
            'dos_data_interpol',
            valid_type=orm.XyData,
            required=False,
        )

        # Here the structure of the workflow is defined
        spec.outline(
            # initialize workflow
            cls.start,
            # validate input
            if_(cls.validate_input)(
                # set DOS contour in parameter node
                cls.set_params_dos,
                # calculate DOS and interpolate result
                # TODO: encapsulate get_dos in restarting mechanism (should be a base class of workflows that start calculations)
                # i.e. use base_restart_calc workchain as parent
                cls.get_dos
            ),
            #  collect results and store DOS output as ArrayData node
            # (dos, lmdos, dos.interpolate, ...)
            cls.return_results
        )

        # definition of exit code in case something goes wrong in this workflow
        spec.exit_code(
            161,
            'ERROR_NO_INPUT_REMOTE_DATA',
            'No remote_data was provided as Input',
        )
        spec.exit_code(
            162,
            'ERROR_KKRCODE_NOT_CORRECT',
            'The code you provided for kkr does not use the plugin kkr.kkr',
        )
        spec.exit_code(
            163,
            'ERROR_CALC_PARAMETERS_INVALID',
            'calc_parameters given are not consistent! Hint: did you give an unknown keyword?',
        )
        spec.exit_code(
            164,
            'ERROR_CALC_PARAMETERS_INCOMPLETE',
            'calc_parameters not complete',
        )
        spec.exit_code(
            165,
            'ERROR_DOS_PARAMS_INVALID',
            'dos_params given in wf_params are not valid',
        )
        spec.exit_code(
            166,
            'ERROR_DOS_CALC_FAILED',
            'KKR dos calculation failed',
        )

    def start(self):
        """
        init context and some parameters
        """
        self.report(f'INFO: started KKR dos workflow version {self._workflowversion}')

        # input para
        wf_dict = self.inputs.wf_parameters.get_dict()
        options_dict = self.inputs.options.get_dict()

        # TODO: check for completeness
        if wf_dict == {}:
            wf_dict = self._wf_default
            self.report('INFO: using default wf parameter')
        if options_dict == {}:
            options_dict = self._options_default
            self.report('INFO: using default options')

        # set values, or defaults
        self.ctx.withmpi = options_dict.get(
            'withmpi',
            self._options_default['withmpi'],
        )
        self.ctx.resources = options_dict.get(
            'resources',
            self._options_default['resources'],
        )
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds',
            self._options_default['max_wallclock_seconds'],
        )
        self.ctx.queue = options_dict.get(
            'queue_name',
            self._options_default['queue_name'],
        )
        self.ctx.custom_scheduler_commands = options_dict.get(
            'custom_scheduler_commands',
            self._options_default['custom_scheduler_commands'],
        )

        self.ctx.dos_params_dict = self.inputs.wf_parameters.get_dict()
        self.ctx.dos_kkrparams = None  # is set in set_params_dos

        self.ctx.description_wf = self.inputs.get(
            'description',
            self._wf_description,
        )
        self.ctx.label_wf = self.inputs.get(
            'label',
            self._wf_label,
        )

        self.report(
            f'INFO: use the following parameter:\n'
            f'withmpi: {self.ctx.withmpi}\n'
            f'Resources: {self.ctx.resources}\n'
            f'Walltime (s): {self.ctx.max_wallclock_seconds}\n'
            f'queue name: {self.ctx.queue}\n'
            f'scheduler command: {self.ctx.custom_scheduler_commands}\n'
            f'description: {self.ctx.description_wf}\n'
            f'label: {self.ctx.label_wf}\n'
            f'dos_params: {self.ctx.dos_params_dict}\n'
        )

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''

    def validate_input(self):
        """
        # validate input and find out which path (1, or 2) to take
        # return True means run voronoi if false run kkr directly
        """
        inputs = self.inputs

        if 'remote_data' in inputs:
            input_ok = True
        else:
            input_ok = False
            return self.exit_codes.ERROR_NO_INPUT_REMOTE_DATA  # pylint: disable=no-member

        # extract correct remote folder of last calculation if input
        # remote_folder node is not from KkrCalculation but kkr_scf_wc workflow
        input_remote = self.inputs.remote_data
        # check if input_remote has single KkrCalculation parent
        parents = input_remote.get_incoming(node_class=orm.CalcJobNode)
        nparents = len(parents.all_link_labels())
        if nparents != 1:
            # extract parent workflow and get uuid of last calc from output node
            parent_workflow = input_remote.inputs.last_RemoteData
            if not isinstance(parent_workflow, orm.WorkChainNode):
                raise InputValidationError(
                    'Input remote_data node neither output of a KKR/voronoi calculation nor of kkr_scf_wc workflow'
                )
            parent_workflow_out = parent_workflow.outputs.output_kkr_scf_wc_ParameterResults
            uuid_last_calc = parent_workflow_out.get_dict().get('last_calc_nodeinfo').get('uuid')
            last_calc = orm.load_node(uuid_last_calc)
            if not isinstance(last_calc, KkrCalculation) and not isinstance(last_calc, VoronoiCalculation):
                raise InputValidationError(
                    'Extracted last_calc node not of type KkrCalculation: check remote_data input node'
                )
            # overwrite remote_data node with extracted remote folder
            output_remote = last_calc.outputs.remote_folder
            self.inputs.remote_data = output_remote

        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                input_ok = False
                return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT  # pylint: disable=no-member

        # set self.ctx.input_params_KKR
        self.ctx.input_params_KKR = get_parent_paranode(self.inputs.remote_data)

        return input_ok

    def set_params_dos(self):
        """
        take input parameter node and change to DOS contour according to input from wf_parameter input
        internally calls the update_params work function to keep track of provenance
        """
        params = self.ctx.input_params_KKR
        input_dict = params.get_dict()
        para_check = kkrparams()

        # step 1 try to fill keywords
        try:
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)
        except:
            return self.exit_codes.ERROR_CALC_PARAMETERS_INVALID  # pylint: disable=no-member

        # step 2: check if all mandatory keys are there
        label = ''
        descr = ''
        missing_list = para_check.get_missing_keys(use_aiida=True)
        if missing_list != []:
            kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
            kkrdefaults_updated = []
            for key_default, val_default in list(kkrdefaults.items()):
                if key_default in missing_list:
                    para_check.set_value(key_default, kkrdefaults.get(key_default), silent=True)
                    kkrdefaults_updated.append(key_default)
                    missing_list.remove(key_default)
            if len(missing_list) > 0:
                self.report(f'ERROR: calc_parameters misses keys: {missing_list}')
                return self.exit_codes.ERROR_CALC_PARAMETERS_INCOMPLETE  # pylint: disable=no-member
            self.report(f'updated KKR parameter node with default values: {kkrdefaults_updated}')
            label = 'add_defaults_'
            descr = 'added missing default keys, '

        # overwrite energy contour to DOS contour no matter what is in input parameter node.
        # Contour parameter given as input to workflow.
        econt_new = self.ctx.dos_params_dict
        # always overwrite NPOL, N1, N3, thus add these to econt_new
        econt_new['NPOL'] = 0
        econt_new['NPT1'] = 0
        econt_new['NPT3'] = 0
        kkr_calc = self.inputs.remote_data.get_incoming().first().node
        ef = kkr_calc.outputs.output_parameters.get_dict()['fermi_energy']  # unit in Ry
        try:
            para_check = set_energy_params(econt_new, ef, para_check)
        except:
            return self.exit_codes.ERROR_DOS_PARAMS_INVALID  # pylint: disable=no-member

        updatenode = orm.Dict(dict=para_check.get_dict())
        updatenode.label = label + 'KKRparam_DOS'
        updatenode.description = descr + \
            'KKR parameter node extracted from parent parameters and wf_parameter input node.'

        paranode_dos = update_params_wf(self.ctx.input_params_KKR, updatenode)
        self.ctx.dos_kkrparams = paranode_dos

    def get_dos(self):
        """
        submit a dos calculation and interpolate result if returns complete
        """

        label = 'KKR DOS calc.'
        dosdict = self.ctx.dos_params_dict
        description = 'dos calc: emin= {}, emax= {}, nepts= {}, tempr={}, kmesh={}'.format(
            dosdict['emin'], dosdict['emax'], dosdict['nepts'], dosdict['tempr'], dosdict['kmesh']
        )
        code = self.inputs.kkr
        remote = self.inputs.remote_data
        params = self.ctx.dos_kkrparams
        options = {
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'resources': self.ctx.resources,
            'queue_name': self.ctx.queue
        }  # ,
        if self.ctx.custom_scheduler_commands:
            options['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        inputs = get_inputs_kkr(
            code,
            remote,
            options,
            label,
            description,
            parameters=params,
            serial=(not self.ctx.withmpi),
        )

        # run the DOS calculation
        self.report('INFO: doing calculation')
        dosrun = self.submit(KkrCalculation, **inputs)

        # for restart workflow:
        self.ctx.last_calc = dosrun

        return ToContext(dosrun=dosrun)

    def return_results(self):
        """
        Collect results, parse DOS output and link output nodes to workflow node
        """

        # check wether or not calculation was taken from cached node
        caching_info = f'INFO: cache_source of dos calc node: {self.ctx.dosrun.get_cache_source()}'
        print(caching_info)
        self.report(caching_info)

        # capture error of unsuccessful DOS run
        if not self.ctx.dosrun.is_finished_ok:
            self.ctx.successful = False
            error = (f'ERROR: DOS calculation failed somehow it is in state {self.ctx.dosrun.process_state}')
            print(error)
            from pprint import pprint
            pprint(self.ctx.dosrun.attributes)
            print('stdout', self.ctx.dosrun.get_scheduler_stdout())
            print('stderr', self.ctx.dosrun.get_scheduler_stderr())
            self.report(error)
            self.ctx.errors.append(error)
            return self.exit_codes.ERROR_DOS_CALC_FAILED  # pylint: disable=no-member

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['withmpi'] = self.ctx.withmpi
        outputnode_dict['resources'] = self.ctx.resources
        outputnode_dict['max_wallclock_seconds'] = self.ctx.max_wallclock_seconds
        outputnode_dict['queue_name'] = self.ctx.queue
        outputnode_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        outputnode_dict['dos_params'] = self.ctx.dos_params_dict
        try:
            outputnode_dict['nspin'] = self.ctx.dosrun.res.nspin
        except:
            error = 'ERROR: nspin not extracted'
            self.report(error)
            self.ctx.successful = False
            self.ctx.errors.append(error)
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['list_of_errors'] = self.ctx.errors

        # create output node with data-provenance
        outputnode = orm.Dict(dict=outputnode_dict)
        outputnode.label = 'kkr_scf_wc_results'
        outputnode.description = ''

        self.report('INFO: create dos results nodes')
        try:
            self.report(f'INFO: create dos results nodes. dos calc retrieved node={self.ctx.dosrun.outputs.retrieved}')
            has_dosrun = True
        except AttributeError as _error:
            self.report('ERROR: no dos calc retrieved node found')
            self.report('Caught AttributeError {}'.format(_error))
            has_dosrun = False

        # interpol dos file and store to XyData nodes
        if has_dosrun:
            dos_retrieved = self.ctx.dosrun.outputs.retrieved
            if 'complex.dos' in dos_retrieved.list_object_names():
                dosXyDatas = parse_dosfiles(dos_retrieved)

        # collect output nodes with link labels
        outdict = {}
        if has_dosrun:
            outdict['dos_data'] = dosXyDatas['dos_data']
            outdict['dos_data_interpol'] = dosXyDatas['dos_data_interpol']
        # create data provenance of results node
        link_nodes = outdict.copy()  # link also dos output nodes
        if has_dosrun:
            link_nodes['doscalc_remote'] = self.ctx.dosrun.outputs.remote_folder
        outdict['results_wf'] = create_out_dict_node(outputnode, **link_nodes)

        # create links to output nodes
        for link_name, node in outdict.items():
            self.out(link_name, node)

        self.report('INFO: done with DOS workflow!\n')


@calcfunction
def parse_dosfiles(dos_retrieved):
    """
    parse dos files to XyData nodes
    """
    from masci_tools.io.common_functions import interpolate_dos
    from masci_tools.io.common_functions import get_Ry2eV

    eVscale = get_Ry2eV()

    with dos_retrieved.open('complex.dos') as dosfolder:
        ef, dos, dos_int = interpolate_dos(dosfolder, return_original=True)

    # convert to eV units
    dos[:, :, 0] = (dos[:, :, 0] - ef) * eVscale
    dos[:, :, 1:] = dos[:, :, 1:] / eVscale
    dos_int[:, :, 0] = (dos_int[:, :, 0] - ef) * eVscale
    dos_int[:, :, 1:] = dos_int[:, :, 1:] / eVscale

    # create output nodes
    dosnode = orm.XyData()
    dosnode.set_x(dos[:, :, 0], 'E-EF', 'eV')
    name = ['tot', 's', 'p', 'd', 'f', 'g']
    name = name[:len(dos[0, 0, 1:]) - 1] + ['ns']
    ylists = [[], [], []]
    for line, _name in enumerate(name):
        ylists[0].append(dos[:, :, 1 + line])
        ylists[1].append('dos {}'.format(_name))
        ylists[2].append('states/eV')
    dosnode.set_y(ylists[0], ylists[1], ylists[2])
    dosnode.label = 'dos_data'
    dosnode.description = 'Array data containing uniterpolated DOS (i.e. dos at finite imaginary part of energy). 3D array with (atoms, energy point, l-channel) dimensions.'

    # now create XyData node for interpolated data
    dosnode2 = orm.XyData()
    dosnode2.set_x(dos_int[:, :, 0], 'E-EF', 'eV')
    ylists = [[], [], []]
    for line, _name in enumerate(name):
        ylists[0].append(dos_int[:, :, 1 + line])
        ylists[1].append('interpolated dos {}'.format(_name))
        ylists[2].append('states/eV')
    dosnode2.set_y(ylists[0], ylists[1], ylists[2])
    dosnode2.label = 'dos_interpol_data'
    dosnode2.description = 'Array data containing interpolated DOS (i.e. dos at real axis). 3D array with (atoms, energy point, l-channel) dimensions.'

    return {'dos_data': dosnode, 'dos_data_interpol': dosnode2}
