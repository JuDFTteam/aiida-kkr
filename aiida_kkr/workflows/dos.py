#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""

from aiida import orm
from masci_tools.io.kkr_params import kkrparams
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
from aiida.common.exceptions import InputValidationError
from aiida_kkr.tools.parse_dos import parse_dosfiles
from aiida_kkr.tools.save_output_nodes import create_out_dict_node
from aiida_kkr.tools.extract_kkrhost_noco_angles import extract_noco_angles
from aiida_kkr.workflows.bs import set_energy_params

__copyright__ = (u'Copyright (c), 2017, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.8.4'
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
        'prepend_text': '',
        'append_text': '',
        'additional_retrieve_list': None
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
        return cls._wf_default.copy()

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
            help='Workflow parameter (see `kkr_dos_wc.get_wf_defaults()`).'
        )
        spec.input(
            'options',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
            help='Computer options used by the workflow.'
        )
        spec.input(
            'remote_data', valid_type=orm.RemoteData, required=True, help='RemoteData node of the parent calculation.'
        )
        spec.input('kkr', valid_type=orm.Code, required=True, help='KKRhost Code node used to run the DOS calculation.')
        spec.input(
            'initial_noco_angles',
            valid_type=orm.Dict,
            required=False,
            help="""Initial non-collinear angles for the magnetic moments. See KkrCalculation for details.
            If this is found in the input potentially extracted nonco angles from the parent calulation are overwritten!"""
        )
        # maybe overwrite some settings from the KKRhost convergence run
        spec.input(
            'params_kkr_overwrite',
            valid_type=orm.Dict,
            required=False,
            help='Overwrite some input parameters of the parent KKR calculation.'
        )
        # expose LDAU input node
        spec.input(
            'settings_LDAU',
            valid_type=orm.Dict,
            required=False,
            help="""
Settings for running a LDA+U calculation. The Dict node should be of the form
    settings_LDAU = Dict(dict={'iatom=0':{
        'L': 3,         # l-block which gets U correction (1: p, 2: d, 3: f-electrons)
        'U': 7.,        # U value in eV
        'J': 0.75,      # J value in eV
        'Eref_EF': 0.,  # reference energy in eV relative to the Fermi energy. This is the energy where the projector wavefunctions are calculated (should be close in energy where the states that are shifted lie (e.g. for Eu use the Fermi energy))
    }})
    Note: you can add multiple entries like the one for iatom==0 in this example. The atom index refers to the corresponding atom in the impurity cluster.
"""
        )

        # define outputs
        spec.output('results_wf', valid_type=orm.Dict, required=True, help='Results collected by the workflow.')
        spec.output('dos_data', valid_type=orm.XyData, required=False, help='XyData node of the parsed DOS output.')
        spec.output(
            'dos_data_interpol',
            valid_type=orm.XyData,
            required=False,
            help='XyData node of the parsed DOS output, interpolated onto the real axis.'
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
        self.ctx.append_text = options_dict.get('append_text', self._options_default['append_text'])
        self.ctx.prepend_text = options_dict.get('prepend_text', self._options_default['prepend_text'])
        self.ctx.additional_retrieve_list = options_dict.get(
            'additional_retrieve_list', self._options_default['additional_retrieve_list']
        )

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
        parent_calc = self.inputs.remote_data.get_incoming(node_class=orm.CalcJobNode).first().node
        if parent_calc.process_label == 'VoronoiCalculation':
            # for the voronoi calculation we need to calculate the Fermi level since it is not in the output parameters directly
            voro_out_para = parent_calc.outputs.output_parameters.get_dict()
            ef = voro_out_para['emin'] - voro_out_para['emin_minus_efermi_Ry']
        else:
            ef = parent_calc.outputs.output_parameters.get_dict()['fermi_energy']  # unit in Ry
        try:
            para_check = set_energy_params(econt_new, ef, para_check)
        except:
            return self.exit_codes.ERROR_DOS_PARAMS_INVALID  # pylint: disable=no-member

        updatenode = orm.Dict(dict=para_check.get_dict())
        updatenode.label = label + 'KKRparam_DOS'
        updatenode.description = descr + \
            'KKR parameter node extracted from parent parameters and wf_parameter input node.'

        paranode_dos = update_params_wf(self.ctx.input_params_KKR, updatenode)

        # maybe overwrite some inputs
        if 'params_kkr_overwrite' in self.inputs:
            self.report(f'found params_kkr_overwrite: {self.inputs.params_kkr_overwrite.get_dict()}')
            updatenode = self.inputs.params_kkr_overwrite
            updatenode.label = 'params overwrite'
            paranode_dos = update_params_wf(paranode_dos, updatenode)

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
        }
        if self.ctx.custom_scheduler_commands:
            options['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        if self.ctx.append_text:
            options['append_text'] = self.ctx.append_text
        if self.ctx.prepend_text:
            options['prepend_text'] = self.ctx.prepend_text
        if self.ctx.additional_retrieve_list:
            options['additional_retrieve_list'] = self.ctx.additional_retrieve_list

        inputs = get_inputs_kkr(
            code, remote, options, label, description, parameters=params, serial=(not self.ctx.withmpi)
        )

        # add nonco angles if found in the parent calculation or in the input
        if 'initial_noco_angles' in self.inputs:
            # overwrite nonco_angles from the input if given
            inputs['initial_noco_angles'] = self.inputs.initial_noco_angles
            self.report('used nonco angles from input to workflow')
        else:
            # extract from the parent calculation
            parent_calc = remote.get_incoming(node_class=orm.CalcJobNode).first().node
            if 'initial_noco_angles' in parent_calc.inputs:
                noco_angles = extract_noco_angles(
                    fix_dir_threshold=orm.Float(1e-6),  # make small enough
                    old_noco_angles=parent_calc.inputs.initial_noco_angles,
                    last_retrieved=parent_calc.outputs.retrieved
                )
                # set nonco angles (either from input or from output if it was updated)
                self.report(noco_angles)
                if noco_angles == {}:
                    noco_angles = parent_calc.inputs.initial_noco_angles
                inputs['initial_noco_angles'] = noco_angles
                self.report(f'extract nonco angles and use from parent ({noco_angles})')

        # LDA+U settings
        if 'settings_LDAU' in self.inputs:
            self.report('Add settings_LDAU input node')
            inputs.settings_LDAU = self.inputs.settings_LDAU

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
            self.report(f'Caught AttributeError {_error}')
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
