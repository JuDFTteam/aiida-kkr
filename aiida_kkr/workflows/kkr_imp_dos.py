#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a impurity DOS calculation and
some helper methods to do so with AiiDA
"""

from aiida.orm import Code, load_node, CalcJobNode, Float, Int, Str, Dict, RemoteData, SinglefileData, XyData, Bool
from aiida.plugins import DataFactory
from aiida.engine import if_, ToContext, WorkChain, calcfunction
from aiida.common import LinkType
from aiida.common.folders import SandboxFolder
from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
from aiida_kkr.workflows.dos import kkr_dos_wc
from aiida_kkr.calculations import KkrimpCalculation
import tarfile
import os
from aiida_kkr.tools.save_output_nodes import create_out_dict_node

__copyright__ = (u'Copyright (c), 2019, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.7.2'
__contributors__ = (u'Fabian Bertoldo', u'Philipp Rüßmann')

# activate verbose output, for debugging only
_VERBOSE_ = True

#TODO: improve workflow output node structure
#TODO: generalise search for imp_info and conv_host from startpot


class kkr_imp_dos_wc(WorkChain):
    """
    Workchain of a DOS calculation for an impurity system starting from a
    converged impurity calculation or workflow

    :param options: (Dict), computer options
    :param wf_parameters: (Dict), specifications for the DOS
    :param kkr: (Code), mandatory: KKR code for gf_writeout step
    :param kkrimp: (Code), mandatory: KKRimp code for DOS calculation
    :param imp_host_pot: (SinglefileData), mandatory: impurity startpotential

    :return workflow_info: (Dict), Information on workflow results
    :return last_calc_output_parameters: (Dict), output parameters of
                                         the last called calculation
    :return last_calc_info: (Dict), information of the last called calculation
    """

    _workflowversion = __version__
    _wf_label = 'kkr_imp_dos_wc'
    _wf_description = 'Workflow for a KKR impurity DOS calculation'

    _options_default = {
        'queue_name': '',  # Queue name to submit jobs too
        'resources': {
            'num_machines': 1
        },  # resources to allocate for the job
        'max_wallclock_seconds': 60 * 60,  # walltime after which the job gets killed (gets parsed to KKR)}
        'custom_scheduler_commands': '',  # some additional scheduler commands
        'withmpi': True
    }  # execute KKR with mpi or without

    _wf_default = {
        'jij_run': False,  # calculate Jij's energy resolved
        'lmdos': False,  # calculate also (l,m) or only l-resolved DOS
        'retrieve_kkrflex': True,  # retrieve kkrflex files to repository or leave on remote computer only
    }
    # add defaults of dos_params since they are passed onto that workflow
    _wf_default['dos_params'] = kkr_dos_wc.get_wf_defaults(silent=True)

    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """

        print(f'Version of workflow: {self._workflowversion}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """

        # Take input of the workflow or use defaults defined above
        super(kkr_imp_dos_wc, cls).define(spec)

        spec.input('kkr', valid_type=Code, required=False, help='KKRhost code, needed if gf_dos_remote is not given.')
        spec.input('kkrimp', valid_type=Code, required=True, help='KKRimp code, always needed.')
        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help='Computer options (resources, quene name, etc.).'
        )
        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._wf_default),
            help='DOS workflow parameters (energy range, etc.).'
        )
        spec.input(
            'host_remote',
            valid_type=RemoteData,
            required=False,
            help='RemoteData node of the (converged) host calculation.'
        )
        spec.input(
            'gf_dos_remote',
            valid_type=RemoteData,
            required=False,
            help='RemoteData node of precomputed host GF for DOS energy contour.'
        )
        spec.input(
            'kkrimp_remote',
            valid_type=RemoteData,
            required=False,
            help='RemoteData node of previous (converged) KKRimp calculation.'
        )
        spec.input(
            'imp_pot_sfd',
            valid_type=SinglefileData,
            required=False,
            help='impurity potential single file data. Needs also impurity_info node.'
        )
        spec.input(
            'impurity_info',
            valid_type=Dict,
            required=False,
            help=
            'impurity info node that specifies the relation between imp_pot_sfd to the host system. Mandatory if imp_pot_sfd is given.'
        )
        spec.input(
            'params_kkr_overwrite',
            valid_type=Dict,
            required=False,
            help='Set some input parameters of the KKR calculation.'
        )
        spec.input(
            'settings_LDAU',
            valid_type=Dict,
            required=False,
            help='Settings for LDA+U run (see KkrimpCalculation for details).'
        )

        spec.expose_inputs(kkr_imp_sub_wc, namespace='BdG', include=('params_overwrite'))
        spec.expose_inputs(kkr_imp_sub_wc, include=('initial_noco_angles', 'rimpshift'))
        spec.expose_inputs(kkr_flex_wc, namespace='gf_writeout', include=('params_kkr_overwrite', 'options'))

        # specify the outputs
        spec.output('workflow_info', valid_type=Dict)
        spec.output('last_calc_output_parameters', valid_type=Dict)
        spec.output('last_calc_info', valid_type=Dict)
        spec.output('dos_data', valid_type=XyData)
        spec.output('dos_data_interpol', valid_type=XyData, required=False)
        spec.output('dos_data_lm', valid_type=XyData, required=False)
        spec.output('dos_data_interpol_lm', valid_type=XyData, required=False)
        spec.output(
            'gf_dos_remote', valid_type=RemoteData, required=False, help='RemoteData node of the computed host GF.'
        )

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,  # start and initialise workflow
            if_(cls.validate_input
                )(  # validate the given input
                    cls.run_gfstep,  # run GF step with DOS energy contour
                    cls.run_imp_dos
                ),  # run DOS for the impurity problem
            cls.return_results  # terminate workflow and return results
        )

        # Define possible exit codes for the workflow
        spec.exit_code(220, 'ERROR_UNKNOWN_PROBLEM', message='Unknown problem detected.')
        spec.exit_code(
            221,
            'ERROR_NO_PARENT_FOUND',
            message='Unable to find the parent remote_data node that led to '
            'the input impurity calculation. You need to specify '
            '`host_remote` and `impurity_info` nodes.'
        )
        spec.exit_code(
            222,
            'ERROR_GF_WRITEOUT_UNSUCCESFUL',
            message='The gf_writeout workflow was not succesful, cannot continue.'
        )
        spec.exit_code(
            223,
            'ERROR_IMP_POT_AND_REMOTE',
            message='The input nodes `imp_pot_sfd` and `kkrimp_remote` are given but are mutually exclusive'
        )
        spec.exit_code(
            224,
            'ERROR_KKR_CODE_MISSING',
            message='KKRhost code node (`inputs.kkr`) is missing if gf_dos_remote is not given.'
        )
        spec.exit_code(
            225, 'ERROR_HOST_REMOTE_MISSING', message='`host_remote` node is missing if gf_dos_remote is not given.'
        )
        spec.exit_code(226, 'ERROR_IMP_SUB_WORKFLOW_FAILURE', message='KKRimp sub-workflow failed.')

    def start(self):
        """
        Initialise context and some parameters
        """

        self.report(f'INFO: started KKR impurity DOS workflow version {self._workflowversion}')
        if _VERBOSE_:
            self.report(f'inputs: {self.inputs}')

        # input both wf and options parameters
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
            if wf_dict == {}:
                wf_dict = self._wf_default
                self.report('INFO: using default wf parameters')

        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
            if options_dict == {}:
                options_dict = self._options_default
                self.report('INFO: using default options parameters')

        # set values, or defaults
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds', self._options_default['max_wallclock_seconds']
        )
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get(
            'custom_scheduler_commands', self._options_default['custom_scheduler_commands']
        )
        self.ctx.options_params_dict = Dict({
            'withmpi': self.ctx.withmpi,
            'resources': self.ctx.resources,
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'queue_name': self.ctx.queue,
            'custom_scheduler_commands': self.ctx.custom_scheduler_commands
        })

        # set workflow parameters for the KKR imputrity calculations
        self.ctx.lmdos = wf_dict.get('lmdos', self._wf_default['lmdos'])
        self.ctx.retrieve_kkrflex = wf_dict.get('retrieve_kkrflex', self._wf_default['retrieve_kkrflex'])

        self.ctx.dos_params_dict = wf_dict.get('dos_params', self._wf_default['dos_params'])
        # fill missing key, value pairs with defaults
        for k, v in self._wf_default['dos_params'].items():
            if k not in self.ctx.dos_params_dict.keys():
                self.ctx.dos_params_dict[k] = v

        # set workflow parameters for the KKR impurity calculation
        self.ctx.jij_run = wf_dict.get('jij_run', self._wf_default['jij_run'])

        # set workflow label and description
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)

        # whether or not to compute the GF writeout step
        self.ctx.skip_gfstep = False

        message = f"""
INFO: use the following parameter:
withmpi: {self.ctx.withmpi}
Resources: {self.ctx.resources}
Walltime (s): {self.ctx.max_wallclock_seconds}
queue name: {self.ctx.queue}
scheduler command: {self.ctx.custom_scheduler_commands}
description: {self.ctx.description_wf}
label: {self.ctx.label_wf}
                  """
        print(message)
        self.report(message)

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''

    def validate_input(self):
        """
        Validate input and catch possible errors
        """

        inputs = self.inputs
        inputs_ok = True
        gf_writeout_calc = None

        if 'imp_pot_sfd' in inputs:
            # check if input potential has incoming return link
            if len(inputs.imp_pot_sfd.get_incoming(link_type=LinkType.RETURN).all()) == 0:
                inputs_ok = self._imp_pot_not_from_wf(inputs_ok)
            else:
                gf_writeout_calc = self._imp_pot_from_wf()

        if 'gf_dos_remote' in self.inputs:
            self.ctx.skip_gfstep = True
        else:
            inputs_ok = self._check_gf_writeout_inputs(inputs_ok, gf_writeout_calc)

        if 'imp_pot_sfd' in self.inputs and 'kkrimp_remote' in self.inputs:
            self.report('[ERROR] both `imp_pot_sfd` and `kkrimp_remote` node in inputs')
            inputs_ok = False
            self.ctx.errors.append(2)  # raises ERROR_IMP_POT_AND_REMOTE
        elif 'imp_pot_sfd' in self.inputs:
            self.report('[INFO] use `imp_pot_sfd` input node')
        elif 'kkrimp_remote' in self.inputs:
            self.report('[INFO] use `kkrimp_remote` input node')
            # extract imp_info node from parent KKRimp calculation
            parent_impcalc = self.inputs.kkrimp_remote.get_incoming(node_class=CalcJobNode).first().node
            self.ctx.imp_info = parent_impcalc.inputs.impurity_info
        else:
            self.report('neither `imp_pot_sfd` nor `kkrimp_remote` node in inputs')
            inputs_ok = False
            self.ctx.errors.append(1)  # raises ERROR_NO_PARENT_FOUND

        message = f'INFO: validated input successfully: {inputs_ok}'
        self.report(message)

        return inputs_ok

    def _imp_pot_not_from_wf(self, inputs_ok):
        """
        input impurity potential is not from a kkrimp workflow
        this means we need also have the `impurity_info` and `host_remote` input nodes
        to be able to run the calculation
        """
        self.report('input potential not from kkrimp workflow: take remote_data folder of host system from input')
        inputs = self.inputs
        if 'impurity_info' in inputs and 'host_remote' in inputs:
            self.ctx.imp_info = inputs.impurity_info
            self.ctx.conv_host_remote = inputs.host_remote
        else:
            message = 'WARNING: startpot has no parent and can not find a converged host RemoteData node'
            self.report(message)
            if 'impurity_info' not in inputs:
                message = '`impurity_info` optional input node not given but needed in this case.'
                self.report(message)
            if 'host_remote' not in inputs:
                message = '`host_remote` optional input node not given but needed in this case.'
                self.report(message)
            inputs_ok = False
            self.ctx.errors.append(1)

        return inputs_ok

    def _imp_pot_from_wf(self):
        """
        we try to extract the input nodes for the GF writeout step from the incoming links to the impurity potential
        """
        # if return link is found get input nodes automatically
        message = 'INFO: get converged host RemoteData node and impurity_info node from database'
        self.report(message)
        self.ctx.kkr_imp_wf = self.inputs.imp_pot_sfd.get_incoming(node_class=kkr_imp_sub_wc).first().node
        message = f'INFO: found underlying kkr impurity workflow (pk: {self.ctx.kkr_imp_wf.pk})'
        self.report(message)
        self.ctx.imp_info = self.ctx.kkr_imp_wf.inputs.impurity_info
        message = f'INFO: found impurity_info node (pk: {self.ctx.imp_info.pk})'
        self.report(message)

        if 'remote_data' in self.ctx.kkr_imp_wf.inputs:
            remote_data_gf_writeout = self.ctx.kkr_imp_wf.inputs.remote_data
            gf_writeout_calc = remote_data_gf_writeout.get_incoming(node_class=CalcJobNode).first().node
            self.ctx.conv_host_remote = gf_writeout_calc.inputs.parent_folder
            message = 'INFO: imported converged_host_remote (pk: {}) and impurity_info from database'.format(
                self.ctx.conv_host_remote.pk
            )
            self.report(message)
        else:
            # follow links in DB upwards to get the remote_folder of the parent calculation
            self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inputs.gf_remote.inputs.remote_folder.inputs.parent_calc_folder.inputs.remote_folder.outputs.remote_folder
            message = 'INFO: imported converged_host_remote (pk: {}) and impurity_info from database'.format(
                self.ctx.conv_host_remote.pk
            )
            self.report(message)

        if 'host_remote' in self.inputs:
            self.ctx.conv_host_remote = self.inputs.host_remote

        return gf_writeout_calc

    def _check_gf_writeout_inputs(self, inputs_ok, gf_writeout_calc):
        """
        check if all necessary inputs are there for the GF writeout step
        """
        if 'kkr' not in self.inputs:
            self.report('[ERROR] `kkr` input node needed if `gf_dos_remote` is not given')
            inputs_ok = False
            self.ctx.errors.append(3)  # raises ERROR_KKR_CODE_MISSING

        if gf_writeout_calc is not None:
            self.report('Use extraced host remote')
        elif 'host_remote' not in self.inputs:
            self.report('[ERROR] `host_remote` input node needed if `gf_dos_remote` is not given')
            inputs_ok = False
            self.ctx.errors.append(4)  # raises ERROR_HOST_REMOTE_MISSING
        else:
            self.ctx.conv_host_remote = self.inputs.host_remote

        return inputs_ok

    def run_gfstep(self):
        """
        Start GF writeout step with DOS energy contour
        """

        if not self.ctx.skip_gfstep:
            options = self.ctx.options_params_dict
            kkrcode = self.inputs.kkr
            converged_host_remote = self.ctx.conv_host_remote
            imp_info = self.ctx.imp_info

            wf_params_gf = {
                'dos_run': True,
                'dos_params': self.ctx.dos_params_dict,
            }
            # maybe set kkrflex_retrieve
            if not self.ctx.retrieve_kkrflex:
                wf_params_gf['retrieve_kkrflex'] = self.ctx.retrieve_kkrflex
            # now convert to AiiDA Dict
            wf_params_gf = Dict(wf_params_gf)

            label_gf = 'GF writeout for imp DOS'
            description_gf = 'GF writeout step with energy contour for impurity DOS'

            builder = kkr_flex_wc.get_builder()
            builder.metadata.label = label_gf  # pylint: disable=no-member
            builder.metadata.description = description_gf  # pylint: disable=no-member
            builder.kkr = kkrcode
            builder.options = options
            builder.wf_parameters = wf_params_gf
            builder.remote_data = converged_host_remote
            builder.impurity_info = imp_info
            if 'params_kkr_overwrite' in self.inputs:
                builder.params_kkr_overwrite = self.inputs.params_kkr_overwrite
            if 'gf_writeout' in self.inputs:
                if 'params_kkr_overwrite' in self.inputs.gf_writeout:
                    builder.params_kkr_overwrite = self.inputs.gf_writeout.params_kkr_overwrite
                if 'options' in self.inputs.gf_writeout:
                    builder.options = self.inputs.gf_writeout.options

            future = self.submit(builder)

            message = f'INFO: running GF writeout (pid: {future.pk})'
            self.report(message)

            return ToContext(gf_writeout=future)

    def run_imp_dos(self):
        """
        Use previous GF step to calculate DOS for the impurity problem
        """

        if not self.ctx.skip_gfstep:
            # use computed gf_writeout
            if not self.ctx.gf_writeout.is_finished_ok:
                return self.exit_codes.ERROR_GF_WRITEOUT_UNSUCCESFUL  # pylint: disable=no-member
            gf_writeout_wf = self.ctx.gf_writeout
            gf_writeout_calc = load_node(self.ctx.gf_writeout.outputs.workflow_info.get_dict().get('pk_flexcalc'))
            gf_writeout_remote = gf_writeout_wf.outputs.GF_host_remote
            self.ctx.pk_flexcalc = self.ctx.gf_writeout.outputs.workflow_info.get_dict().get('pk_flexcalc')
        else:
            # use gf_writeout from input
            gf_writeout_remote = self.inputs.gf_dos_remote
            gf_writeout_calc = gf_writeout_remote.get_incoming(node_class=CalcJobNode).first().node
            self.ctx.pk_flexcalc = gf_writeout_calc.pk

        self.ctx.gf_writeout_remote = gf_writeout_remote

        options = self.ctx.options_params_dict
        kkrimpcode = self.inputs.kkrimp
        if 'imp_pot_sfd' in self.inputs:
            # take impurit potential SingleFileData node
            impurity_pot_or_remote = self.inputs.imp_pot_sfd
        elif 'kkrimp_remote' in self.inputs:
            # or from RemoteData node of previous KKRimp calc
            impurity_pot_or_remote = self.inputs.kkrimp_remote
        imps = self.ctx.imp_info
        nspin = gf_writeout_calc.outputs.output_parameters.get_dict().get('nspin')
        self.ctx.nspin = nspin
        self.report(f'nspin: {nspin}')
        self.ctx.kkrimp_params_dict = Dict({
            'nspin': nspin,
            'nsteps': 1,
            'kkr_runmax': 1,
            'dos_run': True,
            'lmdos': self.ctx.lmdos,
            'jij_run': self.ctx.jij_run,
        })
        kkrimp_params = self.ctx.kkrimp_params_dict
        label_imp = 'KKRimp DOS (GF: {}, imp_pot: {}, Zimp: {}, ilayer_cent: {})'.format(
            gf_writeout_calc.pk,
            impurity_pot_or_remote.pk,  # pylint: disable=possibly-used-before-assignment
            imps.get_dict().get('Zimp'),
            imps.get_dict().get('ilayer_center')
        )
        description_imp = 'KKRimp DOS run (GF: {}, imp_pot: {}, Zimp: {}, ilayer_cent: {}, R_cut: {})'.format(
            gf_writeout_calc.pk, impurity_pot_or_remote.pk,
            imps.get_dict().get('Zimp'),
            imps.get_dict().get('ilayer_center'),
            imps.get_dict().get('Rcut')
        )

        builder = kkr_imp_sub_wc.get_builder()
        builder.metadata.label = label_imp  # pylint: disable=no-member
        builder.metadata.description = description_imp  # pylint: disable=no-member
        builder.kkrimp = kkrimpcode
        builder.options = options
        builder.wf_parameters = kkrimp_params
        builder.remote_data = gf_writeout_remote
        if 'imp_pot_sfd' in self.inputs:
            self.report('Using impurity potential SingelfilData as input')
            builder.host_imp_startpot = impurity_pot_or_remote
        else:
            self.report('Using KKRimp remote folder as input')
            builder.kkrimp_remote = impurity_pot_or_remote
        builder.impurity_info = imps
        # LDA+U settings
        if 'settings_LDAU' in self.inputs:
            self.report('Add settings_LDAU input node')
            builder.settings_LDAU = self.inputs.settings_LDAU

        if 'params_overwrite' in self.inputs.BdG:
            builder.params_overwrite = self.inputs.BdG.params_overwrite
        if 'initial_noco_angles' in self.inputs:
            builder.initial_noco_angles = self.inputs.initial_noco_angles
        if 'rimpshift' in self.inputs:
            builder.rimpshift = self.inputs.rimpshift

        future = self.submit(builder)

        message = f'INFO: running DOS step for impurity system (pk: {future.pk})'
        print(message)
        self.report(message)

        return ToContext(kkrimp_dos=future)

    def return_results(self):
        """
        Return the results and create all of the output nodes
        """

        if self.ctx.errors != []:
            if 1 in self.ctx.errors:
                return self.exit_codes.ERROR_NO_PARENT_FOUND  # pylint: disable=no-member
            elif 2 in self.ctx.errors:
                return self.exit_codes.ERROR_IMP_POT_AND_REMOTE  # pylint: disable=no-member
            elif 3 in self.ctx.errors:
                return self.exit_codes.ERROR_KKR_CODE_MISSING  # pylint: disable=no-member
            elif 4 in self.ctx.errors:
                return self.exit_codes.ERROR_HOST_REMOTE_MISSING  # pylint: disable=no-member
            else:
                return self.exit_codes.ERROR_UNKNOWN_PROBLEM  # pylint: disable=no-member

        message = 'INFO: creating output nodes for the KKR imp DOS workflow ...'
        print(message)
        self.report(message)

        if not self.ctx.kkrimp_dos.is_finished_ok:
            message = 'ERROR: sub workflow for impurity calculation failed'
            print(message)
            self.report(message)
            return self.exit_codes.ERROR_IMP_SUB_WORKFLOW_FAILURE  # pylint: disable=no-member
        else:
            last_calc_pk = self.ctx.kkrimp_dos.outputs.workflow_info.get_dict().get('last_calc_nodeinfo')['pk']
            last_calc = load_node(last_calc_pk)
            last_calc_output_params = last_calc.outputs.output_parameters
            last_calc_info = self.ctx.kkrimp_dos.outputs.workflow_info

            outputnode_dict = {}
            outputnode_dict['impurity_info'] = self.ctx.imp_info.get_dict()
            outputnode_dict['workflow_name'] = self.__class__.__name__
            outputnode_dict['workflow_version'] = self._workflowversion
            if not self.ctx.skip_gfstep:
                outputnode_dict['used_subworkflows'] = {'gf_writeout': self.ctx.gf_writeout.pk}
            else:
                outputnode_dict['used_subworkflows'] = {}
            outputnode_dict['used_subworkflows']['impurity_dos'] = self.ctx.kkrimp_dos.pk

            # interpol dos file and store to XyData nodes
            dos_extracted, dosXyDatas = self.extract_dos_data(last_calc)
            message = f'INFO: extracted DOS data? {dos_extracted}'
            print(message)
            self.report(message)

            # create results node and link rest of results
            link_nodes = dosXyDatas.copy()
            link_nodes['last_calc_output_parameters'] = last_calc_output_params
            link_nodes['last_calc_remote'] = last_calc.outputs.remote_folder
            outputnode_t = create_out_dict_node(Dict(dict=outputnode_dict), **link_nodes)
            outputnode_t.label = 'kkr_imp_dos_wc_inform'
            outputnode_t.description = 'Contains information for workflow'

            if dos_extracted:
                self.out('dos_data', dosXyDatas['dos_data'])
                if 'dos_data_interpol' in dosXyDatas:
                    self.out('dos_data_interpol', dosXyDatas['dos_data_interpol'])
                if self.ctx.lmdos:
                    self.out('dos_data_lm', dosXyDatas['dos_data_lm'])
                    if 'dos_data_interpol_lm' in dosXyDatas:
                        self.out('dos_data_interpol_lm', dosXyDatas['dos_data_interpol_lm'])

            message = f'INFO: workflow_info node: {outputnode_t.uuid}'
            print(message)
            self.report(message)

            self.out('workflow_info', outputnode_t)
            self.out('last_calc_output_parameters', last_calc_output_params)
            self.out('last_calc_info', last_calc_info)
            self.out('gf_dos_remote', self.ctx.gf_writeout_remote)

            message = 'INFO: created output nodes for KKR imp DOS workflow.'
            print(message)
            self.report(message)
            self.report(
                '\n'
                '|------------------------------------------------------------------------------------------------------------------|\n'
                '|-------------------------------------| Done with the KKR imp DOS workflow! |--------------------------------------|\n'
                '|------------------------------------------------------------------------------------------------------------------|'
            )

    def extract_dos_data(self, last_calc):
        """
        Extract DOS data from retrieved folder of KKRimp calculation.
        If output is compressed in tarfile take care of extracting this before parsing `out_ldos*` files.

        The extraction of the DOS data is done in `self.extract_dos_data_from_folder()` calls.

        returns:
          :boolean dos_extracted: logical signalling if extraction of DOS files was successful
          :dictionary dosXyDatas: dictionary containing dos data and interpolated dos data
        """

        # here we look for the dos files or the tarball containing the dos files:
        dos_retrieved = last_calc.outputs.retrieved

        if KkrimpCalculation._FILENAME_TAR in dos_retrieved.list_object_names():
            # deal with packed output files (overwrites dos_retrieved with sandbox into which files are extracted
            # this way the unpacked files are deleted after parsing and only the tarball is kept in the retrieved directory

            # for this we create a Sandbox folder
            with SandboxFolder() as tempfolder:
                # get abs paths of tempfolder and tarfile (needed by extract method of tarfile)

                # get abs path of tempfolder
                with tempfolder.open('.dummy', 'w') as tmpfile:
                    tempfolder_path = os.path.dirname(tmpfile.name)
                # get path of tarfile
                with dos_retrieved.open(KkrimpCalculation._FILENAME_TAR) as tf:
                    tfpath = tf.name

                # extract file from tarfile of retrieved to tempfolder
                with tarfile.open(tfpath) as tf:
                    tar_filenames = [ifile.name for ifile in tf.getmembers()]
                    for filename in tar_filenames:
                        if 'dos' in filename:  # should extract all out_ldos*, out_lmdos* files
                            tf.extract(filename, tempfolder_path)  # extract to tempfolder

                # now files are in tempfolder from where we can extract the dos data
                dos_extracted, dosXyDatas = self.extract_dos_data_from_folder(tempfolder, last_calc)
        else:
            # extract directly from retrieved (no tarball there)
            dos_extracted, dosXyDatas = self.extract_dos_data_from_folder(dos_retrieved, last_calc)

        return dos_extracted, dosXyDatas

    def extract_dos_data_from_folder(self, folder, last_calc):
        """
        Get DOS data and parse files.
        """

        # initialize in case dos data is not extracted
        dosXyDatas = None

        # get list of files in directory (needed since SandboxFolder does not have `list_object_names` method)
        # also extract absolute path of folder (needed by parse_impdosfiles since calcfunction does not work with SandboxFolder as input)
        if isinstance(folder, SandboxFolder):
            folder_abspath = folder.abspath
            filelist = os.listdir(folder_abspath)
        else:
            filelist = folder.list_object_names()

        # check if out_ldos* files are there and parse dos files
        if 'out_ldos.interpol.atom=01_spin1.dat' in filelist or 'out_ldos.interpol.atom=001_spin1.dat' in filelist:
            # extract EF and number of atoms from kkrflex_writeout calculation
            kkrflex_writeout = load_node(self.ctx.pk_flexcalc)
            parent_calc_kkr_converged = kkrflex_writeout.inputs.parent_folder.get_incoming(node_class=CalcJobNode
                                                                                           ).first().node
            ef = parent_calc_kkr_converged.outputs.output_parameters.get_dict().get('fermi_energy')
            last_calc_output_params = last_calc.outputs.output_parameters
            natom = last_calc_output_params.get_dict().get('number_of_atoms_in_unit_cell')
            # parse dosfiles using nspin, EF and Natom inputs
            dosXyDatas = parse_impdosfiles(folder, Int(natom), Int(self.ctx.nspin), Float(ef), Bool(False))
            if self.ctx.lmdos:
                dosXyDatas2 = parse_impdosfiles(folder, Int(natom), Int(self.ctx.nspin), Float(ef), Bool(True))
                dosXyDatas['dos_data_lm'] = dosXyDatas2['dos_data']
                dosXyDatas['dos_data_interpol_lm'] = dosXyDatas2['dos_data_interpol']
            dos_extracted = True
        else:
            dos_extracted = False
            dosXyDatas = None

        return dos_extracted, dosXyDatas


@calcfunction
def parse_impdosfiles(folder, natom, nspin, ef, use_lmdos):
    """
    Read `out_ldos*` files and create XyData node with l-resolved DOS (+node for interpolated DOS if files are found)

    Inputs:
    :param dos_abspath: absolute path to folder where `out_ldos*` files reside (AiiDA Str object)
    :param natom: number of atoms (AiiDA Int object)
    :param nspin: number of spin channels (AiiDA Int object)
    :param ef: Fermi energy in Ry units (AiiDA Float object)

    Returns:
    output dictionary containing
      output = {'dos_data': dosnode, 'dos_data_interpol': dosnode2}
    where `dosnode` and `dosnode2` are AiiDA XyData objects
    """
    from masci_tools.io.common_functions import get_Ry2eV, get_ef_from_potfile
    from numpy import loadtxt, array

    # read dos files
    dos, dos_int = [], []
    for iatom in range(1, natom.value + 1):
        for ispin in range(1, nspin.value + 1):
            name0 = 'out_ldos'
            if use_lmdos.value:
                name0 = 'out_lmdos'
            try:
                # old file names with 2 digits for atom index
                with folder.open(name0 + '.atom=%0.2i_spin%i.dat' % (iatom, ispin)) as dosfile:
                    tmp = loadtxt(dosfile)
                    dos.append(tmp)

                with folder.open(name0 + '.interpol.atom=%0.2i_spin%i.dat' % (iatom, ispin)) as dosfile:
                    tmp = loadtxt(dosfile)
                    if len(tmp) > 0:
                        dos_int.append(tmp)

            except:
                # new file names with 3 digits for atom numbers
                with folder.open(name0 + '.atom=%0.3i_spin%i.dat' % (iatom, ispin)) as dosfile:
                    tmp = loadtxt(dosfile)
                    dos.append(tmp)

                with folder.open(name0 + '.interpol.atom=%0.3i_spin%i.dat' % (iatom, ispin)) as dosfile:
                    try:
                        tmp = loadtxt(dosfile)
                    except:
                        pass
                    # can happen if there are too few points to interpolate on
                    if len(tmp) > 0:
                        dos_int.append(tmp)

    dos = array(dos)
    if len(dos_int) > 0:
        dos_int = array(dos_int)
    else:
        pass

    # convert to eV units
    eVscale = get_Ry2eV()
    dos[:, :, 0] = (dos[:, :, 0] - ef.value) * eVscale
    dos[:, :, 1:] = dos[:, :, 1:] / eVscale
    if len(dos_int) > 0:
        try:
            dos_int[:, :, 0] = (dos_int[:, :, 0] - ef.value) * eVscale
            dos_int[:, :, 1:] = dos_int[:, :, 1:] / eVscale
        except:
            message = 'Not enough data were present in the interpolated dos, due to a lack of energy points to interpolate on'
    else:
        pass

    # create output nodes
    dosnode = XyData()
    dosnode.label = 'dos_data'
    dosnode.description = 'Array data containing uniterpolated DOS (i.e. dos at finite imaginary part of energy). 3D array with (atoms, energy point, l-channel) dimensions.'
    dosnode.set_x(dos[:, :, 0], 'E-EF', 'eV')

    name = ['tot', 's', 'p', 'd', 'f', 'g']
    name = name[:len(dos[0, 0, 1:]) - 1] + ['ns']

    if use_lmdos.value:
        name = [
            's', 'p1', 'p2', 'p3', 'd1', 'd2', 'd3', 'd4', 'd5', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'g1', 'g2',
            'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9'
        ]
        name = name[:len(dos[0, 0, 1:])]

    ylists = [[], [], []]
    for l in range(len(name)):
        ylists[0].append(dos[:, :, 1 + l])
        ylists[1].append('dos ' + name[l])
        ylists[2].append('states/eV')
    dosnode.set_y(ylists[0], ylists[1], ylists[2])

    # node for interpolated DOS
    if len(dos_int) > 0:
        dosnode2 = XyData()
        dosnode2.label = 'dos_interpol_data'
        dosnode2.description = 'Array data containing iterpolated DOS (i.e. dos at finite imaginary part of energy). 3D array with (atoms, energy point, l-channel) dimensions.'
        dosnode2.set_x(dos_int[:, :, 0], 'E-EF', 'eV')
        ylists = [[], [], []]
        for l in range(len(name)):
            ylists[0].append(dos_int[:, :, 1 + l])
            ylists[1].append('interpolated dos ' + name[l])
            ylists[2].append('states/eV')
        dosnode2.set_y(ylists[0], ylists[1], ylists[2])

        output = {'dos_data': dosnode, 'dos_data_interpol': dosnode2}
    else:
        output = {'dos_data': dosnode}

    return output
