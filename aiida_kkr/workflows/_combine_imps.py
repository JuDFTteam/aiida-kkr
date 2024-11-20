# -*- coding: utf-8 -*-
"""
This module contains the workflow which combines pre-converged two single-impurity calculations to a larger impurity calculation
"""

from aiida.engine import WorkChain, if_, ToContext, calcfunction
from aiida.orm import load_node, Dict, WorkChainNode, Int, RemoteData, Bool, ArrayData, CalcJobNode
from aiida_kkr.calculations import KkrCalculation, KkrimpCalculation
from aiida_kkr.workflows import kkr_imp_sub_wc, kkr_flex_wc, kkr_imp_wc
from aiida_kkr.tools.combine_imps import (
    create_combined_imp_info_cf, combine_potentials_cf, get_zimp, get_host_structure, get_nspin_and_pot,
    combine_settings_ldau
)
from aiida_kkr.tools.save_output_nodes import create_out_dict_node
import tarfile
import numpy as np
from masci_tools.io.common_functions import get_Ry2eV

__copyright__ = (u'Copyright (c), 2020, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.3.5'
__contributors__ = (u'Philipp Rüßmann , Rubel Mozumder, David Antognini Silva')

# activate debug writeout
_debug = True


class combine_imps_wc(WorkChain):
    """
    Workchain that combines 2 converged single-impurity calculations to a bigger impurity,
    reusing the preconverged potentials. This is useful, for example, to study co-doping.

   Inputs:
    :param impurity1_output_node:(Dict), required, output node from singel impurity wc, and should be one of the
                                 following:
                                 * for `kkr_imp_wc`: single_imp_worlfow.outputs.workflow_info
                                 * for `kkr_imp_sub_wc`: single_imp_worlfow.outputs.workflow_info
                                 * for `KkrimpCalculation`: single_imp_worlfow.outputs.output_parameters

    :param impurity2_output_node:(Dict), required, output node from second singel impurity wc, and should be one of
                                 the following:
                                 * for `kkr_imp_wc`: single_imp_worlfow.outputs.workflow_info
                                 * for `kkr_imp_sub_wc`: single_imp_worlfow.outputs.workflow_info
                                 * for `KkrimpCalculation`: single_imp_worlfow.outputs.output_parameters
    :offset_imp2:(Dict), required, offset of the second impurity with respect to the first impurity position.
                 e.g. {'index:0 or 1}, the replacement by the second impurity will take place at the same cell
                        or at the next nearest cell respectively.

    :param scf.kkrimp: (Code), mandatory, KKRimp code needed to submit kkr_imp_wc
    :param scf.wf_parameters: (Dict), optional, KKRimp code needed to submit kkr_imp_sub_wc
    :param scf.options: (Dict), optional, computer options for kkr_imp_sub_wc

    :param host_gf.kkr: (Code), optional, KKR code for submit kkr_flex_wc, needed if remote_data_gf is not given
    :param host_gf.options: (Dict), optional, computer options for kkr_flex_wc
    :param host_gf.params_kkr_overwrite: (Dict), optional, needed for kkr calculation for GF writeout

    :param wf_parameters_overwrite: (Dict), optional, specifications for wf_parameters of kkr_imp_sub_wc as well
                                  as well as wf_parameters of kkr_flex_wc.
    :param gf_host_remote: (RemoteData), optional, remote folder of a previous kkrflex writeout step
                          calculations containing the flexfiles and will be used for combine host GF.

   Returns:
    :return workflow_info: (Dict), Information of workflow results
    :return last_calc_output_parameters: (Dict), link to output parameters of the last called calculation of the
                                        scf kkr_imp_sub_wc.
    :return last_potential: (SingleFileData) link to last output potential of scf kkr_imp_sub_wc step.
    :return last_calc_remote: (RemoteData) link to remote data of last called calculation of the scf step.
    :return remote_data_gf: (RemoteData) link to GF_host_remote of outputs of kkr_flex_wc e.g. gf_writeou
                           step (only present of host GF was generated here).
    :return JijData: (ArrayData) Consists magnetic interaction data among the magnetic impurity atoms,
                    such as vector distance(rx, ry, rz) between atoms, spin interaction magnetude J,
                    Dzyaloshinskii-Moriya vector magnitude, and Dzyaloshinskii-Moriya vector component(Dx, Dy, Dz)
    :return JijInfo :(Dict) Consists description about the JijData.

    """

    _workflowversion = __version__
    _wf_default = {
        'jij_run': False,  # Any kind of addition in _wf_default should be updated into the start() as well.
        'allow_unconverged_inputs': False,  # start calculation only if previous impurity calculations are converged?
    }

    @classmethod
    def get_wf_defaults(cls, silent=False):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create
        set of wf_parameters_overwrite.
        returns _wf_defaults
        """
        if not silent:
            print(f'Version of workflow: {cls._workflowversion}')
        return {
            'global': cls._wf_default.copy(),
            'scf': kkr_imp_sub_wc.get_wf_defaults(),
            'ghost_gf': kkr_flex_wc.get_wf_defaults()
        }

    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """

        # take define from AiiDA base class and extend it then
        super(combine_imps_wc, cls).define(spec)

        # expose these inputs from sub-workflows
        spec.expose_inputs(
            kkr_imp_sub_wc,
            namespace='scf',
            include=('kkrimp', 'options', 'wf_parameters', 'settings_LDAU', 'params_overwrite')
        )
        spec.expose_inputs(
            kkr_flex_wc,
            namespace='host_gf',
            include=(
                'kkr',
                'options',
                'params_kkr_overwrite',
                'wf_parameters',
            ),  # expose only those port which are not set automatically
            namespace_options={
                'required': False,
                'populate_defaults': False
            },  # this makes sure the kkr code input is not needed if gf_host_remote is provided and the entire namespace is omitted
        )

        # mandatory inputs
        spec.input(
            'impurity1_output_node',
            required=True,
            valid_type=Dict,  #TODO make validator for input node to make sure it is the output of kkr_imp_wc
            help="""
Output node of a single impurity calculation. This can be the output of either the `kkr_imp_wc`, `kkr_imp_sub_wc`
workflows or of an `KkrimpCalculation`.

Use these output Dict nodes:
  * for `kkr_imp_wc`: single_imp_workfow.outputs.workflow_info
  * for `kkr_imp_sub_wc`: single_imp_workfow.outputs.workflow_info
  * for `KkrimpCalculation`: single_imp_workfow.outputs.output_parameters
"""
        )

        spec.input(
            'impurity2_output_node',
            required=True,
            valid_type=Dict,
            help=
            'Output node of second single impurity calculation. See help string of `impurity1_output_node` for more details.'
        )

        spec.input(
            'offset_imp2',
            valid_type=Dict,
            required=True,
            help="""Offset of the secon impurity with respect to the first impurity.
Can be given either via the 'vector' or the 'index' keys in the dictionary.
The 'vector' option allows to give the offset vector in cartesian units and
the 'index' option allows to five the offset vector in units of the lattice
vectors of the host system's structure."""
        )
        spec.input(
            'wf_parameters_overwrite',
            valid_type=Dict,
            required=False,
            help='To add or edit wf_parameters in scf namespace and add run optioins, if needed'
        )
        spec.input(
            'gf_host_remote',
            valid_type=RemoteData,
            required=
            False,  #TODO Add validator that makes sure this is not given together with the host_gf sub-workflow namespace
            help="""RemoteData node of pre-calculated host Green function (i.e. with kkr_flex_wc).
If given then the writeout step of the host GF is omitted."""
        )

        # structure of the workflow
        spec.outline(
            cls.start,  # initialize workflow (set things in context and some consistency checks)
            cls.combine_single_single,  # check the inputs given node whether from single impurity wc.
            cls.create_big_cluster,  # Create the big cluster for both single-single and multi-single
            cls.update_params,  # update wf_parameters of kkr_imp_sub_wc, kkr_flex_wc and run_options
            if_(cls.need_gf_run)(  # check if GF was given in input and can be reused
                cls.run_gf_writeout
            ),  # write out the host GF
            cls.check_host_gf,  # update the host GF
            cls.create_big_potential,  # combine preconverged potentials to big one
            cls.run_kkrimp_scf,  # run the kkrimp_sub workflow to converge the host-imp startpot
            if_(cls.run_jij)(  # Check Jij step should run or not
                cls.run_jij_step
            ),  # run jij step
            cls.return_results  # check if the calculation was successful and return the result nodes
        )

        # define the possible exit codes
        spec.exit_code(999, 'ERROR_SOMETHING_WENT_WRONG', message='ERROR: take a look')
        spec.exit_code(
            900,
            'ERROR_HOST_STRUCTURES_INCONSISTENT',
            message='Host structures of impurity 1 and impurity 2 are not identical!'
        )
        spec.exit_code(
            800, 'ERROR_INPUT_NOT_SINGLE_IMP_CALC', message='Impurity input is not a single impurity calculation.'
        )
        spec.exit_code(
            850,
            'ERROR_INPLANE_NEIGHBOR_TOO_SMALL',
            message='i_neighbor_inplane needs to be positive and bigger than 0 for in-plane neighbors'
        )
        spec.exit_code(
            950, 'ERROR_INCONSISTENT_NSPIN_VALUES', message='The impurity calculations have different NSPIN values'
        )
        spec.exit_code(700, 'ERROR_HOST_GF_CALC_FAILED', message='The writeout of the host GF failed')

        #TODO to fix this create_combined_imp_info_cf need to take the different layers into account
        # when the difference vector and the neighbors are created

        # define the outputs of the workflow
        spec.output('workflow_info')
        spec.output('last_calc_output_parameters')
        spec.output('last_potential')
        spec.output('last_calc_remote')
        spec.output('remote_data_gf', required=False)
        spec.output('JijData', required=False)
        spec.output('JijInfo', required=False)

    def start(self):  # pylint: disable=inconsistent-return-statements
        """
        prepare context, run_option, wf_parameter, run_option, wf_parameterss and do some consistency checks
        """
        message = f'INFO: started combine_imps_wc workflow version {self._workflowversion}'
        self.report(message)

        self.ctx.run_options = self._wf_default
        self.ctx.allow_unconverged_inputs = self._wf_default.get('allow_unconverged_inputs', False)
        if 'wf_parameters_overwrite' in self.inputs:
            self.ctx.wf_parameters_overwrite = self.inputs.wf_parameters_overwrite.get_dict()
            self.ctx.allow_unconverged_inputs = self.ctx.wf_parameters_overwrite.get('global', {}).get(
                'allow_unconverged_inputs', False
            )

        exit_code, self.ctx.imp1 = self.get_imp_node_from_input(iimp=1)
        if exit_code != 0:
            return exit_code
        exit_code, self.ctx.imp2 = self.get_imp_node_from_input(iimp=2)
        if exit_code != 0:
            return exit_code
        self.ctx.single_vs_single = True  # It is assumed

        # find and compare host structures for the two imps to make sure the impurities are consistent
        host_structure1 = get_host_structure(self.ctx.imp1)
        host_structure2 = get_host_structure(self.ctx.imp2)
        #TODO this can be relaxed to make sure the same structure is used even if it is not the same node
        if host_structure1.uuid != host_structure2.uuid:
            self.report('host structures inconsistent')
            return self.exit_codes.ERROR_HOST_STRUCTURES_INCONSISTENT  # pylint: disable=maybe-no-member

        # save host structure in context
        self.ctx.host_structure = host_structure1

        # settings for offset between imps
        self.ctx.i_neighbor_inplane = self.inputs.offset_imp2['index']

        # preserve the inputs from scf namespace to context
        self.ctx.scf_kkrimp = self.inputs.scf.kkrimp
        if 'options' in self.inputs.scf:
            self.ctx.scf_options = self.inputs.scf.options
        if 'wf_parameters' in self.inputs.scf:
            self.ctx.scf_wf_parameters = self.inputs.scf.wf_parameters
        # Add some run option here
        self.ctx.jij_option = False
        self.ctx.offset_imp2 = self.inputs.offset_imp2

        # TODO: PRESERVE THE INPUTS FROM host_gf NAMESPACE TO CONTEXT
        # TODO: ALSO EDIT THE RUN_GF_WRITEOUT() FOR THIS CORRESPONDING CHANGES

    # To check the inputs both impurity1_output_node impurity2_output_node from the single kkrimp calc or wf, or combine and single calc or wf.
    def combine_single_single(self):
        """
            This function checks whether the impurity1_output_node and the impurity2_output_node are from (kkr_imp_wc, kkr_imp_wc) or (combine_imps_wc, kkr_imp_wc), and always keeps the combine_imps_wc as the first impurity node for combine cluster.
        """
        single_single = True

        single_imp_1 = True
        single_imp_2 = True

        imp_1 = self.ctx.imp1
        imp_2 = self.ctx.imp2
        # check for the impurity1 whether from single kkr_imp_wc or not
        if imp_1.process_class == KkrimpCalculation:
            Zimp_num_1 = imp_1.inputs.impurity_info.get_dict().get('Zimp')
            if isinstance(Zimp_num_1, list):
                if len(Zimp_num_1) > 1:
                    single_imp_1 = False

        elif imp_1.process_class == kkr_imp_sub_wc:
            combine_wc = imp_1.get_incoming(node_class=combine_imps_wc).all()
            Zimp_num_1 = imp_1.inputs.impurity_info.get_dict().get('Zimp')
            if len(Zimp_num_1) > 1:
                single_imp_1 = False
        elif imp_1.process_class == combine_imps_wc:
            single_imp_1 = False

        # check for the impurity2 whether from single kkr_imp_wc or not
        if imp_2.process_class == KkrimpCalculation:
            Zimp_num_2 = imp_2.inputs.impurity_info.get_dict().get('Zimp')
            if isinstance(Zimp_num_2, list):
                if len(Zimp_num_2) > 1:
                    single_imp_2 = False
        elif imp_2.process_class == kkr_imp_sub_wc:
            combine_wc = imp_2.get_incoming(node_class=combine_imps_wc).all()
            Zimp_num_2 = imp_2.inputs.impurity_info.get_dict().get('Zimp')
            if len(Zimp_num_2) > 1:
                single_imp_2 = False
        elif imp_2.process_class == combine_imps_wc:
            single_imp_2 = False

        if not (single_imp_1 and single_imp_2):
            single_single = False
            if single_imp_2 == False:
                if single_imp_1 == False:
                    self.report(
                        f"ERROR: Both 'impurity1_output_node' {self.inputs.impurity1_output_node} and 'impurity2_output_node {self.inputs.impurity2_output_node} are from combine_imps_wc."
                    )
                    return self.exit_codes.ERROR_SOMETHING_WENT_WRONG  # pylint: disable=no-member
                else:
                    self.ctx.imp1 = imp_2
                    self.ctx.imp2 = imp_1
        # TODO: Delete this print line for created for deguging
        self.ctx.single_vs_single = single_single
        self.report(
            f'DEBUG: The self.combine_single_single() is successfully done. This is the single-sigle: {self.ctx.single_vs_single}'
        )

    def extract_imps_info_exact_cluster(self):
        """
                This method collects the all exist impurity info as in the exact crystal rather than in the crystal centering the first impurity at (0,0,0) position. Returns the imps_info_in_exact_cluster dict.
        """
        if self.ctx.single_vs_single:
            #TODO what is self.ctx.imp1 i.e self.ctx.imp1==combine_imps_wc for single
            imps_info_in_exact_cluster = self.imps_info_exact_cluster_2imps(
                self.ctx.imp1, self.ctx.imp2, self.inputs.offset_imp2
            )
            self.report(
                f'DEBUG: The is the imps_info_in_exact_cluster dict for single-single imp calc: {imps_info_in_exact_cluster}\n'
            )
            return 0, imps_info_in_exact_cluster
        else:
            imp1_input = self.ctx.imp1
            # This 'if clause' to extract the imps_info_in_exact_cluster from  workflow_info of the input impurity node
            if imp1_input.process_class == combine_imps_wc:
                parent_combine_wc = imp1_input
                out_workflow_info = parent_combine_wc.outputs.workflow_info

            elif imp1_input.process_class == KkrimpCalculation:
                kkrimp_sub = imp1_input.get_incoming(node_class=kkr_imp_sub_wc).all()[0].node
                parent_combine_wc = kkrimp_sub.get_incoming(node_class=combine_imps_wc).all()[0].node
                out_workflow_info = imp1_input.outputs.output_parameters

            elif imp1_input.process_class == kkr_imp_sub_wc:
                parent_combine_wc = imp1_input.get_incoming(node_class=combine_imps_wc).all()[0].node
                out_workflow_info = imp1_input.outputs.workflow_info

            imp2_impurity_info = self.ctx.imp2.inputs.impurity_info

            try:
                imps_info_in_exact_cluster = out_workflow_info.get_dict()['imps_info_in_exact_cluster']  # pylint: disable=possibly-used-before-assignment
            except KeyError:
                parent_input_imp1 = parent_combine_wc.inputs.impurity1_output_node  # TODO: rename combine_wc to the parent_combine_wc
                parent_input_imp2 = parent_combine_wc.inputs.impurity2_output_node
                parent_input_offset = parent_combine_wc.inputs.offset_imp2

                exit_code, parent_imp1_wc_or_calc = self.get_imp_node_from_input(impurity_output_node=parent_input_imp1)
                if exit_code != 0:
                    return exit_code, None
                exit_code, parent_imp2_wc_or_calc = self.get_imp_node_from_input(impurity_output_node=parent_input_imp2)
                if exit_code != 0:
                    return exit_code, None

                # Here below imps_info_in_exact_cluster is construct from the inputs impurity1_output_node, and impurity2_output_node as the idea was not calculated in the old version attempt of the combine_imps_wc.

                imps_info_in_exact_cluster = self.imps_info_exact_cluster_2imps(
                    parent_imp1_wc_or_calc, parent_imp2_wc_or_calc, parent_input_offset
                )
            # Now to add the input impurity info and off set of the present combine_imps_wc
            imps_info_in_exact_cluster['offset_imps'].append(self.ctx.offset_imp2['index'])

            Zimp_2 = imp2_impurity_info.get_dict()['Zimp']
            if isinstance(Zimp_2, list):
                imps_info_in_exact_cluster['Zimps'].append(Zimp_2[0])
            else:
                imps_info_in_exact_cluster['Zimps'].append(Zimp_2)

            imps_info_in_exact_cluster['ilayers'].append(imp2_impurity_info.get_dict()['ilayer_center'])
            # TODO: Delete the below print line as it is for deburging
            self.report(f'DEBUG: The is the imps_info_in_exact_cluster dict: {imps_info_in_exact_cluster}\n')
            return 0, imps_info_in_exact_cluster  # return also exit code

    def get_impinfo_from_hostGF(self, imp_calc):
        """
        Extract impurity infor node from the incoming host GF folder
        """
        if 'impurity_info' in imp_calc.inputs:
            return imp_calc.inputs.impurity_info

        if imp_calc.process_class == KkrimpCalculation:
            GF_input = imp_calc.inputs.host_Greenfunction_folder
        elif imp_calc.process_class == kkr_imp_sub_wc:
            GF_input = imp_calc.inputs.remote_data
        elif imp_calc.process_class == kkr_imp_wc:
            GF_input = imp_calc.inputs.remote_data_gf
        else:
            raise ValueError('Unable to get impinfo')

        parent_calc = GF_input.get_incoming(node_class=CalcJobNode).first().node
        impinfo = parent_calc.inputs.impurity_info

        return impinfo

    def imps_info_exact_cluster_2imps(self, single_imp1_wc, single_imp2_wc, offset_imp2):
        """
            This construct a python dict keeping info about two single inpurities with respect to the original host structure e.i. before transforming the center to the first impurity position.
        """
        impinfo1 = self.get_impinfo_from_hostGF(single_imp1_wc)
        impinfo2 = self.get_impinfo_from_hostGF(single_imp2_wc)
        # imp_info_in_exact_cluster keeps the exact data before creating the cluster will help to add more imps later.
        imps_info_in_exact_cluster = {
            'Zimps': [],
            'ilayers': [],
            'offset_imps': [0]
        }  # This number is taking the first imp.
        #offset_imp contains offset_index for imps 2nd, 3rd so on
        zimp_1 = impinfo1.get_dict().get('Zimp')
        zimp_2 = impinfo2.get_dict().get('Zimp')
        if isinstance(zimp_1, list):
            zimp_1 = zimp_1[0]
        if isinstance(zimp_2, list):
            zimp_2 = zimp_2[0]
        imps_info_in_exact_cluster['Zimps'].append(zimp_1)
        imps_info_in_exact_cluster['Zimps'].append(zimp_2)

        imps_info_in_exact_cluster['ilayers'].append(impinfo1.get_dict().get('ilayer_center'))
        imps_info_in_exact_cluster['ilayers'].append(impinfo2.get_dict().get('ilayer_center'))
        imps_info_in_exact_cluster['offset_imps'].append(offset_imp2.get_dict().get('index'))

        return imps_info_in_exact_cluster

    def get_imp_node_from_input(self, impurity_output_node=None, iimp=1):
        """
        extract impurty calculation from impurity_output_node as inputs
        """
        if impurity_output_node == None:
            if iimp == 1:
                imp_out = self.inputs.impurity1_output_node
            else:
                imp_out = self.inputs.impurity2_output_node
        else:
            imp_out = impurity_output_node

        kkrimpcalc_parents = imp_out.get_incoming(node_class=KkrimpCalculation).all()
        if len(kkrimpcalc_parents) > 0:
            imp = kkrimpcalc_parents[0].node
        else:
            inc = imp_out.get_incoming(link_label_filter='workflow_info').all()
            if len(inc) != 1:
                self.report(f'input node of imp {iimp} inconsistent')
                return self.exit_codes.ERROR_INPUT_NODE_INCONSISTENT, None  # pylint: disable=maybe-no-member
            imp = inc[0].node

        # consistency checks of input nodes
        # check if input calc was converged etc.
        if not self._check_input_imp(imp):
            self.report(f'something wrong with imp {iimp}: {imp}')
            return self.exit_codes.ERROR_SOMETHING_WENT_WRONG, None  # pylint: disable=maybe-no-member

        # return exit code and imp calculation
        return 0, imp

    def _check_input_imp(self, imp_calc_or_wf):
        """
        check if input calculation is of the right process_class and if it did converge
        """

        if imp_calc_or_wf.process_class != KkrimpCalculation:
            # if imp_calc_or_wf is not a KkrimpCalculation
            # check if imp_calc_or_wf is kkr_imp_wc, kkr_imp_sub_wc or combine_imps_wc workflow
            if not isinstance(imp_calc_or_wf, WorkChainNode):
                self.report(f'impurity_workflow not a WorkChainNode: {imp_calc_or_wf}')
                return False

            if not (
                imp_calc_or_wf.process_class == kkr_imp_wc or imp_calc_or_wf.process_class == kkr_imp_sub_wc or
                imp_calc_or_wf.process_class == self.__class__
            ):
                self.report(f'impurity_workflow class is wrong: {imp_calc_or_wf}')
                return False

        if (not self.ctx.allow_unconverged_inputs):
            # calculation should be converged
            if imp_calc_or_wf.process_class == KkrimpCalculation:
                if not imp_calc_or_wf.outputs.output_parameters['convergence_group']['calculation_converged']:
                    self.report(f'impurity calculation not converged: {imp_calc_or_wf}')
                    return False
            if imp_calc_or_wf.process_class == kkr_imp_wc:
                if not imp_calc_or_wf.outputs.workflow_info.get_dict().get('converged'):
                    self.report('impurity_workflow not converged')
                    return False
            elif (imp_calc_or_wf.process_class == kkr_imp_sub_wc or imp_calc_or_wf.process_class == self.__class__):
                if not imp_calc_or_wf.outputs.workflow_info.get_dict().get('convergence_reached'):
                    self.report('impurity_workflow not converged')
                    return False

        # all checks passed
        return True

    def create_big_cluster(self):  # pylint: disable=inconsistent-return-statements
        """
        combine imp clusters of the two imps
        """
        imp1 = self.ctx.imp1
        imp2 = self.ctx.imp2
        single_single = self.ctx.single_vs_single
        if single_single:
            if _debug:
                print('DEBUG:', list(imp1.inputs))
            impinfo1 = self.get_impinfo_from_hostGF(imp1)
            impinfo2 = self.get_impinfo_from_hostGF(imp2)
            # impinfo1 = imp1.inputs.impurity_info
        else:
            if imp1.process_class == self.__class__:
                imp1 = imp1.get_outgoing(node_class=kkr_imp_sub_wc).all()[0].node
            impinfo1 = imp1.inputs.impurity_info
            self.report(f'DEBUG: impinfo1 : {impinfo1.get_dict()} .')

            if imp2.process_class == self.__class__:
                imp2 = imp2.get_outgoing(node_class=kkr_imp_sub_wc).all()[0].node
            impinfo2 = imp2.inputs.impurity_info
            self.report(f'DEBUG: impinfo2 : {impinfo2.get_dict()} .')

        host_structure = self.ctx.host_structure
        offset_imp2 = self.inputs.offset_imp2

        self.report('create combined imp_info:')
        self.report(f'host structure: {host_structure}')
        self.report(f'imp info 1: {impinfo1}')
        self.report(f'imp info 2: {impinfo2}')

        exit_code, self.ctx.imps_info_in_exact_cluster = self.extract_imps_info_exact_cluster()
        if exit_code != 0:
            return exit_code

        imps_info_in_exact_cluster = self.ctx.imps_info_in_exact_cluster
        if single_single:
            if offset_imp2.get_dict()['index'] < 0:
                return self.exit_codes.ERROR_INPLANE_NEIGHBOR_TOO_SMALL  # pylint: disable=maybe-no-member
            if 'ilayer_center' in impinfo1.get_dict() and 'ilayer_center' in impinfo2.get_dict():
                if impinfo1['ilayer_center'] == impinfo2['ilayer_center'] and self.inputs.offset_imp2['index'] < 1:
                    return self.exit_codes.ERROR_INPLANE_NEIGHBOR_TOO_SMALL  # pylint: disable=maybe-no-member
        else:
            imp_offset_index = offset_imp2['index']
            imp2_ilayer = impinfo2['ilayer_center']

            for offset, ilayer in zip(
                imps_info_in_exact_cluster['offset_imps'][:-1], imps_info_in_exact_cluster['ilayers'][:-1]
            ):
                if (offset, ilayer) == (imp_offset_index, imp2_ilayer):
                    self.report(
                        f"ERROR: The new impurity is overlaping with the existing impurities. Change the 'ilayer_center': {ilayer} or 'offset_index'{offset}."
                    )
                    return self.exit_codes.ERROR_SOMETHING_WENT_WRONG  # pylint: disable=no-member

        # get zimp of imp1
        if single_single:
            _, is_single_imp = self.get_and_check_zimp_list(impinfo1)
            if not is_single_imp:
                return self.exit_codes.ERROR_INPUT_NOT_SINGLE_IMP_CALC  # pylint: disable=maybe-no-member

        # do the same for imp2
        _, is_single_imp = self.get_and_check_zimp_list(impinfo2)
        if not is_single_imp:
            return self.exit_codes.ERROR_INPUT_NOT_SINGLE_IMP_CALC  # pylint: disable=maybe-no-member
        self.report(f'DEBURG: This is the exact imps: {imps_info_in_exact_cluster}')

        # create combined cluster, offset of second imp is extracted from i_neighbor_inplane
        out_dict = create_combined_imp_info_cf(
            host_structure, impinfo1, impinfo2, offset_imp2, Dict(dict=imps_info_in_exact_cluster), Bool(single_single)
        )

        self.ctx.imp_info_combined = out_dict['imp_info_combined']
        self.ctx.kickout_info = out_dict['kickout_info']
        self.report(
            f"DEBUG: imp_info_combined: {out_dict['imp_info_combined']}\n kickout_info: {out_dict['kickout_info']}"
        )

    def get_and_check_zimp_list(self, impurity_info):
        """
        extract zimp from impurity_info node and check if it is consistent (needs to be single impurity calculation)
        """
        is_single_imp = True

        zimp = get_zimp(impurity_info)

        # check if calculation is single imp calculation
        if len(zimp) != 1:
            is_single_imp = False

        return zimp, is_single_imp

    def need_gf_run(self):
        """
        check if GF was given in input and can be reused (then return Falser which means no gf needs to be calculated)
        """
        if 'gf_host_remote' in self.inputs:
            return False

        return True

    def run_gf_writeout(self):
        """
        Write out the host GF
        """

        # create process builder for gf_writeout workflow
        builder = kkr_flex_wc.get_builder()
        builder.impurity_info = self.ctx.imp_info_combined
        builder.kkr = self.inputs.host_gf.kkr

        if 'wf_parameters' in self.inputs.host_gf:
            self.report(f'set wf_parameters {self.inputs.host_gf.wf_parameters.get_dict()}')
            builder.wf_parameters = self.inputs.host_gf.wf_parameters

        if 'options' in self.inputs.host_gf:
            builder.options = self.inputs.host_gf.options

        if 'params_kkr_overwrite' in self.inputs.host_gf:
            self.report(
                'INFO: using params_kkr_overwrite in host_gf step: {}'.format(
                    self.inputs.host_gf.params_kkr_overwrite.get_dict()
                )
            )
            builder.params_kkr_overwrite = self.inputs.host_gf.params_kkr_overwrite

        # find converged_host_remote input (converged potential of host system)
        gf_writeout_calc = None
        if self.ctx.imp1.process_class == KkrimpCalculation:
            #take gf_writeout directly from input to KkrimpCalculation
            gf_writeout_calc = self.ctx.imp1.inputs.host_Greenfunction_folder.get_incoming(node_class=KkrCalculation
                                                                                           ).first().node
        elif self.ctx.imp1.process_class == kkr_imp_sub_wc:
            imp1_sub = self.ctx.imp1
        else:
            if _debug:
                print('DEBUG:', self.ctx.imp1, list(self.ctx.imp1.inputs))
            imp1_sub = self.ctx.imp1.get_outgoing(node_class=kkr_imp_sub_wc).first().node
        if gf_writeout_calc is None:
            gf_writeout_calc = imp1_sub.inputs.remote_data.get_incoming(node_class=KkrCalculation).first().node  # pylint: disable=possibly-used-before-assignment
        builder.remote_data = gf_writeout_calc.inputs.parent_folder

        # set label and description of the calc
        sub_label = 'GF writeout combined imps'
        sub_description = 'GF writeout sub workflow for combine_imps_wc '
        builder.metadata.label = sub_label  # pylint: disable=no-member
        builder.metadata.description = sub_description  # pylint: disable=no-member

        # now submit the workflow
        future = self.submit(builder)

        self.report(f'INFO: running GF writeout (pk: {future.pk})')

        return ToContext(gf_writeout=future)

    def check_host_gf(self):
        """
        Check if host gf is there
        """
        self.ctx.host_gf_ok = True

        if self.need_gf_run():  # check only if the calculation was run
            if not self.ctx.gf_writeout.is_finished_ok:
                self.ctx.host_gf_ok = False

        #TODO check if input host gf remote is consistent

        if not self.ctx.host_gf_ok:
            return self.exit_codes.ERROR_HOST_GF_CALC_FAILED  # pylint: disable=no-member

    def create_big_potential(self):  # pylint: disable=inconsistent-return-statements
        """
        combine preconverged potentials to big one
        """

        # get data from context
        imp1 = self.ctx.imp1
        imp2 = self.ctx.imp2
        kickout_info = self.ctx.kickout_info

        nspin1, pot_imp1 = get_nspin_and_pot(imp1)
        nspin2, pot_imp2 = get_nspin_and_pot(imp2)

        # check consistency of nspin for the two calculations
        if nspin1 != nspin2:
            return self.exit_codes.ERROR_INCONSISTENT_NSPIN_VALUES  # pylint: disable=maybe-no-member

        # now combine potentials
        output_potential_sfd_node = combine_potentials_cf(kickout_info, pot_imp1, pot_imp2, Int(nspin1))

        self.ctx.combined_potentials = output_potential_sfd_node

    # To collate and combine the wf_parameters_overwrite and scf_wf_parameters
    def update_params(self):
        """
        Update the parameters in scf_wf_parameters according to wf_parameters_overwrite if
        any change occur there and also add the run options.
        """

        scf_wf_parameters = self.ctx.scf_wf_parameters.get_dict()

        if 'wf_parameters' in self.inputs.host_gf:
            wf_parameters_flex = self.inputs.host_gf.wf_parameters.get_dict()
        else:
            wf_parameters_flex = {}
        run_options = self.ctx.run_options
        # Update the scf_wf_parameters from the wf_parameters_overwrite
        if 'wf_parameters_overwrite' in self.inputs:
            wf_parameters_overwrite = self.ctx.wf_parameters_overwrite

            for key, val in wf_parameters_overwrite.items():
                # Update the scf_wf_parameters from wf_parameters_overwrite
                if key in list(scf_wf_parameters.keys()):
                    if wf_parameters_overwrite[key] != scf_wf_parameters[key]:
                        scf_old_val = scf_wf_parameters[key]
                        scf_wf_parameters[key] = val
                        msg = f'INFO: Parameter value of {key} set from {scf_old_val} to {val}'
                        self.report(msg)
                else:
                    scf_wf_parameters[key] = val

        # Update the wf_parameters_flex and run_options from the scf_wf_parameters
        key_list = []
        for key, val in scf_wf_parameters.items():
            if key in list(wf_parameters_flex.keys()) or key in list(run_options.keys()):
                # Here preparing the wf parameters for kkr_flex_wc
                if key in list(wf_parameters_flex.keys()):
                    deflt_val = wf_parameters_flex[key]
                    wf_parameters_flex[key] = scf_wf_parameters.get(key, deflt_val)
                # Here preparing the some run option
                if key in list(run_options.keys()):
                    deflt_val = run_options[key]
                    run_options[key] = scf_wf_parameters.get(key, deflt_val)
                    self.report(f'INFO: Probable run option <{key}> is updated here as <{run_options[key]}>')
                key_list.append(key)
        # Here to remove keys from scf_wf_parameters that are needed only for gf_writeout_step and run_option
        val_list = [scf_wf_parameters.pop(key, None) for key in key_list[:]]

        self.report('INFO: The wf_parameters Dict for kkr_imp_sub_wc is ready.')

        self.ctx.run_options = run_options
        self.ctx.wf_parameters_flex = wf_parameters_flex
        self.ctx.scf_wf_parameters = Dict(scf_wf_parameters)

    def run_kkrimp_scf(self):
        """
        run the kkrimp_sub workflow to converge the host-imp startpot
        """

        # construct process builder for kkrimp scf workflow
        builder = kkr_imp_sub_wc.get_builder()
        builder.metadata.label = 'kkrimp scf combined imps'  # pylint: disable=no-member
        builder.metadata.description = f'scf workflow for combined impurities: {self.ctx.imp1.label}, {self.ctx.imp2.label}'  # pylint: disable=no-member

        # add combined impurity-info and startpot
        builder.impurity_info = self.ctx.imp_info_combined
        builder.host_imp_startpot = self.ctx.combined_potentials

        # add host GF (either calculated or form input)
        if 'gf_host_remote' not in self.inputs:
            gf_remote = self.ctx.gf_writeout.outputs.GF_host_remote
        else:
            gf_remote = self.inputs.gf_host_remote
        builder.remote_data = gf_remote

        # settings from scf namespace
        builder.kkrimp = self.inputs.scf.kkrimp
        if 'options' in self.inputs.scf:
            builder.options = self.inputs.scf.options
        if 'wf_parameters' in self.inputs.scf:
            builder.wf_parameters = self.inputs.scf.wf_parameters
        if 'params_overwrite' in self.inputs.scf:
            builder.params_overwrite = self.inputs.scf.params_overwrite

        # take care of LDA+U settings
        add_ldausettings, settings_LDAU_combined = self.get_ldau_combined()
        self.report(f'add LDA+U settings? {add_ldausettings}')
        if add_ldausettings:
            self.report(f'settings_combined: {settings_LDAU_combined.get_dict()}')
            builder.settings_LDAU = settings_LDAU_combined
        if 'settings_LDAU' in self.inputs.scf:
            self.report('Overwriting settings_LDAU_combined with provided settings_LDAU inputs')
            builder.settings_LDAU = self.inputs.scf.settings_LDAU

        # now submit workflow
        future = self.submit(builder)

        self.report(f'INFO: running kkrimp scf workflow for combined impts (uuid= {future.uuid})')

        return ToContext(kkrimp_scf_sub=future)

    def run_jij(self):
        if not self.ctx.kkrimp_scf_sub.is_finished_ok:
            msg = self.exit_codes.ERROR_SOMETHING_WENT_WRONG  # pylint: disable=no-member
            self.report(msg)
            return False
        # TODO : work here from the RUNOPT instead of scf_wf_parameters, because all the RUNOPT and wf_parameters_flex are updated in the update parameter funcion.
        run_options = self.ctx.run_options
        if 'jij_run' in list(run_options.keys()):
            self.ctx.jij_option = run_options['jij_run']

        return self.ctx.jij_option

    # To launch jij step
    def run_jij_step(self):
        """Run the jij calculation with the converged combined_imp_host_pot"""
        if self.ctx.kkrimp_scf_sub.is_finished_ok:
            msg = 'INFO: kkr_imp_sub_wc for combined impurity cluster is succefully done and jij step is getting prepared for launching.'
            self.report(msg)
            combined_scf_node = self.ctx.kkrimp_scf_sub
        else:
            return self.exit_codes.ERROR_SOMETHING_WENT_WRONG  # pylint: disable=no-member

        last_calc = load_node(combined_scf_node.outputs.workflow_info['last_calc_nodeinfo']['uuid'])
        builder = last_calc.get_builder_restart()

        builder.pop('parent_calc_folder')
        builder.impurity_potential = combined_scf_node.outputs.host_imp_pot
        param_dict = {k: v for k, v in builder.parameters.get_dict().items() if v is not None}
        param_dict['CALCJIJMAT'] = 1  # activate Jij calculation, leave the rest as is

        builder.parameters = Dict(param_dict)
        builder.metadata.label = 'KKRimp_Jij (' + last_calc.label.split('=')[1][3:]

        future = self.submit(builder)
        self.report(f'INFO: Submiting Jij calculation (uuid={future.uuid})')

        return ToContext(imp_scf_combined_jij=future)

    def get_ldau_combined(self):
        """
        check if impurity input calculations have LDA+U settings in input and add this here if needed
        """

        if self.ctx.imp1.process_class == KkrimpCalculation or self.ctx.imp1.process_class == kkr_imp_sub_wc:
            imp1_has_ldau = 'settings_LDAU' in self.ctx.imp1.inputs
            if imp1_has_ldau:
                settings_LDAU1 = self.ctx.imp1.inputs.settings_LDAU
                self.report('found LDA+U settings for imp1')
                if self.ctx.imp1.process_class == KkrimpCalculation:
                    retrieved1 = self.ctx.imp1.outputs.retrieved
                elif self.ctx.imp1.process_class == kkr_imp_sub_wc:
                    retrieved1 = self.ctx.imp1.get_outgoing(node_class=KkrimpCalculation).first().node.outputs.retrieved
        elif self.ctx.imp1.process_class == combine_imps_wc:
            imp1_has_ldau = 'settings_LDAU' in self.ctx.imp1.get_outgoing(node_class=kkr_imp_sub_wc).first().node.inputs
            if imp1_has_ldau:
                settings_LDAU1 = self.ctx.imp1.get_outgoing(node_class=kkr_imp_sub_wc).first().node.inputs.settings_LDAU
                self.report('found LDA+U settings for imp1')
                retrieved1 = self.ctx.imp1.get_outgoing(node_class=kkr_imp_sub_wc
                                                        ).first().node.get_outgoing(node_class=KkrimpCalculation
                                                                                    ).first().node.outputs.retrieved

        if self.ctx.imp2.process_class == KkrimpCalculation or self.ctx.imp2.process_class == kkr_imp_sub_wc:
            imp2_has_ldau = 'settings_LDAU' in self.ctx.imp2.inputs
            if imp2_has_ldau:
                settings_LDAU2 = self.ctx.imp2.inputs.settings_LDAU
                self.report('found LDA+U settings for imp2')
                if self.ctx.imp2.process_class == KkrimpCalculation:
                    retrieved2 = self.ctx.imp2.outputs.retrieved
                elif self.ctx.imp2.process_class == kkr_imp_sub_wc:
                    retrieved2 = self.ctx.imp2.get_outgoing(node_class=KkrimpCalculation).first().node.outputs.retrieved
        elif self.ctx.imp2.process_class == combine_imps_wc:
            imp2_has_ldau = 'settings_LDAU' in self.ctx.imp2.get_outgoing(node_class=kkr_imp_sub_wc).first().node.inputs
            if imp2_has_ldau:
                settings_LDAU2 = self.ctx.imp2.get_outgoing(node_class=kkr_imp_sub_wc).first().node.inputs.settings_LDAU
                self.report('found LDA+U settings for imp2')
                retrieved2 = self.ctx.imp2.get_outgoing(node_class=kkr_imp_sub_wc
                                                        ).first().node.get_outgoing(node_class=KkrimpCalculation
                                                                                    ).first().node.outputs.retrieved

        if imp1_has_ldau and imp2_has_ldau:
            # combine LDA+U settings of the two imps
            settings_LDAU_combined = combine_settings_ldau(
                settings_LDAU1=settings_LDAU1,  # pylint: disable=used-before-assignment
                retrieved1=retrieved1,  # pylint: disable=used-before-assignment
                settings_LDAU2=settings_LDAU2,  # pylint: disable=used-before-assignment
                retrieved2=retrieved2,  # pylint: disable=used-before-assignment
                kickout_info=self.ctx.kickout_info
            )
        elif imp1_has_ldau:
            # use only LDA+U settings of imp 1
            settings_LDAU_combined = combine_settings_ldau(
                settings_LDAU1=settings_LDAU1, retrieved1=retrieved1, kickout_info=self.ctx.kickout_info
            )
        elif imp2_has_ldau:
            # add offset to atom index for combined LDA+U settings
            settings_LDAU_combined = combine_settings_ldau(
                settings_LDAU2=settings_LDAU2, retrieved2=retrieved2, kickout_info=self.ctx.kickout_info
            )
        else:
            # return builder unchanged if none of the impurt calculations has LDA+U settings
            return False, {}

        # now add settings_LDAU input to builder
        self.report(
            f'add combined LDAU settings (uuid={settings_LDAU_combined.uuid}): {settings_LDAU_combined.get_dict()}'
        )

        return True, settings_LDAU_combined

    def return_results(self):
        """
        check if the calculation was successful and return the result nodes
        """

        import numpy as np

        self.report('INFO: Return_results:INFO: Return_results:  To collect the WF info and Other results')

        if not self.ctx.kkrimp_scf_sub.is_finished_ok:
            self.report('ERROR: kkrimp convergence step is not finished successfully')
            return self.exit_codes.ERROR_SOMETHING_WENT_WRONG  # pylint: disable=no-member

        # collect results of kkrimp_scf sub-workflow
        kkrimp_scf_sub = self.ctx.kkrimp_scf_sub
        results_kkrimp_sub = kkrimp_scf_sub.outputs.workflow_info
        last_calc = load_node(results_kkrimp_sub['last_calc_nodeinfo']['uuid'])
        output_parameters = last_calc.outputs.output_parameters
        last_remote = last_calc.outputs.remote_folder
        last_pot = kkrimp_scf_sub.outputs.host_imp_pot
        out_dict = {}
        out_dict['workflow_name'] = self.__class__.__name__
        out_dict['workflow_version'] = self._workflowversion

        # collect info of sub workflow
        out_dict['sub_workflows'] = {'kkrimp_scf': {'pk': kkrimp_scf_sub.pk, 'uuid': kkrimp_scf_sub.uuid}}

        # collect some results from scf sub-workflow
        for key in ['successful', 'convergence_value', 'convergence_reached', 'convergence_values_all_steps']:
            out_dict[key] = results_kkrimp_sub[key]

        try:
            magmom_all = np.array(output_parameters['magnetism_group']['spin_moment_per_atom'], dtype=float)[:, -1]
            out_dict['magmoms'] = magmom_all
        except KeyError:
            self.report("No 'magmons' available in 'workflow_info'. Probably, 'mag_init' set to False in the settings.")

        # Parse_jij and collect some info
        is_jij_exist = self.ctx.jij_option
        if is_jij_exist:
            impurity1_output_node = self.inputs.impurity1_output_node
            impurity2_output_node = self.inputs.impurity2_output_node
            jij_calc = self.ctx.imp_scf_combined_jij
            jij_retrieved = jij_calc.outputs.retrieved
            impurity_info = kkrimp_scf_sub.inputs.impurity_info
            out_dict['run_option_info'] = {
                'jij_calc': {
                    'pk': jij_calc.pk,
                    'uuid': jij_calc.uuid,
                    'is_finished_ok': jij_calc.is_finished_ok
                }
            }
            jij_parsed_dict = parse_Jij(jij_retrieved, impurity_info, impurity1_output_node, impurity2_output_node)
        # collect outputs of host_gf sub_workflow if it was done
        if 'gf_host_remote' not in self.inputs:
            gf_writeout = self.ctx.gf_writeout
            gf_sub_remote = gf_writeout.outputs.GF_host_remote

            # add info about sub-workflow to dict output
            out_dict['sub_workflows']['host_gf'] = {'pk': gf_writeout.pk, 'uuid': gf_writeout.uuid}
            # add as output node
            self.out('remote_data_gf', gf_sub_remote)

        # add information on combined cluster and potential
        out_dict['imp_info_combined'] = self.ctx.imp_info_combined.get_dict()
        out_dict['potential_kickout_info'] = self.ctx.kickout_info.get_dict()
        out_dict['imps_info_in_exact_cluster'] = self.ctx.imps_info_in_exact_cluster

        # create results node with input links
        # TODO: Add the imp_scf_combined_jij to link the run_jij_step output on provanace graph
        link_nodes = {'kkrimp_scf_results': results_kkrimp_sub}
        if 'gf_host_remote' in self.inputs:
            link_nodes['GF_host_remote'] = self.inputs.gf_host_remote
        if is_jij_exist:
            link_nodes['Jij_retrieved'] = jij_retrieved
        outputnode = create_out_dict_node(Dict(dict=out_dict), **link_nodes)
        outputnode.label = 'combine_imps_wc_results'

        # add output nodes
        self.out('workflow_info', outputnode)
        self.out('last_potential', last_pot)
        self.out('last_calc_remote', last_remote)
        self.out('last_calc_output_parameters', output_parameters)
        if is_jij_exist:
            self.out('JijData', jij_parsed_dict['Jijdata'])
            self.out('JijInfo', jij_parsed_dict['info'])


@calcfunction
def parse_Jij(retrieved, impurity_info, impurity1_output_node, impurity2_output_node):
    """parser output of Jij calculation and return as ArrayData node"""

    # Collect the zimp for impurity_output_node
    try:
        imp1_z = impurity1_output_node.get_incoming(node_class=kkr_imp_wc
                                                    ).first().node.inputs.impurity_info.get_dict()['Zimp']
    except AttributeError:
        try:
            impurity1_output_node_combine = impurity1_output_node.get_incoming(node_class=combine_imps_wc).first().node
            imp1_z = impurity1_output_node_combine.get_outgoing(node_class=kkr_imp_sub_wc
                                                                ).first().node.inputs.impurity_info.get_dict()['Zimp']
        except AttributeError:
            try:
                imp1_z = impurity1_output_node.get_incoming(node_class=kkr_imp_sub_wc
                                                            ).first().node.inputs.impurity_info.get_dict()['Zimp']
            except AttributeError:
                imp1_z = impurity1_output_node.get_incoming(node_class=KkrimpCalculation
                                                            ).first().node.inputs.impurity_info.get_dict()['Zimp']

    try:
        imp2_z = impurity2_output_node.get_incoming(node_class=kkr_imp_wc
                                                    ).first().node.inputs.impurity_info.get_dict()['Zimp']
    except AttributeError:
        try:
            impurity2_output_node_combine = impurity2_output_node.get_incoming(node_class=combine_imps_wc).first().node
            imp2_z = impurity2_output_node_combine.get_outgoing(node_class=kkr_imp_sub_wc
                                                                ).first().node.inputs.impurity_info.get_dict()['Zimp']
        except AttributeError:
            try:
                imp2_z = impurity2_output_node.get_incoming(node_class=kkr_imp_sub_wc
                                                            ).first().node.inputs.impurity_info.get_dict()['Zimp']
            except AttributeError:
                imp2_z = impurity2_output_node.get_incoming(node_class=KkrimpCalculation
                                                            ).first().node.inputs.impurity_info.get_dict()['Zimp']

    with retrieved.open('out_Jijmatrix') as _f:
        jijdata = np.loadtxt(_f)

    impurity_info = impurity_info.get_dict()
    pos = np.array(impurity_info['imp_cls'])
    z = np.array(impurity_info['imp_cls'])[:, 4]

    if not isinstance(imp1_z, list):
        Vpos = list(np.where(z == imp1_z)[0]) + list(np.where(z == imp2_z)[0])
    else:
        Vpos = [np.where(z == i)[0][0] for i in imp1_z] + list(np.where(z == imp2_z)[0])
    Ry2eV = get_Ry2eV()

    # extract number of atoms
    natom = int(np.sqrt(jijdata.shape[0] / 3 / 3))

    # reshape data
    jij_reshape = jijdata.reshape(3, natom, natom, 3, 3)  # iter, i, j, k, l (Jij_k,l matrix)

    # now combine iterations to get full 3 by 3 Jij matrices for all atom pairs
    jij_combined_iter = np.zeros((natom, natom, 3, 3))
    for iatom in range(natom - 1):
        for jatom in range(natom)[iatom + 1:]:
            for iiter in range(3):
                if iiter == 0:
                    # first iteration with theta, phi = 0, 0
                    # take complete upper block from here since this calculation should be converged best
                    # (rotated moments only one-shot calculations)
                    jij_combined_iter[iatom, jatom, 0, 0] = jij_reshape[iiter, iatom, jatom, 0, 0]
                    jij_combined_iter[iatom, jatom, 0, 1] = jij_reshape[iiter, iatom, jatom, 0, 1]
                    jij_combined_iter[iatom, jatom, 1, 0] = jij_reshape[iiter, iatom, jatom, 1, 0]
                    jij_combined_iter[iatom, jatom, 1, 1] = jij_reshape[iiter, iatom, jatom, 1, 1]
                elif iiter == 1:
                    # second iteraton with theta, phi = 90, 0
                    jij_combined_iter[iatom, jatom, 1, 2] = jij_reshape[iiter, iatom, jatom, 1, 2]
                    jij_combined_iter[iatom, jatom, 2, 1] = jij_reshape[iiter, iatom, jatom, 2, 1]
                    jij_combined_iter[iatom, jatom, 2, 2] = jij_reshape[iiter, iatom, jatom, 2, 2]
                else:
                    # from third iteration with theta, phi = 90, 90
                    jij_combined_iter[iatom, jatom, 0, 2] = jij_reshape[iiter, iatom, jatom, 0, 2]
                    jij_combined_iter[iatom, jatom, 2, 0] = jij_reshape[iiter, iatom, jatom, 2, 0]
                    # add this value to z-z component and average
                    jij_combined_iter[iatom, jatom, 2, 2] += jij_reshape[iiter, iatom, jatom, 2, 2]
                    jij_combined_iter[iatom, jatom, 2, 2] *= 0.5

    # finally convert to meV units (and sign change to have positive number indicate ferromagnetism and negative number antiferromagnetism)
    jij_combined_iter *= -1. * Ry2eV * 1000

    jij_trace = (jij_combined_iter[:, :, 0, 0] + jij_combined_iter[:, :, 1, 1] + jij_combined_iter[:, :, 2, 2]) / 3
    Dij_vec = np.array([(jij_combined_iter[:, :, 1, 2] - jij_combined_iter[:, :, 2, 1]),
                        (jij_combined_iter[:, :, 2, 0] - jij_combined_iter[:, :, 0, 2]),
                        (jij_combined_iter[:, :, 0, 1] - jij_combined_iter[:, :, 1, 0])])

    plotdata = []

    #return jij_combined_iter
    out_txt = 'Output Jij values between V impurities:\n i          j           Jij (meV)           Dij(meV)            D/J             i_zimp          j_zimp \n---------------------------------------------------------------------------\n'
    for iatom in range(natom - 1):
        for jatom in range(natom)[iatom + 1:]:
            if iatom != jatom and iatom in Vpos and jatom in Vpos:
                J = jij_trace[iatom, jatom]
                Dx, Dy, Dz = Dij_vec[0, iatom, jatom], Dij_vec[1, iatom, jatom], Dij_vec[2, iatom, jatom]
                D = np.sqrt(Dx**2 + Dy**2 + Dz**2)
                out_txt += '%3i         %3i         %15.5e          %15.5e          %15.5e          %4i         %4i\n' % (
                    iatom, jatom, J, D, D / J, pos[iatom][4], pos[jatom][4]
                )
                rdiff = pos[jatom] - pos[iatom]
                plotdata.append([rdiff[0], rdiff[1], rdiff[2], J, D, Dx, Dy, Dz])
    plotdata = np.array(plotdata)

    a = ArrayData()
    a.set_array('JijData', plotdata)

    return {'Jijdata': a, 'info': Dict({'text': out_txt})}
