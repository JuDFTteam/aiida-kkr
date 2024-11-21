#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a dos calculation and
some helper methods to do so with AiiDA
"""
import numpy as np
from masci_tools.io.kkr_params import kkrparams
from masci_tools.io.common_functions import get_ef_from_potfile, get_Ry2eV
from masci_tools.io.common_functions import get_alat_from_bravais
from aiida import orm
from aiida.engine import WorkChain, while_, if_, ToContext, calcfunction
from aiida.common.exceptions import NotExistent
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida_kkr.workflows.dos import kkr_dos_wc
from aiida_kkr.tools import find_cluster_radius
from aiida_kkr.tools.common_workfunctions import (
    test_and_get_codenode, update_params, update_params_wf, get_inputs_voronoi
)
from aiida_kkr.tools.save_output_nodes import create_out_dict_node

__copyright__ = (
    u'Copyright (c), 2017-2018, Forschungszentrum Jülich GmbH, '
    'IAS-1/PGI-1, Germany. All rights reserved.'
)
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.13.3'
__contributors__ = u'Philipp Rüßmann'

eV2Ry = 1.0 / get_Ry2eV()


class kkr_startpot_wc(WorkChain):
    """
    Workchain  create starting potential for a KKR calculation by running
    voronoi and getting the starting DOS for first checks on the validity of the input setting.
    Starts from a structure together with a KKR parameter node.

    :param wf_parameters: (Dict), Workchain specifications
    :param options: (Dict), specifications for the computer
    :param structure: (StructureData), aiida structure node to begin
        calculation from (needs to contain vacancies, if KKR needs empty spheres)
    :param kkr: (Code)
    :param voronoi: (Code)
    :param calc_parameters: (Dict), KKR parameter set, passed on to voronoi run.

    :return result_kkr_startpot_wc: (Dict), Information of workflow results
        like Success, last result node, dos array data
    """

    _workflowversion = __version__
    _wf_default = {
        'num_rerun': 4,  # number of times voronoi+starting dos+checks is rerun to ensure non-negative DOS etc
        'fac_cls_increase':
        1.15,  # alat         # factor by which the screening cluster is increased each iteration (up to num_rerun times)
        'natom_in_cls_min': 79,  # minimum number of atoms in screening cluster
        'delta_e_min': 1.,  # eV                  # minimal distance in DOS contour to emin and emax in eV
        'threshold_dos_zero': 10**-2,  #states/eV #
        'check_dos': False,  # logical to determine if DOS is computed and checked
        'ef_set': None  # set Fermi level of starting potential to this value
    }
    _options_default = {
        'queue_name': '',  # Queue name to submit jobs to
        'resources': {
            'num_machines': 1
        },  # resources to allowcate for the job
        'max_wallclock_seconds': 60 * 60,  # walltime after which the job gets killed (gets parsed to KKR)
        'withmpi': True,  # execute KKR with mpi or without
        'custom_scheduler_commands': '',  # some additional scheduler commands
    }
    # add defaults of dos_params since they are passed onto that workflow
    _wf_default['dos_params'] = kkr_dos_wc.get_wf_defaults(silent=True)

    _wf_label = ''
    _wf_description = ''

    _kkr_default_params = kkrparams.get_KKRcalc_parameter_defaults()

    # intended to guide user interactively in setting up a valid wf_params node
    @classmethod
    def get_wf_defaults(cls, silent=False):
        """Print and return _wf_defaults dictionary.

        Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """
        if not silent:
            print(f'Version of workflow: {cls._workflowversion}')
        return cls._wf_default

    @classmethod
    def define(cls, spec):
        """Defines the outline of the workflow."""
        # Take input of the workflow or use defaults defined above
        super(kkr_startpot_wc, cls).define(spec)
        spec.input(
            'wf_parameters',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._wf_default),
            help='Parameters that control the behavior of the workflow'
        )
        spec.input(
            'options',
            valid_type=orm.Dict,
            required=False,
            default=lambda: orm.Dict(dict=cls._options_default),
            help='Computer options passed onto the calculations'
        )
        spec.input(
            'structure',
            valid_type=orm.StructureData,
            required=False,
            help="""
            Structure for which the starting potential should be constructed,
            not needed if parent_KKR is given (typically used to increase the
            lmax but use the output potential of the parent_KKR as starting
            potential).
            """
        )
        spec.input(
            'kkr',
            valid_type=orm.Code,
            required=False,
            help='Kkr code, only needed only if DOS is calculated.',
        )
        spec.input(
            'voronoi',
            valid_type=orm.Code,
            required=True,
            help='Voronoi code',
        )
        spec.input(
            'calc_parameters',
            valid_type=orm.Dict,
            required=False,
            help="""
            KKR-specific parameters passed onto the VoronoiCalculation (lmax etc.).
            """,
        )
        spec.input(
            'startpot_overwrite',
            valid_type=orm.SinglefileData,
            required=False,
            help="""
            Potential which can be used instead of the output potential
            Voronoi constructs.
            """,
        )
        spec.input(
            'parent_KKR',
            valid_type=orm.RemoteData,
            required=False,
            help="""
            RemoteData node of a KKR calculation which is used to overwrite
            the output potential Voronoi constructs
            (typically used to increase lmax). Cannot be used with a different
            structure in the input since the structure is extracted from the
            parent_KKR.
            """,
        )
        # define output nodes
        spec.output(
            'results_vorostart_wc',
            valid_type=orm.Dict,
            required=True,
            help='',
        )
        spec.output(
            'last_doscal_results',
            valid_type=orm.Dict,
            required=False,
            help='',
        )
        spec.output(
            'last_voronoi_results',
            valid_type=orm.Dict,
            required=False,
            help='',
        )
        spec.output(
            'last_voronoi_remote',
            valid_type=orm.RemoteData,
            required=False,
            help='',
        )
        spec.output(
            'last_params_voronoi',
            valid_type=orm.Dict,
            required=False,
            help='',
        )
        spec.output(
            'last_doscal_dosdata',
            valid_type=orm.XyData,
            required=False,
            help='',
        )
        spec.output(
            'last_doscal_dosdata_interpol',
            valid_type=orm.XyData,
            required=False,
            help='',
        )
        # definition of exit codes if the workflow needs to be terminated
        spec.exit_code(
            201,
            'ERROR_INVALID_KKRCODE',
            message='The code you provided for kkr does not use the plugin kkr.kkr',
        )
        spec.exit_code(
            202,
            'ERROR_INVALID_VORONOICODE',
            message='The code you provided for voronoi does not use the plugin kkr.voro',
        )
        spec.exit_code(
            203,
            'ERROR_VORONOI_FAILED',
            message='Voronoi calculation unsuccessful. Check inputs',
        )
        spec.exit_code(
            204,
            'ERROR_VORONOI_PARSING_FAILED',
            message='Voronoi calculation unsuccessful. Check inputs.',
        )
        spec.exit_code(
            205,
            'ERROR_VORONOI_INVALID_RADII',
            message='Voronoi calculation unsuccessful. Structure inconsistent. Maybe you need empty spheres?',
        )
        spec.exit_code(
            206,
            'ERROR_DOSRUN_FAILED',
            message='DOS run unsuccessful. Check inputs.',
        )
        spec.exit_code(
            207,
            'ERROR_BOTH_STRUCTURE_AND_PARENT_KKR_GIVEN',
            message='Can only take either structure or parent_KKR as input.',
        )
        spec.exit_code(
            208,
            'ERROR_NEITHER_STRUCTURE_NOR_PARENT_KKR_GIVEN',
            message='Need either structure or parent_KKR as input.',
        )

        # Here the structure of the workflow is defined
        spec.outline(
            # initialize workflow and check inputs
            cls.start,
            # check if another iteration is done (in case of either voro_ok, doscheck_ok is False)
            while_(cls.do_iteration_check)(
                # run voronoi calculation
                # TODO: encapsulate this in restarting mechanism (should be a base class of workflows that start calculations)
                # i.e. use base_restart_calc workchain as parent
                cls.run_voronoi,
                # check voronoi output (also sets ctx.voro_ok)
                if_(cls.check_voronoi)(
                    # create starting DOS using dos sub-workflow
                    cls.get_dos,
                    # perform some checks and set ctx.doscheck_ok accordingly
                    cls.check_dos
                )
            ),
            # collect results and return
            cls.return_results,
            cls.error_handler
        )

    def start(self):
        """
        init context and some parameters
        """
        self.report(f'INFO: started VoroStart workflow version {self._workflowversion}')

        # check input port consistency and set structure
        if 'structure' in self.inputs:
            has_struc = True
            # save structure in context
            self.ctx.structure = self.inputs.structure
        else:
            has_struc = False
        if 'parent_KKR' in self.inputs:
            has_parent_KKR = True
            # find structure from parent_KKR
            self.ctx.structure, _ = VoronoiCalculation.find_parent_structure(self.inputs.parent_KKR)
        else:
            has_parent_KKR = False

        if has_struc and has_parent_KKR:
            self.report('ERROR: caught exit code ERROR_BOTH_STRUCTURE_AND_PARENT_KKR_GIVEN')
            return self.exit_codes.ERROR_BOTH_STRUCTURE_AND_PARENT_KKR_GIVEN  # pylint: disable=no-member
        if (not has_struc) and (not has_parent_KKR):
            self.report('ERROR: caught exit code ERROR_NEITHER_STRUCTURE_NOR_PARENT_KKR_GIVEN')
            return self.exit_codes.ERROR_NEITHER_STRUCTURE_NOR_PARENT_KKR_GIVEN  # pylint: disable=no-member

        ####### init    #######

        # internal para /control para
        self.ctx.exit_code = None

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

        # set values, or defaults for computer options
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

        # set DOS parameters
        self.ctx.dos_params_dict = wf_dict.get(
            'dos_params',
            self._wf_default['dos_params'],
        )

        # set label and description of the workflow
        self.ctx.description_wf = self.inputs.get(
            '_description',
            self._wf_description,
        )
        self.ctx.label_wf = self.inputs.get(
            '_label',
            self._wf_label,
        )

        # iterative rerunning parameters
        self.ctx.iter = 0
        self.ctx.Nrerun = wf_dict.get(
            'num_rerun',
            self._wf_default['num_rerun'],
        )

        # initialize checking booleans
        self.ctx.is_starting_iter = True
        self.ctx.doscheck_ok = True
        self.ctx.voro_ok = False
        self.ctx.check_dos = wf_dict.get(
            'check_dos',
            self._wf_default['check_dos'],
        )
        self.ctx.dos_check_fail_reason = None

        # some physical parameters that are reused
        self.ctx.nclsmin = wf_dict.get(
            'natom_in_cls_min',
            self._wf_default['natom_in_cls_min'],
        )
        self.ctx.fac_clsincrease = wf_dict.get(
            'fac_cls_increase',
            self._wf_default['fac_cls_increase'],
        )
        self.ctx.efermi = None

        # find starting cluster radius
        # find cluster radius (in alat units)
        self.ctx.r_cls = self.find_cluster_radius_alat()

        # difference in eV to emin (e_fermi) if emin (emax) are
        # larger (smaller) than emin (e_fermi)
        self.ctx.delta_e = wf_dict.get(
            'delta_e_min',
            self._wf_default['delta_e_min'],
        )
        # threshold for dos comparison (comparison of dos at emin)
        self.ctx.threshold_dos_zero = wf_dict.get(
            'threshold_dos_zero',
            self._wf_default['threshold_dos_zero'],
        )

        # set Fermi level with input value
        self.ctx.ef_set = wf_dict.get(
            'ef_set',
            self._wf_default['ef_set'],
        )

        # TODO add missing info
        # print the inputs
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
            f'Max. number of voronoi reruns: {self.ctx.Nrerun}\n'
            f'factor cluster increase: {self.ctx.fac_clsincrease}\n'
            f'default cluster radius (in alat): {self.ctx.r_cls}\n'
            f'min. number of atoms in screening cls: {self.ctx.nclsmin}\n'
            f'min. dist in DOS contour to emin/emax: {self.ctx.delta_e} eV\n'
            f'threshold where DOS is zero: {self.ctx.threshold_dos_zero} states/eV\n'
        )

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''

        # get kkr and voronoi codes from input
        if self.ctx.check_dos:
            try:
                test_and_get_codenode(self.inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_INVALID_KKRCODE  # pylint: disable=no-member
        try:
            test_and_get_codenode(self.inputs.voronoi, 'kkr.voro', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_INVALID_VORONOICODE  # pylint: disable=no-member

    def run_voronoi(self):
        """Run voronoi calculation with parameters from input."""
        # incerement iteration counter
        self.ctx.iter += 1

        # increase some parameters
        # check if cluster size is actually the reason for failure
        if self.ctx.iter > 1 and self.ctx.dos_check_fail_reason not in [
            'EMIN too high',
        ]:
            self.ctx.r_cls = self.ctx.r_cls * self.ctx.fac_clsincrease

        structure = self.ctx.structure
        self.ctx.formula = structure.get_formula()
        label = f'voronoi calculation step {self.ctx.iter}'
        description = f'{self.ctx.description_wf} voronoi on {self.ctx.formula}'

        voronoicode = self.inputs.voronoi

        # get valid KKR parameters
        if self.ctx.iter > 1:
            # take value from last run to continue
            params = self.ctx.last_params
            first_iter = False
        else:
            # used input or defaults in first iteration
            first_iter = True
            if 'calc_parameters' in self.inputs:
                params = self.inputs.calc_parameters
            else:
                kkrparams_default = kkrparams()
                para_version = self._kkr_default_params[1]
                for key, val in self._kkr_default_params[0].items():
                    kkrparams_default.set_value(key, val, silent=True)
                # create Dict node
                params = orm.Dict(dict=kkrparams_default.get_dict())
                params.label = 'Defaults for KKR parameter node'
                params.description = f'defaults as defined in kkrparams of version {para_version}'
            #  set last_params accordingly (used below for provenance tracking)
            self.ctx.last_params = params

        self.report(f'INFO: input params: {params}')

        # check if RCLUSTZ is set and use setting from wf_parameters instead
        # (calls update_params_wf to keep track of provenance)
        updated_params = False
        update_list = []
        kkr_para = kkrparams()
        for key, val in params.get_dict().items():
            kkr_para.set_value(key, val, silent=True)
        set_vals = kkr_para.get_set_values()
        set_vals = [keyvalpair[0] for keyvalpair in set_vals]
        default_values = kkrparams.get_KKRcalc_parameter_defaults()[0]
        if 'RCLUSTZ' in set_vals:
            rcls_input = params.get_dict()['RCLUSTZ']
            # set r_cls by default or from input in first iteration
            if self.ctx.r_cls < rcls_input and first_iter:
                self.ctx.r_cls = rcls_input
                updated_params = True
                update_list.append('RCLUSTZ')
            elif self.ctx.nclsmin > 0 and self.ctx.r_cls > rcls_input:
                # change rcls with iterations, do this only if nclsmin is >0
                # (0 or negative numbers trigger use of RCLUSTZ)
                updated_params = True
                update_list.append('RCLUSTZ')
        else:
            updated_params = True
            update_list.append('RCLUSTZ')
        # in case of dos check verify that RMAX, GMAX are set and use setting
        # from wf_parameters otherwise
        if 'RMAX' in set_vals:
            update_list.append('RMAX')
            rmax_input = params.get_dict()['RMAX']
        elif self.ctx.check_dos:  # add only if doscheck is done
            updated_params = True
            update_list.append('RMAX')
            rmax_input = default_values.get('RMAX')
        if 'GMAX' in set_vals:
            update_list.append('GMAX')
            gmax_input = params.get_dict()['GMAX']
        elif self.ctx.check_dos:  # add only if doscheck is done
            updated_params = True
            update_list.append('GMAX')
            gmax_input = default_values.get('GMAX')
        # check if any mandatory keys are not set and set them with the default
        # values if missing in input parameters
        for key, value in default_values.items():
            if key not in update_list and key not in set_vals:
                self.report(f'INFO: setting {key} to default value {value}')
                kkr_para.set_value(key, value)

        # check if Fermi lavel should be set with input value
        if self.ctx.ef_set is not None:
            update_list.append('ef_set')

        # check if emin should be changed:
        if self.ctx.iter > 1 and (self.ctx.dos_check_fail_reason == 'EMIN too high'):
            # decrease emin by self.ctx.delta_e for DOS run only
            emin_old = self.ctx.emin
            emin_new = emin_old - self.ctx.delta_e
            self.ctx.dos_params_dict['emin'] = emin_new
            # change emin_new to internal Ry units
            emin_new = (emin_new - self.ctx.efermi) / eV2Ry
            update_list.append('EMIN')

        # store updated nodes (also used via last_params in kkr_scf_wc)
        if updated_params:
            # set values that are updated
            if 'RCLUSTZ' in update_list:
                kkr_para.set_value('RCLUSTZ', self.ctx.r_cls)
                self.report(f'INFO: setting RCLUSTZ to {self.ctx.r_cls}')
            if 'RMAX' in update_list:
                kkr_para.set_value('RMAX', rmax_input)  # pylint: disable=possibly-used-before-assignment
                self.report(f'INFO: setting RMAX to {rmax_input} (needed for DOS check with KKRcode)')
            if 'GMAX' in update_list:
                kkr_para.set_value('GMAX', gmax_input)  # pylint: disable=possibly-used-before-assignment
                self.report(f'INFO: setting GMAX to {gmax_input} (needed for DOS check with KKRcode)')
            if 'EMIN' in update_list:
                kkr_para.set_value('EMIN', emin_new)
                self.report(f'INFO: setting EMIN to {emin_new} (DOS check showed that EMIN is too high)')
            if 'ef_set' in update_list:
                kkr_para.set_value('EF_SET', self.ctx.ef_set)
                self.report(f'INFO: setting Fermi level of stating potential to {self.ctx.ef_set}')

            updatenode = orm.Dict(dict=kkr_para.get_dict())
            updatenode.description = f'changed values: {update_list}'
            if first_iter:
                updatenode.label = 'initial params from wf input'
                # used workfunction for provenance tracking if parameters have been changed
                params = update_params_wf(self.ctx.last_params, updatenode)
                self.ctx.last_params = params
            else:
                updatenode.label = f'updated params: {update_list}'
                # also keep track of last voronoi output if that has been used
                voro_out = self.ctx.voro_calc.outputs.output_parameters
                params = update_voro_input(self.ctx.last_params, updatenode, voro_out)
                self.ctx.last_params = params

        # run voronoi step
        options = {
            'queue_name': self.ctx.queue,
            'resources': self.ctx.resources,
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'custom_scheduler_commands': self.ctx.custom_scheduler_commands
        }

        if 'structure' in self.inputs:
            # normal mode starting from structure
            builder = get_inputs_voronoi(
                voronoicode,
                structure,
                options,
                label,
                description,
                params=params,
            )
        else:
            # parent_KKR mode
            builder = get_inputs_voronoi(
                voronoicode,
                None,
                options,
                label,
                description,
                params=params,
                parent_KKR=self.inputs.parent_KKR,
            )
        if 'startpot_overwrite' in self.inputs:
            builder.potential_overwrite = self.inputs.startpot_overwrite
        self.report(f'INFO: run voronoi step {self.ctx.iter}')
        future = self.submit(builder)

        # return remote_voro (passed to dos calculation as input)
        return ToContext(voro_calc=future)

    def check_voronoi(self):
        """Check voronoi output.

        return True/False if voronoi output is ok/problematic
        if output is problematic try to increase some parameters
        (e.g. cluster radius) and rerun up tp N_rerun_max times
        initializes with returning True
        """

        # check wether or not calculation was taken from cached node
        caching_info = f'INFO: cache_source of voro calc node: {self.ctx.voro_calc.get_cache_source()}'
        print(caching_info)
        self.report(caching_info)

        # do some checks with the voronoi output (finally sets self.ctx.voro_ok)
        self.ctx.voro_ok = True

        # TODO: print some info if voro_calc does not finish ok (need to find the reason for it)

        print('check_voro_calc', self.ctx.voro_calc, self.ctx.voro_calc.is_finished_ok)
        from pprint import pprint
        pprint(self.ctx.voro_calc.attributes)

        # check calculation state (calculation must be completed)
        if not self.ctx.voro_calc.is_finished_ok:
            self.report('ERROR: Voronoi calculation not in FINISHED state')
            self.ctx.voro_ok = False
            print('abort with', self.exit_codes.ERROR_VORONOI_FAILED)  # pylint: disable=no-member
            self.ctx.exit_code = self.exit_codes.ERROR_VORONOI_FAILED  # pylint: disable=no-member
            return False

        # check if parser returned some error
        voro_parser_errors = self.ctx.voro_calc.res.parser_errors
        if voro_parser_errors != []:
            self.report(f'ERROR: Voronoi Parser returned Error(s): {voro_parser_errors}')
            self.ctx.voro_ok = False
            self.ctx.exit_code = self.exit_codes.ERROR_VORONOI_PARSING_FAILED  # pylint: disable=no-member
            return False

        # check self.ctx.nclsmin condition
        clsinfo = self.ctx.voro_calc.res.cluster_info_group
        ncls = clsinfo.pop('number_of_clusters')

        nclsmin_last_calc = 1000
        for icls in range(len(clsinfo['cluster_info_atoms'])):
            tmp_ncls = clsinfo['cluster_info_atoms'][icls]['sites']
            if tmp_ncls < nclsmin_last_calc:
                nclsmin_last_calc = tmp_ncls
        self.report(f'INFO: number of atoms in smallest cluster: {nclsmin_last_calc}')

        if self.ctx.nclsmin > 0 and self.ctx.nclsmin > nclsmin_last_calc or ncls < 1:
            self.report(f'WARNING: minimal cluster smaller than threshold of {self.ctx.nclsmin}')
            self.ctx.voro_ok = False

        # check radii condition
        radii = self.ctx.voro_calc.res.radii_atoms_group
        r_ratio1 = radii[0]['rout_over_dist_nn']
        r_ratio2 = radii[0]['rmt0_over_rout']
        if r_ratio1 >= 100. or r_ratio2 >= 100.:
            self.report(f'ERROR: radii information inconsistent: Rout/dis_NN={r_ratio1}, RMT0/Rout={r_ratio2}')
            self.ctx.voro_ok = False
            self.ctx.exit_code = self.exit_codes.ERROR_VORONOI_INVALID_RADII  # pylint: disable=no-member
            return False

        # fix emin/emax
        # remember: efermi, emin and emax are in internal units (Ry) but delta_e is in eV!
        ret = self.ctx.voro_calc.outputs.retrieved
        if 'potential_overwrite' in self.ctx.voro_calc.inputs:
            potfile_overwrite = self.ctx.voro_calc.inputs.potential_overwrite
            with potfile_overwrite.open(potfile_overwrite.filename) as _file:
                self.ctx.efermi = get_ef_from_potfile(_file)
        else:
            potfile_name = VoronoiCalculation._OUT_POTENTIAL_voronoi
            if 'parent_KKR' in self.ctx.voro_calc.inputs:
                potfile_name = VoronoiCalculation._POTENTIAL_IN_OVERWRITE
            print(ret.list_object_names(), potfile_name)
            with ret.open(potfile_name) as _file:
                self.ctx.efermi = get_ef_from_potfile(_file)

        emin_dos = self.ctx.dos_params_dict['emin'] * eV2Ry + self.ctx.efermi
        if self.ctx.iter == 1:
            self.ctx.emin = self.ctx.voro_calc.res.emin
        emin_out = self.ctx.emin
        self.report(f'INFO: emin dos input: {emin_dos}, emin voronoi output: {emin_out}')
        if emin_out - self.ctx.delta_e * eV2Ry < emin_dos:
            emin_new = emin_out - self.ctx.delta_e * eV2Ry
            emin_new = (emin_new - self.ctx.efermi) / eV2Ry  # convert to relative eV units
            self.ctx.dos_params_dict['emin'] = emin_new
            self.report(f'INFO: set emin for dos automatically'
                        f' to {emin_out - self.ctx.delta_e * eV2Ry} Ry')

        emax = self.ctx.dos_params_dict['emax']
        self.report(f'INFO: emax dos input: {emax}, efermi voronoi output: {self.ctx.efermi}')
        if emax < self.ctx.delta_e:
            self.ctx.dos_params_dict['emax'] = self.ctx.delta_e
            self.report(
                f'INFO: self.ctx.efermi ({self.ctx.efermi} Ry) +'
                f' delta_e ({self.ctx.delta_e * eV2Ry} Ry) larger than'
                f' emax ({emax} Ry).'
                f' Setting automatically to {self.ctx.efermi + self.ctx.delta_e * eV2Ry} Ry'
            )

        #TODO implement other checks?

        self.report(f'INFO: Voronoi check finished with result: {self.ctx.voro_ok}')

        # finally return result of check
        return self.ctx.voro_ok

    def do_iteration_check(self):
        """Check if another iteration should be done."""
        if self.ctx.exit_code is not None:
            # break out of while loop if an error was caught
            return False
        if self.ctx.is_starting_iter:
            # initial iteration (at least one has to be done)
            # reset starting iter flag
            self.ctx.is_starting_iter = False
            return True
        if self.ctx.iter >= self.ctx.Nrerun:
            # check if maximal number of iterations is reached
            return False
        if self.ctx.voro_ok and self.ctx.doscheck_ok:
            # if both steps succeed we are done
            return False
        return True

    def get_dos(self):
        """Call dos sub workflow and pass the input and submit the calculation."""
        if self.ctx.check_dos:
            self.report(f'INFO: Doing DOS calculation in iteration {self.ctx.iter}')
            # take subset of input and prepare parameter node for dos workflow
            options_dict = {
                'queue_name': self.ctx.queue,
                'resources': self.ctx.resources,
                'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
                'withmpi': self.ctx.withmpi,
                'custom_scheduler_commands': self.ctx.custom_scheduler_commands
            }
            options_node = orm.Dict(dict=options_dict)
            options_node.label = 'options'
            wfdospara_node = orm.Dict(dict=self.ctx.dos_params_dict)
            wfdospara_node.label = 'DOS params'
            wfdospara_node.description = 'DOS parameters passed from kkr_startpot_wc input to DOS sub-workflow'

            code = self.inputs.kkr
            print('VORO_CALC:', self.ctx.voro_calc)
            print('VORO_CALC ok?', self.ctx.voro_calc.is_finished_ok)
            print('OUTGOINGS:', self.ctx.voro_calc.get_outgoing().all())
            remote = self.ctx.voro_calc.outputs.remote_folder
            wf_label = 'DOS calculation'
            wf_desc = 'subworkflow of a DOS calculation that perform a singe-shot KKR calc.'

            builder = kkr_dos_wc.get_builder()
            builder.metadata.description = wf_desc  # pylint: disable=no-member
            builder.metadata.label = wf_label  # pylint: disable=no-member
            builder.kkr = code
            builder.wf_parameters = wfdospara_node
            builder.options = options_node
            builder.remote_data = remote

            future = self.submit(builder)

            return ToContext(doscal=future)

    def check_dos(self):
        """Checks if dos of starting potential is ok."""
        dos_ok = True
        self.ctx.dos_check_fail_reason = None

        if self.ctx.check_dos:
            # check parser output
            doscal = self.ctx.doscal
            print('check doscal', doscal)
            from pprint import pprint
            pprint(doscal.attributes)
            try:
                dos_outdict = doscal.outputs.results_wf.get_dict()
                print(dos_outdict)
                if not dos_outdict['successful']:
                    self.report('ERROR: DOS workflow unsuccessful')
                    self.ctx.doscheck_ok = False
                    self.ctx.exit_code = self.exit_codes.ERROR_DOSRUN_FAILED  # pylint: disable=no-member

                if dos_outdict['list_of_errors'] != []:
                    self.report(f'ERROR: DOS wf output contains errors: {dos_outdict["list_of_errors"]}')
                    self.ctx.doscheck_ok = False
                    self.ctx.exit_code = self.exit_codes.ERROR_DOSRUN_FAILED  # pylint: disable=no-member
            except:
                self.ctx.doscheck_ok = False
                self.ctx.exit_code = self.exit_codes.ERROR_DOSRUN_FAILED  # pylint: disable=no-member

            # needed for checks
            emin = self.ctx.voro_calc.res.emin

            # check for negative DOS
            try:
                dosdata = doscal.outputs.dos_data
                natom = len(self.ctx.voro_calc.res.shapes)
                nspin = dos_outdict['nspin']

                ener = dosdata.get_x()[1]  # shape= natom*nspin, nept
                totdos = dosdata.get_y()[0][1]  # shape= natom*nspin, nept

                if len(ener) != nspin * natom:
                    self.report(
                        f'ERROR: DOS output shape does not fit nspin,'
                        f' natom information: len(energies)={len(ener)},'
                        f' natom={natom}, nspin={nspin}'
                    )
                    self.ctx.doscheck_ok = False
                    self.ctx.exit_code = self.exit_codes.ERROR_DOSRUN_FAILED  # pylint: disable=no-member

                # deal with snpin==1 or 2 cases and check negtive DOS
                for iatom in range(natom // nspin):
                    for ispin in range(nspin):
                        x, y = ener[iatom * nspin + ispin], totdos[iatom * nspin + ispin]
                        if nspin == 2 and ispin == 0:
                            y = -y
                        if y.min() < 0:
                            self.report(
                                f'INFO: negative DOS value found in'
                                f' (atom, spin)=({iatom},{ispin}) at'
                                f' iteration {self.ctx.iter}'
                            )
                            dos_ok = False
                            self.ctx.dos_check_fail_reason = 'DOS negative'

                # check starting EMIN
                dosdata_interpol = doscal.outputs.dos_data_interpol

                ener = dosdata_interpol.get_x()[1]  # shape= natom*nspin, nept
                totdos = dosdata_interpol.get_y()[0][1]  # shape= natom*nspin, nept [0] for total DOS
                Ry2eV = get_Ry2eV()

                for iatom in range(natom // nspin):
                    for ispin in range(nspin):
                        x, y = ener[iatom * nspin + ispin], totdos[iatom * nspin + ispin]
                        xrel = abs(x - (emin - self.ctx.efermi) * Ry2eV)
                        mask_emin = np.where(xrel == xrel.min())
                        ymin = abs(y[mask_emin])
                        if ymin > self.ctx.threshold_dos_zero:
                            self.report(f'INFO: DOS at emin not zero! {ymin}>{self.ctx.threshold_dos_zero}')
                            dos_ok = False
                            self.ctx.dos_check_fail_reason = 'EMIN too high'
            except AttributeError:
                dos_ok = False

            # TODO check for semi-core-states

            # TODO check rest of dos_output node if something seems important to check

        # finally set the value in context (needed in do_iteration_check)
        if dos_ok:
            self.ctx.doscheck_ok = True
        else:
            self.ctx.doscheck_ok = False

        return None

    def return_results(self):
        """Return the results of the dos calculations

        This should run through and produce output nodes even if everything failed,
        therefore it only uses results from context.
        """

        # finalchecks before write out
        if not self.ctx.voro_ok or not self.ctx.doscheck_ok:
            self.ctx.successful = False

        self.report('INFO: create vorostart results nodes.')

        # voronoi outputs
        try:
            voro_pk = self.ctx.voro_calc.pk
        except AttributeError:
            voro_pk = None
        try:
            voro_calc = self.ctx.voro_calc.outputs.output_parameters
        except NotExistent:
            self.report(f'ERROR: Results ParameterNode of voronoi (pk={voro_pk}) not found')
            voro_calc = None
        try:
            voro_remote = self.ctx.voro_calc.outputs.remote_folder
        except NotExistent:
            self.report(f'ERROR: RemoteFolderNode of voronoi (pk={voro_pk}) not found')
            voro_remote = None
        try:
            last_params = self.ctx.last_params
        except AttributeError:
            self.report(f'ERROR: Input ParameterNode of voronoi (pk={voro_pk}) not found')
            last_params = None

        # dos calculation outputs
        try:
            doscal = self.ctx.doscal.outputs.results_wf
        except AttributeError:
            self.report('WARNING: Results ParameterNode of DOS calc not found')
            doscal = None
        try:
            dosdata = self.ctx.doscal.outputs.dos_data
        except AttributeError:
            self.report('WARNING: DOS data of DOS calc not found')
            dosdata = None
        try:
            dosdata_interpol = self.ctx.doscal.outputs.dos_data_interpol
        except AttributeError:
            self.report('WARNING: interpolated DOS data of DOS calc not found')
            dosdata_interpol = None

        self.report(f'INFO: last_voro_calc={self.ctx.voro_calc}')
        self.report(f'INFO: voro_results={voro_calc}')
        self.report(f'INFO: voro_remote={voro_remote}')
        self.report(f'INFO: last_params={last_params}')
        try:
            self.report(f'INFO: last doscal={self.ctx.doscal}')
            self.report(f'INFO: doscal_results={doscal}')
            self.report(f'INFO: dosdata={dosdata}')
            self.report(f'INFO: dosdata_interpol={dosdata_interpol}')
        except:
            self.report('INFO: no doscal data')

        # collect output nodes
        link_nodes = {}

        if voro_calc is not None and voro_remote is not None and last_params is not None:
            link_nodes['last_voronoi_results'] = voro_calc
            link_nodes['last_voronoi_remote'] = voro_remote
            link_nodes['last_params_voronoi'] = last_params
        if doscal is not None and dosdata is not None and dosdata_interpol is not None:
            link_nodes['last_doscal_results'] = doscal
            link_nodes['last_doscal_dosdata'] = dosdata
            link_nodes['last_doscal_dosdata_interpol'] = dosdata_interpol

        # create dict to store results of workflow output
        res_node_dict = {}
        res_node_dict['workflow_name'] = self.__class__.__name__
        res_node_dict['workflow_version'] = self._workflowversion
        res_node_dict['successful'] = self.ctx.successful
        res_node_dict['list_of_errors'] = self.ctx.errors
        res_node_dict['withmpi'] = self.ctx.withmpi
        res_node_dict['resources'] = self.ctx.resources
        res_node_dict['max_wallclock_seconds'] = self.ctx.max_wallclock_seconds
        res_node_dict['queue_name'] = self.ctx.queue
        res_node_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        res_node_dict['dos_params'] = self.ctx.dos_params_dict
        res_node_dict['description'] = self.ctx.description_wf
        res_node_dict['label'] = self.ctx.label_wf
        res_node_dict['last_rclustz'] = self.ctx.r_cls
        res_node_dict['min_num_atoms_in_cluster'] = self.ctx.nclsmin
        res_node_dict['factor_rcls_increase'] = self.ctx.fac_clsincrease
        res_node_dict['last_iteration'] = self.ctx.iter
        res_node_dict['max_iterations'] = self.ctx.Nrerun
        res_node_dict['last_voro_ok'] = self.ctx.voro_ok
        res_node_dict['last_dos_ok'] = self.ctx.doscheck_ok
        res_node_dict['starting_fermi_energy'] = self.ctx.efermi

        # create output Dict node and connect to the rest of the output nodes
        res_node = create_out_dict_node(orm.Dict(dict=res_node_dict), **link_nodes)
        res_node.label = 'vorostart_wc_results'
        res_node.description = ''

        # fill output_nodes dict with
        self.out('results_vorostart_wc', res_node)
        for link_label, node in link_nodes.items():
            self.out(link_label, node)

        self.report('INFO: done with kkr_startpot workflow!\n')

    def find_cluster_radius_alat(self):
        """Find an estimate for the cluster radius.

        The radius is chosen such that it comes close to having nclsmin atoms
        in the cluster.
        """

        # some settings
        ncls_target = self.ctx.nclsmin
        structure = self.ctx.structure

        # first find cluster radius in Ang. units
        r_cls_ang = None
        if self.ctx.nclsmin > 0:
            r_cls_ang = find_cluster_radius(structure, ncls_target)[0]
        else:
            if 'calc_parameters' in self.inputs:
                r_cls_ang = self.inputs.calc_parameters.get_dict().get('RCLUSTZ', None)
        if r_cls_ang is None:
            raise ValueError('nclsmin<=0 but RCLUSTZ not in calc_parameters!')

        # get lattic constant
        cell = np.array(structure.cell)
        alat = get_alat_from_bravais(cell, structure.pbc[2])

        # now check whether input parameter contain use_input_alat key
        # and overwrite alat with used input alat if given
        if 'calc_parameters' in self.inputs:
            params = self.inputs.calc_parameters
            use_input_alat = params.get_dict().get('use_input_alat', False)
            alat_input = params.get_dict().get('ALATBASIS', None)
            if use_input_alat and alat_input is not None:
                alat = alat_input

        # now get r_cls in proper alat units
        r_cls_alat = r_cls_ang / alat

        # return cluster radius
        return r_cls_alat

    def error_handler(self):
        """Capture errors raised in validate_input"""
        if self.ctx.exit_code is not None:
            return self.ctx.exit_code


@calcfunction
def update_voro_input(params_old, updatenode, voro_output):
    """Pseudo wf used to keep track of updated parameters in voronoi calculation.

    voro_output only enters as dummy argument for correct connection but logic using this value is done somewhere else.
    """
    dummy = voro_output
    # voro_output is only dummy input to draw connection in graph

    updatenode_dict = updatenode.get_dict()
    new_parameternode = update_params(params_old, nodename=None, nodedesc=None, **updatenode_dict)
    return new_parameternode
