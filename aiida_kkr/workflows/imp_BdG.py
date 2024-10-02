# Workflow for impurity BdG calculation from converged normal state impurity portential and BdG host calculation

from aiida.engine import WorkChain, ToContext, if_
from aiida.orm import Dict, RemoteData, Code, CalcJobNode, WorkChainNode, Float, Bool, XyData, SinglefileData
from aiida_kkr.workflows import kkr_imp_wc, kkr_imp_dos_wc
from aiida_kkr.tools.find_parent import get_calc_from_remote
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.2'
__contributors__ = (u'David Antognini Silva, Philipp Rüßmann')

# TODO: add _wf_default parameters and activate get_wf_defaults method
# TODO: add input interdependencies in the workchain description


class kkrimp_BdG_wc(WorkChain):
    """
    Workchain for one-shot BdG impurity DOS calculation from a converged normal state host calculation.
    1) Normal state impurity potential scf
    2) One-shot BdG DOS calcluation
        a) Writing out of the BdG kkrflex DOS files
        b) One-shot BdG impurity DOS
    The intermediate steps (1 & 2a) can be skipped by providing the corresponding nodes as inputs to the workflow.

    inputs::
        :param options: (Dict), computer options
        :param impurity_info: (Dict), information about the impurity cluster
        :param voronoi: (Code), Voronoi code for creating the impurity starting potential
        :param kkr: (Code), KKR host code for the writing out kkrflex files
        :param kkrimp: (Code), KKR impurity code for the normal state impurity scf and BdG impurity DOS calculation
        :param BdG_settings: (Dict), set BdG parameters
        :param imp_scf.startpot: (SinglefileData), pre-converged impurity potential used to start the impurity scf calculation
        :param imp_scf.converged_potential: (SinglefileData), converged impurity potential, skips the impurity scf calculation if provided
        :param imp_scf.wf_parameters: (Dict), parameters for the kkr impurity scf
        :param imp_scf.gf_writeout.params_kkr_overwrite: (Dict), set some input parameters of the KKR calculation for the GF writeout step of impurity scf workchain
        :param imp_scf.gf_writeout.options: (Dict), computer settings
        :param imp_scf.scf.params_overwrite: (Dict), set some input parameters of the KKR impurity scf
        :param imp_scf.options: (Dict), computer settings
        :param imp_scf.remote_data_host: (RemoteData), parent folder of converged host normal state KkrCalculation

        :param dos.wf_parameters: (Dict), parameters for the DOS calculation
        :param dos.gf_dos_remote: (RemoteData), node of precomputed host GF for DOS energy contour
        :param dos.gf_writeout.params_kkr_overwrite: (Dict), set some input parameters of the KKR calculation for the GF writeout step of imßpurity dos workchain
        :param dos.gf_writeout.host_remote: (RemoteData), parent folder of kkrflex writeout step for DOS calculation
        :param dos.gf_writeout.kkr: (Code), KKR code for writing out of kkrflex files for impurity DOS calculation
        :param dos.gf_writeout.options: (Dict), computer settings
        :param dos.options: (Dict), computer settings

    returns::
        :return workflow_info: (Dict), Information on workflow results
        :return output_parameters: (Dict), output parameters of the workflow
        :return dos_data: (XyData), impurity DOS data output node
        :return dos_data_interpol: (XyData), interpolated DOS data output node
        :return impurity_potential: (SinglefileData), converged normal state impurity potential node
        :return gf_host_BdG: (RemoteData), kkrflex writeout step files of DOS calculation
    """

    _wf_version = __version__

    # _wf_default = {
    #     'default': None,  #Put in here default input parameters
    # }
    _options_default = {
        'max_wallclock_seconds': 36000,
        'resources': {
            'num_machines': 1
        },
        'withmpi': True,
        'queue_name': ''
    }

    # @classmethod
    # def get_wf_defaults(self, silent=False):
    #     """
    #     Return the default values of the workflow parameters (wf_parameters input node)
    #     """
    #     if not silent:
    #         print(f'Version of the kkrimp_BdG_wc workflow: {self._wf_version}')
    #     return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(kkrimp_BdG_wc, cls).define(spec)

        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help=
            'Computer options (walltime etc.) passed onto KkrCalculation, fall back to settings from parent calculation if not given'
        )

        spec.input(
            'impurity_info',
            valid_type=Dict,
            required=False,
            help='Information of the impurity like position in the unit cell, screening cluster, atom type.'
        )

        spec.input('kkr', valid_type=Code, required=False, help='KKRhost code, needed to run the KkrCalculation')

        spec.input(
            'kkrimp', valid_type=Code, required=True, help='KKRimp code used to converge the impurity calculation'
        )

        spec.input(
            'voronoi',
            valid_type=Code,
            required=False,
            help='Voronoi code used to create the impurity starting potential.'
        )

        spec.input(
            'calc_DOS',
            valid_type=Bool,
            required=False,
            default=lambda: Bool(False),
            help='Set this to TRUE to calculate DOS'
        )

        spec.input(
            'BdG_settings',
            valid_type=Dict,
            default=lambda: Dict(
                dict={
                    # 'DELTA_BDG': 1e-4, # 1.36 meV # could be set in the future?
                    'USE_BdG': True,
                    'USE_E_SYMM_BDG': True,
                    # 'FIX_NONCO_ANGLES': True, # could be set in the future?
                }
            ),
            help='Define BdG parameters'
        )

        # expose inputs for impurity normal state scf
        spec.expose_inputs(
            kkr_imp_wc,
            namespace='imp_scf',
            include=(
                'startpot', 'wf_parameters', 'gf_writeout', 'scf.params_overwrite', 'scf.initial_noco_angles',
                'scf.settings_LDAU'
            )
        )
        spec.inputs['imp_scf']['gf_writeout']['kkr'].required = False
        spec.input('imp_scf.options', required=False, help='computer options for impurity scf step')
        spec.input(
            'imp_scf.converged_potential',
            valid_type=SinglefileData,
            required=False,
            help='Converged impurity potential. Skips impurity scf step.'
        )

        spec.input(
            'imp_scf.remote_data_host',
            valid_type=RemoteData,
            required=False,
            help='Parent folder of previously converged host normal state KkrCalculation'
        )

        spec.input(
            'imp_scf.remote_data_gf',
            valid_type=RemoteData,
            required=False,
            help='RemoteData node of precomputed host Green function'
        )

        # inputs for impurity BdG scf
        spec.expose_inputs(
            kkr_imp_wc,
            namespace='BdG_scf',
            include=('startpot', 'remote_data_gf', 'gf_writeout', 'scf.initial_noco_angles')
        )
        spec.inputs['BdG_scf']['gf_writeout']['kkr'].required = False
        spec.input('BdG_scf.options', required=False, help='computer options for BdG impurity scf step')

        spec.input(
            'BdG_scf.remote_data_host', required=False, help='Parent folder of previously converged BdG KkrCalculation'
        )

        # inputs for impurity dos
        spec.expose_inputs(
            kkr_imp_dos_wc,
            namespace='dos',
            include=('wf_parameters', 'gf_dos_remote', 'gf_writeout', 'initial_noco_angles')
        )

        spec.input(
            'dos.gf_writeout.host_remote',
            valid_type=RemoteData,
            required=False,
            help='Parent calculation from where the GF writeout starts. '
        )

        spec.input(
            'dos.gf_writeout.kkr',
            valid_type=Code,
            required=False,
            help='KKRhost code used to create DOS kkrflex files'
        )

        spec.input('dos.options', valid_type=Dict, required=False, help='Computer options for DOS step')

        # Here outputs are defined

        spec.output('workflow_info', valid_type=Dict, required=False)
        spec.output('output_parameters', valid_type=Dict, required=False)
        spec.output('dos_data', required=False, valid_type=XyData)
        spec.output('dos_data_interpol', required=False, valid_type=XyData)
        spec.output('dos_data_lm', required=False, valid_type=XyData)
        spec.output('dos_data_lm_interpol', required=False, valid_type=XyData)
        spec.output('impurity_potential', valid_type=SinglefileData)
        spec.output('gf_host_BdG', valid_type=RemoteData, required=False)

        # Here outlines are being specified
        spec.outline(
            # For initializing workflow
            cls.start,
            cls.validate_input,
            if_(cls.do_imp_pot_calc)(cls.imp_pot_calc),
            if_(cls.do_BdG_scf)(cls.imp_BdG_calc),
            if_(cls.do_calc_DOS)(cls.DOS_calc),
            cls.results
        )

        # Define all possible error messages
        spec.exit_code(
            100, 'ERROR_KKRCODE_NOT_CORRECT', 'The code you provided for kkr does not use the plugin kkr.kkr'
        )
        spec.exit_code(
            101, 'ERROR_KKRIMPCODE_NOT_CORRECT', 'The code you provided for kkrimp does not use the plugin kkr.kkrimp'
        )
        spec.exit_code(
            102, 'ERROR_VORONOICODE_NOT_CORRECT',
            'The code you provided for voronoi does not use the plugin kkr.voronoi'
        )

        spec.exit_code(200, 'ERROR_INVALID_PARENT', 'Parent calculation is not valid')

        #Now here we define all the functions from the outline

    def start(self):
        """
        Set up context of the workflow
        """
        self.report(f'INFO: started KKR BdG impurity version {self._wf_version}')

    def validate_input(self):
        """
        validate inputs
        """

        # validate for kkr code
        if 'kkr' in self.inputs:
            try:
                test_and_get_codenode(self.inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT  # pylint: disable=no-member

        # validate for kkrimp code, always done because this is a required input
        try:
            test_and_get_codenode(self.inputs.kkrimp, 'kkr.kkrimp', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_KKRIMPCODE_NOT_CORRECT  # pylint: disable=no-member

        # validate for voronoi code
        if 'voronoi' in self.inputs:
            try:
                test_and_get_codenode(self.inputs.voronoi, 'kkr.voro', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_VORONOICODE_NOT_CORRECT  # pylint: disable=no-member

    def do_imp_pot_calc(self):
        """
        run impurity potential calculation only if converged impurity potential not provided already
        """
        if not 'converged_potential' in self.inputs.imp_scf:
            return True

    def imp_pot_calc(self):
        """
        run normal state impurity scf calculation
        """

        builder = kkr_imp_wc.get_builder()

        builder.impurity_info = self.inputs.impurity_info
        builder.voronoi = self.inputs.voronoi
        builder.kkr = self.inputs.kkr
        builder.kkrimp = self.inputs.kkrimp
        if 'remote_data_gf' in self.inputs.imp_scf:
            builder.remote_data_gf = self.inputs.imp_scf.remote_data_gf
        else:
            builder.remote_data_host = self.inputs.imp_scf.remote_data_host
        builder.wf_parameters = self.inputs.imp_scf.wf_parameters
        if 'options' in self.inputs.imp_scf:
            builder.options = self.inputs.imp_scf.options
        else:
            builder.options = self.inputs.options
        if 'initial_noco_angles' in self.inputs.imp_scf:
            builder.scf.initial_noco_angles = self.inputs.imp_scf.initial_noco_angles  # pylint: disable=no-member

        if 'gf_writeout' in self.inputs.imp_scf:
            if 'options' in self.inputs.imp_scf.gf_writeout:
                builder.gf_writeout.options = self.inputs.imp_scf.gf_writeout.options  # pylint: disable=no-member
            if 'kkr' in self.inputs.imp_scf.gf_writeout:
                builder.kkr = self.inputs.imp_scf.gf_writeout.kkr
            if 'params_kkr_overwrite' in self.inputs.imp_scf.gf_writeout:
                builder.params_kkr_overwrite = self.inputs.imp_scf.gf_writeout.params_kkr_overwrite
        builder.gf_writeout.kkr = builder.kkr  # pylint: disable=no-member

        if 'params_overwrite' in self.inputs.imp_scf.scf:
            builder.scf.params_overwrite = self.inputs.imp_scf.scf.params_overwrite  # pylint: disable=no-member

        if 'startpot' in self.inputs.imp_scf:
            builder.startpot = self.inputs.imp_scf.startpot

        if 'settings_LDAU' in self.inputs.imp_scf.scf:
            builder.scf.settings_LDAU = self.inputs.imp_scf.scf.settings_LDAU  # pylint: disable=no-member

        imp_calc = self.submit(builder)

        return ToContext(last_imp_calc=imp_calc)

    def imp_BdG_calc(self):
        """
        run BdG one-shot impurity calculation
        """

        builder = kkr_imp_wc.get_builder()
        if 'impurity_info' in self.inputs:
            builder.impurity_info = self.inputs.impurity_info
        builder.voronoi = self.inputs.voronoi
        builder.kkr = self.inputs.kkr
        builder.kkrimp = self.inputs.kkrimp

        if 'startpot' in self.inputs.imp_scf:
            builder.startpot = self.inputs.imp_scf.startpot
        else:
            builder.startpot = self.ctx.last_imp_calc.outputs.converged_potential

        if 'remote_data_gf' in self.inputs.BdG_scf:
            builder.remote_data_gf = self.inputs.BdG_scf.remote_data_gf

        if 'gf_writeout' in self.inputs.BdG_scf:
            if 'kkr' in self.inputs.BdG_scf.gf_writeout:
                builder.kkr = self.inputs.BdG_scf.gf_writeout.kkr
            if 'params_kkr_overwrite' in self.inputs.BdG_scf.gf_writeout:
                builder.params_kkr_overwrite = self.inputs.BdG_scf.gf_writeout.params_kkr_overwrite
        if 'kkr' in self.inputs:
            builder.gf_writeout.kkr = builder.kkr  # pylint: disable=no-member
        if 'initial_noco_angles' in self.inputs.BdG_scf:
            builder.scf.initial_noco_angles = self.inputs.BdG_scf.initial_noco_angles  # pylint: disable=no-member

        builder.remote_data_host = self.inputs.BdG_scf.remote_data_host

        if 'options' in self.inputs.BdG_scf:
            builder.options = self.inputs.BdG_scf.options
        else:
            if 'options' in self.inputs.imp_scf:
                builder.options = self.inputs.imp_scf
            else:
                builder.options = self.inputs.options

        settings = kkr_imp_wc.get_wf_defaults()[1]
        settings['strmix'] = 0.01
        settings['aggrmix'] = 0.01
        settings['nsteps'] = 1
        settings['kkr_runmax'] = 1
        settings['mag_init'] = True
        builder.wf_parameters = Dict(dict=settings)
        BdG_params = self.inputs.BdG_settings.get_dict()
        builder.scf.params_overwrite = Dict(dict=BdG_params)  # pylint: disable=no-member

        imp_calc_BdG = self.submit(builder)

        return ToContext(last_imp_calc_BdG=imp_calc_BdG)

    def do_BdG_scf(self):
        """
        run BdG scf step only if BdG impurity potential not provided and DOS calculation is not planned
        """
        # BdG scf not implemented in the KKRimp code yet
        return False
        # if (not 'startpot' in self.inputs.BdG_scf) and (not self.inputs.calc_DOS):
        # return True

    def do_calc_DOS(self):
        """
        check if DOS calculation should be done
        """
        if self.inputs.calc_DOS:
            return True

    def DOS_calc(self):
        """
        run DOS calculation
        """
        from aiida_kkr.workflows import kkr_imp_dos_wc

        builder = kkr_imp_dos_wc.get_builder()

        if 'kkr' in self.inputs:
            builder.kkr = self.inputs.kkr
        builder.kkrimp = self.inputs.kkrimp

        #define computer options
        if 'options' in self.inputs.dos:
            builder.options = self.inputs.dos.options
        else:
            if 'options' in self.inputs.BdG_scf:
                builder.options = self.inputs.BdG_scf.options
            else:
                if 'options' in self.inputs.imp_scf:
                    builder.options = self.inputs.imp_scf.options
                else:
                    builder.options = self.inputs.options

        # set nonco angles
        if 'initial_noco_angles' in self.inputs.dos:
            builder.initial_noco_angles = self.inputs.dos.initial_noco_angles

        # skip BdG step and just use the starting potential instead?
        # faster and same accuracy?!
        if 'startpot' in self.inputs.BdG_scf:
            builder.imp_pot_sfd = self.inputs.BdG_scf.startpot
        else:
            if not self.inputs.calc_DOS:
                builder.imp_pot_sfd = self.ctx.last_imp_calc_BdG.outputs.converged_potential
            else:
                if 'converged_potential' in self.inputs.imp_scf:
                    builder.imp_pot_sfd = self.inputs.imp_scf.converged_potential
                else:
                    builder.imp_pot_sfd = self.ctx.last_imp_calc.outputs.converged_potential

        if 'impurity_info' in self.inputs:
            builder.impurity_info = self.inputs.impurity_info

        builder.wf_parameters = Dict(dict=kkr_imp_dos_wc.get_wf_defaults())

        if 'wf_parameters' in self.inputs.dos:
            builder.wf_parameters = self.inputs.dos.wf_parameters

        BdG_params = self.inputs.BdG_settings.get_dict()
        builder.BdG.params_overwrite = Dict(dict=BdG_params)  # pylint: disable=no-member

        if 'gf_writeout' in self.inputs.dos:
            # set optional inputs
            if 'kkr' in self.inputs.dos.gf_writeout:
                builder.kkr = self.inputs.dos.gf_writeout.kkr
            if 'params_kkr_overwrite' in self.inputs.dos.gf_writeout:
                builder.gf_writeout.params_kkr_overwrite = self.inputs.dos.gf_writeout.params_kkr_overwrite  # pylint: disable=no-member
            if 'host_remote' in self.inputs.dos.gf_writeout:
                builder.host_remote = self.inputs.dos.gf_writeout.host_remote
            if 'options' in self.inputs.dos.gf_writeout:
                builder.gf_writeout.options = self.inputs.dos.gf_writeout.options  # pylint: disable=no-member
        if 'kkr' in self.inputs:
            builder.gf_writeout.kkr = builder.kkr  # pylint: disable=no-member
        if 'gf_dos_remote' in self.inputs.dos:
            builder.gf_dos_remote = self.inputs.dos.gf_dos_remote

        DOS_calc = self.submit(builder)

        return ToContext(DOS_node=DOS_calc)

    def results(self):
        """
        return the results nodes of the workchain
        """
        if self.inputs.calc_DOS:
            if 'dos_data' in self.ctx.DOS_node.outputs.dos_data:
                self.out('dos_data', self.ctx.DOS_node.outputs.dos_data)
            if 'dos_data_interpol' in self.ctx.DOS_node.outputs.dos_data_interpol:
                self.out('dos_data_interpol', self.ctx.DOS_node.outputs.dos_data_interpol)
            if 'dos_data_lm' in self.ctx.DOS_node.outputs:
                self.out('dos_data_lm', self.ctx.DOS_node.outputs.dos_data_lm)
            if 'dos_data_interpol_lm' in self.ctx.DOS_node.outputs:
                self.out('dos_data_interpol_lm', self.ctx.DOS_node.outputs.dos_data_interpol_lm)

        if self.do_BdG_scf():
            if 'startpot' not in self.inputs.BdG_scf:
                self.out('workflow_info', self.ctx.last_imp_calc_BdG.outputs.workflow_info)
                self.out('output_parameters', self.ctx.last_imp_calc_BdG.outputs.last_calc_output_parameters)
                if 'remote_data_gf' not in self.inputs.BdG_scf:
                    self.out('gf_host_BdG', self.ctx.last_imp_calc_BdG.outputs.remote_data_gf)
                else:
                    self.out('gf_host_BdG', self.inputs.BdG_scf.remote_data_gf)

        if 'converged_potential' not in self.inputs.imp_scf:
            self.out('impurity_potential', self.ctx.last_imp_calc.outputs.converged_potential)
        else:
            self.out('impurity_potential', self.inputs.imp_scf.converged_potential)
