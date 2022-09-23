# Workflow for impurity BdG calculation from converged normal state impurity portential and BdG host calculation

from aiida.engine import WorkChain, ToContext
from aiida.orm import Dict, RemoteData, Code, CalcJobNode, WorkChainNode, Float
from aiida_kkr.workflows import kkr_imp_wc
from aiida_kkr.tools.find_parent import get_calc_from_remote
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = (u'David Antognini Silva')


class kkrimp_BdG_wc(WorkChain):
    """
    Workchain for blablabla
    inputs::
        :blabla: blabla
    returns::
        :blabla : blabla
    """

    _wf_version = __version__
    _wf_default = {
        'blabla': None,  #Put in here default input parameters
    }
    _options_default = {
        'max_wallclock_seconds': 36000,
        'resources': {
            'num_machines': 1
        },
        'withmpi': True,
        'queue_name': ''
    }

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Return the default values of the workflow parameters (wf_parameters input node)
        """
        if not silent:
            print(f'Version of the kkrimp_BdG_wc workflow: {self._wf_version}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(kkrimp_BdG_wc, cls).define(spec)

        # here inputs are defined
        #spec.input(
        #    'wf_parameters',
        #    valid_type=Dict,
        #    required=False,
        #    default=lambda: Dict(dict=cls._wf_default),
        #    help='Parameters of the BdG impurity workflow (see output of kkrimp_BdG_wc.get_wf_default() for more details).'
        #)

        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help=
            'Computer options (walltime etc.) passed onto KkrCalculation, fall back to settings from parent calculation if not given'
        )

        spec.input(
            'remote_data_host',
            valid_type=RemoteData,
            required=True,
            help='Parent folder of previously converged host normal state KkrCalculation'
        )

        spec.input(
            'remote_data_host_BdG',
            valid_type=RemoteData,
            required=True,
            help='Parent folder of previously converged BdG KkrCalculation'
        )

        spec.input(
            'impurity_info',
            valid_type=Dict,
            required=True,
            help='Information of the impurity like position in the unit cell, screening cluster, atom type.'
        )

        spec.input('kkr', valid_type=Code, required=True, help='KKRhost code, needed to run the KkrCalculation')

        spec.input(
            'kkrimp', valid_type=Code, required=True, help='KKRimp code used to converge the impurity calculation'
        )

        spec.input(
            'voronoi',
            valid_type=Code,
            required=True,
            help='Voronoi code used to create the impurity starting potential.'
        )

        spec.expose_inputs(kkr_imp_wc, namespace='imp_scf', include=('startpot', 'wf_parameters'))
        spec.expose_inputs(kkr_imp_wc, namespace='BdG_scf', include=('startpot', 'remote_data_gf'))

        # Here outputs are defined

        #spec.output('results_wf', valid_type=WorkChainNode)
        #spec.output('total_energy')
        spec.output('workflow_info', valid_type=Dict)
        spec.output('output_parameters', valid_type=Dict)

        # Here outlines are being specified
        spec.outline(
            # For initialiging workflow
            cls.start,
            cls.validate_input,
            cls.imp_pot_calc,
            cls.imp_BdG_calc,
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
        try:
            test_and_get_codenode(self.inputs.kkr, 'kkr.kkr', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT  # pylint: disable=no-member

        # validate for kkrimp code
        try:
            test_and_get_codenode(self.inputs.kkrimp, 'kkr.kkrimp', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_KKRIMPCODE_NOT_CORRECT

        # validate for voronoi code
        try:
            test_and_get_codenode(self.inputs.voronoi, 'kkr.voro', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_VORONOICODE_NOT_CORRECT

        # save parent calculation
        input_remote = self.inputs.remote_data_host
        parents = input_remote.get_incoming(node_class=CalcJobNode).all()
        if len(parents) != 1:
            # check if parent is unique
            return self.exit_codes.ERROR_INVALID_PARENT  # pylint: disable=no-member
        self.ctx.parent_calc = get_calc_from_remote(input_remote)

    def imp_pot_calc(self):
        """
        calculate normal state impurity potential
        """
        if 'startpot' not in self.inputs.BdG_scf:

            builder = kkr_imp_wc.get_builder()
            builder.impurity_info = self.inputs.impurity_info
            builder.voronoi = self.inputs.voronoi
            builder.kkr = self.inputs.kkr
            builder.kkrimp = self.inputs.kkrimp
            builder.options = self.inputs.options
            builder.remote_data_host = self.inputs.remote_data_host

            settings = kkr_imp_wc.get_wf_defaults()[1]
            settings['strmix'] = 0.01
            settings['aggrmix'] = 0.01
            #settings['nsteps'] = 200
            settings['nsteps'] = 2
            settings['mag_init'] = True

            #remove the following line that is there for testing purposes
            settings['kkr_runmax'] = 1

            builder.wf_parameters = Dict(dict=settings)
            # builder.wf_parameters = self.inputs.imp_scf.wf_parameters

            imp_calc = self.submit(builder)
            #self.ctx.imp_calc = imp_calc

            return ToContext(last_imp_calc=imp_calc)

    def imp_BdG_calc(self):
        """
        BdG one-shot impurity calculation
        """

        builder = kkr_imp_wc.get_builder()
        builder.impurity_info = self.inputs.impurity_info
        builder.voronoi = self.inputs.voronoi
        builder.kkr = self.inputs.kkr
        builder.kkrimp = self.inputs.kkrimp

        if 'startpot' in self.inputs.BdG_scf:
            builder.startpot = self.inputs.BdG_scf.startpot
        else:
            builder.startpot = self.ctx.last_imp_calc.outputs.converged_potential

        if 'remote_data_gf' in self.inputs.BdG_scf:
            builder.remote_data_gf = self.inputs.BdG_scf.remote_data_gf

        builder.remote_data_host = self.inputs.remote_data_host_BdG
        builder.options = self.inputs.options

        settings = kkr_imp_wc.get_wf_defaults()[1]
        settings['strmix'] = 0.01
        settings['aggrmix'] = 0.01
        settings['nsteps'] = 1
        settings['kkr_runmax'] = 1
        settings['mag_init'] = True
        builder.wf_parameters = Dict(dict=settings)
        builder.scf.params_overwrite = Dict(dict={'USE_BdG': True, 'USE_E_SYMM_BdG': True})

        imp_calc_BdG = self.submit(builder)

        return ToContext(last_imp_calc_BdG=imp_calc_BdG)

    def results(self):
<<<<<<< HEAD
        #self.out('results_wf', self.ctx.last_imp_calc_BdG)
        self.out('workflow_info',self.ctx.last_imp_calc_BdG.outputs.workflow_info)
        self.out('output_parameters', self.ctx.last_imp_calc_BdG.outputs.last_calc_output_parameters)
        #tot_energy = self.ctx.last_imp_calc_BdG.outputs.last_calc_output_parameters.get_attribute('energy')
        #self.out('total_energy', tot_energy)
=======
        result = self.ctx.last_imp_calc_BdG
        self.out('results_wf', result)
>>>>>>> 33683b5137286b16ea9a77bcd0c4b1022cc2dc89
