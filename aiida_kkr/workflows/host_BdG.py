# It will calculate host BdG calculatin

from aiida.engine import Worchain, ToContext
from aiida_kkr.workflows import kkr_scf_wc, kkr_dos_wc
from aiida.orm import  StructureData

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = (u'Mohammad Hemmati')


class kkrhost_BdG_wc(WorkChain):
    """
    Workchain for blabla
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
            print(f'Version of the kkrhost_BdG_wc workflow: {self._wf_version}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(kkrhost_BdG_wc, cls).define(spec)

        # here inputs are defined
        #spec.input(
        #    'wf_parameters',
        #    valid_type=Dict,
        #    required=False,
        #    default=lambda: Dict(dict=cls._wf_default),
        #    help='Parameters of the BdG impurity workflow (see output of kkrimp_BdG_wc.get_wf_default() for more details).'
        #)
        
        spec.expose_inputs(
            kkr_scf_wc, namespace='kkr_scf', include=('wf_parameters', 'calc_parameters')
        )

        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help=
            'Computer options (walltime etc.) passed onto KkrCalculation, fall back to settings from parent calculation if not given'
        )
        spec.input(
            'structure',
            valid_type=StructureData,
            required=True,
            help=
            'This is the structure we want to use'
        )


        spec.input(
            'kkr', valid_type=Code,
            required=True,
            help='KKRhost code, needed to run the KkrCalculation')

       

        spec.input(
            'voronoi',
            valid_type=Code,
            required=True,
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

 
        spec.input(
            'dos.gf_writeout.host_remote',
            valid_type=RemoteData,
            required=False,
            help='Parent calculation from where the GF writeout starts. '
        )


        # Here outputs are defined

        spec.output('output_parameters', valid_type=Dict, required=False)
 

        # Here outlines are being specified
        spec.outline(
            # For initialiging workflow
            cls.start,
            cls.validate_input,
            cls.kkr_scf_wc,
            
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

   
    def kkr_scf_wf(self):
#///////////////////
#///////////////////
        
        
        
        
        builder = kkr_scf_wc.get_builder()
        
        builder.structure = self.inputs.structure
        
        
        settings = kkr_scf_wc.get_wf_defaults()[0]
        
        # I need to double check this, do I need to be in the defaults
        settings.convergence_setting_coarse = settings.convergence_setting_fine
        builder.wf_parameters = Dict(dict=settings)
        
        if 'wf_parameters' in 'kkr_scf':
            builder.wf_parameters = self.inputs.kkr_scf.wf_parameters

        
        if 'calc_parameters' in 'kkr_scf':
            builder.clac_paramters = self.inputs.kkr_scf.calc_parameters
        
        
        builder.kkr = self.inputs.kkr
        
        builder.voronoi = self.inputs.voronoi
        
        if 'options' in self.inputs:
            builder.options = self.inputs.options
        
        scf_start = self.submit(builder)
        
        return ToContext(kkr_scf_start=scf_start)

#///////////////////

    def results(self):
        """
        return the results nodes of the workchain
        """

