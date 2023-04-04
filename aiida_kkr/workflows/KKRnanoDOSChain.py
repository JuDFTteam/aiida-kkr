from aiida.orm import Code, CalcJobNode, Dict, Bool, RemoteData, StructureData, load_node
from aiida.engine import calcfunction, WorkChain, ToContext, if_
from aiida.plugins.factories import CalculationFactory
from aiida_kkr.calculations.kkrnano import KKRnanoCalculation
from aiida_kkr.data.strucwithpot import StrucWithPotData

__copyright__ = (u'Copyright (c), 2021, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.1'
__contributors__ = ('Markus Struckmann')


class KKRnanoDOSChain(WorkChain):
    """WorkChain to multiply two numbers and add a third, for testing and demonstration purposes."""
    #TODO: Add functionality for NOCO mode
    #TODO: change scfsteps to realistic numbers, e.g. 20, 100, 100
    _DEFAULT_KKRNANO_WC_PARA = {
        'coarse_mixing': {
            'scfsteps': {
                'value': 1
            },
            'imix': {
                'value': 0
            },
            'mixing': {
                'value': 0.05
            },
            'target_rms': {
                'value': 1e-02
            }
        },
        'colinear': {
            'scfsteps': {
                'value': 1
            },
            'imix': {
                'value': 6
            },
            'mixing': {
                'value': 0.01
            },
            'target_rms': {
                'value': 1e-08
            }
        },
        'SOC': {
            'scfsteps': {
                'value': 1
            },
            'imix': {
                'value': 6
            },
            'mixing': {
                'value': 0.01
            },
            'target_rms': {
                'value': 1e-05
            }
        },
        'DOS': {
            'imix': {
                'value': 6
            },
            'mixing': {
                'value': 0.01
            },
            'target_rms': {
                'value': 1e-08
            },
            'emin': {
                'value': -0.35
            },
            'emax': {
                'value': 0.35
            },
            'npnt1': {
                'value': 0
            },
            'npnt2': {
                'value': 70
            },
            'npnt3': {
                'value': 0
            },
            'tempr': {
                'value': 300
            }
        }
    }

    _OPTIONS_DEFAULT = {
        'queue_name': 'th1-2020-32',  # Queue name to submit jobs to
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 15
        },  # ressources to allocate for the job
        'withmpi': True,  # execute KKR with mpi or without
        'custom_scheduler_commands': ''  # some additional scheduler commands
    }
    _OPTIONS_ADVANCED = {
        'queue_name': 'th1-2020-32',
        'max_wallclock_seconds': 72000,  # 20 h
        'resources': {
            'tot_num_mpiprocs': 15,
            'num_machines': 1
        },
        'withmpi': True
    }

    def str_mixing(self):
        builder = KKRnanoCalculation.get_builder()
        params = self.inputs.calc_parameters.get_dict()
        strparams = self.inputs.wf_parameters.get_dict()['coarse_mixing']

        #Start from parent or from strucwithpot
        if hasattr(self.inputs, 'strucwithpot') and not hasattr(self.inputs, 'parent_folder'):
            builder.strucwithpot = self.inputs.strucwithpot
        else:
            parent_folder = self.inputs.parent_folder
            builder.parent_folder = parent_folder

        for key in strparams:
            params[key] = strparams[key]


#             params['mixing']=strparams["mixing"]
#             params["imix"]=strparams["imix"]
#             params['target_rms']={"value": 1e-02}
#             params["scfsteps"]=strparams["scfsteps"]

#KKRnanoCalculation._check_input_dict(load_node('58f35be4-3ddc-4b24-947b-c3ddbb8959ac'),params)

#preparing the approptiate parameters
        params['soc'] = {'value': False}
        params['KORBIT'] = {'value': 0}

        #builder = KKRnanoCalculation.get_builder()

        builder.metadata.label = 'WC_coarse_mixing'  # pylint: disable=no-member

        builder.parameters = Dict(params)

        #builder.metadata.dry_run = True
        #machine=#
        builder.metadata.options = self.inputs.options.get_dict()  # pylint: disable=no-member

        builder.code = self.inputs.code  # orm.Code.get_from_string("KKRnanoVersuch2@iffslurm")
        str_calc = self.submit(builder)

        return ToContext(str_result=str_calc)

    def col(self):
        builder = KKRnanoCalculation.get_builder()
        params = self.inputs.calc_parameters.get_dict()
        colparams = self.inputs.wf_parameters.get_dict()['colinear']
        parent_folder = self.ctx.str_result.outputs.remote_folder

        for key in colparams:
            params[key] = colparams[key]

        #preparing the approptiate parameters
        params['soc'] = {'value': False}
        params['KORBIT'] = {'value': 0}

        #KKRnanoCalculation._check_input_dict(load_node('58f35be4-3ddc-4b24-947b-c3ddbb8959ac'),params)

        #builder = KKRnanoCalculation.get_builder()

        builder.metadata.label = 'WC_colinear'  # pylint: disable=no-member

        builder.parameters = Dict(params)

        #builder.metadata.dry_run = True
        #machine=#
        builder.metadata.options = self.inputs.options.get_dict()  # pylint: disable=no-member
        builder.parent_folder = parent_folder

        builder.code = self.inputs.code  # orm.Code.get_from_string("KKRnanoVersuch2@iffslurm")
        col_calc = self.submit(builder)

        return ToContext(col_result=col_calc)

    def SOC(self):
        builder = KKRnanoCalculation.get_builder()
        params = self.inputs.calc_parameters.get_dict()
        socparams = self.inputs.wf_parameters.get_dict()['SOC']

        if self.inputs.start_from_col.value:
            #Start from parent or from strucwithpot
            if hasattr(self.inputs, 'strucwithpot') and not hasattr(self.inputs, 'parent_folder'):
                builder.strucwithpot = self.inputs.strucwithpot
            else:
                parent_folder = self.inputs.parent_folder
                builder.parent_folder = parent_folder

        else:
            parent_folder = self.ctx.col_result.outputs.remote_folder
            builder.parent_folder = parent_folder

        for key in socparams:
            params[key] = socparams[key]

        #preparing the approptiate parameters
        params['soc'] = {'value': True}
        params['KORBIT'] = {'value': 1}

        #KKRnanoCalculation._check_input_dict(load_node('58f35be4-3ddc-4b24-947b-c3ddbb8959ac'),params)

        #builder = KKRnanoCalculation.get_builder()

        builder.metadata.label = 'WC_SOC'  # pylint: disable=no-member

        builder.parameters = Dict(params)

        #builder.metadata.dry_run = True

        builder.metadata.options = self.inputs.options.get_dict()  # pylint: disable=no-member

        builder.code = self.inputs.code  # orm.Code.get_from_string("KKRnanoVersuch2@iffslurm")
        soc_calc = self.submit(builder)

        return ToContext(soc_result=soc_calc)

    def DOS(self):

        builder = KKRnanoCalculation.get_builder()
        params = self.inputs.calc_parameters.get_dict()
        dosparams = self.inputs.wf_parameters.get_dict()['DOS']
        parent_folder = self.ctx.soc_result.outputs.remote_folder

        for key in dosparams:
            params[key] = dosparams[key]

        #preparing the approptiate parameters
        params['soc'] = {'value': True}
        params['KORBIT'] = {'value': 1}
        params['scfsteps'] = {'value': 1}
        params['npol'] = {'value': 0}

        #KKRnanoCalculation._check_input_dict(load_node('58f35be4-3ddc-4b24-947b-c3ddbb8959ac'),params)

        builder = KKRnanoCalculation.get_builder()

        builder.metadata.label = 'WC_DOS'  # pylint: disable=no-member

        builder.parameters = Dict(params)

        #builder.metadata.dry_run = True

        builder.metadata.options = self.inputs.options.get_dict()  # pylint: disable=no-member
        builder.parent_folder = parent_folder

        builder.code = self.inputs.code  # orm.Code.get_from_string("KKRnanoVersuch2@iffslurm")
        dos_calc = self.submit(builder)

        return ToContext(dos_result=dos_calc)

    #@calcfunction
    #def create_output_node(outself):

    def prepare_output(self):
        #         output_dictionary=self.ctx.soc_result.outputs.output_parameters.get_dict()
        #         if self.inputs.calculate_DOS.value:
        #             dos_dict=self.ctx.dos_result.outputs.output_parameters.get_dict()['DOS']
        #             output_dictionary['DOS']=dos_dict
        #         self.out('result', Dict(dict=output_dictionary))
        if self.inputs.calculate_DOS.value:
            output_dictionary = self.ctx.dos_result.outputs.output_parameters
        else:
            output_dictionary = self.ctx.soc_result.outputs.output_parameters
        self.out('result', output_dictionary)

    def check_input(self):
        #TODO add inputdict check like in calculation module
        inputdict = self.inputs.calc_parameters.get_dict()
        try:
            if not inputdict['NSPIND']['value'] == 2:
                return self.exit_codes.ERROR_INPUT_SPIN  # pylint: disable=no-member
        except:
            return self.exit_codes.ERROR_INPUT_SPIN  # pylint: disable=no-member

    def check_if_DOS_to_run(self):
        if self.inputs.calculate_DOS.value:
            return True
        else:
            return False

    def check_if_str_to_run(self):
        if self.inputs.start_from_col:
            return False
        else:
            return True

    def check_if_col_to_run(self):
        if self.inputs.start_from_col:
            return False
        else:
            return True

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)

        # define input nodes (optional ones have required=False)
        spec.input(
            'calc_parameters',
            valid_type=Dict,
            required=True,
            help='Dict node that specifies the input calc_parameters for KKRnano (k-point density etc.)'
        )

        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._DEFAULT_KKRNANO_WC_PARA),
            help='Dict node that specifies the WorkChain calc_parameters for KKRnano (thresholds, mixings, etc.)'
        )
        spec.input(
            'parent_folder',
            valid_type=RemoteData,
            required=False,
            help='Use a node that specifies a parent KKRnano or voronoi calculation'
        )
        spec.input(
            'strucwithpot',
            valid_type=StrucWithPotData,
            required=False,
            help=
            'Use a StrucWithPotData node, that contains the AiiDA structure and the shape and potential files of the individual atoms.'
        )

        spec.input('options', valid_type=Dict, required=False, default=lambda: Dict(dict=cls._OPTIONS_DEFAULT))
        spec.input(
            'calculate_DOS',
            valid_type=Bool,
            required=False,
            default=lambda: Bool(False),
            help='Trigger DOS calculation.'
        )
        spec.input(
            'start_from_col',
            valid_type=Bool,
            required=False,
            default=lambda: Bool(False),
            help='Start with SOC calculation.'
        )

        spec.input('code', valid_type=Code)

        spec.outline(
            cls.check_input,
            if_(cls.check_if_str_to_run)(cls.str_mixing),
            if_(cls.check_if_col_to_run)(cls.col), cls.SOC,
            if_(cls.check_if_DOS_to_run)(cls.DOS), cls.prepare_output
        )
        spec.output('result', valid_type=Dict)
        spec.exit_code(400, 'ERROR_INPUT_SPIN', message='Number of spin directions has to be set to 2 for SOC!')
