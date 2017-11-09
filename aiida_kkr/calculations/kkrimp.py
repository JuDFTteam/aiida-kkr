import json

from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.orm import DataFactory

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

class KkrimpCalculation(JobCalculation):
    """
    AiiDA calculation plugin for a KKR calculation
    .
    """

    def _init_internal_params(self):
        """
        Init internal parameters at class load time
        """
        # reuse base class function
        super(KkrCalculation, self)._init_internal_params()

        # List of mandatory input files
        self._INPUTCARD = 'inputcard'
        self._POTENTIAL = 'potential'

        # List of optional input files (may be mandatory for some setting in inputputcard)
        self._SHAPEFUN = 'shapefun'
        self._SCOEF = 'scoef'
        self._NONCO_ANGLES = 'nonco_angles.dat'

	
	# List of output files that should always be present
        self._OUT_POTENTIAL = 'out_potential'
        self._OUTPUT_0_INIT = 'output.0.txt'
        self._OUTPUT_000 = 'output.000.txt'
        self._OUT_TIMING_000 = 'out_timing.000.txt'


        # template.product entry point defined in setup.json
        self._default_parser = 'kkr.kkrimpParser'

    @classproperty
    def _use_methods(cls):
        """
        Add use_* methods for calculations.
        
        Code below enables the usage
        my_calculation.use_parameters(my_parameters)
        """
        use_dict = JobCalculation._use_methods
        use_dict.update({
            "parameters": {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'parameters',
                'docstring':
                ("Use a node that specifies the input parameters ")
            },
            "structure": {
                'valid_types': StructureData,
                'additional_parameter': None,
                'linkname': 'parameters',
                'docstring':
                ("Use a node that specifies the input crystal structure ")
                },
            })
        return use_dict

    def _prepare_for_submission(self, tempfolder, inputdict):
        """
        Create input files.

            :param tempfolder: aiida.common.folders.Folder subclass where
                the plugin should put all its files.
            :param inputdict: dictionary of the input nodes as they would
                be returned by get_inputs_dict
        """
        # Check inputdict
        try:
            parameters = inputdict.pop(self.get_linkname('parameters'))
        except KeyError:
            raise InputValidationError("No parameters specified for this "
                                       "calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters not of type "
                                       "ParameterData")
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this "
                                       "calculation")
        if inputdict:
            raise ValidationError("Unknown inputs")

        # In this example, the input file is simply a json dict.
        # Adapt for your particular code!
        input_dict = parameters.get_dict()

        # Write input to file
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        with open(input_filename, 'w') as infile:
            json.dump(input_dict, infile)

        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = []
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = [self._OUTPUT_FILE_NAME]

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = [
            self._INPUT_FILE_NAME, self._OUTPUT_FILE_NAME
        ]
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo
