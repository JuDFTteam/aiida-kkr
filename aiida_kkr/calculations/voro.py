# -*- coding: utf-8 -*-
"""
Input plug-in for a voronoi calculation.
"""

from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.orm import DataFactory
from aiida_kkr.tools.common_functions import generate_inputcard_from_structure, check_2Dinput

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

class VoronoiCalculation(JobCalculation):
    """
    AiiDA calculation plugin for a voronoi calculation (creation of starting potential and shapefun)
    .
    """

    def _init_internal_params(self):
        """
        Init internal parameters at class load time
        """
        # reuse base class (i.e. JobCalculation) functions
        super(VoronoiCalculation, self)._init_internal_params()

        # Default input and output files
        self._DEFAULT_INPUT_FILE = 'inputcard' # will be shown with inputcat
        self._DEFAULT_OUTPUT_FILE = 'out_voronoi' #'shell output will be shown with outputca

        # List of mandatory input files
        self._INPUT_FILE_NAME = 'inputcard'
        #self._INPUTCARD = 'inputcard'
	
	   # List of output files that should always be present
        self._OUTPUT_FILE_NAME = 'out_voronoi'
       
        # template.product entry point defined in setup.json
        self._default_parser = 'kkr.voroparser'
        
        # File names
        self._ATOMINFO = 'atominfo.txt'
        self._RADII = 'radii.dat'
        self._SHAPEFUN = 'shapefun'
        self._VERTICES = 'vertices.dat'
        self._OUT_POTENTIAL_voronoi = 'output.pot'

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
                'linkname': 'structure',
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
            print parameters
        except KeyError:
            raise InputValidationError("No parameters specified for this "
                                       "calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters not of type "
                                       "ParameterData")
        try:
            structure = inputdict.pop(self.get_linkname('structure'))
        except KeyError:
            raise InputValidationError("No structure specified for this "
                                       "calculation")
        if not isinstance(structure, StructureData):
            raise InputValidationError("structure not of type "
                                        "StructureData")
        
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this "
                                       "calculation")
        if inputdict:
                raise ValidationError("Unknown inputs: {}".format(inputdict))


        ###################################
        # Check for 2D case
        twoDimcheck, msg = check_2Dinput(structure, parameters)
        if not twoDimcheck:
            raise InputValidationError(msg)
        
        # Prepare inputcard from Structure and input parameter data
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        generate_inputcard_from_structure(parameters, structure, input_filename)


        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = []
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = [self._OUTPUT_FILE_NAME, self._ATOMINFO, 
                                  self._RADII, self._SHAPEFUN, self._VERTICES, 
                                  self._OUT_POTENTIAL_voronoi, 
                                  self._INPUT_FILE_NAME]

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo
