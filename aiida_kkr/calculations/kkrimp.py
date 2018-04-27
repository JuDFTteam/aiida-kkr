import json

from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError, UniquenessError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.orm import DataFactory
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.calculations.kkr import KkrCalculation
import os

ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
SinglefileData = DataFactory('singlefile')

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
        super(KkrimpCalculation, self)._init_internal_params()

        # List of mandatory input files
        self._CONFIG = 'config.cfg'
        self._POTENTIAL = 'potential'
        self._KKRFLEX_GREEN = KkrCalculation()._KKRFLEX_GREEN
        self._KKRFLEX_TMAT = KkrCalculation()._KKRFLEX_TMAT
        self._KKRFLEX_ATOMINFO = KkrCalculation()._KKRFLEX_ATOMINFO
        self._KKRFLEX_INTERCELL_REF = KkrCalculation()._KKRFLEX_INTERCELL_REF
        self._KKRFLEX_INTERCELL_CMOMS = KkrCalculation()._KKRFLEX_INTERCELL_CMOMS

        # List of optional input files (may be mandatory for some setting in inputputcard)
        self._SHAPEFUN = 'shapefun'
        self._KKRFLEX_ANGLE = 'kkrflex_angle'
        self._KKRFLEX_LLYFAC = 'kkrflex_llyfac'
        
        # full list of kkrflex files
        self._ALL_KKRFLEX_FILES = KkrCalculation()._ALL_KKRFLEX_FILES.append(self._KKRFLEX_ANGLE, self._KKRFLEX_LLYFAC)

        # List of output files that should always be present
        self._OUT_POTENTIAL = 'out_potential'
        self._OUTPUT_000 = 'output.000.txt'
        self._OUT_TIMING_000 = 'out_timing.000.txt'


        # template.product entry point defined in setup.json
        self._default_parser = 'kkr.kkrimpparser'

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
                ("Use a node that specifies the input parameters (calculation settings)")
                },
            "impurity_potential": {
                'valid_types': SinglefileData,
                'additional_parameter': None,
                'linkname': 'potential',
                'docstring':
                ("Use a node contains the input potential")
                },
            "parent_calc_folder": {
                'valid_types': RemoteData,
                'additional_parameter': None,
                'linkname': 'parent_calc_folder',
                'docstring':
                ("Use a node that specifies a parent KKRimp calculation")
                },
            "host_Greenfunction_folder": {
                'valid_types': RemoteData,
                'additional_parameter': None,
                'linkname': 'GFhost_folder',
                'docstring':
                ("Use a node that specifies the host KKR calculation contaning the host Green function and tmatrix (KkrCalculation with impurity_info input)")
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
            raise InputValidationError("No parameters specified for this calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters not of type ParameterData")
            
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this calculation")
            
        imp_info, kkrflex_file_paths, shapfun_path = self._get_and_verify_hostfiles()
        
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
        
        
        def _get_and_verify_hostfiles(self, inputdict):
            """
            Check inputdict for host_Greenfunction_folder and extract impurity_info, paths to kkrflex-files and path of shapefun file
            
            :param inputdict: input dictionary containing all input nodes to KkrimpCalculation
            :returns: 
                * imp_info: ParameterData node containing impurity information like position, Z_imp, cluster size, etc.
                * kkrflex_file_paths: dict of absolute file paths for the kkrflex files
                * shapefun_path: absolute path of the shapefunction file in the host calculation (needed to construct shapefun_imp)
                * shapes: mapping array of atoms to shapes (<SHAPE> input)
            :note: shapefun_path is None if host_Greenfunction calculation was not full-potential
            :raises: 
                * InputValidationError, if inputdict does not contain 'host_Greenfunction'
                * 
            """
            # get host_parent node and check consistency
            try:
                host_parent = inputdict.pop(self.get_linkname('host_Greenfunction_folder'))
            except KeyError:
                raise InputValidationError("No host_Greenfunction_folder specified for this calculation")
            # check host parent type
            if not isinstance(host_parent, RemoteData):
                raise InputValidationError("host_Greenfunction_folder not of type RemoteData")
            # extract parent calculation
            parent_calcs = host_parent.get_inputs(node_type=JobCalculation)
            n_parents = len(parent_calcs)
            if n_parents != 1:
                raise UniquenessError(
                        "Input RemoteData is child of {} "
                        "calculation{}, while it should have a single parent"
                        "".format(n_parents, "" if n_parents == 0 else "s"))
            else:
                parent_calc = parent_calcs[0]
                
            # extract impurity_info
            imp_info = parent_calc.get_inputs_dict().get('impurity_info', None)
            if imp_info is None:
                raise InputValidationError("host_Greenfunction calculation does not have an input node impurity_info")
                
            # check if host parent was KKRFLEX calculation
            hostfolderpath = parent_calc.out.retrieved.folder.abspath
            input_file = os.path.join(hostfolderpath, 'path', KkrCalculation()._DEFAULT_INPUT_FILE)
            params_host_calc = kkrparams().read_keywords_from_inputcard(input_file)
            
            if 'RUNOPT' not in params_host_calc.get_dict().keys():
                host_ok = False
            elif 'KKRFLEX' not in params_host_calc.get_dict().get('RUNOPT', []):
                host_ok = False
            else:
                host_ok = True
                
            if not host_ok:
                raise InputValidationError("host_Greenfunction calculation was not a KKRFLEX run")
            
            # extract absolute paths of kkrflex_* files
            if 
            
            # extract absolute path of host shapefun and shapes array
                
            return imp_info, kkrflex_file_paths, shapefun_path, shapes

