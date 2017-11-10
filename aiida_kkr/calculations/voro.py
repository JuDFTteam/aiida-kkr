# -*- coding: utf-8 -*-
"""
Input plug-in for a voronoi calculation.
"""

from scipy import array
from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.common.constants import elements as PeriodicTableElements
from aiida.orm import DataFactory
from aiida_kkr.tools.kkrcontrol import write_kkr_inputcard_template, fill_keywords_to_inputcard, create_keyword_default_values

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
a_to_bohr = 1.8897261254578281

class VoronoiCalculation(JobCalculation):
    """
    AiiDA calculation plugin for a voronoi calculation (creation of starting potential and shapefun)
    .
    """

    def _init_internal_params(self):
        """
        Init internal parameters at class load time
        """
        # reuse base class function
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
        self._ATOMINFO = 'atominfo.dat'
        self._RADII = 'radii.dat'
        self._SHAPEFUN = 'shapefun'
        self._VERTIVES = 'vertices.dat'
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
        # Prepare Structure

        # Get the connection between coordination number and element symbol
        # maybe do in a differnt way
        
        _atomic_numbers = {data['symbol']: num for num,
                        data in PeriodicTableElements.iteritems()}
        
        # KKr wants units in bohr and relativ coordinates
        bravais = array(structure.cell)*a_to_bohr
        alat = bravais.max()
        bravais = bravais/alat
        sites = structure.sites
        natyp = len(sites)
        positions = []
        charges = []
        for site in sites:
            pos = site.position 
            #TODO convert to rel pos and make sure that type is rigth for script (array or tuple)
            relpos = array(pos) 
            positions.append(relpos)
            sitekind = structure.get_kind(site.kind_name)
            site_symbol = sitekind.symbol
            charges.append(_atomic_numbers[site_symbol])
            # TODO does not work for Charged atoms, find out how...
        # TODO get empty spheres
        positions = array(positions)
        
        ######################################
        # Prepare keywords for kkr
        input_dict = parameters.get_dict()
        keywords = create_keyword_default_values()
        for key, val in input_dict.iteritems():
            keywords[key] = val 
            # TODO IF the input node scheme is changed from [val, format] to val this needs to be changed

        # we always overwride these keys:
        keywords['NATYP'][0] = natyp
        keywords['ALATBASIS'][0] = alat

        # Write input to file
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        write_kkr_inputcard_template(bravais, natyp, positions, charges, outfile=input_filename)#'inputcard.tmpl')
        print input_filename
        fill_keywords_to_inputcard(keywords, runops=[], testops=[], CPAconc=[], template=input_filename, output=input_filename)


        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = []
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = [self._OUTPUT_FILE_NAME, self._ATOMINFO, self._RADII,
                                        self._SHAPEFUN, self._VERTIVES, self._OUT_POTENTIAL_voronoi]

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo
