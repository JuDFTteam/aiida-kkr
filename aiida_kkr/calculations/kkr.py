# -*- coding: utf-8 -*-
"""
Input plug-in for a KKR calculation.
"""
import os
from scipy import array

from aiida.orm.calculation.job import JobCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.utils import classproperty
from aiida.common.constants import elements as PeriodicTableElements
from aiida.common.exceptions import (InputValidationError, ValidationError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.orm import DataFactory
from aiida.common.exceptions import UniquenessError
from aiida_kkr.tools.kkrcontrol import write_kkr_inputcard_template, fill_keywords_to_inputcard, create_keyword_default_values

RemoteData = DataFactory('remote')
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
a_to_bohr = 1.8897261254578281


class KkrCalculation(JobCalculation):
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
       
       
       # Default input and output files
        self._DEFAULT_INPUT_FILE = 'inputcard' # will be shown with inputcat
        self._DEFAULT_OUTPUT_FILE = 'out_kkr' #'shell output will be shown with outputca
        
        self._OUTPUT_FILE_NAME = 'out_kkr'

        # List of mandatory input files
        self._INPUT_FILE_NAME = 'inputcard'
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
        self._default_parser = 'kkr.kkrparser'
        
        # files that will be copied from local computer if parent was KKR calc
        self._copy_filelist_kkr = [self._SHAPEFUN, self._OUT_POTENTIAL]

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
            "parent_folder": {
                'valid_types': RemoteData,
                'additional_parameter': None,
                'linkname': 'parent_calc_folder',
                'docstring': (
                    "Use a remote or local repository folder as parent folder "
                    "(also for restarts and similar). It should contain all the "
                    "needed files for a KKR calc, only edited files should be "
                    "uploaded from the repository.")
            },
            })
            #"structure": {
            #    'valid_types': StructureData,
            #    'additional_parameter': None,
            #    'linkname': 'structure',
            #    'docstring':
            #    ("Use a node that specifies the input crystal structure ")
            #},
        return use_dict

    def _prepare_for_submission(self, tempfolder, inputdict):
        """
        Create input files.

            :param tempfolder: aiida.common.folders.Folder subclass where
                the plugin should put all its files.
            :param inputdict: dictionary of the input nodes as they would
                be returned by get_inputs_dict
        """
        
        has_parent = False        
        local_copy_list = []
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
        try:
            parent_calc_folder = inputdict.pop(self.get_linkname('parent_folder'))
        except KeyError:
            raise InputValidationError("Voronoi files needed for KKR calculation, "
                                       "you need to provide a Parent Folder/RemoteData node.")
                                       
        if not isinstance(parent_calc_folder, RemoteData):
            raise InputValidationError("parent_calc_folder, if specified,"
                                           "must be of type RemoteData")

        # extract parent calculation
        parent_calcs = parent_calc_folder.get_inputs(node_type=JobCalculation)
        n_parents = len(parent_calcs)
        if n_parents != 1:
            raise UniquenessError(
                    "Input RemoteData is child of {} "
                    "calculation{}, while it should have a single parent"
                    "".format(n_parents, "" if n_parents == 0 else "s"))
            parent_calc = parent_calcs[0]
            has_parent = True
        if n_parents:
            parent_calc = parent_calcs[0]
            has_parent = True          
        
        # check that it is a valid parent
        #self._check_valid_parent(parent_calc)


        # if voronoi calc do
        # check if folder from db given, or get folder from rep.
        # Parent calc does not has to be on the same computer.
        #TODO so far we copy every thing from local computer ggf if kkr we want to copy remotely


                
        # get StructureData node from Parent if Voronoi
        structure = None        
        if isinstance(parent_calc, VoronoiCalculation):
            try:            
                structure = parent_calc.get_inputs_dict()['structure']    
            except KeyError:
                # raise InputvaluationError # TODO raise some error
                pass
                print('Could not get structure from parent.')
            
        if inputdict:
            raise ValidationError("Unknown inputs")


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
        # Prepare keywords for kkr from input

        # TODO: policy, override from user input, or from previous calculation?
        # my opinion, user input should override..., because otherwise in db
        # you have input was EMIN = x und used was EMIN=Y
        input_dict = parameters.get_dict()
        keywords = create_keyword_default_values()
        for key, val in input_dict.iteritems():
            keywords[key] = val 
            # TODO IF the input node scheme is changed from [val, format] to val this needs to be changed

        # we always overwride these keys, and since we start from the defaults they are present
        keywords['NATYP'][0] = natyp
        keywords['ALATBASIS'][0] = alat
        
        # Set some keywords from voronoi output
        emin = parent_calc.res.EMIN
        keywords['EMIN'][0] = emin

        # Write input to file
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        write_kkr_inputcard_template(bravais, natyp, positions, charges, outfile=input_filename)
        fill_keywords_to_inputcard(keywords, runops=[], testops=[], CPAconc=[], template=input_filename, output=input_filename)

        #################
        # Decide what files to copy
        if has_parent:
            # copy the right files #TODO check first if file, exists and throw
            # warning, now this will throw an error
            outfolderpath = parent_calc.out.retrieved.folder.abspath
            self.logger.info("out folder path {}".format(outfolderpath))
            
            copylist = []
            if isinstance(parent_calc, KkrCalculation):
                copylist = self._copy_filelist_kkr
                # TODO ggf copy remotely...
            if isinstance(parent_calc, VoronoiCalculation):
                copylist = [parent_calc._SHAPEFUN, 
                            parent_calc._OUT_POTENTIAL_voronoi]              
            
            for file1 in copylist:
                filename = file1
                if file1 == 'output.pot' or file1 == self._OUT_POTENTIAL:
                    filename = self._POTENTIAL
                local_copy_list.append((
                        os.path.join(outfolderpath, 'path', file1),
                        os.path.join(filename)))
            # TODO different copy lists, depending on the keywors input


        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = []
        
        # TODO retrieve list needs some logic, retrieve certain files, 
        # only if certain input keys are specified....
        calcinfo.retrieve_list = [self._DEFAULT_OUTPUT_FILE, 
                                  self._INPUT_FILE_NAME,
                                  self._POTENTIAL,
                                  self._SHAPEFUN,
                                  self._SCOEF,
                                  self._NONCO_ANGLES,
                                  self._OUT_POTENTIAL,
                                  self._OUTPUT_0_INIT,
                                  self._OUTPUT_000,
                                  self._OUT_TIMING_000]

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.code_uuid = code.uuid
        codeinfo.stdout_name = self._DEFAULT_OUTPUT_FILE
        calcinfo.codes_info = [codeinfo]

        return calcinfo


    def _check_valid_parent(self, calc):
        """
        Check that calc is a valid parent for a FleurCalculation.
        It can be a VoronoiCalculation, KKRCalculation
        """

        try:
            if (((not isinstance(calc, VoronoiCalculation)))
                            and (not isinstance(calc, KkrCalculation))):
                raise ValueError("Parent calculation must be a VoronoiCalculation, a "
                                 "KkrCalculation or a CopyonlyCalculation")
        except ImportError:
            if ((not isinstance(calc, VoronoiCalculation))
                            and (not isinstance(calc, KkrCalculation)) ):
                raise ValueError("Parent calculation must be a VoronoiCalculation or "
                                 "a KkrCalculation")



    def use_parent_calculation(self, calc):
        """
        Set the parent calculation of KKR,
        from which it will inherit the outputsubfolder.
        The link will be created from parent RemoteData to KkrCalculation
        """
        from aiida.common.exceptions import NotExistent

        self._check_valid_parent(calc)

        remotedatas = calc.get_outputs(type=RemoteData)
        if not remotedatas:
            raise NotExistent("No output remotedata found in "
                                  "the parent")
        if len(remotedatas) != 1:
            raise UniquenessError("More than one output remotedata found in "
                                  "the parent")
        remotedata = remotedatas[0]

        self._set_parent_remotedata(remotedata)


    def _set_parent_remotedata(self,remotedata):
        """
        Used to set a parent remotefolder in the restart of fleur.
        """
        if not isinstance(remotedata,RemoteData):
            raise ValueError('remotedata must be a RemoteData')

        # complain if another remotedata is already found
        input_remote = self.get_inputs(node_type=RemoteData)
        if input_remote:
            raise ValidationError("Cannot set several parent calculation to a "
                                  "KKR calculation")

        self.use_parent_folder(remotedata)
