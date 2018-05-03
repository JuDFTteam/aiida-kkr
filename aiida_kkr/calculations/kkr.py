# -*- coding: utf-8 -*-
"""
Input plug-in for a KKR calculation.
"""

import os

from aiida.orm.calculation.job import JobCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.orm import DataFactory
from aiida.common.exceptions import UniquenessError
from aiida_kkr.tools.common_workfunctions import (generate_inputcard_from_structure,
                                                  check_2Dinput_consistency, update_params_wf)

#define aiida structures from DataFactory of aiida
RemoteData = DataFactory('remote')
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.4"
__contributors__ = ("Jens Broeder", "Philipp Rüßmann")



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
        
        # calculation plugin version
        self._CALCULATION_PLUGIN_VERSION = __version__
       
        # Default input and output files
        self._DEFAULT_INPUT_FILE = 'inputcard' # will be shown with inputcat
        self._DEFAULT_OUTPUT_FILE = 'out_kkr'  # verdi shell output will be shown with outputcat
        
        # same as _DEFAULT_OUTPUT_FILE: piped output of kkr execution to this file
        self._OUTPUT_FILE_NAME = self._DEFAULT_OUTPUT_FILE

        # List of mandatory input files
        self._INPUT_FILE_NAME = self._DEFAULT_INPUT_FILE
        self._POTENTIAL = 'potential'

        # List of optional input files (may be mandatory for some settings in inputcard)
        self._SHAPEFUN = 'shapefun' # mandatory if nonspherical calculation
        self._SCOEF = 'scoef' # mandatory for KKRFLEX calculation and some functionalities
        self._NONCO_ANGLES = 'nonco_angles.dat' # mandatory if noncollinear directions are used that are not (theta, phi)= (0,0) for all atoms
        self._NONCO_ANGLES_IMP = 'nonco_angles_imp.dat' # mandatory for GREENIMP option (scattering code)
        self._SHAPEFUN_IMP = 'shapefun_imp' # mandatory for GREENIMP option (scattering code)
        self._POTENTIAL_IMP = 'potential_imp' # mandatory for GREENIMP option (scattering code)
	
	   # List of output files that should always be present
        self._OUT_POTENTIAL = 'out_potential'
        self._OUTPUT_0_INIT = 'output.0.txt'
        self._OUTPUT_000 = 'output.000.txt'
        self._OUTPUT_2 = 'output.2.txt'
        self._OUT_TIMING_000 = 'out_timing.000.txt'
        self._NONCO_ANGLES_OUT = 'nonco_angles_out.dat'
        
        # special files (some runs)
        # DOS files
        self._COMPLEXDOS = 'complex.dos'
        self._DOS_ATOM = 'dos.atom%i'
        self._LMDOS = 'lmdos.%2i.%i.dat'
        # kkrflex files for impurity calculation
        self._KKRFLEX_GREEN = 'kkrflex_green'
        self._KKRFLEX_TMAT = 'kkrflex_tmat'
        self._KKRFLEX_ATOMINFO = 'kkrflex_atominfo'
        self._KKRFLEX_INTERCELL_REF = 'kkrflex_intercell_ref'
        self._KKRFLEX_INTERCELL_CMOMS = 'kkrflex_intercell_cmoms'
        self._ALL_KKRFLEX_FILES = [self._KKRFLEX_GREEN, self._KKRFLEX_TMAT, self._KKRFLEX_ATOMINFO, self._KKRFLEX_INTERCELL_REF, self._KKRFLEX_INTERCELL_CMOMS]
        
        # template.product entry point defined in setup.json
        self._default_parser = 'kkr.kkrparser'
        
        # files that will be copied from local computer if parent was KKR calc
        self._copy_filelist_kkr = [self._SHAPEFUN, self._OUT_POTENTIAL]
        
        # list of keywords that are not allowed to be modified (new calculation 
        # starting from structure and voronoi run is needed instead):
        self._do_never_modify = ['ALATBASIS', 'BRAVAIS', 'NAEZ', '<RBASIS>', 'CARTESIAN', 
                                 'INTERFACE', '<NLBASIS>', '<RBLEFT>', 'ZPERIODL', 
                                 '<NRBASIS>', '<RBRIGHT>', 'ZPERIODR', 'KSHAPE', '<SHAPE>', 
                                 '<ZATOM>', 'NATYP', '<SITE>', '<CPA-CONC>', '<KAOEZL>', '<KAOEZR>']
        #TODO implement workfunction to modify structure (e.g. to use VCA)

        
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
            "impurity_info": {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'impurity_info',
                'docstring': ("Use a ParameterNode that specifies properties "
                              "for a follwoing impurity calculation (e.g. setting "
                              "of impurity cluster in scoef file that is "
                              "automatically created).")
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
        
        has_parent = False        
        local_copy_list = []
        # Check inputdict
        try:
            parameters = inputdict.pop(self.get_linkname('parameters'))
        except KeyError:
            raise InputValidationError("No parameters specified for this calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters not of type ParameterData")
            
        try:
            imp_info = inputdict.pop(self.get_linkname('impurity_info'))
            found_imp_info = True
        except KeyError:
            imp_info = None
            found_imp_info = False
        if found_imp_info and not isinstance(imp_info, ParameterData):
            raise InputValidationError("impurity_info not of type ParameterData")
            
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this calculation")
            
        try:
            parent_calc_folder = inputdict.pop(self.get_linkname('parent_folder'))
        except KeyError:
            raise InputValidationError("Voronoi or previous KKR files needed for KKR calculation, "
                                       "you need to provide a Parent Folder/RemoteData node.")
                                  
        #TODO deal with data from folder data if calculation is continued on a different machine
        if not isinstance(parent_calc_folder, RemoteData):
            raise InputValidationError("parent_calc_folder must be of type RemoteData")

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
        if n_parents == 1:
            parent_calc = parent_calcs[0]
            has_parent = True         
        
        # check if parent is either Voronoi or previous KKR calculation
        self._check_valid_parent(parent_calc)
        
        # extract parent input parameter dict for following check
        try:
            parent_inp_dict = parent_calc.inp.parameters.get_dict()
        except:
            self.logger.error("Failed trying to find input parameter of parent {}".format(parent_calc))
            raise InputValidationError("No parameter node found of parent calculation.")
            
        # check if no keys are illegally overwritten (i.e. compare with keys in self._do_never_modify)
        for key in parameters.get_dict().keys():
            value = parameters.get_dict()[key]
            #self.logger.info("Checking {} {}".format(key, value))
            if not value is None:
                if key in self._do_never_modify:
                    oldvalue = parent_inp_dict[key]
                    if value != oldvalue:
                        self.logger.error("You are trying to set keyword {} = {} but this is not allowed since the structure would be modified. Please use a suitable workfunction instead.".format(key, value))
                        raise InputValidationError("You are trying to modify a keyword that is not allowed to be changed!")


        #TODO check for remote folder (starting from folder data not implemented yet)
        # if voronoi calc check if folder from db given, or get folder from rep.
        # Parent calc does not has to be on the same computer.
        # so far we copy every thing from local computer ggf if kkr we want to copy remotely

                
        # get StructureData node from Parent if Voronoi
        structure = None        
        self.logger.info("KkrCalculation: Get structure node from voronoi parent")
        if isinstance(parent_calc, VoronoiCalculation):
            self.logger.info("KkrCalculation: Parent is Voronoi calculation")
            try:            
                structure, voro_parent = VoronoiCalculation.find_parent_structure(parent_calc) 
            except:
                self.logger.error('KkrCalculation: Could not get structure from Voronoi parent.')
                raise ValidationError("Cound not find structure node")
        elif isinstance(parent_calc, KkrCalculation):
            self.logger.info("KkrCalculation: Parent is KKR calculation")
            try:            
                self.logger.info('KkrCalculation: extract structure from KKR parent')
                structure, voro_parent = VoronoiCalculation.find_parent_structure(parent_calc) 
            except:
                self.logger.error('Could not get structure from parent.')
                raise ValidationError('Cound not find structure node starting from parent {}'.format(parent_calc))
        else:
            self.logger.error("KkrCalculation: Parent is neither Voronoi nor KKR calculation!")
            raise ValidationError('Cound not find structure node')
            
        if inputdict:
            self.logger.error('KkrCalculation: Unknown inputs for structure lookup')
            raise ValidationError("Unknown inputs")


        ###################################
        
        # Check for 2D case
        twoDimcheck, msg = check_2Dinput_consistency(structure, parameters)
        if not twoDimcheck:
            raise InputValidationError(msg)
            
        # set shapes array either from parent voronoi run or read from inputcard in kkrimporter calculation
        if parent_calc.get_parser_name() != 'kkr.kkrimporterparser':
            # get shapes array from voronoi parent
            shapes = voro_parent.res.shapes
        else:
            # extract shapes from input parameters node constructed by kkrimporter calculation
            shapes = voro_parent.inp.parameters.get_dict().get('<SHAPE>')
        
        # Prepare inputcard from Structure and input parameter data
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        natom, nspin, newsosol = generate_inputcard_from_structure(parameters, structure, input_filename, parent_calc, shapes=shapes)

        # prepare scoef file if impurity_info was given
        write_scoef = False
        runopt = parameters.get_dict().get('RUNOPT', None)
        kkrflex_opt = False
        if runopt is not None:
            if 'KKRFLEX' in runopt:
                kkrflex_opt = True
        if kkrflex_opt:
            write_scoef = True
        elif found_imp_info:
            self.logger.info('Found impurity_info in inputs of the calculation, automatically add runopt KKRFLEX')
            write_scoef = True
            runopt = parameters.get_dict().get('RUNOPT', [])
            runopt.append('KKRFLEX')
            parameters = update_params_wf(parameters, ParameterData(dict={'RUNOPT':runopt, 'nodename': '', 'nodedesc':''}))
        if found_imp_info and write_scoef:
            # TODO check completeness of impurity info!
            # TODO implement this!
            # placeholder: take Fabian's functions later on
            scoef = ['      13\n',
                     '  0.0000000000000000000E+00  0.0000000000000000000E+00  0.0000000000000000000E+00    1 29.0  0.000000000E+00\n',
                     '  0.2220446049250313081E-14 -0.7071067811865453523E+00 -0.7071067811865453523E+00    1 29.0  0.100000000E+01\n',
                     ' -0.7071067811865453523E+00  0.2220446049250313081E-14 -0.7071067811865453523E+00    1 29.0  0.100000000E+01\n',
                     '  0.7071067811865497932E+00  0.2220446049250313081E-14 -0.7071067811865453523E+00    1 29.0  0.100000000E+01\n',
                     '  0.2220446049250313081E-14  0.7071067811865497932E+00 -0.7071067811865453523E+00    1 29.0  0.100000000E+01\n',
                     ' -0.7071067811865453523E+00 -0.7071067811865453523E+00  0.2220446049250313081E-14    1 29.0  0.100000000E+01\n',
                     '  0.7071067811865497932E+00 -0.7071067811865453523E+00  0.2220446049250313081E-14    1 29.0  0.100000000E+01\n',
                     ' -0.7071067811865453523E+00  0.7071067811865497932E+00  0.2220446049250313081E-14    1 29.0  0.100000000E+01\n',
                     '  0.7071067811865497932E+00  0.7071067811865497932E+00  0.2220446049250313081E-14    1 29.0  0.100000000E+01\n',
                     '  0.2220446049250313081E-14 -0.7071067811865453523E+00  0.7071067811865497932E+00    1 29.0  0.100000000E+01\n',
                     ' -0.7071067811865453523E+00  0.2220446049250313081E-14  0.7071067811865497932E+00    1 29.0  0.100000000E+01\n',
                     '  0.7071067811865497932E+00  0.2220446049250313081E-14  0.7071067811865497932E+00    1 29.0  0.100000000E+01\n',
                     '  0.2220446049250313081E-14  0.7071067811865497932E+00  0.7071067811865497932E+00    1 29.0  0.100000000E+01\n']
            scoef_filename = os.path.join(tempfolder.get_abs_path(''), self._SCOEF)
            self.logger.info('Writing scoef file {}'.format(scoef_filename))
            with open(scoef_filename, 'w') as file:
                file.writelines(scoef)
        elif write_scoef:
            self.logger.info('Need to write scoef file but no impurity_info given!')
            raise ValidationError('Found RUNOPT KKRFLEX but no impurity_info in inputs')

            

        #################
        # Decide what files to copy based on settings to the code (e.g. KKRFLEX option needs scoef)
        if has_parent:
            # copy the right files #TODO check first if file, exists and throw
            # warning, now this will throw an error
            outfolderpath = parent_calc.out.retrieved.folder.abspath
            outfolderpath = os.path.join(outfolderpath, 'path')
            self.logger.info("out folder path {}".format(outfolderpath))
            
            copylist = []
            if isinstance(parent_calc, KkrCalculation):
                copylist = self._copy_filelist_kkr
                # TODO ggf copy remotely...
                
            if isinstance(parent_calc, VoronoiCalculation):
                copylist = [parent_calc._SHAPEFUN]
                # copy either overwrite potential or voronoi output potential 
                # (voronoi caclualtion retreives only one of the two)
                if parent_calc._POTENTIAL_IN_OVERWRITE in os.listdir(outfolderpath):
                    copylist.append(parent_calc._POTENTIAL_IN_OVERWRITE)
                else:
                    copylist.append(parent_calc._OUT_POTENTIAL_voronoi)
                    
            #change copylist in case the calculation starts from an imported calculation
            if parent_calc.get_parser_name() == 'kkr.kkrimporterparser':
                copylist = []
                if not os.path.exists(os.path.join(outfolderpath, self._OUT_POTENTIAL)):
                    copylist.append(self._POTENTIAL)
                else:
                    copylist.append(self._OUT_POTENTIAL)
                if os.path.exists(os.path.join(outfolderpath, self._SHAPEFUN)):
                    copylist.append(self._SHAPEFUN)
        
            # create local_copy_list from copylist and change some names automatically
            for file1 in copylist:
                filename = file1
                if (file1 == 'output.pot' or file1 == self._OUT_POTENTIAL or 
                   (isinstance(parent_calc, VoronoiCalculation) and file1 == parent_calc._POTENTIAL_IN_OVERWRITE)):
                    filename = self._POTENTIAL
                local_copy_list.append((
                        os.path.join(outfolderpath, file1),
                        os.path.join(filename)))
            # TODO different copy lists, depending on the keywors input
            self.logger.info('local copy list: {}'.format(local_copy_list))


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
                                  self._NONCO_ANGLES_OUT,
                                  self._OUT_POTENTIAL,
                                  self._OUTPUT_0_INIT,
                                  self._OUTPUT_000,
                                  self._OUTPUT_2,
                                  self._OUT_TIMING_000]
        
        # for special cases add files to retireve list:
            
        # 1. dos calculation, add *dos* files if NPOL==0
        retrieve_dos_files = False
        print('NPOL in parameter input:', parameters.get_dict()['NPOL'])
        if 'NPOL' in  parameters.get_dict().keys():
            if parameters.get_dict()['NPOL'] == 0:
                retrieve_dos_files = True
        if 'TESTOPT' in  parameters.get_dict().keys():
            testopts = parameters.get_dict()['TESTOPT']
            if testopts is not None :
                stripped_test_opts = [i.strip() for i in testopts]
                if 'DOS' in stripped_test_opts:
                    retrieve_dos_files = True
        if retrieve_dos_files:
            print('adding files for dos output', self._COMPLEXDOS, self._DOS_ATOM, self._LMDOS)
            add_files = [self._COMPLEXDOS]
            for iatom in range(natom):
                add_files.append(self._DOS_ATOM%(iatom+1))
                for ispin in range(nspin):
                    add_files.append((self._LMDOS%(iatom+1, ispin+1)).replace(' ','0'))
            print(add_files)
            calcinfo.retrieve_list += add_files
            
        # 2. KKRFLEX calculation
        retrieve_kkrflex_files = False
        if 'RUNOPT' in  parameters.get_dict().keys():
            runopts = parameters.get_dict()['RUNOPT']
            if runopts is not None :
                stripped_run_opts = [i.strip() for i in runopts]
                if 'KKRFLEX' in stripped_run_opts:
                    retrieve_kkrflex_files = True
        if retrieve_kkrflex_files:
            add_files = self._ALL_KKRFLEX_FILES
            print('adding files for KKRFLEX output', add_files)
            calcinfo.retrieve_list += add_files
            

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
                raise ValueError("Parent calculation must be a VoronoiCalculation or a KkrCalculation")
        except ImportError:
            if ((not isinstance(calc, KkrCalculation)) ):
                raise ValueError("Parent calculation must be a VoronoiCalculation or a KkrCalculation")


    def _set_parent_remotedata(self, remotedata):
        """
        Used to set a parent remotefolder in the restart of fleur.
        """
        if not isinstance(remotedata,RemoteData):
            raise ValueError('remotedata must be a RemoteData')

        # complain if another remotedata is already found
        input_remote = self.get_inputs(node_type=RemoteData)
        if input_remote:
            raise ValidationError("Cannot set several parent calculation to a KKR calculation")

        self.use_parent_folder(remotedata)

        
