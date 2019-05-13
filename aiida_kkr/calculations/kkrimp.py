# -*- coding: utf-8 -*-
"""
Input plug-in for a KKRimp calculation.
"""
from __future__ import print_function

from builtins import range
from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.common.exceptions import (InputValidationError, ValidationError, UniquenessError)
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.orm import DataFactory
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.tools.tools_kkrimp import modify_potential
from aiida_kkr.tools.tools_kkrimp import make_scoef
from masci_tools.io.common_functions import search_string
from aiida_kkr.calculations.voro import VoronoiCalculation
import os
from numpy import array, sqrt, sum, where

ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
SinglefileData = DataFactory('singlefile')


__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.3"
__contributors__ = ("Philipp Rüßmann")

#TODO: implement 'ilayer_center' consistency check

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
        
        # calculation plugin version
        self._CALCULATION_PLUGIN_VERSION = __version__

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
        self._ALL_KKRFLEX_FILES = KkrCalculation()._ALL_KKRFLEX_FILES
        self._ALL_KKRFLEX_FILES.append(self._KKRFLEX_ANGLE)
        self._ALL_KKRFLEX_FILES.append(self._KKRFLEX_LLYFAC)

        # List of output files that should always be present
        self._OUT_POTENTIAL = 'out_potential'
        self._OUTPUT_000 = 'out_log.000.txt'
        self._OUT_TIMING_000 = 'out_timing.000.txt'
        self._OUT_ENERGYSP = 'out_energysp_eV'
        self._OUT_ENERGYSP_PER_ATOM = 'out_energysp_per_atom_eV'
        self._OUT_ENERGYTOT = 'out_energytotal_eV'
        self._OUT_ENERGYTOT_PER_ATOM =  'out_energytotal_per_atom_eV'
        
        # Lift of output files that are retrieved if special conditions are fulfilled
        self._OUT_JIJDIJ = 'out_JijDij'
        self._OUT_JIJDIJ_LOCAL = 'out_JijDij_local'
        self._OUT_JIJMAT = 'out_Jijmatrix'
        self._OUT_JIJMAT_LOCAL = 'out_Jijmatrix_local'
        self._OUT_LDOS_BASE = 'out_ldos.atom=%2i_spin%i.dat'
        self._OUT_LDOS_INTERPOL_BASE = 'out_ldos.interpol.atom=%2i_spin%i.dat'
        self._OUT_LMDOS_BASE = 'out_lmdos.atom=%2i_spin%i.dat'
        self._OUT_LMDOS_INTERPOL_BASE = 'out_lmdos.interpol.atom=%2i_spin%i.dat'
        self._OUT_MAGNETICMOMENTS = 'out_magneticmoments'
        self._OUT_ORBITALMOMENTS = 'out_orbitalmoments'

        # template.product entry point defined in setup.json
        self._default_parser = 'kkr.kkrimpparser'
        
        # default names used within aiida (verdi calculation out(in)putcat)
        self._OUTPUT_FILE_NAME = 'out_kkrimp'
        self._INPUT_FILE_NAME = self._CONFIG
        self._DEFAULT_OUTPUT_FILE = self._OUTPUT_FILE_NAME # will be shown with outputcat
        self._DEFAULT_INPUT_FILE = self._INPUT_FILE_NAME # will be shown with inputcat

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
                ("Use a node that specifies the input parameters (calculation settings).")
                },
            "host_Greenfunction_folder": {
                'valid_types': RemoteData,
                'additional_parameter': None,
                'linkname': 'GFhost_folder',
                'docstring':
                ("Use a node that specifies the host KKR calculation contaning "
                 "the host Green function and tmatrix (KkrCalculation with "
                 "impurity_info input).")
                },
            "impurity_potential": {
                'valid_types': SinglefileData,
                'additional_parameter': None,
                'linkname': 'potential',
                'docstring':
                ("Use a node contains the input potential.")
                },
            "parent_calc_folder": {
                'valid_types': RemoteData,
                'additional_parameter': None,
                'linkname': 'parent_calc_folder',
                'docstring':
                ("Use a node that specifies a parent KKRimp calculation.")
                },
            "impurity_info": {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'impurity_info',
                'docstring':
                ("Use a Parameter node that specifies properties "
                 "for a impurity calculation.")
                }
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
        # Check inputdict and extrace nodes
        tmp = self._check_and_extract_input_nodes(inputdict)
        parameters = tmp[0]
        code = tmp[1]
        imp_info = tmp[2]
        kkrflex_file_paths = tmp[3]
        shapefun_path = tmp[4]
        shapes = tmp[5]
        host_parent_calc = tmp[6]
        params_host = tmp[7]
        impurity_potential = tmp[8]
        parent_calc_folder = tmp[9]
        structure = tmp[10]
        
        # Prepare input files for KKRimp calculation
        # 1. fill kkr params for KKRimp, write config file and eventually also kkrflex_llyfac file
        host_file_path = kkrflex_file_paths[self._KKRFLEX_GREEN].split(self._KKRFLEX_GREEN)[0]
        self._extract_and_write_config(parent_calc_folder, params_host, parameters, tempfolder, host_file_path)
        # 2. write shapefun from impurity info and host shapefun and copy imp. potential
        potfile = self._get_pot_and_shape(imp_info, shapefun_path, shapes, impurity_potential, parent_calc_folder, tempfolder, structure)
        # 3. change kkrflex_atominfo to match impurity case
        self._change_atominfo(imp_info, kkrflex_file_paths, tempfolder)
        
        # prepare copy and retrieve lists
        local_copy_list = [(potfile, self._POTENTIAL)]
        for key in list(kkrflex_file_paths.keys()):
            if key!=self._KKRFLEX_ATOMINFO:
                src_path = kkrflex_file_paths[key]
                dst_path = key
                local_copy_list.append((src_path, dst_path))
                
        retrieve_list = [self._OUTPUT_FILE_NAME, self._CONFIG, self._KKRFLEX_ATOMINFO, 
                         self._SHAPEFUN, self._KKRFLEX_ANGLE, self._KKRFLEX_LLYFAC, 
                         self._OUT_POTENTIAL, self._OUTPUT_000, self._OUT_TIMING_000,
                         self._OUT_ENERGYSP_PER_ATOM, self._OUT_ENERGYTOT_PER_ATOM]

        # retrieve l(m)dos files
        runopts = parameters.get_dict().get('RUNFLAG')
        if runopts is None:
            runopts = []
        testopts = parameters.get_dict().get('TESTFLAG')
        if testopts is None:
            testopts = []
        allopts = runopts+testopts
        if 'lmdos' in allopts or 'ldos' in allopts:
            file = open(tempfolder.get_abs_path(self._CONFIG))
            config = file.readlines()
            file.close()
            itmp = search_string('NSPIN', config)
            if itmp>=0:
                nspin = int(config[itmp].split()[-1])
            else:
                raise ValueError("Could not extract NSPIN value from config.cfg")
            file = open(tempfolder.get_abs_path(self._KKRFLEX_ATOMINFO))
            atominfo = file.readlines()
            file.close()
            itmp = search_string('NATOM', atominfo)
            if itmp>=0:
                natom = int(atominfo[itmp+1].split()[0])
            else:
                raise ValueError("Could not extract NATOM value from kkrflex_atominfo")
            for iatom in range(1,natom+1):
                for ispin in range(1,nspin+1):
                    retrieve_list.append((self._OUT_LDOS_BASE%(iatom, ispin)).replace(' ', '0'))
                    retrieve_list.append((self._OUT_LDOS_INTERPOL_BASE%(iatom, ispin)).replace(' ', '0'))
                    retrieve_list.append((self._OUT_LMDOS_BASE%(iatom, ispin)).replace(' ', '0'))
                    retrieve_list.append((self._OUT_LMDOS_INTERPOL_BASE%(iatom, ispin)).replace(' ', '0'))

        file = open(tempfolder.get_abs_path(self._CONFIG))
        config = file.readlines()
        file.close()
        itmp = search_string('NSPIN', config)
        if itmp>=0:
            nspin = int(config[itmp].split()[-1])
        else:
            raise ValueError("Could not extract NSPIN value from config.cfg")
        if 'tmatnew' in allopts and nspin>1:
            retrieve_list.append(self._OUT_MAGNETICMOMENTS) 
            file = open(tempfolder.get_abs_path(self._CONFIG))
            outorb = file.readlines()
            file.close()
            itmp = search_string('CALCORBITALMOMENT', outorb)
            if itmp>=0:
                calcorb = int(outorb[itmp].split()[-1])
            else:
                calcorb = 0
            if calcorb==1:
                retrieve_list.append(self._OUT_ORBITALMOMENTS)
        
        
        # Prepare CalcInfo to be returned to aiida (e.g. retreive_list etc.)
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = retrieve_list

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo
        
    #####################################################
    # helper functions
    
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
            * InputValidationError, if host_Greenfunction_folder not of right type
            * UniquenessError, if host_Greenfunction_folder does not have exactly one parent
            * InputValidationError, if host_Greenfunction does not have an input node impurity_info
            * InputValidationError, if host_Greenfunction was not a KKRFLEX calculation
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
        try:
            imp_info_inputnode = inputdict.pop(self.get_linkname('impurity_info'))
            if not isinstance(imp_info_inputnode, ParameterData):
                raise InputValidationError("impurity_info not of type ParameterData")
            imp_info = parent_calc.get_inputs_dict().get('impurity_info', None)
            if imp_info is None:
                raise InputValidationError("host_Greenfunction calculation does not have an input node impurity_info")
            found_impurity_inputnode = True
            found_host_parent = True
        except KeyError:
            imp_info = parent_calc.get_inputs_dict().get('impurity_info', None)
            if imp_info is None:
                raise InputValidationError("host_Greenfunction calculation does not have an input node impurity_info")
            found_impurity_inputnode = False
        # if impurity input is seperate input, check if it is the same as 
        # the one from the parent calc (except for 'Zimp'). If that's not the 
        # case, raise an error
        if found_impurity_inputnode and found_host_parent:
            #TODO: implement also 'ilayer_center' check
            if imp_info_inputnode.get_attr('Rcut') == imp_info.get_attr('Rcut'):
                check_consistency_imp_info = True
                try:
                    if (imp_info_inputnode.get_attr('hcut') == imp_info.get_attr('hcut')
                        and imp_info_inputnode.get_attr('cylinder_orient') == imp_info.get_attr('cylinder_orient')
                        and imp_info_inputnode.get_attr('Rimp_rel') == imp_info.get_attr('Rimp_rel')
                        and imp_info_inputnode.get_attr('imp_cls') == imp_info.get_attr('imp_cls')):
                        print('impurity_info node from input and from previous GF calculation are compatible')
                        check_consistency_imp_info = True
                    else:
                        print('impurity_info node from input and from previous GF calculation are NOT compatible!. ' 
                              'Please check your impurity_info nodes for consistency.')
                        check_consistency_imp_info = False
                except AttributeError:
                    print("Non default values of the impurity_info node from input and from previous "
                          "GF calculation are compatible. Default values haven't been checked")
                    check_consistency_imp_info = True
            else:
                print('impurity_info node from input and from previous GF calculation are NOT compatible!. ' 
                      'Please check your impurity_info nodes for consistency.')
                check_consistency_imp_info = False           
            if check_consistency_imp_info:
                imp_info = imp_info_inputnode
            else:
                raise InputValidationError("impurity_info nodes (input and GF calc) are not compatible") 
                    
        # check if host parent was KKRFLEX calculation
        hostfolderpath = parent_calc.out.retrieved.folder.abspath
        hostfolderpath = os.path.join(hostfolderpath, 'path')
        input_file = os.path.join(hostfolderpath, KkrCalculation()._DEFAULT_INPUT_FILE)
        params_host_calc = kkrparams(params_type='kkr') # initialize kkrparams instance to use read_keywords_from_inputcard
        params_host_calc.read_keywords_from_inputcard(inputcard=input_file)
        
        if 'RUNOPT' not in list(params_host_calc.get_dict().keys()):
            host_ok = False
        elif 'KKRFLEX' not in params_host_calc.get_dict().get('RUNOPT', []):
            host_ok = False
        else:
            host_ok = True
            
        if not host_ok:
            raise InputValidationError("host_Greenfunction calculation was not a KKRFLEX run")
        
        # extract absolute paths of kkrflex_* files
        kkrflex_file_paths = {}
        for file in self._ALL_KKRFLEX_FILES:
            file_abspath = os.path.join(hostfolderpath, file)
            if os.path.exists(file_abspath):
                kkrflex_file_paths[file] = file_abspath
                
        # extract absolute path of host shapefun
        file_abspath = os.path.join(hostfolderpath, KkrCalculation()._SHAPEFUN)
        if os.path.exists(file_abspath):
            shapefun_path = file_abspath
        else:
            shapefun_path = None
            
        # extract shapes array from parameters read from inputcard
        shapes = params_host_calc.get_dict().get('<SHAPE>', None)
        if type(shapes)==int:
            shapes = [shapes]
            
        # extract input structure
        try:
            structure, voro_parent = VoronoiCalculation.find_parent_structure(parent_calc)
        except:
            structure, voro_parent = None, None
        if structure is None:
            raiseInputValidationError("No structure node found from host GF parent")
            
        return imp_info, kkrflex_file_paths, shapefun_path, shapes, parent_calc, params_host_calc, structure
        
        
    def _check_and_extract_input_nodes(self, inputdict):
        """
        Extract input nodes from inputdict and check consitency of input nodes
        :param inputdict: dict of inputnodes
        :returns: 
            * parameters (aiida_kkr.tools.kkr_params.kkrparams), optional: parameters of KKRimp that end up in config.cfg
            * code (KKRimpCodeNode): code of KKRimp on some machine
            * imp_info (ParameterDataNode): parameter node of the impurity information, extracted from host_parent_calc
            * kkrflex_file_paths (dict): dictionary of {filenames: absolute_path_to_file} for the kkrflex-files
            * shapfun_path (str): absolute path of the shapefunction of the host parent calculation
            * host_parent_calc (KkrCalculation): node of the parent host calculation where the kkrflex-files were created
            * impurity_potential (SinglefileData): single file data node containing the starting potential for the impurity calculation
            * parent_calc_folder (RemoteData): remote directory of a parent KKRimp calculation
        """
        # 1. extract parameters, optional (defaults extracted from host calculation)
        try:
            parameters = inputdict.pop(self.get_linkname('parameters'))
            if not isinstance(parameters, ParameterData):
                raise InputValidationError("parameters not of type ParameterData")
        except KeyError:
            self.logger.info("No parameters specified for this calculation, use defaults")
            parameters = None
        if parameters is not None: # convert to kkrparams instance
            parameters = kkrparams(params_type='kkrimp', **parameters.get_dict())
        # 2. extract code
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this calculation")
        # 3. get hostfiles
        imp_info, kkrflex_file_paths, shapfun_path, shapes, host_parent_calc, params_host, structure = self._get_and_verify_hostfiles(inputdict)
        
        # 4. check impurity potential or parent calc input
        # imp potential
        try:
            impurity_potential = inputdict.pop(self.get_linkname('impurity_potential'))
            if not isinstance(impurity_potential, SinglefileData):
                raise InputValidationError("impurity_potential not of type SinglefileData")
            found_imp_pot = True
        except KeyError:
            impurity_potential = None
            found_imp_pot = False
        # parent calc folder
        try:
            parent_calc_folder = inputdict.pop(self.get_linkname('parent_calc_folder'))
            if not isinstance(parent_calc_folder, RemoteData):
                raise InputValidationError("parent_calc_folder not of type RemoteData")
            found_parent_calc = True
        except KeyError:
            parent_calc_folder = None
            found_parent_calc = False
        # consistency checks
        if not found_parent_calc and not found_imp_pot:
            raise InputValidationError("Neither impurity_potential nor parent_calc_folder specified for this calculation.\n"
                                       "Please provide either impurity_potential or parent_calc_folder.")
        elif found_parent_calc and found_imp_pot:
            raise InputValidationError("Both impurity_potential and parent_calc_folder specified for this calculation.\n"
                                       "Please provide one one, i.e. either impurity_potential or parent_calc_folder.")
        # 5. anything left that is unknown?
        if inputdict:
            raise ValidationError("Unknown inputs: {}".format(inputdict))
        # Done checking inputdict, returning ...
        return parameters, code, imp_info, kkrflex_file_paths, shapfun_path, shapes, host_parent_calc, params_host, impurity_potential, parent_calc_folder, structure


    def _extract_and_write_config(self, parent_calc_folder, params_host, parameters, tempfolder, kkrflex_file_path0):
        """
        fill kkr params for KKRimp and write config file
        also writes kkrflex_llyfac file if Lloyd is used in the host system
        """
        # initialize kkrimp parameter set with default values
        params_kkrimp = kkrparams(params_type='kkrimp', NPAN_LOGPANELFAC=2, 
                                  RADIUS_MIN=-1, NCOLL=0, SPINORBIT=0, SCFSTEPS=1,
                                  IMIX=0, MIXFAC=0.05, ITDBRY=20, BRYMIX=0.05,
                                  QBOUND=10**-7, RUNFLAG=[], TESTFLAG=[], HFIELD=[0.0, 0],
                                  CALCFORCE=0, CALCJIJMAT=0, CALCORBITALMOMENT=0, ICST=2)
        
        # keys that are being overwritten from host calculation settings
        keys_overwrite = ['NSPIN', 'KVREL', 'XC', 'INS', 'ICST', 'RADIUS_LOGPANELS',
                          'NPAN_EQ', 'NPAN_LOG', 'NCHEB', 'QBOUND']
        for key in keys_overwrite:
            if key=='XC':
                key0 = 'KEXCOR'
            elif key=='RADIUS_LOGPANELS':
                key0 = 'R_LOG'
            elif key=='MIXFAC':
                key0 = 'STRMIX'
            else:
                key0 = key
            val = params_host.get_value(key0)
            if val is not None:
                params_kkrimp.set_value(key, val)
        # settings for SOC solver
        runopts = params_host.get_value('RUNOPT')
        if 'NEWSOSOL' in runopts:
            params_kkrimp.set_multiple_values(NCOLL=1, SPINORBIT=1, CALCORBITALMOMENT=1, TESTFLAG=['tmatnew'])
        else:
            params_kkrimp.set_multiple_values(NCOLL=0, SPINORBIT=0, CALCORBITALMOMENT=0, TESTFLAG=[])
        # special settings
        runopts = params_host.get_value('RUNOPT')
        if 'SIMULASA' in runopts or (params_kkrimp.get_value('NCOLL')>0 and params_kkrimp.get_value('INS')==0):
            runflag = ['SIMULASA']
        else:
            runflag = []
        # take care of LLYsimple (i.e. Lloyd in host system)
        if 'LLOYD' in runopts:
            # add runflag for imp code
            runflag.append('LLYsimple')
            # also extract renormalization factor and create kkrflex_llyfac file (contains one value only)
            outfile = os.path.join(kkrflex_file_path0, 'output.000.txt')
            llyfac_file = tempfolder.get_abs_path(self._KKRFLEX_LLYFAC)
            with open(outfile) as f:
                txt = f.readlines()
                iline = search_string('RENORM_LLY: Renormalization factor of total charge', txt)
                if iline>=0:
                    llyfac = txt[iline].split()[-1]
                    # now write kkrflex_llyfac to tempfolder where later on config file is also written
                    with open(llyfac_file, 'w') as f2:
                        f2.writelines([llyfac])

        # now set runflags
        params_kkrimp.set_value('RUNFLAG', runflag)
            
        # overwrite keys if found in parent_calc
        if parent_calc_folder is not None:
            params_parent = parent_calc_folder.get_inputs_dict().get('parameters', None)
        else:
            params_parent = None
        if params_parent is not None:
            params_parent = kkrparams(params_type='kkr', **params_parent.get_dict())
            for (key, val) in params_parent.get_set_values():
                self._check_key_setting_consistency(params_kkrimp, key, val)
                params_kkrimp.set_value(key, val)
        
        # finally overwrite from input parameters
        if parameters is not None:
            for (key, val) in parameters.get_set_values():
                self._check_key_setting_consistency(params_kkrimp, key, val)
                params_kkrimp.set_value(key, val)
                
        # write config.cfg
        config_filename = tempfolder.get_abs_path(self._CONFIG)
        params_kkrimp.fill_keywords_to_inputfile(output=config_filename)
    

    def _change_atominfo(self, imp_info, kkrflex_file_paths, tempfolder):
        """
        change kkrflex_atominfo to match impurity case
        """
        # read in atominfo file as it is written out
        with open(kkrflex_file_paths.get(KkrCalculation()._KKRFLEX_ATOMINFO)) as file:
            atominfo = file.readlines()
            
        #TODO implement logic to extract this info from imp_info
        replace_zatom_imp = []
        
        # read scoef for comparison with Rimp_rel
        scoef = []
        with open(tempfolder.get_abs_path(KkrCalculation()._SCOEF), 'r') as file:
            Nscoef = int(file.readline().split()[0])
            for iline in range(Nscoef):
                tmpline = file.readline().split()
                scoef.append([float(i) for i in tmpline[:3]])
        scoef = array(scoef)
        
        # find replaceZimp list from Zimp and Rimp_rel
        imp_info_dict = imp_info.get_dict()
        Zimp_list = imp_info_dict.get('Zimp')
        Rimp_rel_list = imp_info_dict.get('Rimp_rel', [[0,0,0]])
        for iatom in range(len(Zimp_list)):
            rtmp = Rimp_rel_list[iatom]
            diff = sqrt(sum((rtmp-scoef)**2, axis=1))
            Zimp = Zimp_list[iatom]
            ipos_replace = where(diff==diff.min())[0][0]
            replace_zatom_imp.append([ipos_replace, Zimp])
        
        for (iatom, zimp) in replace_zatom_imp:
            tmp = atominfo[iatom+4].split()
            x, y, z = float(tmp[0]), float(tmp[1]), float(tmp[2])
            zatom = float(tmp[3])
            virt, remove, lmax = int(tmp[4]), int(tmp[5]), int(tmp[6])
            zatom = zimp
            tmp = ' %24.16f %24.16f %24.16f %5.2f %4i %4i %4i\n'%(x, y, z, zatom, virt, remove, lmax)
            atominfo[iatom+4] = tmp
            
        # write atominfo file
        atominfo_path = tempfolder.get_abs_path(self._KKRFLEX_ATOMINFO)
        with open(atominfo_path, 'w') as file:
            file.writelines(atominfo)
    

    def _get_pot_and_shape(self, imp_info, shapefun_path, shapes, impurity_potential, parent_calc_folder, tempfolder, structure):
        """
        write shapefun from impurity info and host shapefun and copy imp. potential
        """

        scoef_filename = os.path.join(tempfolder.get_abs_path(''), KkrCalculation()._SCOEF)
        imp_info_dict = imp_info.get_dict()
        Rcut = imp_info_dict.get('Rcut', None)
        hcut = imp_info_dict.get('hcut', -1.)
        cylinder_orient = imp_info_dict.get('cylinder_orient', [0., 0., 1.])
        ilayer_center = imp_info_dict.get('ilayer_center', 0)
        make_scoef(structure, Rcut, scoef_filename, hcut, cylinder_orient, ilayer_center)
        
        # create impurity shapefun
        shapefun_new_path = tempfolder.get_abs_path(self._SHAPEFUN)
        modify_potential().shapefun_from_scoef(scoef_filename, shapefun_path, shapes, shapefun_new_path)
            
        # find path to input potential
        if impurity_potential is not None:
            potfile = impurity_potential.get_file_abs_path()
        elif parent_calc_folder is not None:
            self.logger.info('parent_calc_folder {} {}'.format(parent_calc_folder, parent_calc_folder.get_inputs_dict()))
            retrieved = parent_calc_folder.get_inputs(node_type=JobCalculation)[0].get_outputs_dict().get('retrieved', None)
            self.logger.info('potfile {} {}'.format(retrieved, self._OUT_POTENTIAL))
            potfile = retrieved.get_abs_path(self._OUT_POTENTIAL)
        else:
            raise InputValidationError('ERROR in _get_pot_and_shape: neither impurity potential nor parent_calc_folder given!')
            
        # return path of input potential (added to copy_list)
        return potfile
        
    
    def _check_key_setting_consistency(self, params_kkrimp, key, val):
        """
        Check if key/value pair that is supposed to be set is not in conflict with previous settings of parameters in params_kkrimp
        """
        param_ok = True
        
        #TODO implement checks
        
        if not param_ok:
            raise ValueError('Trying to set key "{}" with value "{}" which is in conflict to previous settings!'.format(key, val))
        
