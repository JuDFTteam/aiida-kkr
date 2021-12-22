# -*- coding: utf-8 -*-
"""
Input plug-in for a KKRnano calculation.
"""
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode, Dict, Bool, RemoteData, StructureData
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.common.exceptions import UniquenessError
from aiida_kkr.tools import find_parent_structure
from masci_tools.io.common_functions import get_Ang2aBohr
from aiida.common.exceptions import InputValidationError

import numpy as np

__copyright__ = (u'Copyright (c), 2021, Forschungszentrum Jülich GmbH, ' 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.1'
__contributors__ = ('Markus Struckmann', 'Philipp Rüßmann')


class KKRnanoCalculation(CalcJob):
    """
    AiiDA calculation plugin for a KKRnano calculation
    """

    ####################
    # File names etc.
    ####################
    # calculation plugin version
    _CALCULATION_PLUGIN_VERSION = __version__
    # Default input and output files
    _DEFAULT_INPUT_FILE = 'input.conf'  # will be shown with inputcat
    _DEFAULT_EFERMI_FILE = 'EFERMI'
    _DEFAULT_NOCO_INPUT_FILE = 'nonco_angle.dat'
    _DEFAULT_OUTPUT_FILE = 'out'  #'shell output will be shown with outputcat
    _DEFAULT_OUTPUT_PREP_FILE = 'output.0.txt'
    _DEFAULT_NOCO_OUTPUT_FILE = 'nonco_angle_out.dat'
    # template.product entry point defined in setup.json
    _DEFAULT_PARSER = 'kkr.kkrnanoparser'
    # File names
    _SHAPEFUN = 'shapefun'
    _POTENTIAL = 'potential'
    _OUT_POTENTIAL = 'out_potential'
    _RBASIS = 'rbasis.xyz'

    
    

    
    _DEFAULT_KKRNANO_PARA = {
      "bzdivide": {"value": [10,10,10],"required": True,"description": "number of k-points in each direction"},
      "emin": {"value": -1.2, "unit":"Rydberg", "required": True, "description":  "lower energy of contour"},
      "emax": {"value": 1.0 , "unit":"Rydberg",  "required": True,\
               "description":  "upper energy of contour (relevant only for DOS calculation)"},
      "npnt1": {"value": 3, "unit":"Rydberg",  "required": True,\
                "description": "number of points starting at emin, parallel to imaginary axis"},
      "npnt2": {"value":  20,  "required": True,\
                "description": "number of points parallel to real axis starting from emin + imag part."},
      "npnt3": {"value":  3,  "required": True,\
                "description": "number of points parallel to real axis in interval (E_F - 30*k*T + imag, E_F + imag)"},
      "npol": {"value":  7,  "required": True,\
               "description": "Number of Matsubara poles, npol=0 triggers DOS calculation"} ,
      "tempr": {"value":  800,  "unit":"Kelvin",\
                "required": False, "description": "artificial temperature (Kelvin) for energy broadening,\
                determines together with npol the distance from real axis"},
      "scfsteps": {"value":  1,"required": True, "description": "Number of scf steps"},
      "imix": {"value":  6,"required": True,\
               "description": "mixing method: imix = 0 -> straight mixing; imix = 1 -> straight mixing;\
               imix = 4 -> Broyden's 2nd method; imix = 5 -> gen. Anderson mixing;\
               imix = 6 -> Broyden's 2nd method with support for >1 atom per process"}, 

      "mixing": {"value":  0.01, "required": True,\
                 "description": "straight mixing parameter"},
      "rmax": {"value":  8.0, "required": True,\
               "description": "Ewald sum cutoff in real space in units of lattice constants"} ,
      "gmax": {"value":  48.0, "required": True,\
               "description": "Ewald sum cutoff in reciprocal space in units of 2*Pi/alat"},
      "nsra": {"value":  2, "required": True, "description": "1=non-scalar-relativistic 2=scalar-relativistic"},
      "kte": {"value":  1, "required": True,\
              "description": "1=calculate energies, -1 = total energy only, less I/O"} ,
      "rclust_voronoi": {"value":  2.00, "required": True,\
                         "description": "radius of cluster used for Voronoi (Should be irrelevant for use with aiida)"}#,
      #"soc":{"value": True}
    }

    
    @classmethod
    def define(cls, spec):
        """
        define internals and inputs / outputs of calculation
        """
        # reuse base class (i.e. CalcJob) functions
        super(KKRnanoCalculation, cls).define(spec)
        # now define input files and parser
        spec.inputs['metadata']['options']['parser_name'].default = cls._DEFAULT_PARSER
        spec.inputs['metadata']['options']['input_filename'].default = cls._DEFAULT_INPUT_FILE
        spec.inputs['metadata']['options']['output_filename'].default = cls._DEFAULT_OUTPUT_FILE
        #spec.inputs['metadata']['options']['output_prep_filename'].default = cls._DEFAULT_OUTPUT_PREP_FILE ->does not work
        
        # define input nodes (optional ones have required=False)
        spec.input('parameters', valid_type=Dict, required=False,
                   default= lambda: Dict(dict=cls._DEFAULT_KKRNANO_PARA),
                   help='Dict node that specifies the input parameters for KKRnano (k-point density etc.)')
        spec.input('nocoangles', valid_type=Dict, required=False, default= lambda: Dict(dict={}),
                   help='Dict node that specifies the starting angles for non-colinear calculations\
                   (only needed in conjunction with non-colinear calculations, i. e. KORBIT=1\
                   (which is also necessary for SOC calculations!))' )
        spec.input(
            'parent_folder',
            valid_type=RemoteData,
            required=True,
            help='Use a node that specifies a parent KKRnano or voronoi calculation'
        )
        
        spec.input('convert', valid_type=Bool, required=False, default=lambda: Bool(False),
                  help='Activate to use together with set up convert code in order to retrieve potential files.')
        
        #spec.input('structure', valid_type=StructureData, required=False, default:)
        # define outputs
        spec.output('output_parameters', valid_type=Dict, required=True, help='results of the calculation')
        spec.default_output_node = 'output_parameters'

        # define exit codes, also used in parser
        spec.exit_code(301, 'ERROR_NO_OUTPUT_FILE', message='KKRnano output file not found')
        spec.exit_code(302, 'ERROR_PARSING_FAILED', message='KKRnano parser retuned an error')


    def prepare_for_submission(self, tempfolder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param tempfolder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        #prepare inputs
        parameters = self.inputs.parameters.get_dict()
        nonco_angles=self.inputs.nocoangles.get_dict()
        parent_calc=self.inputs.parent_folder#.get_incoming(node_class=CalcJobNode).first().node
        parent_calc_calc_node=parent_calc.get_incoming(node_class=CalcJobNode).first().node
        structure = find_parent_structure(parent_calc_calc_node).get_pymatgen_structure()
        parent_outfolder = parent_calc_calc_node.outputs.retrieved
        code = self.inputs.code

        # Local copy list for continuing KKRnano calculations
        # if necessary, change some names in order to start from a preconverged calculation
        # (cf. jukkr/source/KKRnano/scripts/prepare.sh)
        local_copy_list_for_continued=[(parent_outfolder.uuid,"bin.vpotnew", "bin.vpotnew.0"),
                           (parent_outfolder.uuid,"bin.vpotnew.idx", "bin.vpotnew.0.idx"),
                           (parent_outfolder.uuid,"bin.meshes","bin.meshes.0"  ),
                           (parent_outfolder.uuid,"bin.meshes.idx","bin.meshes.0.idx" ),
                           (parent_outfolder.uuid,"bin.energy_mesh","bin.energy_mesh.0" ),
                           (parent_outfolder.uuid,"bin.atoms","bin.atoms")]
  
        
                           #(parent_outfolder.uuid,"nonco_angle_out.dat","bin.atoms"  )]
        #TODO: Parse noco_angle_out.dat and write a new file from that
        
        
        #Check if convert mode has been activated
        convert=self.inputs.convert.value
        if convert:
            if not parent_calc_calc_node.process_label=='KKRnanoCalculation':
                raise InputValidationError("A convert process can only be started from a KKRnano calculation. Check also that convert code is used!")
            
            #mark files from previous calc as such
            local_copy_list_for_continued=[(parent_outfolder.uuid,"bin.vpotnew", "bin.vpotnew"),
                   (parent_outfolder.uuid,"bin.vpotnew.idx", "bin.vpotnew.idx"),
                   (parent_outfolder.uuid,"bin.meshes","bin.meshes"  ),
                   (parent_outfolder.uuid,"bin.meshes.idx","bin.meshes.idx" ),
                   (parent_outfolder.uuid,"bin.energy_mesh","bin.energy_mesh" ),
                   (parent_outfolder.uuid,"bin.atoms","bin.atoms"),
                   (parent_outfolder.uuid,self._DEFAULT_OUTPUT_PREP_FILE,self._DEFAULT_OUTPUT_PREP_FILE),
                   (parent_outfolder.uuid,self._DEFAULT_OUTPUT_FILE,self._DEFAULT_OUTPUT_FILE)]

        
        # Check inputdict
        self._check_input_dict(parameters)
        #Check if parent is a KKRnano calculation
        write_efermi=False
        self._check_valid_parent(parent_outfolder)
        
        
        if parent_calc_calc_node.process_label=='KKRnanoCalculation':
            fermi=parent_calc_calc_node.outputs.output_parameters.get_dict()['fermi_energy_in_ryd'][-1]
            write_efermi=True

        # Check if non-colinear calculation mode is activated
        noco=False
        if "KORBIT" in parameters:
            if parameters["KORBIT"]['value']==1:
                noco=True
        if "soc" in parameters:
            print("SOC in parameters")
            print("NOCO0",noco)
            if parameters["soc"]["value"]==True:
                print("NOCO1",noco)
                if noco==False:
                    print("NOCO2",noco)
                    parameters["KORBIT"]['value']=True
                    self.logger.warn('KORBIT was set to 1 -> SOC is implemented for the NOCO-Chebyshev solver, only! This is however not changed in the input node and might lead to inconsistencies if this feature has been added in KKRnano!')
                    noco=True
        
        
        print("NOCO",noco)
        if noco:
            with tempfolder.open(self._DEFAULT_NOCO_INPUT_FILE, u'w') as nonco_angles_handle:
                self._write_nonco_angles(nonco_angles_handle,nonco_angles,structure)
        
        
        # Prepare rbasis.xyz and input.conf from Structure and input parameter data
        with tempfolder.open(self._DEFAULT_INPUT_FILE, u'w') as input_file_handle:
            self._write_input_file(input_file_handle, parameters, structure)
        with tempfolder.open(self._RBASIS, u'w') as rbasis_handle:
            self._write_rbasis(rbasis_handle, structure)
        if write_efermi:
            with tempfolder.open(self._DEFAULT_EFERMI_FILE, u'w') as efermi_file_handle:
                self._write_efermi_file(efermi_file_handle, fermi)
            
        
        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list= self._get_local_copy_list(parent_calc,local_copy_list_for_continued)
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = [self._DEFAULT_OUTPUT_PREP_FILE,
                                 self._DEFAULT_OUTPUT_FILE,
                                 self._DEFAULT_NOCO_OUTPUT_FILE,
                                 "vpot*","DOS*"]
        
    #TODO: Add fancy functionality to start from potential files "built" from vpot-files 
        
        #add the binary files that are necessary to restart a calculation to the retrieved list
        for j in range(len(local_copy_list_for_continued)):
               calcinfo.retrieve_list.append(local_copy_list_for_continued[j][1])
     

        codeinfo = CodeInfo()
        if convert:
            codeinfo.cmdline_params = ['--convert']
        else:
            codeinfo.cmdline_params = []
        codeinfo.stdout_name = self._DEFAULT_OUTPUT_FILE
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo
    
    
    def _list2string(self,list0):
        string=''
        for item in list0:
            string+=str(item) + ' '
        return string

    def _array2string(self,array):
        if array.shape==():
            return str(array)
        else:
            return self._list2string(array) 
    
    def _getParametersEntry(self,key, value):
        """writing out various entry types from a dict into the format suitable for KKRnano"""
        
        if type(value)==str or type(value)==int:
            content="{} = {}".format(key, value)
        elif type(value)==float or type(value)==np.float_:
            content="{} = {}".format(key, str.replace(str(value),"e", "D"))
        elif type(value)==np.ndarray:
            content="{} = {}".format(key,self._array2string(value))
        elif type(value)==list:
             content="{} = {}".format(key,self._list2string(value))
        elif type(value) is bool:
            value2write = "f"
            if value:
                value2write="t"
            self.logger.info("TEST123")
            content="{} = {}".format(key,value2write)
        else:
            print('WARNING: Unknown datatype of entry "{}". \
                  Assume array and proceed'.format(value))
            try:
                value=np.array(value)
                content="{} = {}".format(key,self._array2string(value))
                print("{} = {}".format(key,self._array2string(value)))
                return content
            except:
                self.logger.error('ERROR: Datatype cannot be used as array!')
        content=content+"\n"
        return content

    def _write_input_file(self, input_file_handle, parameters, structure):
        """write the input.conf file for KKRnano"""
        
        lattice_param_angs=max(structure.lattice.abc)
        write_list=[]
        
        #header of input.conf
        write_list.append("# LATTICE \n \n# lattice parameter in units of the Bohr radius \n")


        
        write_list.append('alat = {}'.format(lattice_param_angs*get_Ang2aBohr()))

        write_list.append("\n# scale basis coordinates by these factors \nbasisscale = 1.0  1.0  1.0 \n \n")
        latt_dict=structure.lattice.as_dict()

        write_list.append('#BRAVAIS\n')
        rbrav=np.array(latt_dict['matrix'])/lattice_param_angs

        directions=['a','b','c']
        for i in range(3):
            write_list.append('bravais_{}=  {}  {}  {}\n'.format(directions[i],rbrav[i][0],rbrav[i][1],rbrav[i][2]))
        write_list.append("\ncartesian = t\n")
        #writing other parameters from parameter dict
        #params=parameters.get_dict()
        params=parameters
        
        for key in params:
            write_list.append(self._getParametersEntry(key,params[key]["value"]))
        
        write_list.append("# CHANGED SCRIPT")
        
        
        input_file_handle.writelines(write_list)

        
    def _write_nonco_angles(self,nonco_angles_handle,nonco_angles,structure):
        """
        write nonco_angles.dat file for KKRnano
        created from dictionary with structure dict={'atom':{1:{'theta':0,'phi':0, 'fix_angle_mode':1},...}
        The angles 'theta' (polar angle going from z- to x-direction) and 'phi' (azimuthal angle) are given in deg.
        'fix_angle_mode' takes values 0,1,2,3. 0 is for relaxation of the spin-direciton, 1 is for fixing it; 2 and 3
        are for constraining fields calculations.
        """
        n_atoms=len(structure.atomic_numbers)
        
        #TODO: Adapt to KKRhost style
        #Create dictionary if none was given in the input
        if nonco_angles=={}:
            nonco_angles={'atom':{}}
            for i in range(n_atoms):
                nonco_angles['atom'][i+1]={'theta': 0.0, 'phi': 0.0, 'fix_angle_mode': 1}
        print(n_atoms)
        print(nonco_angles)
        print(len(nonco_angles['atom'].keys()))
        if n_atoms != len(nonco_angles['atom'].keys()):
            raise InputValidationError("The number of atoms in the structure must\
            match the number of atoms specified for the non-colinear angles.")
        
        write_list=[]
        for key in nonco_angles['atom']:
            write_list.append("{} {} {}\n".format(nonco_angles['atom'][key]['theta'],\
                                                  nonco_angles['atom'][key]['phi'],\
                                                  nonco_angles['atom'][key]['fix_angle_mode']))
        print(write_list)                                         
        nonco_angles_handle.writelines(write_list)
        
                              
                              
    def _write_rbasis(self, rbasis_handle, structure):
        """write the rbasis.xyz file for KKRnano"""
        #PSE element symbols 
        symbols = ('__',
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
        'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
        'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',
        'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
        'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
        'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo'
    )
                
        
        
        write_list=[str(len(structure.atomic_numbers)),'\n',\
                    str(structure.composition), ", ", str(structure.get_space_group_info()),'\n']
        
        lattice_param_angs=max(structure.lattice.abc)
        
        for i in range(np.shape(structure.cart_coords)[0]):
               write_list.append('{}   {}   {}   {}\n'.format(symbols[structure.atomic_numbers[i]], \
                  structure.cart_coords[i,0]/lattice_param_angs, \
                  structure.cart_coords[i,1]/lattice_param_angs,\
                  structure.cart_coords[i,2]/lattice_param_angs))
        
        rbasis_handle.writelines(write_list)
        
    def _write_efermi_file(self,efermi_handle, efermi):
        """ write file EFERMI necessary to restart a calculation. """
        efermi_handle.writelines([str(efermi)])
                                 
    def _get_local_copy_list(self, parent_calc_remote,local_copy_list_for_continued):
        """find files to copy from the parent_folder's retrieved to the iniput of a KKRnano calculation"""

        parent_calc = parent_calc_remote.get_incoming(node_class=CalcJobNode).first().node
        retrieved = parent_calc.outputs.retrieved

        local_copy_list = []
        if not self._is_KkrnanoCalc(parent_calc):
            # copy input potential from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._OUT_POTENTIAL_voronoi, self._POTENTIAL)]

            # copy shapefun from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._SHAPEFUN, self._SHAPEFUN)]

        elif parent_calc.process_label=='KKRnanoCalculation':
            local_copy_list=local_copy_list_for_continued
        return local_copy_list


    def _check_valid_parent(self, parent_calc_folder):
        """
        Check that calc is a valid parent for a KKRnano.
        It can be a VoronoiCalculation, KKRnanoCalculation
        """

        # extract parent calculation
        parent_calcs = parent_calc_folder.get_incoming(node_class=CalcJobNode)
        n_parents = len(parent_calcs.all_link_labels())
        if n_parents != 1:
            raise UniquenessError(
                'Input RemoteData is child of {} '
                'calculation{}, while it should have a single parent'
                ''.format(n_parents, '' if n_parents == 0 else 's')
            )
        else:
            parent_calc = parent_calcs.first().node
            overwrite_pot = True

#         if ((not self._is_KkrnanoCalc(parent_calc))):
#             raise ValueError('Parent calculation must be a KkrCalculation.')
        if (not parent_calc.process_label=='KKRnanoCalculation')\
    and (not parent_calc.process_label=='VoronoiCalculation'):
            raise ValueError("Parent Calculation has to be another KKRnano or a Voronoi calculation.")
    
    def _check_input_dict(self,inputdict):
        """ checks if all essential keys are contained in the inputdict and if it has the right format """
        #Check if all essential keys are present
        all_present=True
        missing_keys=[]
        for key in self._DEFAULT_KKRNANO_PARA:
            if not key in inputdict:
                all_present=False
                missing_keys.append(key)
        if not all_present:
            raise InputValidationError('Not all essential keys were given in the input dictionary. \
            At least the following keys are missing: \
            {}'.format(missing_keys))
        #Check for format
        keys_in_wrong_format=[]
        correct_format=True
        for key in inputdict:
            try:
                if not "value" in inputdict[key]:
                    correct_format=False
                    keys_in_wrong_format.append(key)
            except:
                    correct_format=False
                    keys_in_wrong_format.append(key)
        if not correct_format:
            raise InputValidationError('Some keys were provided in an incorrect format. Dict entries must have entry `value`\
(see e. g. KKRnano-Plugin/aiida-kkr/aiida_kkr/calculations/kkrnano.py). At least the following keys are in an incorrect format: \
{}'.format(keys_in_wrong_format))    
                
    def _is_KkrnanoCalc(self, calc):
        """
        check if calc contains the file out_potential
        """
        is_KKR = False
        if calc.process_type == 'aiida.calculations:kkr.kkrnano':
            retrieved_node = calc.get_retrieved_node()
            if True:#'out_potential' in retrieved_node.list_object_names():
                is_KKR = True

        return is_KKR