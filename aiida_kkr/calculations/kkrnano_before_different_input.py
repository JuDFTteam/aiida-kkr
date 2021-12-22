# -*- coding: utf-8 -*-
"""
Input plug-in for a KKRnano calculation.
"""
from aiida.engine import CalcJob
from aiida.orm import CalcJobNode, Dict, RemoteData
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.common.exceptions import UniquenessError
from aiida_kkr.tools import find_parent_structure
from masci_tools.io.common_functions import get_Ang2aBohr
    
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
    _DEFAULT_OUTPUT_FILE = 'out'  #'shell output will be shown with outputcat
    _DEFAULT_OUTPUT_PREP_FILE = 'output.0.txt'
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
                         "description": "radius of cluster used for Voronoi (Should be irrelevant for use with aiida)"},
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
        spec.input(
            'parent_calc',
            valid_type=RemoteData,
            required=True,
            help='Use a node that specifies a parent KKRnano or voronoi calculation'
        )

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
        # Check inputdict
        parameters = self.inputs.parameters
        structure = find_parent_structure(self.inputs.parent_calc).get_pymatgen_structure()
        code = self.inputs.code

        # Prepare inputcard from Structure and input parameter data
        with tempfolder.open(self._DEFAULT_INPUT_FILE, u'w') as input_file_handle:
            self._write_input_file(input_file_handle, parameters, structure)
        with tempfolder.open(self._RBASIS, u'w') as rbasis_handle:
            self._write_rbasis(rbasis_handle, structure)
        
        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = self._get_local_copy_list(self.inputs.parent_calc)
        calcinfo.remote_copy_list = []
        calcinfo.retrieve_list = ["output.0.txt","out"] #self._DEFAULT_OUTPUT_PREP_FILE,
                                 #self._DEFAULT_OUTPUT_FILE]#["output.0.txt"]
            # TODO fill retrieve list with binary output file
     

        codeinfo = CodeInfo()
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
        else:
            print('WARNING: Unknown datatype of entry "{}". \
                  Assume array and proceed'.format(value))
            try:
                value=np.array(value)
                content="{} = {}".format(key,self._array2string(value))
                print("{} = {}".format(key,self._array2string(value)))
                return content
            except:
                print('ERROR: Datatype cannot be used as array!')
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
        params=parameters.get_dict()
        for key in params:
            write_list.append(self._getParametersEntry(key,params[key]["value"]))
        
        input_file_handle.writelines(write_list)

        
    
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
        

    def _get_local_copy_list(self, parent_calc_remote):
        """find files to copy from the parent_folder's retrieved to the iniput of a KKRnano calculation"""

        parent_calc = parent_calc_remote.get_incoming(node_class=CalcJobNode).first().node
        retrieved = parent_calc.outputs.retrieved

        local_copy_list = []
        if not self._is_KkrnanoCalc(parent_calc):
            # copy input potential from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._OUT_POTENTIAL_voronoi, self._POTENTIAL)]

            # copy shapefun from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._SHAPEFUN, self._SHAPEFUN)]

        else:  #if parent_calc.process_class == KkrImporterCalculation:
            #TODO add functionality to start from KKRnano parent
            raise NotImplementedError('starting from KKRnano parent not implemented yet.') 

        return local_copy_list


    def _check_valid_parent(self, parent_calc_folder):
        """
        Check that calc is a valid parent for a KKRnano.
        It can be a VoronoiCalculation, KKRnanoCalculation
        """
        overwrite_pot = False

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

        if ((not self._is_KkrnanoCalc(parent_calc))):
            raise ValueError('Parent calculation must be a KkrCalculation')

        return overwrite_pot, parent_calc

    def _is_KkrnanoCalc(self, calc):
        """
        check if calc contains the file out_potential
        """
        is_KKR = False
        if calc.process_type == 'aiida.calculations:kkr.kkrnano':
            retrieved_node = calc.get_retrieved_node()
            if 'out_potential' in retrieved_node.list_object_names():
                is_KKR = True

        return is_KKR
