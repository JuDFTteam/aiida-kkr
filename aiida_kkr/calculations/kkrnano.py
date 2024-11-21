# -*- coding: utf-8 -*-
"""
Input plug-in for a KKRnano calculation.
"""

import numpy as np
from masci_tools.io.common_functions import get_Ang2aBohr
from aiida.orm import CalcJobNode, Dict, Bool, Float, RemoteData, StructureData
from aiida.engine import CalcJob
from aiida.common import NotExistent
from aiida.common.datastructures import (CalcInfo, CodeInfo)
from aiida.common.exceptions import InputValidationError, UniquenessError
from aiida_kkr.tools.find_parent import get_remote, get_parent
from aiida_kkr.data.strucwithpot import StrucWithPotData

__copyright__ = (u'Copyright (c), 2021, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.2'
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
    _CONVERT_OUTPUT_FILE = 'convert_out'
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
      'bzdivide': {'value': [10,10,10],'required': True,'description': 'number of k-points in each direction'},
      'emin': {'value': -1.2, 'unit':'Rydberg', 'required': True, 'description':  'lower energy of contour'},
      'emax': {'value': 1.0 , 'unit':'Rydberg',  'required': True,\
               'description':  'upper energy of contour (relevant only for DOS calculation)'},
      'npnt1': {'value': 3, 'unit':'Rydberg',  'required': True,\
                'description': 'number of points starting at emin, parallel to imaginary axis'},
      'npnt2': {'value':  20,  'required': True,\
                'description': 'number of points parallel to real axis starting from emin + imag part.'},
      'npnt3': {'value':  3,  'required': True,\
                'description': 'number of points parallel to real axis in interval (E_F - 30*k*T + imag, E_F + imag)'},
      'npol': {'value':  7,  'required': True,\
               'description': 'Number of Matsubara poles, npol=0 triggers DOS calculation'} ,

      'scfsteps': {'value':  1,'required': True, 'description': 'Number of scf steps'},
      'imix': {'value':  6,'required': True,\
               'description': "mixing method: imix = 0 -> straight mixing; imix = 1 -> straight mixing;\
               imix = 4 -> Broyden's 2nd method; imix = 5 -> gen. Anderson mixing;\
               imix = 6 -> Broyden's 2nd method with support for >1 atom per process"
              },

      'mixing': {'value':  0.01, 'required': True,\
                 'description': 'straight mixing parameter'},
      'rmax': {'value':  8.0, 'required': True,\
               'description': 'Ewald sum cutoff in real space in units of lattice constants'} ,
      'gmax': {'value':  48.0, 'required': True,\
               'description': 'Ewald sum cutoff in reciprocal space in units of 2*Pi/alat'},
      'nsra': {'value':  2, 'required': True, 'description': '1=non-scalar-relativistic 2=scalar-relativistic'},
      'kte': {'value':  1, 'required': True,\
              'description': '1=calculate energies, -1 = total energy only, less I/O'} ,
      'rclust_voronoi': {'value':  2.00, 'required': True,\
                         'description': 'radius of cluster used for Voronoi (Should be irrelevant for use with aiida)'}#,
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
        spec.input(
            'parameters',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._DEFAULT_KKRNANO_PARA),
            help='Dict node that specifies the input parameters for KKRnano (k-point density etc.)'
        )
        spec.input(
            'nocoangles',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict={}),
            help='Dict node that specifies the starting angles for non-colinear calculations\
                   (only needed in conjunction with non-colinear calculations, i. e. KORBIT=1\
                   (which is also necessary for SOC calculations!))'
        )
        spec.input(
            'parent_folder',
            valid_type=RemoteData,
            required=False,
            help='Use a node that specifies a parent KKRnano or voronoi calculation'
        )

        spec.input(
            'convert',
            valid_type=Bool,
            required=False,
            default=lambda: Bool(False),
            help='Activate to use together with set up convert code in order to retrieve potential files.'
        )

        spec.input(
            'passed_lattice_param_angs',
            valid_type=Float,
            required=False,
            default=lambda: Float(-10000.0),
            help='Use a prespecified lattice constant in Angstrom as input for KKRnano, i. e. in the input.conf file. \
                   Default is the length of the longest Bravais vector in the structure object used for the voronoi calculation. \
                   This can be useful in the context of treating supercells.'
        )

        spec.input('strucwithpot', valid_type=StrucWithPotData, required=False)

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
        #passed_lattice_const = self.inputs.passed_lattice_param_angs.value
        nonco_angles = self.inputs.nocoangles.get_dict()
        #parent_calc=self.inputs.parent_folder#.get_incoming(node_class=CalcJobNode).first().node
        #parent_calc_calc_node=parent_calc.get_incoming(node_class=CalcJobNode).first().node
        #structure = self._find_parent_struc_from_voro_or_stwpd(parent_calc_calc_node)[0].get_pymatgen_structure()
        #parent_outfolder = parent_calc_calc_node.outputs.retrieved
        code = self.inputs.code

        #determine number of MPI processes
        try:
            num_mpi_procs = self.metadata.options.resources['tot_num_mpiprocs']
        except:
            try:
                num_mpi_procs=self.metadata.options.resources['num_machines']*\
                self.metadata.options.resources['num_mpiprocs_per_machine']
            except:
                raise InputValidationError(
                    "The total number of MPI processes could not be determined. In case of doubt: Specify `builder.metadata.options.resources['tot_num_mpiprocs']=?` Do not forget number of machines"
                )

        #Check if convert mode has been activated
        convert = self.inputs.convert.value

        #StrucWithPot object as starting point -> contains passed_lattice_constant, shapefun, potential, and structure
        if hasattr(self.inputs, 'strucwithpot') and not hasattr(self.inputs, 'parent_folder'):
            use_strucwithpot = True
            strucwithpot = self.inputs.strucwithpot
            structure = strucwithpot.structure.get_pymatgen_structure()

            if self.inputs.passed_lattice_param_angs.value == -10000.0:  #check if default value is used
                if hasattr(strucwithpot, 'specified_lattice_constant'):
                    passed_lattice_const = strucwithpot.specified_lattice_constant.value
                else:
                    passed_lattice_const = self.inputs.passed_lattice_param_angs.value
            else:
                passed_lattice_const = self.inputs.passed_lattice_param_angs.value

            parent_outfolder_uuid = None  #use dummy so that copylists are still available
        elif not hasattr(self.inputs, 'strucwithpot') and hasattr(self.inputs, 'parent_folder'):
            use_strucwithpot = False
            passed_lattice_const = self.inputs.passed_lattice_param_angs.value
            parent_calc = self.inputs.parent_folder  #.get_incoming(node_class=CalcJobNode).first().node
            parent_calc_calc_node = parent_calc.get_incoming(node_class=CalcJobNode).first().node

            #Find structure unless convert mode is activated, as for that no structure is needed
            if not convert:
                structure = self.find_parent_struc_from_voro_or_stwpd(parent_calc_calc_node)[0].get_pymatgen_structure()
            parent_outfolder = parent_calc_calc_node.outputs.retrieved

            parent_outfolder_uuid = parent_outfolder.uuid

        else:
            raise InputValidationError(
                'Either `strucwithpot` or a `parent_folder` has to be provided. If necessary remove one of the inputs.'
            )

        print('passed lattice constant=', passed_lattice_const)

        # Local copy list for continuing KKRnano calculations (Note: parent_outfolder_uuid is None if use_strucwithpot)
        # if necessary, change some names in order to start from a preconverged calculation
        # (cf. jukkr/source/KKRnano/scripts/prepare.sh)
        local_copy_list_for_continued = [(parent_outfolder_uuid, 'bin.vpotnew', 'bin.vpotnew.0'),
                                         (parent_outfolder_uuid, 'bin.vpotnew.idx', 'bin.vpotnew.0.idx'),
                                         (parent_outfolder_uuid, 'bin.meshes', 'bin.meshes.0'),
                                         (parent_outfolder_uuid, 'bin.meshes.idx', 'bin.meshes.0.idx'),
                                         (parent_outfolder_uuid, 'bin.energy_mesh', 'bin.energy_mesh.0'),
                                         (parent_outfolder_uuid, 'bin.atoms', 'bin.atoms')]

        #(parent_outfolder_uuid,"nonco_angle_out.dat","bin.atoms"  )]
        #TODO: Parse noco_angle_out.dat and write a new file from that

        #Convert mode is not available for call using strucwithpot
        if use_strucwithpot and convert:
            raise InputValidationError('Convert mode cannot be used for call using a strucwithpot object.')
        if convert:
            if not parent_calc_calc_node.process_label == 'KKRnanoCalculation':
                raise InputValidationError(
                    'A convert process can only be started from a KKRnano calculation. Check also that convert code is used!'
                )

            #mark files from previous calc as such (Note: parent_outfolder_uuid is None if use_strucwithpot)
            local_copy_list_for_continued = [
                (parent_outfolder_uuid, 'bin.vpotnew', 'bin.vpotnew'),
                (parent_outfolder_uuid, 'bin.vpotnew.idx', 'bin.vpotnew.idx'),
                (parent_outfolder_uuid, 'bin.meshes', 'bin.meshes'),
                (parent_outfolder_uuid, 'bin.meshes.idx', 'bin.meshes.idx'),
                (parent_outfolder_uuid, 'bin.energy_mesh', 'bin.energy_mesh'),
                (parent_outfolder_uuid, 'bin.atoms', 'bin.atoms'), (parent_outfolder_uuid, 'bin.dims', 'bin.dims'),
                (parent_outfolder_uuid, self._DEFAULT_OUTPUT_PREP_FILE, self._DEFAULT_OUTPUT_PREP_FILE),
                (parent_outfolder_uuid, self._DEFAULT_OUTPUT_FILE, self._DEFAULT_OUTPUT_FILE)
            ]

        # Check inputdict
        parameters = self._check_input_dict(parameters, num_mpi_procs, convert)

        #Check if parent is a KKRnano calculation
        write_efermi = False
        if not use_strucwithpot:
            self._check_valid_parent(parent_outfolder)
            if parent_calc_calc_node.process_label == 'KKRnanoCalculation':
                fermi = parent_calc_calc_node.outputs.output_parameters.get_dict()['fermi_energy_in_ryd'][-1]
                write_efermi = True

        # Check if non-colinear calculation mode is activated
        noco = False
        if 'KORBIT' in parameters:
            if parameters['KORBIT']['value'] == 1:
                noco = True
        if 'soc' in parameters:
            if parameters['soc']['value'] == True:
                if noco == False:
                    parameters['KORBIT']['value'] = 1
                    self.logger.warn(
                        'KORBIT was set to 1 -> SOC is implemented for the NOCO-Chebyshev solver, only! This is however not changed in the input node and might lead to inconsistencies if this feature has been added in KKRnano!'
                    )
                    noco = True

        if noco:
            with tempfolder.open(self._DEFAULT_NOCO_INPUT_FILE, u'w') as nonco_angles_handle:
                self._write_nonco_angles(nonco_angles_handle, nonco_angles, structure)  # pylint: disable=possibly-used-before-assignment

        # Prepare rbasis.xyz and input.conf from Structure and input parameter data unless convert mode
        if not convert:
            with tempfolder.open(self._DEFAULT_INPUT_FILE, u'w') as input_file_handle:
                self._write_input_file(input_file_handle, parameters, structure, passed_lattice_const)
            with tempfolder.open(self._RBASIS, u'w') as rbasis_handle:
                self._write_rbasis(rbasis_handle, structure, passed_lattice_const)
        if write_efermi:
            with tempfolder.open(self._DEFAULT_EFERMI_FILE, u'w') as efermi_file_handle:
                self._write_efermi_file(efermi_file_handle, fermi)  # pylint: disable=used-before-assignment

        # Prepare potential and shapefun file from strucwithpot, if necessary
        if use_strucwithpot:
            with tempfolder.open(self._POTENTIAL, u'w') as potential_file_handle:
                self._write_potential_file(potential_file_handle, strucwithpot)
            with tempfolder.open(self._SHAPEFUN, u'w') as shapefun_file_handle:
                self._write_shapefun_file(shapefun_file_handle, strucwithpot)

        # Prepare CalcInfo to be returned to aiida
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid

        if not use_strucwithpot:
            calcinfo.local_copy_list = self._get_local_copy_list(parent_calc, local_copy_list_for_continued)
        if use_strucwithpot:
            calcinfo.local_copy_list = []
        calcinfo.remote_copy_list = []

        #retrieve list
        #TODO add NOCO
        calcinfo.retrieve_list = [
            self._DEFAULT_OUTPUT_PREP_FILE, self._DEFAULT_OUTPUT_FILE, self._DEFAULT_NOCO_OUTPUT_FILE, 'bin.dims',
            self._CONVERT_OUTPUT_FILE, 'vpot*', 'DOS*', 'nonco_angle_out.dat'
        ]  # bin.dims is added here, as this should be retrieved for generating vpot-files (negligable in size anyway)

        for j in range(
            len(local_copy_list_for_continued)
        ):  #add the binary files that are necessary to restart a calculation to the retrieved list
            calcinfo.retrieve_list.append(local_copy_list_for_continued[j][1])

        codeinfo = CodeInfo()
        if convert:
            codeinfo.cmdline_params = ['--convert']
            codeinfo.stdout_name = self._CONVERT_OUTPUT_FILE
        else:
            codeinfo.cmdline_params = []
            codeinfo.stdout_name = self._DEFAULT_OUTPUT_FILE
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo

    def _list2string(self, list0):
        string = ''
        for item in list0:
            string += str(item) + ' '
        return string

    def _array2string(self, array):
        if array.shape == ():
            return str(array)
        else:
            return self._list2string(array)

    def _getParametersEntry(self, key, value):
        """writing out various entry types from a dict into the format suitable for KKRnano"""

        if type(value) == str or type(value) == int:
            content = f'{key} = {value}'
        elif type(value) == float or type(value) == np.float_:
            content = f"{key} = {str.replace(str(value), 'e', 'D')}"
        elif type(value) == np.ndarray:
            content = f'{key} = {self._array2string(value)}'
        elif type(value) == list:
            content = f'{key} = {self._list2string(value)}'
        elif type(value) is bool:
            value2write = 'f'
            if value:
                value2write = 't'
            self.logger.info('TEST123')
            content = f'{key} = {value2write}'
        else:
            print(f'WARNING: Unknown datatype of entry "{value}".                   Assume array and proceed')
            try:
                value = np.array(value)
                content = f'{key} = {self._array2string(value)}'
                print(f'{key} = {self._array2string(value)}')
                return content
            except:
                self.logger.error('ERROR: Datatype cannot be used as array!')
        content = content + '\n'
        return content

    def _get_lattice_constant(self, structure, passed_lattice_const):
        """determine whether a passed lattice constant should be used"""
        # default value, see above. (Probably) No one would use this large a lattice constant to replace the other
        if passed_lattice_const == -10000.0:
            lattice_param_angs = max(structure.lattice.abc)
        elif passed_lattice_const > 0.0:
            lattice_param_angs = passed_lattice_const
            self.logger.info('Using passed lattice constant. The vectors are scaled accordingly!')
        else:
            raise InputValidationError('The passed lattice constant must not be negative!')
        return lattice_param_angs

    def _write_input_file(self, input_file_handle, parameters, structure, passed_lattice_const):
        """write the input.conf file for KKRnano"""

        #lattice_param_angs=max(structure.lattice.abc)

        lattice_param_angs = self._get_lattice_constant(structure, passed_lattice_const)

        write_list = []

        #header of input.conf
        write_list.append('# LATTICE \n \n# lattice parameter in units of the Bohr radius \n')

        write_list.append(f'alat = {lattice_param_angs * get_Ang2aBohr()}')

        write_list.append('\n# scale basis coordinates by these factors \nbasisscale = 1.0  1.0  1.0 \n \n')
        latt_dict = structure.lattice.as_dict()

        write_list.append('#BRAVAIS\n')
        rbrav = np.array(latt_dict['matrix']) / lattice_param_angs

        directions = ['a', 'b', 'c']
        for i in range(3):
            write_list.append(f'bravais_{directions[i]}=  {rbrav[i][0]:.15f}  {rbrav[i][1]:.15f}  {rbrav[i][2]:.15f}\n')
        write_list.append('\ncartesian = t\n')
        #writing other parameters from parameter dict
        #params=parameters.get_dict()
        params = parameters

        for key in params:
            write_list.append(self._getParametersEntry(key, params[key]['value']))

        write_list.append('# CHANGED SCRIPT')

        input_file_handle.writelines(write_list)

    def _write_nonco_angles(self, nonco_angles_handle, nonco_angles, structure):
        """
        write nonco_angles.dat file for KKRnano
        created from dictionary with structure dict={'atom':{1:{'theta':0.0,'phi':00.0, 'fix_angle_mode':1},...}
        The angles 'theta' (polar angle going from z- to x-direction) and 'phi' (azimuthal angle) are given in deg.
        'fix_angle_mode' takes values 0,1,2,3. 0 is for relaxation of the spin-direciton, 1 is for fixing it; 2 and 3
        are for constraining fields calculations.
        """
        n_atoms = len(structure.atomic_numbers)

        #TODO: Adapt to KKRhost style
        #Create dictionary if none was given in the input
        if nonco_angles == {}:
            nonco_angles = {'atom': {}}
            for i in range(n_atoms):
                nonco_angles['atom'][i + 1] = {'theta': 0.0, 'phi': 0.0, 'fix_angle_mode': 1}
        if n_atoms != len(nonco_angles['atom'].keys()):
            raise InputValidationError(
                'The number of atoms in the structure must\
            match the number of atoms specified for the non-colinear angles.'
            )

        write_list = []
        for key in nonco_angles['atom']:
            write_list.append('{} {} {}\n'.format(nonco_angles['atom'][key]['theta'],\
                                                  nonco_angles['atom'][key]['phi'],\
                                                  nonco_angles['atom'][key]['fix_angle_mode']))

        nonco_angles_handle.writelines(write_list)

    def _get_rbasis_atom_symbol(self, atomic_number):
        """returns either an element symbol or a vacuum symbol"""
        #PSE element symbols
        symbols = (
            '__', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
            'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
            'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo'
        )
        if abs(atomic_number) > len(
            symbols
        ) + 30:  #pymatgen stores X sites with very high random `atomic numbers` that can also be negative
            return symbols[0]
        else:
            return symbols[atomic_number]

    def _write_rbasis(self, rbasis_handle, structure, passed_lattice_const):
        """write the rbasis.xyz file for KKRnano"""

        try:
            space_group_info = str(structure.get_space_group_info())
        except:
            space_group_info = 'No space group could be identifed'

        write_list=[str(len(structure.atomic_numbers)),'\n',\
                    str(structure.composition), ', ', space_group_info,'\n']

        lattice_param_angs = self._get_lattice_constant(structure, passed_lattice_const)

        for i in range(np.shape(structure.cart_coords)[0]):
            write_list.append('{}   {:.15f}   {:.15f}   {:.15f}\n'.format(self._get_rbasis_atom_symbol(structure.atomic_numbers[i]), \
               structure.cart_coords[i,0]/lattice_param_angs, \
               structure.cart_coords[i,1]/lattice_param_angs,\
               structure.cart_coords[i,2]/lattice_param_angs))

        rbasis_handle.writelines(write_list)

    def _write_efermi_file(self, efermi_handle, efermi):
        """ write file EFERMI necessary to restart a calculation. """
        efermi_handle.writelines([str(efermi)])

    def _write_potential_file(self, potential_handle, strucwithpot):
        """ write file potential from strucwithpot object's potential files. """
        for vpot in strucwithpot.potentials:
            potential_handle.writelines(vpot.get_content())  #.split('\n'))

    def _write_shapefun_file(self, shapefun_handle, strucwithpot):
        """ write file potential from strucwithpot object's potential files. """

        shapefun_handle.writelines(strucwithpot.make_shapefun())  #.split('\n'))

    def _get_local_copy_list(self, parent_calc_remote, local_copy_list_for_continued):
        """find files to copy from the parent_folder's retrieved to the input of a KKRnano calculation"""
        #Note: the check if a non-empty copylist is needed, i. e. if started
        # from a parent calculation, is done in prepare_for_submission
        parent_calc = parent_calc_remote.get_incoming(node_class=CalcJobNode).first().node
        retrieved = parent_calc.outputs.retrieved

        local_copy_list = []
        if not self._is_KkrnanoCalc(parent_calc):
            # copy input potential from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._OUT_POTENTIAL_voronoi, self._POTENTIAL)]

            # copy shapefun from voronoi output
            local_copy_list += [(retrieved.uuid, parent_calc.process_class._SHAPEFUN, self._SHAPEFUN)]

        elif parent_calc.process_label == 'KKRnanoCalculation':
            local_copy_list = local_copy_list_for_continued
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


        if (not parent_calc.process_label=='KKRnanoCalculation')\
    and (not parent_calc.process_label=='VoronoiCalculation'):
            raise ValueError('Parent Calculation has to be another KKRnano or a Voronoi calculation.')

    def _check_input_dict(self, inputdict, num_mpi_procs, convert):
        """ checks if all essential keys are contained in the inputdict and if it has the right format """
        #Check if all essential keys are present
        all_present = True
        missing_keys = []
        for key in self._DEFAULT_KKRNANO_PARA:
            if not key in inputdict:
                all_present = False
                missing_keys.append(key)
        if not all_present:
            raise InputValidationError(
                'Not all essential keys were given in the input dictionary. \
            At least the following keys are missing: \
            {}'.format(missing_keys)
            )
        #Check for format
        keys_in_wrong_format = []
        correct_format = True
        for key in inputdict:
            try:
                if not 'value' in inputdict[key]:
                    correct_format = False
                    keys_in_wrong_format.append(key)
            except:
                correct_format = False
                keys_in_wrong_format.append(key)
        if not correct_format:
            self.logger.warn(
                'Some keys were provided in an incorrect format. Dict entries must have entry `value`\
(see e. g. KKRnano-Plugin/aiida-kkr/aiida_kkr/calculations/kkrnano.py). At least the following keys are in an incorrect format: \
{} \nMoving on by assuming just the key `value` was missing and adapt accordingly'.format(keys_in_wrong_format)
            )
            #try to correct incorrect entries, most likely they do not have the 'value' feature
            try:
                for key in keys_in_wrong_format:
                    likely_value = inputdict[key]
                    inputdict[key] = {'value': likely_value}
            except:
                raise InputValidationError(
                    'Some keys were provided in an incorrect format. Dict entries must have entry `value`\
(see e. g. KKRnano-Plugin/aiida-kkr/aiida_kkr/calculations/kkrnano.py). At least the following keys are in an incorrect format: \
{}.'.format(keys_in_wrong_format)
                )

        #Check if the values are provided in the correct data type
        all_correct_dtype = True
        keys_with_incorrect_dtypes = []
        values_with_incorrect_dtypes = []
        correct_dtypes_for_incorrect_keys = []
        for key in self._DEFAULT_KKRNANO_PARA:
            default_value = self._DEFAULT_KKRNANO_PARA[key]['value']
            dtype_required = type(default_value)

            entered_value = inputdict[key]['value']
            dtype_entered = type(entered_value)

            if dtype_required == list:
                try:
                    for item_index in range(len(default_value)):
                        default_entry = default_value[item_index]
                        entered_entry = entered_value[item_index]
                        if not type(default_entry) == type(entered_entry):
                            all_correct_dtype = False
                            keys_with_incorrect_dtypes.append(key)
                            values_with_incorrect_dtypes.append(str(entered_value))
                            correct_dtypes_for_incorrect_keys.append(
                                f'array/list of {str(type(default_entry))} (wrong entry type)'
                            )
                except IndexError:
                    all_correct_dtype = False
                    keys_with_incorrect_dtypes.append(key)
                    values_with_incorrect_dtypes.append(str(entered_value))
                    correct_dtypes_for_incorrect_keys.append(
                        'array/list of {} (wrong length, \
length should be {} )'.format(str(type(default_entry)), str(len(default_value)))
                    )

            elif not dtype_required == dtype_entered:
                all_correct_dtype = False
                keys_with_incorrect_dtypes.append(key)
                values_with_incorrect_dtypes.append(str(entered_value))
                correct_dtypes_for_incorrect_keys.append(f'{str(type(default_entry))}')

            if not all_correct_dtype:
                #Generate helpful error message
                dtype_error_message = 'Some keys were provided with an incorrect datatype. These values are \n'
                for j in range(len(keys_with_incorrect_dtypes)):
                    dtype_error_message+='key: {}, entered value: {}, needed datatype: {}\n'.format(keys_with_incorrect_dtypes[j],\
                                                values_with_incorrect_dtypes[j],correct_dtypes_for_incorrect_keys[j])
                raise InputValidationError(dtype_error_message)

        #Check if the number of mpi processes, specified with the builder matches the settings in the inputdict
        mpikeys = ['NTHRDS', 'EMPID', 'SMPID', 'num_atom_procs']
        needed_num_mpi_procs = 1
        for key in mpikeys:
            if key in inputdict:
                needed_num_mpi_procs *= inputdict[key]['value']
        if num_mpi_procs % needed_num_mpi_procs != 0 and not convert:
            raise (
                InputValidationError(
                    "Number of MPI processes does not match! In builder.metadata.options.resources['num_mpi_procs']= {}, however according to the builder.parameters node, the needed number of MPI process is {}"
                    .format(num_mpi_procs, needed_num_mpi_procs)
                )
            )
        elif convert and num_mpi_procs != 1:
            #TODO check if convert step works with more cores
            raise (
                InputValidationError(
                    "Number of MPI processes does not match! In builder.metadata.options.resources['num_mpi_procs']= {}, however for the convert step the needed number of MPI processes is 1."
                    .format(num_mpi_procs)
                )
            )

        return inputdict

    def _is_KkrnanoCalc(self, calc):
        """
        check if calc contains the file out_potential
        """
        is_KKR = False
        if calc.process_type == 'aiida.calculations:kkr.kkrnano':
            retrieved_node = calc.get_retrieved_node()
            if True:  #'out_potential' in retrieved_node.list_object_names():
                is_KKR = True

        return is_KKR

    @classmethod
    def _get_struc(self, parent_calc):
        """
        Get strucwithpot from a parent_folder (result of a calculation, typically a remote folder)
        """
        try:
            return parent_calc.inputs.strucwithpot.structure
        except:
            return parent_calc.inputs.structure

    @classmethod
    def _has_struc(self, parent_folder):
        """
        Check if parent_folder has strucwithpot information in its input
        """
        success = True
        link_labels = parent_folder.get_incoming().all_link_labels()
        if 'strucwithpot' not in link_labels and 'structure' not in link_labels:
            success = False
        return success

    @classmethod
    def find_parent_struc_from_voro_or_stwpd(self, parent_folder):
        """
        Find the StructureData node recuresively in chain of parent calculations (structure node is input to voronoi calculation or in input StrucWithPotData node of KKRnano calculation)

        returns structure, parent_folder
        """
        iiter = 0
        Nmaxiter = 1000
        parent_folder_tmp = get_remote(parent_folder)
        while not self._has_struc(parent_folder_tmp) and iiter < Nmaxiter:
            parent_folder_tmp = get_remote(get_parent(parent_folder_tmp))
            iiter += 1
            if iiter % 200 == 0:
                print(
                    'Warning: find_parent_structure takes quite long (already searched {} ancestors). Stop after {}'.
                    format(iiter, Nmaxiter)
                )
        if self._has_struc(parent_folder_tmp):
            struc = self._get_struc(parent_folder_tmp)
            return struc, parent_folder_tmp
        else:
            raise ValueError(f'structure not found {parent_folder_tmp}')
