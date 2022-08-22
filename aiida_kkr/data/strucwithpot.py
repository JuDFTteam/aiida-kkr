from aiida.orm import CalcJobNode, Dict, StructureData, SinglefileData, load_node, Data
from aiida.common.exceptions import InputValidationError
from aiida.common import NotExistent
from aiida_kkr.tools.find_parent import find_parent_structure, get_parent, get_remote

import os

__copyright__ = (u'Copyright (c), 2022, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.1'
__contributors__ = ('Markus Struckmann')


class StrucWithPotData(Data):
    def __init__(self, passedStructure=None, list_of_shapes=None, list_of_pots=None, specified_lattice_constant=False,KKRnanoCalcNode=None,\
                 VoronoiCalcNode=None, KKRhostCalcNode=None, **kwargs):
        """
        Can be created from an AiiDA structure object and a list of shapes and potentials, corresponding to this object. Also possible to create from KKRnanoCalcNode or a VoronoiCalcNode. Instead of a Voronoi calculation also a KKRhost calculation can be read in using VoronoiCalcNode
        """

        super().__init__(**kwargs)

        if not KKRnanoCalcNode == None:
            passedStructure, list_of_shapes, list_of_pots = self.get_strucwithpot_from_KKRnanoConvert(KKRnanoCalcNode)
        if not VoronoiCalcNode == None:
            passedStructure, list_of_shapes, list_of_pots = self.get_strucwithpot_from_Voronoi(VoronoiCalcNode)
        if not KKRhostCalcNode == None:
            passedStructure, list_of_shapes, list_of_pots = self.get_strucwithpot_from_Voronoi(KKRhostCalcNode)

        if not passedStructure == None and not list_of_shapes == None and not list_of_pots == None:
            self.structure = passedStructure
            self.shapes = list_of_shapes
            self.potentials = list_of_pots
        else:
            raise InputValidationError(
                'Please check input. Either a KKRnanoCalcNode or a StructureData object,\
and lists of SingleFileData potential and shape files have to be provided.'
            )

        #Possibly check if input nodes were store to the database, using e.g. node.is_stored
        if not specified_lattice_constant == False:
            self.specified_lattice_constant = specified_lattice_constant

    # Getter and Setter of the properties

    @property
    def structure(self):
        """
        The AiiDA original StructureData object
        """
        return load_node(self.get_attribute('structure'))

    @structure.setter
    def structure(self, struc):
        """
        Set the AiiDA original StructureData object.
        (Note that only the corresponding UUID is stored in the DataBase)

        :raise ValueError:
        """
        struc.store()
        self.set_attribute('structure', struc.uuid)

    @property
    def shapes(self):
        """
        The shape list in exactly the same order as stored in the StructureData object
        """
        return [load_node(item) for item in self.get_attribute('shapes')]

    @shapes.setter
    def shapes(self, shape_list):
        """
        Set the the shape list
        (Note that only the corresponding UUIDs in a list are stored in the DataBase)
        :raise ValueError:
        """
        for item in shape_list:
            item.store()
        shape_list_uuids = [item.uuid for item in shape_list]
        self.set_attribute('shapes', shape_list_uuids)

    @property
    def potentials(self):
        """
        The potential list in exactly the same order as stored in the StructureData object
        """
        return [load_node(item) for item in self.get_attribute('potentials')]

    @potentials.setter
    def potentials(self, potential_list):
        """
        Set the the potential list
        (Note that only the corresponding UUIDs in a list are stored in the DataBase)
        :raise ValueError:
        """
        for item in potential_list:
            item.store()
        potential_list_uuids = [item.uuid for item in potential_list]
        self.set_attribute('potentials', potential_list_uuids)

    @property
    def specified_lattice_constant(self):
        """
        The specified_lattice_constant to use in the input file of KKR calculations in angs
        """
        return load_node(self.get_attribute('specified_lattice_constant'))

    @specified_lattice_constant.setter
    def specified_lattice_constant(self, specified_lattice_constant):
        """
        Set the the specified_lattice_constant in angs

        :raise ValueError:
        """
        specified_lattice_constant.store()
        self.set_attribute('specified_lattice_constant', specified_lattice_constant.uuid)

    def make_shapefun(self):
        """
        Generate a shapefun file for the structure in string
        """
        #number of shapefuns in file
        shapefun = f'    {len(self.shapes)}\n'

        #header (this part is ignored)
        if len(self.shapes) % 4 == 0:
            headerlength = int(len(self.shapes) / 4)
        else:
            headerlength = int(len(self.shapes) / 4) + 1

        for k in range(headerlength):
            shapefun += '  0.100000000000E+01  0.100000000000E+01  0.100000000000E+01  0.100000000000E+01\n'

        #actual shapefuns
        for shape in self.shapes:
            shapefun += shape.get_content()

        return shapefun

    def make_potential(self):
        """
        Generate a potential file for the structure in string
        """
        potential = ''
        for vpot in self.potentials:
            potential += vpot.get_content()
        return potential

    def _check_for_CPA(self, structure, potentiallist):
        """
        returns a list of sites, which have more than one element and a potential list with tuples for theses sites
        """
        CPA_sites = []
        potential_no = 0
        new_potlist = []
        for site in structure.get_pymatgen().sites:
            no_species = len(site.species.elements)
            if no_species > 1:
                CPA_sites.append(site)
                new_potlist.append(tuple(potentiallist[potential_no:potential_no + no_species]))
                potential_no += no_species
            else:
                new_potlist.append(potentiallist[potential_no])
                potential_no += 1
        return CPA_sites, new_potlist

    def get_strucwithpot_from_Voronoi(self, calcnode):
        """
        Generate a strucwithpot object from a KKRnano calculation
        """
        #TODO: Check if input is valid

        structure = find_parent_structure(calcnode)[0]
        cwd = os.getcwd()

        #shapefun
        shapelist = []
        try:
            with calcnode.outputs.retrieved.open('shapefun') as _f:
                shapes = _f.readlines()
        except:
            calcnode.outputs.remote_folder.getfile('shapefun', f'{cwd}/shapefun')
            with open('shapefun', 'r') as reader:
                shapes = reader.readlines()
        lines = []
        for line in range(len(shapes)):
            if shapes[line].find('Shape') > 0:
                lines.append(line)
        lines.append(len(shapes))
        print('Printing out paths of shape files \n__________')
        for j in range(len(lines) - 1):
            shape_string = ''
            for k in range(lines[j], lines[j + 1]):
                shape_string += shapes[k]

            shape_no_filename = 'shape.' + str(j + 1).zfill(7)
            with open(shape_no_filename, 'w') as file:
                file.write(shape_string)
                path = os.path.realpath(file.name)

            print(path)
            abs_path = f'{cwd}/{shape_no_filename}'
            self.put_object_from_file(abs_path, shape_no_filename)  #Problem has to be called via instance
            self.set_attribute(shape_no_filename.replace('.', ''), shape_no_filename)
            with self.open(shape_no_filename, 'r') as _f:
                shapelist.append(SinglefileData(_f.name))
        print('__________')
        #potential
        potentiallist = []
        try:
            with calcnode.outputs.retrieved.open('output.pot') as _f:
                potentials = _f.readlines()
        except:
            filename = 'out_potential'
            calcnode.outputs.remote_folder.getfile(filename, f'{cwd}/{filename}')
            with open(filename, 'r') as reader:
                potentials = reader.readlines()
        lines = []
        #finding lines where the potentials start/stop
        for line in range(len(potentials)):
            if potentials[line].find('exc') > 0:
                if potentials[line].find('SPIN') > 0:
                    if potentials[line].find('SPIN DOWN') > 0:
                        lines.append(line)
                        print('Added spin-resolved potential')
                else:
                    lines.append(line)
        lines.append(len(potentials))
        print('__________')
        #creating new files with the individual potentials
        print('Printing out paths of potential files\n__________')
        for j in range(len(lines) - 1):
            potential_string = ''
            for k in range(lines[j], lines[j + 1]):
                potential_string += potentials[k]
                #potential_string+="\n"
            potential_no_filename = 'vpot.' + str(j + 1).zfill(7)
            with open(potential_no_filename, 'w') as file:
                file.write(potential_string)
                path = os.path.realpath(file.name)
            print('local path: ', path)
            abs_path = f'{cwd}/{potential_no_filename}'
            self.put_object_from_file(abs_path, potential_no_filename)  #Problem has to be called via instance
            self.set_attribute(potential_no_filename.replace('.', ''), potential_no_filename)
            with self.open(potential_no_filename, 'r') as _f:
                potentiallist.append(SinglefileData(_f.name))
                print('repository path: ', _f)
        print('__________')

        #Check for CPA
        if not len(structure.sites) == len(potentiallist):
            CPA_sites, _ = self._check_for_CPA(structure, potentiallist)

            if len(CPA_sites) == 0:
                raise InputValidationError(
                    'The number of sites in the found parent structure does not match the number of obtained potentials.'
                )
            else:
                print(
                    'WARNING: The number of sites in the found parent structure does not match the number of obtained potentials. This might be due to some CPA input which cannot be processed by KKRnano. The following sites are not (single) chemical elements: \n',
                    CPA_sites
                )

        return structure, shapelist, potentiallist

    def get_strucwithpot_from_KKRnanoConvert(self, calcnode):
        """
        Generate a strucwithpot object from a KKRnano calculation
        """
        #Check if input is valid
        if calcnode.inputs.convert.value == False:
            raise InputValidationError(
                'Only the convert step output can be processed! If this has not been done, yet, the parameter `builder.convert=Bool(False)` can be used and the process be run with 1 MPI, to obtain ASCII-potential files'
            )
        structure = find_parent_structure(calcnode)[0]
        pot_list = []
        for item in calcnode.outputs.retrieved.list_object_names():
            if item.find('vpot') == 0:
                with calcnode.outputs.retrieved.open(item) as _f:
                    path = _f.name
                pot_list.append(SinglefileData(path))

        if not len(structure.sites) == len(pot_list):
            raise InputValidationError(
                'The number of sites in the found parent structure does not match the number of obtained potentials. Found {} sites in structure, but {} potential files.'
                .format(len(structure.sites), len(pot_list))
            )

        #using the shapefun of either the input voronoi calculation or input strucwithpot object
        shape_list = []
        # best approach probably not to piggyback on the found structure and find the corresponding structure based on the first voronoi calculation that appears, or strucwithpot input
        #parent_calc=self.inputs.parent_folder#.get_incoming(node_class=CalcJobNode).first().node

        #find shapefunction
        shape_list = self.find_parent_shapefun(calcnode)

        return structure, shape_list, pot_list

    def find_parent_shapefun(self, parent_folder):
        """
        Find the shape files recursively in chain of parent calculations, either to be extracted from "shapefun" file or "shapes" files
        """
        iiter = 0
        Nmaxiter = 1000

        parent_folder_tmp = get_parent(parent_folder)

        print(parent_folder_tmp)
        parent_folder_tmp_listdir = parent_folder_tmp.listdir(
        )  #requires remote ssh connection, therefore much quicker this way
        while not 'shape.0000001' in parent_folder_tmp_listdir and not 'shapefun' in parent_folder_tmp_listdir and iiter < Nmaxiter:
            parent_folder_tmp = get_remote(get_parent(parent_folder_tmp))  #at this point the result is a CalcNode!
            parent_folder_tmp_listdir = parent_folder_tmp.outputs.remote_folder.listdir()
            iiter += 1
            print(parent_folder_tmp)
            if iiter % 200 == 0:
                print(
                    'Warning: finding shapes takes quite long (already searched {} ancestors). Stop after {}'.format(
                        iiter, Nmaxiter
                    )
                )
        if 'shape.0000001' in parent_folder_tmp_listdir or 'shapefun' in parent_folder_tmp_listdir:
            try:
                shapelist = self.extract_shapefuns(parent_folder_tmp.outputs.remote_folder)
            except AttributeError:
                shapelist = self.extract_shapefuns(parent_folder_tmp)
            return shapelist
        else:
            raise ValueError('shape files not found ')

    def extract_shapefuns(self, parent_folder):
        """
        Process shapefun files from remote folder data in the way needed for StrucWithPot
        """
        shapelist = []
        cwd = os.getcwd()
        print('Extracting shapes')
        has_shapefun_file = False
        found_shapes = False  #to see if KKRnano produced some shapes as output (better than using shapefun)
        parent_folder_listdir = parent_folder.listdir()  #quicker this way as based on SSH tunneling
        for filename in parent_folder_listdir:
            if filename.find('shape.') >= 0:
                if has_shapefun_file:
                    shapelist = []
                abs_path = f'{cwd}/{filename}'
                parent_folder.getfile(filename, f'{cwd}/{filename}')
                self.put_object_from_file(abs_path, filename)
                self.set_attribute(filename.replace('.', ''), filename)
                os.remove(filename)
                with self.open(filename, 'r') as _f:
                    shapelist.append(SinglefileData(_f.name))
                    print('Found shape in repsitory:')
                    print(_f.name)
                has_shapefun_file = False
                found_shapes = True
        if 'shapefun' in parent_folder_listdir and not found_shapes:
            filename = 'shapefun'
            print('Shapefun in dir, this part of the program might need more testing')
            abs_path = f'{cwd}/{filename}'
            parent_folder.getfile(filename, f'{cwd}/{filename}')

            with open(filename, 'r') as reader:
                shapes = reader.readlines()


#                 print(os.path.realpath(reader.name))
            lines = []
            for line in range(len(shapes)):
                if shapes[line].find('Shape') > 0:
                    lines.append(line)
            lines.append(len(shapes))
            for j in range(len(lines) - 1):
                shape_string = ''
                for k in range(lines[j], lines[j + 1]):
                    shape_string += shapes[k]

                shape_no_filename = 'shape.' + str(j + 1).zfill(7)
                with open(shape_no_filename, 'w') as file:
                    file.write(shape_string)
                    path = os.path.realpath(file.name)
                print(path)
                abs_path = f'{cwd}/{shape_no_filename}'
                self.put_object_from_file(abs_path, shape_no_filename)  #Problem has to be called via instance
                self.set_attribute(shape_no_filename.replace('.', ''), shape_no_filename)
                with self.open(shape_no_filename, 'r') as _f:
                    shapelist.append(SinglefileData(_f.name))
                os.remove(shape_no_filename)
            has_shapefun_file = True
        if has_shapefun_file:
            print(
                'WARNING: Only a shapefun from some Voronoi input was found, it is possible that the potential does not match the shapefun parameters, unless they are set this way explicitly in the respective input file! It is advisable to use the `write_shapes=1` command in input.conf'
            )
        print('Found shapelist:')
        print(shapelist)
        return shapelist

    def sites(self):
        """
        Gives back a list of sites. Particularly useful for processing CPA input.
        (AiiDA structure site, potential file(s), shapefun)
        """
        struc = self.structure
        pots = self.potentials
        shapes = self.shapes

        _, new_pots = self._check_for_CPA(struc, pots)
        site_list = []
        for j in range(len(struc.sites)):
            site_list.append({
                'StructureDataSite': struc.sites[j],
                'potential': new_pots[j],
                'shapefun': shapes[j % len(shapes)]
            })
        return site_list
