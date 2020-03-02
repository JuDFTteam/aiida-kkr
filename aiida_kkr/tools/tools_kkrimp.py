#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tools for the impurity caluclation plugin and its workflows
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals
from builtins import object, str
from six.moves import range
from six.moves import input
from masci_tools.io.common_functions import open_general
from masci_tools.io.common_functions import search_string


__copyright__ = (u"Copyright (c), 2018, Forschungszentrum Jülich GmbH,"
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.6.1"
__contributors__ = u"Philipp Rüßmann"


class modify_potential(object):
    """
    Class for old modify potential script, ported from modify_potential script, initially by D. Bauer
    """

    def _check_potstart(self, str1, mode='pot', shape_ver='new'):
        if mode=='shape':
            if shape_ver=='new':
                check1='Shape number' in str1
            else:
                check1= (len(str1)==11)
        else:
            check1='exc:' in str1
        return check1

    def _read_input(self, filepath):
        with open_general(filepath) as f:
            data = f.readlines()
            filepathname = f.name
        
        if 'shapefun' in filepathname:
            mode = 'shape'
        else:
            mode = 'pot'

        # read file
        index1=[];index2=[]
        for i in range(len(data)):
          if self._check_potstart(data[i], mode=mode):
            index1.append(i)
            if len(index1)>1: index2.append(i-1)
        index2.append(i)

        # read shapefun if old style is used
        if mode=='shape' and len(index1)<1:
            index1=[];index2=[]
            for i in range(len(data)):
                if self._check_potstart(data[i], mode=mode, shape_ver='old'):
                    index1.append(i)
                if len(index1)>1: index2.append(i-1)
                index2.append(i)

        """
        print(index1)
        print(index2)

        print('Potential file read')
        print('found %i potentials in file'%len(index1))
        print('')
        """

        return index1, index2, data

    def shapefun_from_scoef(self, scoefpath, shapefun_path, atom2shapes, shapefun_new):
        """
        Read shapefun and create impurity shapefun using scoef info and shapes array

        :param scoefpath: absolute path to scoef file
        :param shapefun_path: absolute path to input shapefun file
        :param shapes: shapes array for mapping between atom index and shapefunction index
        :param shapefun_new: absolute path to output shapefun file to which the new shapefunction will be written
        """
        index1, index2, data = self._read_input(shapefun_path)

        order=list(range(len(index1)))

        natomtemp = int(open_general(scoefpath).readlines()[0])
        filedata=open_general(scoefpath).readlines()[1:natomtemp+1]
        listnew=[]
        for line in filedata:
            if (len(line.split())>1):
                listnew.append(atom2shapes[int(line.split()[3])-1]-1)
        order = listnew

        datanew=[]
        for i in range(len(order)):
          for ii in range(index1[order[i]], index2[order[i]]+1  ):
            datanew.append(data[ii])

        # add header to shapefun_new
        tmp = datanew
        datanew = []
        datanew.append('   %i\n' %(len(order)))
        datanew.append('  1.000000000000E+00\n')
        datanew += tmp
        with open_general(shapefun_new,'w') as f:
            f.writelines(datanew)

    def neworder_potential(self, potfile_in, potfile_out, neworder, potfile_2=None, replace_from_pot2=None):
        """
        Read potential file and new potential using a list describing the order of the new potential.
        If a second potential is given as input together with an index list, then the corresponding of
        the output potential are overwritten with positions from the second input potential.

        :param potfile_in: absolute path to input potential
        :type potfile_in: str
        :param potfile_out: absolute path to output potential
        :type potfile_out: str
        :param neworder: list after which output potential is constructed from input potential
        :type neworder: list
        :param potfile_2: optional, absolute path to second potential file if
            positions in new list of potentials shall be replaced by positions of
            second potential, requires *replace_from_pot* to be given as well
        :type potfile_2: str
        :param replace_from_pot: optional, list containing tuples of (position
            in newlist that is to be replaced, position in pot2 with which position
            is replaced)
        :type replace_from_pot: list

        :usage:
            1. modify_potential().neworder_potential(<path_to_input_pot>, <path_to_output_pot>, [])
        """
        from numpy import array, shape

        index1, index2, data = self._read_input(potfile_in)

        if potfile_2 is not None:
            index12, index22, data2 = self._read_input(potfile_2)
            # check if also replace_from_pot2 is given correctly
            if replace_from_pot2 is None:
                raise ValueError('replace_from_pot2 not given')
            else:
                replace_from_pot2 = array(replace_from_pot2)
                if shape(replace_from_pot2)[1]!=2:
                    raise ValueError('replace_from_pot2 needs to be a 2D array!')
        else:
            if replace_from_pot2 is not None:
                raise ValueError('replace_from_pot2 given but potfile_2 not given')

        # set order in which potential file is written
        # ensure that numbers are integers:
        order = [int(i) for i in neworder]

        datanew=[]
        for i in range(len(order)):
            # check if new position is replaced with position from old pot
            if replace_from_pot2 is not None and i in replace_from_pot2[:,0]:
                    replace_index = replace_from_pot2[replace_from_pot2[:,0]==i][0][1]
                    for ii in range(index12[replace_index], index22[replace_index]+1 ):
                        datanew.append(data2[ii])
            else: # otherwise take new potntial according to input list
                    for ii in range(index1[order[i]], index2[order[i]]+1 ):
                        datanew.append(data[ii])

        # write out new potential
        with open_general(potfile_out,'w') as f:
            f.writelines(datanew)


####################################################################################
# for create scoef functions
import numpy as np

def get_structure_data(structure):
    """
    Function to take data from AiiDA's StructureData type and store it into a single numpy array of the following form:
    a = [[x-Position 1st atom, y-Position 1st atom, z-Position 1st atom, index 1st atom, charge 1st atom, 0.],
         [x-Position 2nd atom, y-Position 2nd atom, z-Position 2nd atom, index 2nd atom, charge 1st atom, 0.],
         [..., ..., ..., ..., ..., ...],
         ...
         ]

    :param structure: input structure of the type StructureData

    :return: numpy array a[# of atoms in the unit cell][5] containing the structure related data (positions in units
             of the unit cell length)

    :note:
    """

    #import packages
    from aiida.common.constants import elements as PeriodicTableElements
    import numpy as np

    #get the connection between coordination number and element symbol
    _atomic_numbers = {data['symbol']:num for num,
                data in PeriodicTableElements.items()}

    #initialize the array that will be returned later (it will be a (# of atoms in the cell) x 6-matrix)
    a = np.zeros((len(structure.sites),6))
    k = 0 #running index for filling up the array with data correctly
    charges = [] #will be needed to return the charge number for the different atoms later

    #loop to fill up the a-array with positions, index, charge and a 0. for every atom in the cell
    sites = structure.sites
    n = len(structure.sites) + 1 #needed to do the indexing of atoms
    m = len(structure.sites) #needed to do the indexing of atoms
    for site in sites:
        for j in range(3):
            a[k][j] = site.position[j]
        sitekind = structure.get_kind(site.kind_name)
        naez = n - m
        m = m - 1
        #convert the element symbol from StructureData to the charge number
        for ikind in range(len(sitekind.symbols)):
            site_symbol = sitekind.symbols[ikind]
            charges.append(_atomic_numbers[site_symbol])
        i = len(charges) - 1
        a[k][3] = int(naez)
        a[k][4] = float(charges[i])
        k = k + 1

    return a

def select_reference(structure_array, i):
    """
    Function that references all of the atoms in the cell to one particular atom i in the cell and calculates
    the distance from the different atoms to atom i. New numpy array will have the form:
    x = [[x-Position 1st atom, y-Position 1st atom, z-Position 1st atom, index 1st atom, charge 1st atom,
            distance 1st atom to atom i],
         [x-Position 2nd atom, y-Position 2nd atom, z-Position 2nd atom, index 2nd atom, charge 1st atom,
            distance 1st atom to atom i],
         [..., ..., ..., ..., ..., ...],
         ...
         ]

    :param structure_array: input array of the cell containing all the atoms (obtained from get_structure_data)
    :param i: index of the atom which should be the new reference

    :return: new structure array with the origin at the selected atom i (for KKRimp: impurity atom)

    :note: the first atom in the structure_array is labelled with 0, the second with 1, ...
    """

    #import packages
    import numpy as np

    #initialize new array for data centered around atom i
    x = np.zeros((len(structure_array),6))

    #take the positions of atom i as new origin
    x_ref = np.array([structure_array[i][0], structure_array[i][1], structure_array[i][2], 0, 0, 0])

    #calculate the positions and distances for all atoms in the cell with respect to the chosen atom i
    for k in range(len(structure_array)):
        x[k][5] = get_distance(structure_array, i, k)
        for j in range(5):
            x[k][j] = structure_array[k][j] - x_ref[j]

    return x

def get_distance(structure_array, i, j):
    """
    Calculates and returns the distances between to atoms i and j in the given structure_array

    :param structure_array: input numpy array of the cell containing all the atoms ((# of atoms) x 6-matrix)
    :params i, j: indices of the atoms for which the distance should be calculated (indices again start at 0)

    :return: distance between atoms i and j in units of alat

    :note:
    """

    #import math package for square root calculation
    import math

    #calculate x-, y-, z-components of distance of i and j
    del_x = structure_array[i][0] - structure_array[j][0]
    del_y = structure_array[i][1] - structure_array[j][1]
    del_z = structure_array[i][2] - structure_array[j][2]

    #return absolute value of the distance of atom i and j
    return math.sqrt(del_x*del_x + del_y*del_y + del_z*del_z)

def rotate_onto_z(structure, structure_array, vector):
    """
    Rotates all positions of a structure array of orientation 'orient' onto the z-axis. Needed to implement the
    cylindrical cutoff shape.

    :param structure: input structure of the type StructureData
    :param structure_array: input structure array, obtained by select_reference for the referenced system.
    :param vector: reference vector that has to be mapped onto the z-axis.

    :return: rotated system, now the 'orient'-axis is aligned with the z-axis
    """

    from masci_tools.io.common_functions import vec_to_angles
    import math
    import numpy as np

    #get angles, from vector
    angles = vec_to_angles(vector)
    theta = angles[1]
    phi = angles[2]

    #initialize needed arrays
    x_res = np.delete(structure_array, np.s_[3:6], 1)
    x_temp_1 = np.delete(structure_array, np.s_[3:6], 1)
    x_temp_2 = np.delete(structure_array, np.s_[3:6], 1)

    #define rotation matrices
    #========================
    #rotation around z-axis with angle phi
    R_z = np.array([[math.cos(-phi), -math.sin(-phi), 0.],
                    [math.sin(-phi), math.cos(-phi), 0.],
                    [0., 0., 1]])
    #rotation around y-axis with angle theta
    R_y = np.array([[math.cos(-theta), 0, math.sin(-theta)],
                    [0., 1., 0.],
                    [-math.sin(-theta), 0., math.cos(-theta)]])

    #first rotate around z-axis
    for i in range(len(structure_array)):
        x_temp_1[i] = np.dot(R_z, x_res[i])
        x_temp_2[i] = np.dot(R_y, x_temp_1[i])

    return x_temp_2

def find_neighbors(structure, structure_array, i, radius, clust_shape='spherical', h=0., vector=[0., 0., 1.]):
    """
    Applies periodic boundary conditions and obtains the distances between the selected atom i in the cell and
    all other atoms that lie within a cutoff radius r_cut. Afterwards an numpy array with all those atoms including
    atom i (x_res) will be returned.

    :param structure: input parameter of the StructureData type containing the three bravais lattice cell vectors
    :param structure_array: input numpy structure array containing all the structure related data
    :param i: centered atom at which the origin lies (same one as in select_reference)
    :param radius: Specifies the radius of the cylinder or of the sphere, depending on clust_shape.
                   Input in units of the lattice constant.
    :param clust_shape: specifies the shape of the cluster that is used to determine the neighbors for the 'scoef' file.
                        Default value is 'spherical'. Other possible forms are 'cylindrical' ('h' and 'orient'
                        needed), ... .
    :param h: needed for a cylindrical cluster shape. Specifies the height of the cylinder. Default=0.
              Input in units of the lattice constant.
    :param vector: needed for a cylindrical cluster shape. Specifies the orientation vector of the cylinder. Default:
                   z-direction.

    :return: array with all the atoms within the cutoff (x_res)

    :ToDo: - dynamical box construction (r_cut determines which values n1, n2, n3 have)
           - different cluster forms (spherical, cylinder, ...), add default parameters, better solution for 'orient'
    """

    #import packages
    import numpy as np
    import math

    # make sure we take the correct cell length for 2D structures
    if not structure.pbc[2]:
        allpos = np.array([isite.position for isite in structure.sites])
        maxpos = allpos[allpos[:,2]==max(allpos[:,2])]
        c3 = maxpos[0]
        c3[2]+=10 # make sure there is no overlap
        sl3 = np.sum(np.sqrt(c3**2))
    else:
        sl3 = structure.cell_lengths[2]
        c3 = structure.cell[2]
    

    #obtain cutoff distance from radius and h
    dist_cut = max(radius, h)

    #initialize arrays and reference the system
    x = select_reference(structure_array, i)
    x_temp = np.array([x[i]])

    #calculate needed amount of boxes in all three directions
    #========================================================
    #spherical approach (same distance in all three directions)
    if clust_shape == 'spherical':
        box_1 = int(radius/structure.cell_lengths[0] + 1)
        box_2 = int(radius/structure.cell_lengths[1] + 1)
        box_3 = int(radius/sl3 + 1)
    #cylindrical shape (different distances for the different directions)
    elif clust_shape == 'cylindrical':
        maxval = max(h/2., radius)
        box_1 = int(maxval/structure.cell_lengths[0] + 1)
        box_2 = int(maxval/structure.cell_lengths[1] + 1)
        box_3 = int(maxval/sl3 + 1)
    #================================================================================================================

    #create array of all the atoms in an expanded system
    box = max(box_1, box_2, box_3)
    cell = np.array(structure.cell)
    cell[2] = c3
    for j in range(len(x)):
        for n in range(-box, box + 1):
            for m in range(-box, box + 1):
                for l in range(-box, box + 1):
                    x_temp = np.append(x_temp, [[x[j][0] + (n*cell[0][0] + m*cell[1][0] + l*cell[2][0]),
                                                 x[j][1] + (n*cell[0][1] + m*cell[1][1] + l*cell[2][1]),
                                                 x[j][2] + (n*cell[0][2] + m*cell[1][2] + l*cell[2][2]),
                                                 x[j][3], x[j][4], 0.]], axis = 0)

    #x_temp now contains all the atoms and their positions regardless if they are bigger or smaller than the cutoff
    x_new = x_temp

    #calculate the distances between all the atoms and the center atom i
    for j in range(len(x_temp)):
        x_new[j][5] = get_distance(x_temp, 0, j)

    #initialize result array
    x_res = np.array([x[i]])

    #only take atoms into account whose distance to atom i is smaller than the cutoff radius
    #dist_cut = dist_cut
    if clust_shape == 'spherical':
        for j in range(len(x_temp)):
            if x_new[j][5] <= dist_cut and x_new[j][5] > 0.:
                x_res = np.append(x_res, [[x_temp[j][0],
                                             x_temp[j][1],
                                             x_temp[j][2],
                                             x_temp[j][3], x_temp[j][4], x_new[j][5]]], axis=0)
    elif clust_shape == 'cylindrical':
        for j in range(len(x_temp)):
            #rotate system into help system that is aligned with the z-axis
            x_help = rotate_onto_z(structure, x_temp, vector)

            #calculate in plane distance and vertical distance
            vert_dist = np.absolute(x_help[j][2])
            inplane_dist = math.sqrt(x_help[j][0]**2 + x_help[j][1]**2)
            #print(vert_dist, inplane_dist)
            if vert_dist <= h/2. and inplane_dist <= radius and x_new[j][5] > 0.:
                x_res = np.append(x_res, [[x_temp[j][0],
                                             x_temp[j][1],
                                             x_temp[j][2],
                                             x_temp[j][3], x_temp[j][4], x_new[j][5]]], axis=0)

    #return an unordered array of all the atoms which are within the cutoff distance with respect to atom i
    return x_res

def write_scoef(x_res, path):
    """
    Sorts the data from find_neighbors with respect to the distance to the selected atom and writes the data
    correctly formatted into the file 'scoef'. Additionally the total number of atoms in the list is written out
    in the first line of the file.

    :param x_res: array of atoms within the cutoff radius obtained by find_neighbors (unsorted)

    :output: returns scoef file with the total number of atoms in the first line, then with the formatted positions,
             indices, charges and distances in the subsequent lines.
    """

    #sort the data from x_res with respect to distance to the centered atom
    m = x_res[:,-1].argsort()
    x_res = x_res[m]

    #write data of x_res into the 'scoef'-file
    with open_general(path, 'w') as file:
        file.write(str("{0:4d}".format(len(x_res))))
        file.write("\n")
        for i in range(len(x_res)):
            file.write(str("{0:26.19e}".format(x_res[i][0])))
            file.write(" ")
            file.write(str("{0:26.19e}".format(x_res[i][1])))
            file.write(" ")
            file.write(str("{0:26.19e}".format(x_res[i][2])))
            file.write(" ")
            file.write(str("{0:4d}".format(int(x_res[i][3]))))
            file.write(" ")
            file.write(str("{0:4.1f}".format(x_res[i][4])))
            file.write(" ")
            file.write(str("{0:26.19e}".format(x_res[i][5])))
            file.write("\n")

def make_scoef(structure, radius, path, h=-1., vector=[0., 0., 1.], i=0, alat_input=None):
    """
    Creates the 'scoef' file for a certain structure. Needed to conduct an impurity KKR calculation.

    :param structure: input structure of the StructureData type.
    :param radius: input cutoff radius in Ang. units.
    :param h: height of the cutoff cylinder (negative for spherical cluster shape). For negative values, clust_shape
              will be automatically assumed as 'spherical'. If there will be given a h > 0, the clust_shape
              will be 'cylindrical'.
    :param vector: orientation vector of the cylinder (just for clust_shape='cylindrical').
    :param i: atom index around which the cluster should be centered. Default: 0 (first atom in the structure).
    :param alat_input: input lattice constant in Ang. If `None` use the lattice constant that is automatically found. Otherwise rescale everything.
    """
    from masci_tools.io.common_functions import get_alat_from_bravais, get_aBohr2Ang
    from numpy import array

    #shape of the cluster is specified
    if h < 0.:
        clust_shape = 'spherical'
    else:
        clust_shape = 'cylindrical'

    #store data from StructureData type in an numpy array
    structure_array = get_structure_data(structure)

    #creates an array with all the atoms within a certain cluster shape and size with respect to atom i
    c = find_neighbors(structure, structure_array, i, radius, clust_shape, h, vector)

    # convert to internal units (units of the lattice constant)
    alat = get_alat_from_bravais(array(structure.cell), structure.pbc[2])
    if alat_input is not None:
        # use input lattice constant instead of automatically found alat
        alat = alat_input
    # now take out alat factor
    c[:,:3] = c[:,:3] / alat # rescale atom positions
    c[:,-1] = c[:,-1] / alat # rescale distances


    #writes out the 'scoef'-file
    write_scoef(c, path)
    return c
