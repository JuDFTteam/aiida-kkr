# -*- coding: utf-8 -*-
"""
This contains code snippets for making a supercell in one direction. It can be run repeatedly in various directions, in order to generate a larger supercell.
"""

#Check if accuracy is high enough!
from pymatgen.core import Lattice, Structure
import numpy as np


def make_super_cell(structure, scale_bravais_c, bravais_axis=2):
    matrix0 = np.array(np.asarray(structure.lattice.matrix), dtype=np.float128)
    coords = np.array(np.asarray(structure.cart_coords), dtype=np.float128)
    atomic_numbers = list(structure.atomic_numbers)
    new_atomic_numbers = list(structure.atomic_numbers)
    new_coords = np.copy(coords)
    for k in np.arange(1, scale_bravais_c):
        # print(np.shape(coords))
        new_coords_block = np.zeros(np.shape(coords))
        for j in np.arange(np.shape(coords)[0]):
            new_coords_block[j, :] = coords[j, :] + k * matrix0[bravais_axis, :]

        new_coords = np.concatenate((new_coords, new_coords_block), axis=0)
        new_atomic_numbers = new_atomic_numbers + atomic_numbers

    # print(new_coords)
    # print(new_atomic_numbers)
    matrix = np.copy(matrix0)
    matrix[2, :] = scale_bravais_c * matrix0[2, :]
    latt = Lattice(matrix)
    struc = Structure(latt, new_atomic_numbers, new_coords, coords_are_cartesian=True)
    return struc
