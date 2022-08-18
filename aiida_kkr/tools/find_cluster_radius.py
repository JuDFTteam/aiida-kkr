# -*- coding: utf-8 -*-
"""
Helper functions to find the cluster radius
"""

import numpy as np
from aiida.orm import StructureData
from masci_tools.io.common_functions import get_alat_from_bravais


def find_cluster_radius_old(structure, nclsmin, n_max_box=50, nbins=100):
    """
    Takes structure information (cell and site positions) and computes the minimal cluster radius needed
    such that all clusters around all atoms contain more than `nclsmin` atoms.

    :note: Here we assume spherical clusters around the atoms!

    :param structure: input structure for which the clusters are analyzed
    :param nclsmin: minimal number of atoms in the screening cluster
    :param n_max_box: maximal number of supercells in 3D volume
    :param nbins: number of bins in which the cluster number is analyzed

    :returns: minimal cluster radius needed in Angstroem
    :returns: minimal cluster radius needed in units of the lattice constant
    """

    # extract values needed from structure
    cell = np.array(structure.cell)
    pos = np.array([site.position for site in structure.sites])

    # settings for supercell box
    box = int((n_max_box / len(pos))**(1 / 3.) + 0.5)
    # print('maximal number of atoms in box (time number of atoms in unit cell):', (box*2+1)**3)

    # find all positions in the supercell
    all_pos_box = np.zeros_like(pos)
    for i in range(-box, box + 1):
        for j in range(-box, box + 1):
            for k in range(-box, box + 1):
                for site in pos:
                    tmppos = site + i * cell[0] + j * cell[1] + k * cell[2]
                    all_pos_box = np.append(all_pos_box, [tmppos], axis=0)
    all_pos_box = all_pos_box[len(pos):]

    # computer number of atoms in the clusters
    # Attention: assumes spherical clusters!
    rclsmax_ang = -1
    for site in pos:
        tmpdiff = np.sort(np.sqrt(np.sum((all_pos_box - site)**2, axis=1)))[1:]
        rmax = tmpdiff.max()
        clssizes = [len(tmpdiff[tmpdiff < i * rmax]) for i in np.linspace(0, 1, nbins)]
        rclsmax_atom = (np.linspace(0, 1, nbins) * rmax)[np.where(np.array(clssizes) < nclsmin)[0].max() + 1]
        if rclsmax_atom > rclsmax_ang:
            rclsmax_ang = rclsmax_atom

    # convert also to alat units
    rclsmax_alat = rclsmax_ang / get_alat_from_bravais(cell, structure.pbc[2])

    # now the minimal cluster radius needed to get the spherical screening clusters around the atoms larger than
    # nclsmin atoms is found and can be returned
    return rclsmax_ang, rclsmax_alat


def find_cluster_radius(structure, nclsmin=15, Rclsmax=10.):
    """
    Find cluster radius to have >=nclsmin atoms in each cluster.

    :param structure: input structure
    :param nclsmin: required minimal cluster size
    :param Rclsmax: max radius used in screening for cluster size

    :return R0: cluster radius in Ang. units
    :return Ncls_all: cluster sizes, list with the lenth of the sites
    """

    # create auxiliary structure (makes sure this also works with 2D structures)
    saux = StructureData(cell=structure.cell)
    for site in structure.sites:
        kind = structure.get_kind(site.kind_name)
        saux.append_atom(position=site.position, symbols=kind.symbols, weights=kind.weights)
    ps = saux.get_pymatgen()

    # find cluster radius that has at least the minimal number of atoms
    R0 = -1.
    neighbors_all = ps.get_all_neighbors(max(Rclsmax, np.linalg.norm(structure.cell, axis=1).max()))
    for neighbors in neighbors_all:
        dist_all = np.sort([np.linalg.norm(n.coords) for n in neighbors])
        # -2 because we take the (nclsmin-1)-th position
        # in the array (position at (0,0,0) is added automatically)
        R0 = max(R0, dist_all[nclsmin - 2])

    # consistency check
    if R0 <= 0.:
        raise ValueError('Cluster radius not found properly')

    # count number of atoms in the clusters
    neighbors_all = ps.get_all_neighbors(R0)
    Ncls_all = [1 + len(n) for n in neighbors_all]

    return R0, Ncls_all
