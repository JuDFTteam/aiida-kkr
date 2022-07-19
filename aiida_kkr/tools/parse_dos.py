#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is the calcfunction that parses the KKRhost DOS output, include interpolation of complex.dos files
"""

from aiida import orm
from aiida.engine import calcfunction
from masci_tools.io.common_functions import (interpolate_dos, get_Ry2eV)


@calcfunction
def parse_dosfiles(dos_retrieved):
    """
    parse dos files to XyData nodes
    take files from the retrieved folder of a KKRhost DOS run
    if BdG files are found we also parse and return those results
    """

    if 'complex.dos' in dos_retrieved.list_object_names():
        with dos_retrieved.open('complex.dos') as dosfolder:
            ef, dos, dos_int = interpolate_dos(dosfolder, return_original=True)
    else:
        raise ValueError('`complex.dos` not found in retrieved')

    dosnode, dosnode2 = make_XyData_nodes(dos, dos_int, ef)
    return_nodes = {'dos_data': dosnode, 'dos_data_interpol': dosnode2}

    for labeladd in ['_eh', '_he', '_hole']:
        # also try tp parse BdG output
        if 'complex.dos' + labeladd in dos_retrieved.list_object_names():
            with dos_retrieved.open('complex.dos' + labeladd) as dosfolder:
                ef, dos, dos_int = interpolate_dos(dosfolder, return_original=True)
            dosnode, dosnode2 = make_XyData_nodes(dos, dos_int, ef, labeladd)
            return_nodes['dos_data' + labeladd] = dosnode
            return_nodes['dos_data_interpol' + labeladd] = dosnode2

    return return_nodes


def make_XyData_nodes(dos, dos_int, ef, labeladd=''):
    """
    create XyData nodes from the dos and dos_int arrays, also transform to eV units before saving
    """
    eVscale = get_Ry2eV()
    # convert to eV units
    dos[:, :, 0] = (dos[:, :, 0] - ef) * eVscale
    dos[:, :, 1:] = dos[:, :, 1:] / eVscale
    dos_int[:, :, 0] = (dos_int[:, :, 0] - ef) * eVscale
    dos_int[:, :, 1:] = dos_int[:, :, 1:] / eVscale

    # create output nodes
    dosnode = orm.XyData()
    dosnode.set_x(dos[:, :, 0], 'E-EF', 'eV')
    name = ['tot', 's', 'p', 'd', 'f', 'g']
    name = name[:len(dos[0, 0, 1:]) - 1] + ['ns']
    ylists = [[], [], []]
    for line, _name in enumerate(name):
        ylists[0].append(dos[:, :, 1 + line])
        ylists[1].append(f'dos {_name}')
        ylists[2].append('states/eV')
    dosnode.set_y(ylists[0], ylists[1], ylists[2])
    dosnode.label = 'dos_data' + labeladd
    dosnode.description = 'Array data containing uniterpolated DOS (i.e. dos at finite imaginary part of energy). 3D array with (atoms, energy point, l-channel) dimensions.'

    # now create XyData node for interpolated data
    dosnode2 = orm.XyData()
    dosnode2.set_x(dos_int[:, :, 0], 'E-EF', 'eV')
    ylists = [[], [], []]
    for line, _name in enumerate(name):
        ylists[0].append(dos_int[:, :, 1 + line])
        ylists[1].append(f'interpolated dos {_name}')
        ylists[2].append('states/eV')
    dosnode2.set_y(ylists[0], ylists[1], ylists[2])
    dosnode2.label = 'dos_interpol_data' + labeladd
    dosnode2.description = 'Array data containing interpolated DOS (i.e. dos at real axis). 3D array with (atoms, energy point, l-channel) dimensions.'
    return dosnode, dosnode2
