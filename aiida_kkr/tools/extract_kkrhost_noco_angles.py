# -*- coding: utf-8 -*-
"""
Calcfucntion that extracts the nonco angles from the output of a KkrCalculation
"""

from aiida.orm import Dict
from aiida.engine import calcfunction
from aiida_kkr.calculations import KkrCalculation
from numpy import sqrt, loadtxt


@calcfunction
def extract_noco_angles(**kwargs):
    """
    Extract noco angles from retrieved nonco_angles_out.dat files and save as Dict node which can be used as initial values for the next KkrCalculation.
    New angles are compared to old angles and if they are closer thanfix_dir_threshold they are not allowed to change anymore
    """
    # for comparison read previous theta and phi values and the threshold after which the moments are kept fixed
    noco_angles_old = kwargs['old_noco_angles'].get_dict()
    natom = len(noco_angles_old['phi'])
    noco_angles_old = [
        [noco_angles_old['theta'][i], noco_angles_old['phi'][i], noco_angles_old['fix_dir'][i]] for i in range(natom)
    ]
    fix_dir_threshold = kwargs['fix_dir_threshold'].value

    last_retrieved = kwargs['last_retrieved']
    noco_out_name = KkrCalculation._NONCO_ANGLES_OUT
    if noco_out_name in last_retrieved.list_object_names():
        with last_retrieved.open(noco_out_name) as noco_file:
            noco_angles_new = loadtxt(noco_file, usecols=[0, 1])
            # check if theta and phi change less than fix_dir_threshold
            fix_dir = [
                sqrt((noco_angles_new[i][0] - noco_angles_old[i][0])**2 +
                     (noco_angles_new[i][1] - noco_angles_old[i][1])**2) < fix_dir_threshold for i in range(natom)
            ]
            new_initial_noco_angles = Dict({
                'theta': list(noco_angles_new[:, 0]),
                'phi': list(noco_angles_new[:, 1]),
                # convert from numpy.bool_ to standard python bool, otherwise this is not json serializable and cannot be stored
                'fix_dir': [bool(i) for i in fix_dir]
            })
        return new_initial_noco_angles
