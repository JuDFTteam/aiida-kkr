# -*- coding: utf-8 -*-
"""
Kick out core states from the potential that are above a certain threshold
"""

from aiida.engine import calcfunction
from masci_tools.io.common_functions import open_general


def kick_out_corestates(potfile, potfile_out, emin):
    """
    Read potential file and kick out all core states that lie higher than emin.
    If no core state lies higher than emin then the output potential will be the same as the input potential
    :param potfile: input potential
    :param potfile_out: output potential where some core states are kicked out
    :param emin: minimal energy above which all core states are kicked out from potential
    :returns: number of lines that have been deleted
    """
    from masci_tools.io.common_functions import get_corestates_from_potential
    from numpy import where, array

    # read core states
    nstates, energies, lmoments = get_corestates_from_potential(potfile)

    # read potential file
    with open_general(potfile) as f:
        txt = f.readlines()

    # get start of each potential part
    istarts = [iline for iline in range(len(txt)) if 'POTENTIAL' in txt[iline]]
    all_lines = list(range(len(txt)))  # index array

    # change list of core states
    for ipot in range(len(nstates)):
        if nstates[ipot] > 0:
            m = where(energies[ipot] > emin)
            if len(m[0]) > 0:
                istart = istarts[ipot]
                # change number of core states in potential
                # print(txt[istart+6])
                txt[istart + 6] = '%i 1\n' % (nstates[ipot] - len(m[0]))
                # now remove energy line accordingly
                for ie_out in m[0][::-1]:
                    m_out = where(array(all_lines) == istart + 6 + ie_out + 1)[0][0]
                    e_out = all_lines.pop(m_out)

    # find number of deleted lines
    num_deleted = len(txt) - len(all_lines)

    if num_deleted > 0:
        # write output potential
        with open_general(potfile_out, u'w') as f2:
            txt_new = []
            for iline in all_lines:
                txt_new.append(str(txt[iline]))
            f2.writelines(txt_new)

    # return number of lines that were deleted
    return num_deleted


@calcfunction
def kick_out_corestates_wf(potential_sfd, emin):
    """
    Workfunction that kicks out all core states from single file data potential that are higher than emin.
    :param potential_sfd: SinglefileData type of potential
    :param emin: Energy threshold above which all core states are removed from potential (Float)
    :returns: potential without core states higher than emin (SinglefileData)
    """
    from aiida.common.folders import SandboxFolder
    from aiida.orm import SinglefileData

    with SandboxFolder() as tmpdir:
        with tmpdir.open('potential_deleted_core_states', 'w') as potfile_out:
            with potential_sfd.open(potential_sfd.filename) as potfile_in:
                num_deleted = kick_out_corestates(potfile_in, potfile_out, emin)
        # store new potential as single file data object
        if num_deleted > 0:
            with tmpdir.open('potential_deleted_core_states', 'rb') as potfile_out:
                potential_nocore_sfd = SinglefileData(file=potfile_out)

    # return potential
    if num_deleted > 0:
        return potential_nocore_sfd
    else:
        return potential_sfd.clone()
