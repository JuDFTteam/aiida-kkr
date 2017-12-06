# -*- coding: utf-8 -*-
"""
Parser for the KKR Code.
The parser should never fail, but it should catch 
all errors and warnings and show them to the user.
"""

if __name__=='__main__':
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_kkr.tools.common_functions import (search_string, get_version_info, get_Ry2eV,
                                              get_corestates_from_potential, get_highest_core_state)
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida.common.exceptions import InputValidationError


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = ("Jens Broeder", "Philipp Rüßmann")


class KkrParser(Parser):
    """
    Parser class for parsing output of KKR code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of KkrParser
        """
        # check for valid input
        if not isinstance(calc, KkrCalculation):
            raise InputValidationError("Input calc must be a KkrCalculation")
        
        self._ParserVersion = __version__

        #reuse init of base class
        super(KkrParser, self).__init__(calc)
        

    # pylint: disable=protected-access
    def parse_with_retrieved(self, retrieved):
        """
        Parse output data folder, store results in database.

        :param retrieved: a dictionary of retrieved nodes, where
          the key is the link name
        :returns: a tuple with two values ``(bool, node_list)``, 
          where:

          * ``bool``: variable to tell if the parsing succeeded
          * ``node_list``: list of new nodes to be stored in the db
            (as a list of tuples ``(link_name, node)``)
        """
        success = False
        node_list = ()

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error("No retrieved folder found")
            return success, node_list

        # check what is inside the folder
        list_of_files = out_folder.get_folder_list()

        # we need at least the output file name as defined in calcs.py
        if self._calc._OUTPUT_FILE_NAME not in list_of_files:
            self.logger.error("Output file not found")
            return success, node_list

        # Parse output files of KKR calculation
        outfile = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
        outfile_0init = out_folder.get_abs_path(self._calc._OUTPUT_0_INIT)
        outfile_000 = out_folder.get_abs_path(self._calc._OUTPUT_000)
        #potfile_in = out_folder.get_abs_path(self._calc._POTENTIAL)
        potfile_out = out_folder.get_abs_path(self._calc._OUT_POTENTIAL)
        timing_file = out_folder.get_abs_path(self._calc._OUT_TIMING_000)
        #scoef_file = out_folder.get_abs_path(self._calc._SCOEF)
        #nonco_out_file = out_folder.get_abs_path(self._calc._NONCO_ANGLES_OUT)
        
        
        out_dict = {'parser_version': self._ParserVersion}
        #TODO job title, compound description
        success, msg_list, out_dict = parse_kkr_outputfile(out_dict, outfile, 
                                                           outfile_0init, 
                                                           outfile_000, 
                                                           timing_file, 
                                                           potfile_out)
        out_dict['parser_warnings'] = msg_list
        
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        
        return success, node_list
    
    
##########################################################
# here follow the parser functions:


def parse_array_float(outfile, searchstring, splitinfo, replacepair=None):
    from numpy import array
    f = open(outfile)
    tmptxt = f.readlines()
    f.close()
    itmp = 0
    res = []
    while itmp>=0:
        itmp = search_string(searchstring, tmptxt)
        if itmp>=0:
            tmpval = tmptxt.pop(itmp)
            if replacepair is not None:
                tmpval = tmpval.replace(replacepair[0], replacepair[1])
            if splitinfo[0]==1:
                tmpval = float(tmpval.split(splitinfo[1])[splitinfo[2]])
            elif splitinfo[0]==2:
                tmpval = float(tmpval.split(splitinfo[1])[splitinfo[2]].split()[splitinfo[3]])
            else:
                raise ValueError("splitinfo[0] has to be either 1 or 2")
            res.append(tmpval)
    res = array(res)
    return res


def get_rms(outfile, outfile2):
    res = parse_array_float(outfile, 'average rms-error', [1, '=', 1], ['D', 'E'])
    res2 = parse_array_float(outfile2, 'rms-error for atom', [1, '=', 1], ['D', 'E'])
    niter = len(res) # number of iterations
    natoms = len(res2)/niter # number of atoms in system, needed to take only atom resolved rms of last iteration
    return res, res2[-natoms:]


def get_neutr(outfile):
    res = parse_array_float(outfile, 'charge neutrality in unit cell', [1, '=', 1])
    return res


def get_magtot(outfile):
    res = parse_array_float(outfile, 'TOTAL mag. moment in unit cell', [1, '=', 1])
    return res


def get_EF(outfile):
    res = parse_array_float(outfile, 'E FERMI', [2, 'FERMI', 1, 0])
    return res


def get_DOS_EF(outfile):
    res = parse_array_float(outfile, 'DOS(E_F)', [1, '=', 1])
    return res


def get_Etot(outfile):
    res = parse_array_float(outfile, 'TOTAL ENERGY in ryd.', [1, ':', 1])
    return res


def find_warnings(outfile):
    from numpy import array
    f = open(outfile)
    tmptxt = f.readlines()
    tmptxt_caps = [txt.upper() for txt in tmptxt]
    f.close()
    itmp = 0
    res = []
    while itmp>=0:
        itmp = search_string('WARNING', tmptxt_caps)
        if itmp>=0:
            tmpval = tmptxt_caps.pop(itmp)
            tmpval = tmptxt.pop(itmp)
            res.append(tmpval.strip())
    return array(res)


def extract_timings(outfile):
    from numpy import array
    f = open(outfile)
    tmptxt = f.readlines()
    f.close()
    itmp = 0
    res = []
    search_keys = ['main0               ', 'main1a - tbref      ', 
                   'main1a              ', 'main1b - calctref13 ', 
                   'main1b              ', 'main1c - serial part', 
                   'main1c              ', 'main2               ', 
                   'Time in Iteration   ']
    while itmp>=0:
        tmpvals = []
        for isearch in search_keys:
            itmp = search_string(isearch, tmptxt)
            if itmp>=0:
                tmpval = float(tmptxt.pop(itmp).split()[-1])
                tmpvals.append(tmpval)
        if len(tmpvals)>0:
            res.append(tmpvals)
    res = array(res)[-1]
    return [[search_keys[i].strip(), res[i]] for i in range(len(res))]


def get_charges_per_atom(outfile_000):
    res1 = parse_array_float(outfile_000, 'charge in wigner seitz sphere', [1, '=', 1])
    res2 = parse_array_float(outfile_000, 'nuclear charge', [2, 'nuclear charge', 1, 0])
    res3 = parse_array_float(outfile_000, 'core charge', [1, '=', 1])
    return res1, res2, res3


def get_single_particle_energies(outfile_000):
    """
    extracts single particle energies from outfile_000 (output.000.txt)
    returns the valence contribution of the single particle energies
    """
    from numpy import array
    f = open(outfile_000)
    tmptxt = f.readlines()
    f.close()
    itmp = 0
    res = []
    while itmp>=0:
        itmp = search_string('band energy per atom', tmptxt)
        if itmp>=0:
            tmpval = float(tmptxt.pop(itmp).split()[-1])
            res.append(tmpval)
    return array(res)


def get_econt_info(outfile_0init):
    f = open(outfile_0init)
    tmptxt = f.readlines()
    f.close()
    
    itmp = search_string('E min', tmptxt)
    emin = float(tmptxt[itmp].split('min')[1].split('=')[1].split()[0])
    
    itmp = search_string('Temperature', tmptxt)
    tempr = float(tmptxt[itmp].split('Temperature')[1].split('=')[1].split()[0])
    
    itmp = search_string('Number of energy points', tmptxt)
    Nepts = int(tmptxt[itmp].split(':')[1].split()[0])
    
    itmp = search_string('poles =', tmptxt)
    Npol = int(tmptxt[itmp].split('=')[1].split()[0])
    
    itmp = search_string('contour:', tmptxt)
    tmp = tmptxt[itmp].replace(',','').split(':')[1].split()
    N1 = int(tmp[2])
    N2 = int(tmp[5])
    N3 = int(tmp[8])
    
    return emin, tempr, Nepts, Npol, N1, N2, N3


def get_core_states(potfile):
    from numpy import array
    ncore, energies, lmoments = get_corestates_from_potential(potfile=potfile)
    emax, lmax, descr_max = [], [], []
    for ipot in range(len(ncore)):
        if ncore[ipot] > 0:
            lvalmax, energy_max, descr = get_highest_core_state(ncore[ipot], energies[ipot], lmoments[ipot])
        else:
            lvalmax, energy_max, descr = None, None, 'no core states'
        emax.append(energy_max)
        lmax.append(lvalmax)
        descr_max.append(descr)
    return array(ncore), array(emax), array(lmax), array(descr_max)


def get_alatinfo(outfile_0init):
    f = open(outfile_0init)
    tmptxt = f.readlines()
    f.close()
    itmp = search_string('Lattice constants :', tmptxt)
    alat = float(tmptxt[itmp].split(':')[1].split('=')[1].split()[0])
    twopialat = float(tmptxt[itmp].split(':')[1].split('=')[2].split()[0])
    return alat, twopialat


def get_scfinfo(outfile_0init, outfile_000, outfile):
    f = open(outfile_000)
    tmptxt = f.readlines()
    f.close()
    
    itmp = search_string('ITERATION :', tmptxt)
    tmpval = tmptxt[itmp].split(':')[1].split()
    niter = int(tmpval[0])
    nitermax = int(tmpval[3])
    
    f = open(outfile)
    tmptxt = f.readlines()
    f.close()
    itmp1 = search_string('SCF ITERATION CONVERGED', tmptxt)
    itmp2 = search_string('NUMBER OF SCF STEPS EXHAUSTED', tmptxt)
    if itmp1>=0:
        converged = True
    else:
        converged = False
    if itmp2>=0:
        nmax_reached = True
    else:
        nmax_reached = False
       
    f = open(outfile_0init)
    tmptxt = f.readlines()
    f.close()
    itmp = search_string('STRMIX        FCM       QBOUND', tmptxt)
    tmpval = tmptxt[itmp+1].split()
    strmix = float(tmpval[0])
    fcm = float(tmpval[1])
    qbound = float(tmpval[2])
    tmpval = tmptxt[itmp+4].split()
    brymix = float(tmpval[0])
    itmp = search_string('IMIX    IGF    ICC', tmptxt)
    imix = int(tmptxt[itmp+1].split()[0])
    idtbry = int(tmptxt[itmp+4].split()[0])
    
    mixinfo = [imix, strmix, qbound, fcm, idtbry, brymix]
        
    return niter, nitermax, converged, nmax_reached, mixinfo


def get_kmeshinfo(outfile_0init, outfile_000):
    """
    Extract kmesh info from output.0.txt and output.000.txt
    """
    # first get info from output.0.txt
    f = open(outfile_0init)
    tmptxt = f.readlines()
    f.close()
    nkmesh = []
    itmp = search_string('number of different k-meshes', tmptxt)
    nkmesh.append( int(tmptxt[itmp].split(':')[1].split()[0]) )
    itmp = search_string('k-mesh NofKs', tmptxt)
    nofks, nkx, nky, nkz = [],[],[],[]
    for ik in range(nkmesh[0]):
        tmpval = tmptxt[itmp+2+ik].split()
        nofks.append(int(tmpval[1]))
        nkx.append(int(tmpval[2]))
        nky.append(int(tmpval[3]))
        nkz.append(int(tmpval[4]))
    
    tmpdict = {'number_of_kpts':nofks, 'n_kx':nkx, 'n_ky':nky, 'n_kz':nkz}
    nkmesh.append(tmpdict)
    
    #next get kmesh_ie from output.000.txt
    f = open(outfile_000)
    tmptxt = f.readlines()
    f.close()
    kmesh_ie = []
    itmp = 0
    while itmp>=0:
        itmp = search_string('KMESH =', tmptxt)
        print itmp
        if itmp>=0:
            tmpval = int(tmptxt.pop(itmp).split()[-1])
            kmesh_ie.append(tmpval)
    
    return nkmesh, kmesh_ie
    

def get_symmetries(outfile_0init):
    f = open(outfile_0init)
    tmptxt = f.readlines()
    f.close()
    itmp = search_string('symmetries found for this lattice:', tmptxt)
    nsym = int(tmptxt[itmp].split(':')[1].split()[0])
    itmp = search_string('<SYMTAUMAT>', tmptxt)
    tmpdict = {}
    for isym in range(nsym):
        tmpval = tmptxt[itmp+5+isym].split()
        desc = tmpval[1]
        inversion = int(tmpval[2])
        euler = [float(tmpval[3]), float(tmpval[4]), float(tmpval[5])]
        unitary = int(tmpval[6].replace('T', '1').replace('F', '0'))
        tmpdict[desc] = {'has_inversion':inversion, 'is_unitary':unitary, 'euler_angles':euler} 
    desc = tmpdict
    return nsym, desc
    

def get_ewald(outfile_0init):
    f = open(outfile_0init)
    tmptxt = f.readlines()
    f.close()
    itmp = search_string('< LATTICE3D >', tmptxt)
    tmpval = tmptxt[itmp+7].split()[2:]
    rsum = float(tmpval[2]), int(tmpval[0]), int(tmpval[1])
    tmpval = tmptxt[itmp+8].split()[2:]
    gsum = float(tmpval[2]), int(tmpval[0]), int(tmpval[1])
    itmp = search_string('setting bulk Madelung coefficients', tmptxt)
    if itmp>=0:
        info = '3D'
    else:
        info = '2D'
    return rsum, gsum, info

def get_spinmom_per_atom(outfile):
    return result, vec, angles

def get_orbmom(outfile):
    return result, vec, angles


def parse_kkr_outputfile(out_dict, outfile, outfile_0init, outfile_000, timing_file, potfile_out):
    """
    Parser method for the kkr outfile. It returns a dictionary with results
    """
    # scaling factors
    Ry2eV = get_Ry2eV()
    
    # collection of parsing error messages
    msg_list = []
        
    try:
        code_version, compile_options, serial_number = get_version_info(outfile)
        tmp_dict = {}
        tmp_dict['code_version'] = code_version
        tmp_dict['compile_options'] = compile_options
        tmp_dict['calculation_serial_number'] = serial_number
        out_dict['code_info_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: Version Info"
        msg_list.append(msg)
    
    tmp_dict = {} # used to group convergence info (rms, rms per atom, charge neutrality)
    # also initialize convegence_group where all info stored for all iterations is kept
    out_dict['convergence_group'] = tmp_dict
    try:
        result, result_atoms_last = get_rms(outfile, outfile_000)
        tmp_dict['rms'] = result[-1]
        tmp_dict['rms_all_iterations'] = result
        tmp_dict['rms_per_atom'] = result_atoms_last
        tmp_dict['rms_unit'] = 'unitless'
        out_dict['convergence_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: rms-error"
        msg_list.append(msg)
    
    try:
        result = get_neutr(outfile)
        tmp_dict['charge_neutrality'] = result[-1]
        out_dict['convergence_group']['charge_neutrality_all_iterations'] = result
        tmp_dict['charge_neutrality_unit'] = 'electrons'
        out_dict['convergence_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: charge neutrality"
        msg_list.append(msg)
       
    tmp_dict = {} # used to group magnetism info (spin and orbital moments)
    try:
        result = get_magtot(outfile)
        if len(result)>0:
            tmp_dict['total_spin_moment'] = result[-1]
            out_dict['convergence_group']['total_spin_moment_all_iterations'] = result
            tmp_dict['total_spin_moment_unit'] = 'mu_Bohr'
            out_dict['magnetism_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: total magnetic moment"
        msg_list.append(msg)
    
    try:
        result, vec, angles = get_spinmom_per_atom(outfile)
        if len(result)>0:
            tmp_dict['spin_moment_per_atom'] = result[-1,:]
            tmp_dict['spin_moment_vector_per_atom'] = vec[-1,:]
            tmp_dict['spin_moment_angles_per_atom'] = angles[-1,:]
            out_dict['convergence_group']['spin_moment_per_atom_all_iterations'] = result[:,:]
            out_dict['convergence_group']['spin_moment_angles_per_atom_all_iterations'] = angles[:,:]
            tmp_dict['spin_moment_unit'] = 'mu_Bohr'
            tmp_dict['spin_moment_angles_per_atom_unit'] = 'degree'
            out_dict['magnetism_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: spin moment per atom"
        msg_list.append(msg)
    
    try:
        result, vec, angles = get_orbmom(outfile)
        if len(result)>0:
            tmp_dict['total_orbital_moment'] = sum(result[-1,:])
            tmp_dict['orbital_moment_per_atom'] = result[-1,:]
            tmp_dict['orbital_moment_vector_per_atom'] = vec[-1,:]
            tmp_dict['orbital_moment_angles_per_atom'] = angles[-1,:]
            out_dict['convergence_group']['orbital_moment_per_atom_all_iterations'] = result[:,:]
            out_dict['convergence_group']['orbital_moment_angles_per_atom_all_iterations'] = angles[:,:]
            tmp_dict['orbital_moment_unit'] = 'mu_Bohr'
            tmp_dict['orbital_moment_angles_per_atom_unit'] = 'degree'
            out_dict['magnetism_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: orbital moment"
        msg_list.append(msg)

    try:
        result = get_EF(outfile)
        out_dict['fermi_energy'] = result[-1]
        out_dict['convergence_group']['fermi_energy_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: EF"
        msg_list.append(msg)

    try:
        result = get_DOS_EF(outfile)
        out_dict['dos_at_fermi_energy'] = result[-1]
        out_dict['convergence_group']['dos_at_fermi_energy_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: DOS@EF"
        msg_list.append(msg)

    try:
        result = get_Etot(outfile)
        out_dict['energy'] = result[-1]*Ry2eV
        out_dict['energy_unit'] = 'eV'
        out_dict['total_energy_Ry'] = result[-1]
        out_dict['total_energy_Ry_unit'] = 'Rydberg'
        out_dict['convergence_group']['total_energy_Ry_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: total energy"
        msg_list.append(msg)

    try:
        result = find_warnings(outfile)
        tmp_dict = {}
        tmp_dict['number_of_warnings'] = len(result)
        tmp_dict['warnings_list'] = result
        out_dict['warnings_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: search for warnings"
        msg_list.append(msg)

    try:
        result = extract_timings(timing_file)
        out_dict['timings_group'] = result
        out_dict['timings_unit'] = 'seconds'
    except:
        msg = "Error parsing output of KKR: timings"
        msg_list.append(msg)
    
    try:
        result = get_single_particle_energies(outfile_000)
        out_dict['single_particle_energies'] = result*Ry2eV
        out_dict['single_particle_energies_unit'] = 'eV'
    except:
        msg = "Error parsing output of KKR: single particle energies"
        msg_list.append(msg)
    
    try:
        result_WS, result_tot, result_C = get_charges_per_atom(outfile_000)
        niter = len(out_dict['convergence_group']['rms_all_iterations'])
        natyp = len(result_tot)/niter
        out_dict['nuclear_charge_per_atom'] = result_tot[-natyp:]
        out_dict['charge_core_states_per_atom'] = result_C[-natyp:]
        out_dict['charge_valence_states_per_atom'] = result_WS[-natyp:]-result_C[-natyp:]
        out_dict['nuclear_charge_per_atom_unit'] = 'electron charge'
        out_dict['charge_core_states_per_atom_unit'] = 'electron charge'
        out_dict['charge_valence_states_per_atom_unit'] = 'electron charge'
    except:
        msg = "Error parsing output of KKR: charges"
        msg_list.append(msg)
    
    try:
        emin, tempr, Nepts, Npol, N1, N2, N3 = get_econt_info(outfile_0init)
        tmp_dict = {}
        tmp_dict['emin'] = emin
        tmp_dict['emin_unit'] = 'Rydberg'
        tmp_dict['number_of_energy_points'] = Nepts
        tmp_dict['temperature'] = tempr
        tmp_dict['temperature_unit'] = 'Kelvin'
        tmp_dict['npol'] = Npol
        tmp_dict['n1'] = N1
        tmp_dict['n2'] = N2
        tmp_dict['n3'] = N3
        out_dict['energy_contour_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: energy contour"
        msg_list.append(msg)
    
    try:
        ncore, emax, lmax, descr_max = get_core_states(potfile_out)
        tmp_dict = {}
        tmp_dict['number_of_core_states_per_atom'] = ncore
        tmp_dict['energy_highest_lying_core_state_per_atom'] = emax
        tmp_dict['energy_highest_lying_core_state_per_atom_unit'] = 'Rydberg'
        tmp_dict['descr_highest_lying_core_state_per_atom'] = descr_max
        out_dict['core_states_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: core_states"
        msg_list.append(msg)
        
    try:
        alat, twopioveralat = get_alatinfo(outfile_0init)
        out_dict['alat_internal'] = alat
        out_dict['two_pi_over_alat_internal'] = twopioveralat
        out_dict['alat_internal_unit'] = 'a_Bohr'
        out_dict['two_pi_over_alat_internal_unit'] = '1/a_Bohr'
    except:
        msg = "Error parsing output of KKR: alat, 2*pi/alat"
        msg_list.append(msg)
        
    try:
        niter, nitermax, converged, nmax_reached, mixinfo = get_scfinfo(outfile_0init, outfile_000, outfile)
        out_dict['convergence_group']['number_of_iterations'] = niter
        out_dict['convergence_group']['number_of_iterations_max'] = nitermax
        out_dict['convergence_group']['calculation_converged'] = converged
        out_dict['convergence_group']['nsteps_exhausted'] = nmax_reached
        out_dict['convergence_group']['imix'] = mixinfo[0]
        out_dict['convergence_group']['strmix'] = mixinfo[1]
        out_dict['convergence_group']['qbound'] = mixinfo[2]
        out_dict['convergence_group']['fcm'] = mixinfo[3]
        out_dict['convergence_group']['idtbry'] = mixinfo[4]
        out_dict['convergence_group']['brymix'] = mixinfo[5]
    except:
        msg = "Error parsing output of KKR: scfinfo"
        msg_list.append(msg)
    
    try:
        nkmesh, kmesh_ie = get_kmeshinfo(outfile_0init, outfile_000)
        tmp_dict = {}
        tmp_dict['number_different_kmeshes'] = nkmesh[0]
        tmp_dict['number_kpoints_per_kmesh'] = nkmesh[1]
        tmp_dict['kmesh_energypoint'] = kmesh_ie
        out_dict['kmesh_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: kmesh"
        msg_list.append(msg)
    
    try:
        nsym, desc = get_symmetries(outfile_0init)
        tmp_dict = {}
        tmp_dict['number_of_symmetries'] = nsym
        tmp_dict['symmetry_description'] = desc
        out_dict['symmetries_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: symmetries"
        msg_list.append(msg)
        
    #TODO Ewald summation
    try:
        rsum, gsum, info = get_ewald(outfile_0init)
        tmp_dict = {}
        tmp_dict['ewald_summation_mode'] = info
        tmp_dict['rsum_cutoff'] = rsum[0]
        tmp_dict['rsum_number_of_vectors'] = rsum[1]
        tmp_dict['rsum_number_of_shells'] = rsum[2]
        tmp_dict['rsum_cutoff_unit'] = 'a_Bohr'
        tmp_dict['gsum_cutoff'] = gsum[0]
        tmp_dict['gsum_number_of_vectors'] = gsum[1]
        tmp_dict['gsum_number_of_shells'] = gsum[2]
        tmp_dict['gsum_cutoff_unit'] = '1/a_Bohr'
        out_dict['ewald_sum_group'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: ewald summation for madelung poterntial"
        msg_list.append(msg)
        
        
    #convert arrays to lists
    from numpy import ndarray
    for key in out_dict.keys():
        if type(out_dict[key])==ndarray:
            out_dict[key] = list(out_dict[key])
        elif type(out_dict[key])==dict:
            for subkey in out_dict[key].keys():
                if type(out_dict[key][subkey])==ndarray:
                    out_dict[key][subkey] = list(out_dict[key][subkey])
                    
                    
    # return output with error messages if there are any
    if len(msg_list)>0:
        return False, msg_list, out_dict
    else:
        return True, [], out_dict
      

#"""
# for testing
if __name__=='__main__':
    from pprint import pprint
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()
    
    #from aiida_kkr.parsers.kkr import parse_kkr_outputfile
    path0 = '/Users/ruess/sourcecodes/aiida/repositories/scratch_local_machine/73/a4/2752-db84-4073-a3fc-489863d74720/'
    outfile = path0+'out_kkr'
    outfile_0init = path0+'output.0.txt'
    outfile_000 = path0+'output.000.txt'
    timing_file = path0+'out_timing.000.txt'
    potfile_out = path0+'out_potential'
    out_dict = {}
    success, msg_list, out_dict = parse_kkr_outputfile(out_dict, outfile, outfile_0init, outfile_000, timing_file, potfile_out)
    
    pprint(out_dict)
    if not success:
        print 'Number of parser_errors', len(msg_list)
        for msg in msg_list:
            print msg
#"""