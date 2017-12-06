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
        
        
        out_dict = {'ParserVersion': self._ParserVersion}
        success, msg, out_dict = parse_kkr_outputfile(out_dict, outfile, outfile_0init, outfile_000, timing_file, potfile_out)

        return success, self._get_nodelist(out_dict)
    
    # here follow the parser functions:
    
    def _get_nodelist(self, out_dict):
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        return node_list


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
            res.append(tmpval)
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
    return [[search_keys[i], res[i]] for i in range(len(res))]


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
    return ncore, emax, lmax, descr_max


def parse_kkr_outputfile(out_dict, outfile, outfile_0init, outfile_000, timing_file, potfile_out):
    """
    Parser method for the kkr outfile. It returns a dictionary with results
    """
    # TODO This still needs to be overworked.
    # If we want to store the arrays we have to store them as array data
    # If only certain results are imported they have to be distilled and stored in the output node
    
    # scaling factors
    Ry2eV = get_Ry2eV()
        
    try:
        code_version, compile_options, serial_number = get_version_info(outfile)
        out_dict['Code_version'] = code_version
        out_dict['Compile_options'] = compile_options
        out_dict['Calculation_serial_number'] = serial_number
    except:
        msg = "Error parsing output of KKR: Version Info"
        return False, msg, out_dict
    
    tmp_dict = {} # used to group convergence info (rms, rms per atom, charge neutrality)
    try:
        result, result_atoms_last = get_rms(outfile, outfile_000)
        tmp_dict['rms'] = result[-1]
        tmp_dict['rms_all_iterations'] = result
        tmp_dict['rms_per_atom'] = result_atoms_last
        out_dict['convergence'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: rms-error"
        return False, msg, out_dict
    
    try:
        result = get_neutr(outfile)
        tmp_dict['charge_neutrality'] = result[-1]
        tmp_dict['charge_neutrality_all_iterations'] = result
        out_dict['convergence'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: charge neutrality"
        return False, msg, out_dict
    
    try:
        result = get_magtot(outfile)
        if len(result)>0:
            out_dict['total_magnetic_moment'] = result[-1]
            out_dict['total_magnetic_moment_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: total magnetic moment"
        return False, msg, out_dict

    try:
        result = get_EF(outfile)
        out_dict['EF'] = result[-1]
        out_dict['EF_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: EF"
        return False, msg, out_dict

    try:
        result = get_DOS_EF(outfile)
        out_dict['DOS_EF'] = result[-1]
        out_dict['DOS_EF_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: DOS@EF"
        return False, msg, out_dict

    try:
        result = get_Etot(outfile)
        out_dict['total_energy'] = result[-1]*Ry2eV
        out_dict['units_total_energy'] = 'eV'
        out_dict['total_energy_Ry'] = result[-1]
        out_dict['total_energy_Ry_all_iterations'] = result
    except:
        msg = "Error parsing output of KKR: total energy"
        return False, msg, out_dict

    try:
        result = find_warnings(outfile)
        out_dict['number_of_warnings'] = len(result)
        out_dict['warnings_list'] = result
    except:
        msg = "Error parsing output of KKR: search for warnings"
        return False, msg, out_dict

    try:
        result = extract_timings(timing_file)
        out_dict['timings'] = result
        out_dict['units_timings'] = 'seconds'
    except:
        msg = "Error parsing output of KKR: timings"
        return False, msg, out_dict
    
    try:
        result = get_single_particle_energies(outfile_000)
        out_dict['single_particle_energies'] = result*Ry2eV
        out_dict['units_single_particle_energies'] = 'eV'
    except:
        msg = "Error parsing output of KKR: single particle energies"
        return False, msg, out_dict
    
    try:
        result_WS, result_tot, result_C = get_charges_per_atom(outfile_000)
        niter = len(out_dict['rms_all_iterations'])
        natyp = len(result_tot)/niter
        out_dict['nuclear_charge_per_atom'] = result_tot[-natyp:]
        out_dict['charge_core_states_per_atom'] = result_C[-natyp:]
        out_dict['charge_valence_states_per_atom'] = result_WS[-natyp:]-result_C[-natyp:]
    except:
        msg = "Error parsing output of KKR: single particle energies"
        return False, msg, out_dict
    
    try:
        emin, tempr, Nepts, Npol, N1, N2, N3 = get_econt_info(outfile_0init)
        tmp_dict = {}
        tmp_dict['EMIN'] = emin
        tmp_dict['Nepts'] = Nepts
        tmp_dict['Temperature'] = tempr
        tmp_dict['Npol'] = Npol
        tmp_dict['N1'] = N1
        tmp_dict['N2'] = N2
        tmp_dict['N3'] = N3
        out_dict['energy_contour'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: energy contour"
        return False, msg, out_dict
    
    try:
        ncore, emax, lmax, descr_max = get_core_states(potfile_out)
        tmp_dict = {}
        tmp_dict['number_of_core_states_per_atom'] = ncore
        tmp_dict['energy_highest_lying_core_state_per_atom'] = emax
        tmp_dict['descr_highest_lying_core_state_per_atom'] = descr_max
        out_dict['core_states'] = tmp_dict
    except:
        msg = "Error parsing output of KKR: core_states"
        return False, msg, out_dict

    
    return True, "Completed parsing of KKR output successfully.", out_dict
      

#"""
# for testing
if __name__=='__main__':
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
    success, msg, out_dict = parse_kkr_outputfile(out_dict, outfile, outfile_0init, outfile_000, timing_file, potfile_out)
    if not success:
        print 'Error-msg:', msg
    print(out_dict)
#"""