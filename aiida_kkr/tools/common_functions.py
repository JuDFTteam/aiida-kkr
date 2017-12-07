#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 13:25:35 2017

@author: ruess
"""


#helper functions used in calculation, parser etc.
def get_alat_from_bravais(bravais, is3D=True):
    from numpy import sqrt, sum
    bravais_tmp = bravais
    if not is3D:
        #take only in-plane lattice to find maximum as alat
        bravais_tmp = bravais[:2,:2]
    return sqrt(sum(bravais_tmp**2, axis=1)).max()
    
def get_Ang2aBohr():
    return 1.8897261254578281
    
def get_aBohr2Ang():
    return 1/get_Ang2aBohr()

def get_Ry2eV():
    return 13.605693009
    
def search_string(searchkey, txt):
    iline = 0
    for line in txt:
        if searchkey in line:
            return iline
        iline+=1
    return -1


def get_version_info(outfile):
    f = open(outfile)
    tmptxt = f.readlines()
    f.close()
    itmp = search_string('Code version:', tmptxt)
    code_version = tmptxt.pop(itmp).split(':')[1].strip()
    itmp = search_string('Compile options:', tmptxt)
    compile_options = tmptxt.pop(itmp).split(':')[1].strip()
    itmp = search_string('serial number for files:', tmptxt)
    serial_number = tmptxt.pop(itmp).split(':')[1].strip()
    return code_version, compile_options, serial_number


def get_corestates_from_potential(potfile='potential'):
    """Read core states from potential file"""
    from scipy import zeros
    txt = open(potfile).readlines()

    #get start of each potential part
    istarts = [iline for iline in range(len(txt)) if 'POTENTIAL' in txt[iline]]

    n_core_states = [] #number of core states per potential
    e_core_states = [] #energies of core states
    l_core_states = [] #angular momentum index, i.e. 0=s, 1=p etc...
    for ipot in range(len(istarts)):
        line = txt[istarts[ipot]+6]
        n = int(line.split()[0])
        n_core_states.append(n)
        elevels = zeros(n) #temp array for energies
        langmom = zeros(n, dtype=int) #temp array for angular momentum index
        for icore in range(n):
            line = txt[istarts[ipot]+7+icore].split()
            langmom[icore] = int(line[0])
            elevels[icore] = float(line[1].replace('D', 'E'))
        e_core_states.append(elevels)
        l_core_states.append(langmom)

    return n_core_states, e_core_states, l_core_states


def get_highest_core_state(nstates, energies, lmoments):
    """Find highest lying core state from list of core states, needed to find and check energy contour"""
    idx = energies.argmax()
    lval = lmoments[idx]
    nquant = sum(lmoments == lval) + lval
    level_descr = '%i%s'%(nquant, 'spdfgh'[lval])

    return lval, energies[idx], level_descr


def generate_inputcard_from_structure(parameters, structure, input_filename, parent_calc=None, shapes=None):
    """
    Takes information from parameter and structure data and writes input file 'input_filename'
    
    :param parameters: input parameters node containing KKR-related input parameter
    :param structure: input structure node containing lattice information
    :param input_filename: input filename, typically called 'inputcard'
    
    
    optional arguments
    :param parent_calc: input parent calculation node used to determine if EMIN 
                        parameter is automatically overwritten (from voronoi output)
                        or not
    :param shapes: input shapes array (set automatically by 
                   aiida_kkr.calculations.Kkrcaluation and shall not be overwritten)
    
    
    :note: assumes valid structure and parameters, i.e. for 2D case all necessary 
           information has to be given. This is checked with function 
           'check_2D_input' called in aiida_kkr.calculations.Kkrcaluation
    """
    
    from aiida.common.constants import elements as PeriodicTableElements
    from numpy import array
    from aiida_kkr.calculations.voro import VoronoiCalculation
    from aiida_kkr.tools.kkr_params import kkrparams
    
    #list of globally used constants
    a_to_bohr = get_Ang2aBohr()


    # Get the connection between coordination number and element symbol
    # maybe do in a differnt way
    
    _atomic_numbers = {data['symbol']: num for num,
                    data in PeriodicTableElements.iteritems()}
    
    # KKR wants units in bohr
    bravais = array(structure.cell)*a_to_bohr
    alat = get_alat_from_bravais(bravais, is3D=structure.pbc[2])
    bravais = bravais/alat
    
    sites = structure.sites
    naez = len(sites)
    positions = []
    charges = []
    weights = [] # for CPA
    isitelist = [] # counter sites array for CPA
    isite = 0
    for site in sites:
        pos = site.position 
        #TODO maybe convert to rel pos and make sure that type is right for script (array or tuple)
        abspos = array(pos)*a_to_bohr/alat # also in units of alat
        positions.append(abspos)
        isite += 1
        sitekind = structure.get_kind(site.kind_name)
        for ikind in range(len(sitekind.symbols)):
            site_symbol = sitekind.symbols[ikind]
            if not sitekind.has_vacancies():
                charges.append(_atomic_numbers[site_symbol])
            else:
                charges.append(0.0)
            #TODO deal with VCA case
            if sitekind.is_alloy():
                weights.append(sitekind.weights[ikind])
            else:
                weights.append(1.)
            
            isitelist.append(isite)
    
    weights = array(weights)
    isitelist = array(isitelist)
    charges = array(charges)
    positions = array(positions)
        

    ######################################
    # Prepare keywords for kkr from input structure
    
    # get parameter dictionary
    input_dict = parameters.get_dict()
    
    # empty kkrparams instance (contains formatting info etc.)
    params = kkrparams()
    
    # for KKR calculation set EMIN automatically from parent_calc (ausways in res.emin of voronoi and kkr)
    if ('EMIN' not in input_dict.keys() or input_dict['EMIN'] is None) and parent_calc is not None:
        print('Overwriting EMIN with value from parent calculation')
        emin = parent_calc.res.emin
        print('Setting emin:',emin, 'is emin None?',emin is None)
        params.set_value('EMIN', emin)
        
    # overwrite keywords with input parameter
    for key in input_dict.keys():
        params.set_value(key, input_dict[key])

    # Write input to file (the parameters that are set here are not allowed to be modfied externally)
    params.set_multiple_values(BRAVAIS=bravais, ALATBASIS=alat, NAEZ=naez, 
                               ZATOM=charges, RBASIS=positions, CARTESIAN=True)
    # for CPA case:
    if len(weights)>naez:
        natyp = len(weights)
        params.set_value('NATYP', natyp)
        params.set_value('<CPA-CONC>', weights)
        params.set_value('<SITE>', isitelist)
        
    # write shapes (extracted from voronoi parent automatically in kkr calculation plugin)
    if shapes is not None:
        params.set_value('<SHAPE>', shapes)
    
    # write inputfile
    params.fill_keywords_to_inputfile(output=input_filename)
    
    
def check_2Dinput(structure, parameters):
    """
    Check if structure and parameter data are complete and matching.
    
    :param input: structure, needs to be a valid aiida StructureData node
    :param input: parameters, needs to be valid aiida ParameterData node
    
    returns (False, errormessage) if an inconsistency has been found, otherwise return (True, '2D consistency check complete')
    """
    # default is bulk, get 2D info from structure.pbc info (periodic boundary contitions)
    is2D = False
    if not all(structure.pbc):
        # check periodicity, assumes finite size in z-direction
        if structure.pbc != (True, True, False):
            return (False, "Structure.pbc is neither (True, True, True) for bulk nor (True, True, False) for surface calculation!")
        is2D = True
    
    # check for necessary info in 2D case
    inp_dict = parameters.get_dict()
    set_keys = [i for i in inp_dict.keys() if inp_dict[i] is not None]
    has2Dinfo = True
    for icheck in ['INTERFACE', '<NRBASIS>', '<RBLEFT>', '<RBRIGHT>', 'ZPERIODL', 'ZPERIODR', '<NLBASIS>']:
        if icheck not in set_keys:
            has2Dinfo = False
    if has2Dinfo and not inp_dict['INTERFACE'] and is2D:
        return (False, "'INTERFACE' parameter set to False but structure is 2D")
        
    if has2Dinfo!=is2D:
        return (False, "2D info given in parameters but structure is 3D")
    
    # if everything is ok:
    return (True, "2D consistency check complete")
    


"""
# testing ...

def submission(calc, submit_test=False):
    import os
    if submit_test:
        subfolder, script_filename = calc.submit_test()
        print "Test_submit for calculation (uuid='{}')".format(
            calc.uuid)
        print "Submit file in {}".format(os.path.join(
            os.path.relpath(subfolder.abspath),
            script_filename
        ))
        return -1
    else:
        calc.store_all()
        print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
            calc.uuid, calc.dbnode.pk)
        calc.submit()
        print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
            calc.uuid, calc.dbnode.pk)
        return calc.dbnode.pk
    
    
if __name__=='__main__':
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()
    from aiida.orm import Code, load_node
    from aiida_kkr.calculations.kkr import KkrCalculation
    resubmit = True
    # Load aiida structure  and initial parameter node node 
    Cu = load_node(577)
    ParaNode = load_node(583)
    
    # step0: voronoi calculation
    VoroCalc = load_node(596)
    ResultedParameterVoro = VoroCalc
    ParaNode4 = load_node(645)
    print ParaNode4.get_dict()
    remote_voro = ResultedParameterVoro.get_outputs_dict()['remote_folder']
    
    # setup new kkr calculation
    KKRcode = Code.get_from_string('KKRcode@my_mac')
    KKRcalc = KkrCalculation()
    KKRcalc.label = 'KKR calculation'
    KKRcalc.set_withmpi(False)
    KKRcalc.set_resources({"num_machines" : 1})
    KKRcalc.set_max_wallclock_seconds(300)
    KKRcalc.set_computer('my_mac')
    KKRcalc.use_code(KKRcode)
    KKRcalc.use_parameters(ParaNode4)
    KKRcalc.use_parent_folder(remote_voro)
    
    from subprocess import call
    call('rm ../*/*pyc', shell=True)
    
    if resubmit:
        sid = submission(KKRcalc, submit_test=False)
    else:
        sid = 211
       
    KKRcalc = load_node(sid)
    
    print KKRcalc.has_finished_ok()
    print KKRcalc.has_finished()
    
    #!verdi calculation logshow 195
    
"""