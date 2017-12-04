#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 13:25:35 2017

@author: ruess
"""


#helper functions used in calculation, parser etc.
def get_alat_from_bravais(bravais):
    from numpy import sqrt, sum
    #return bravais.max()
    return sqrt(sum(bravais**2, axis=1)).max()
    
def get_Ang2aBohr():
    return 1.8897261254578281
    
def get_aBohr2Ang():
    return 1/get_Ang2aBohr()
    
def search_string(searchkey, txt):
    iline = 0
    for line in txt:
        if searchkey in line:
            return iline
        iline+=1
    return -1

def generate_inputcard_from_structure(parameters, structure, input_filename, parent_calc=None):
    """
    Takes information from parameter and structure data and writes input file 'input_filename'
    
    :param input: parameters, parameter node containing KKR-related input parameter
    :param input: structure, structure node containing lattice information
    :param input: input_filename, filename created, typically called 'inputcard'
    
    optional argument
    :param input: parent_calc, parent calculation node used to determine if EMIN parameter is automatically overwritten (from voronoi output) or not
    """
    
    from aiida_kkr.tools.common_functions import get_alat_from_bravais, get_Ang2aBohr
    from aiida.common.constants import elements as PeriodicTableElements
    from numpy import array
    from aiida_kkr.calculations.voro import VoronoiCalculation
    from aiida_kkr.tools.kkr_params import kkrparams
    
    
    # check if calculation is voronoi calculatino or not
    if isinstance(parent_calc, VoronoiCalculation):
        parentisvoronoi = True
    else:
        parentisvoronoi = False
    
    #list of globally used constants
    a_to_bohr = get_Ang2aBohr()


    # Get the connection between coordination number and element symbol
    # maybe do in a differnt way
    
    _atomic_numbers = {data['symbol']: num for num,
                    data in PeriodicTableElements.iteritems()}
    
    # KKR wants units in bohr
    bravais = array(structure.cell)*a_to_bohr
    alat = get_alat_from_bravais(bravais)
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
    #TODO default is bulk, get 2D from structure.pbc info (periodic boundary contitions)

    ######################################
    # Prepare keywords for kkr from input structure
    
    # get parameter dictionary
    input_dict = parameters.get_dict()
    
    # empty kkrparams instance (contains formatting info etc.)
    if parentisvoronoi:
        params = kkrparams(params_type='voronoi')
    else:
        params = kkrparams()
    
    # in case of starting from Voronoi set EMIN automatically
    if parentisvoronoi:
        print('Overwriting EMIN with value from voronoi output')
        emin = parent_calc.res.EMIN
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
        
    # write inputfile
    params.fill_keywords_to_inputfile(output=input_filename)
    
    
# tests etc.

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
    
