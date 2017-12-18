#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Here commonly used functions that do not need aiida-stuff (i.e. can be tested 
without a database) are collected.
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

def angles_to_vec(magnitude, theta, phi):
    """
    convert (magnitude, theta, phi) to (x,y,z)
    
    theta/phi need to be in radians!
    
    Input can be single number, list of numpy.ndarray data
    Returns x,y,z vector 
    """
    from numpy import ndarray, array, cos, sin
    # correct data type if necessary
    if type(magnitude) != ndarray:
        magnitude = array(magnitude)
    if type(theta) != ndarray:
        theta = array(theta)
    if type(phi) != ndarray:
        phi = array(phi)
    vec = []
    for ivec in range(len(magnitude)):
        r_inplane = magnitude[ivec]*sin(theta[ivec])
        x = r_inplane*cos(phi[ivec])
        y = r_inplane*sin(phi[ivec])
        z = cos(theta[ivec])*magnitude[ivec]
        vec.append([x,y,z])
    return array(vec)

def vec_to_angles(vec):
    """
    converts vector (x,y,z) to (magnitude, theta, phi)
    """
    from numpy import array, arctan2, sqrt, shape
    magnitude, theta, phi = [], [], []
    if len(vec)<=3 and len(shape(vec))<2:
        vec = [vec]
    for ivec in range(len(vec)):
        phi.append(arctan2(vec[ivec, 1], vec[ivec, 0]))
        r_inplane = sqrt(vec[ivec, 0]**2+vec[ivec, 1]**2)
        theta.append(arctan2(r_inplane, vec[ivec, 2]))
        magnitude.append(sqrt(r_inplane**2+vec[ivec, 2]**2))
    return array(magnitude), array(theta), array(phi)
    


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


def generate_inputcard_from_structure(parameters, structure, input_filename, parent_calc=None, shapes=None, isvoronoi=False):
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
    if not isvoronoi:
        params = kkrparams()
    else:
        params = kkrparams(params_type='voronoi')
    
    # for KKR calculation set EMIN automatically from parent_calc (ausways in res.emin of voronoi and kkr)
    if ('EMIN' not in input_dict.keys() or input_dict['EMIN'] is None) and parent_calc is not None:
        print('Overwriting EMIN with value from parent calculation')
        emin = parent_calc.res.emin
        print('Setting emin:',emin, 'is emin None?',emin is None)
        params.set_value('EMIN', emin)
        
    # overwrite keywords with input parameter
    for key in input_dict.keys():
        params.set_value(key, input_dict[key], silent=True)

    # Write input to file (the parameters that are set here are not allowed to be modfied externally)
    params.set_multiple_values(BRAVAIS=bravais, ALATBASIS=alat, NAEZ=naez, 
                               ZATOM=charges, RBASIS=positions, CARTESIAN=True)
    # for CPA case:
    if len(weights)>naez:
        natyp = len(weights)
        params.set_value('NATYP', natyp)
        params.set_value('<CPA-CONC>', weights)
        params.set_value('<SITE>', isitelist)
    else:
        natyp = naez
        
    # write shapes (extracted from voronoi parent automatically in kkr calculation plugin)
    if shapes is not None:
        params.set_value('<SHAPE>', shapes)
    
    # write inputfile
    params.fill_keywords_to_inputfile(output=input_filename)
    
    nspin = params.get_value('NSPIN')
    
    newsosol = False
    if 'NEWSOSOL' in params.get_value('RUNOPT'):
        newsosol = True
    
    return natyp, nspin, newsosol
    
    
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


def interpolate_dos(dospath, return_original=False, ):
    """
    interpolation function copied from complexdos3 fortran code
    
    Principle of DOS here: Two-point contour integration
    for DOS in the middle of the two points. The input DOS
    and energy must be complex. Parameter deltae should be
    of the order of magnitude of eim. 
    
        
          <-2*deltae->   _
               /\        |     DOS=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
              /  \       |
            (1)  (2)   2*i*eim=2*i*pi*Kb*Tk
            /      \     |
           /        \    |
    ------------------------ (Real E axis)
    
    :param input: dospath, path where 'complex.dos' file can be found
    
    :returns: E_Fermi, numpy array of interpolated dos 
    
    :note: output units are in Ry!
    """
    from numpy import array, real, imag
    
    f = open(dospath+'/complex.dos', 'r')
    text = f.readline() # dummy readin of header, may be replaced later
    npot = int(f.readline().split()[0])
    iemax = int(f.readline().split()[0])
    lmax = int(f.readline().split()[0])
    
    dosnew_all_atoms = []
    dos_all_atoms = []
    
    for i1 in range(npot):
        #print('Reading potential',i1)
        # Read header (not used)
        for iheader in range(3):
            text = f.readline()
        
        # extract EF
        ef = float(f.readline().split()[7])
        
        # some more dummy lines
        for iheader in range(5,9+1):
            text = f.readline()
        
        # now header is done. start reading DOS
        # Read dos: (total dos stored at DOS(LMAX+1,IE))
        dos_l_cmplx = []
        for ie in range(iemax):
            tmpline = f.readline().replace('(','').replace(')','').replace(',','').split()
            ez = float(tmpline[0])+1j*float(tmpline[1])
            dostmp_complex = [[tmpline[len(tmpline)-2], tmpline[len(tmpline)-1]]]
            dostmp_complex += [[tmpline[iline], tmpline[iline+1]] for iline in range(2,len(tmpline)-2,2)]
            dostmp = [ez]+[float(ds[0])+1j*float(ds[1]) for ds in dostmp_complex]
            dos_l_cmplx.append(dostmp)
        dos_l_cmplx = array(dos_l_cmplx)
        dos_l = imag(dos_l_cmplx.copy())
        dos_l[:,0] = real(dos_l_cmplx.copy()[:,0])
        dos_all_atoms.append(dos_l)
        
        # Compute and write out corrected dos at new (middle) energy points:
        dosnew = []
        ez = dos_l_cmplx[:,0]
        for ie in range(1, iemax-1):
            deltae = real(ez[ie+1] - ez[ie])
            eim = imag(ez[ie])
            enew = real(ez[ie]) # Real quantity
        
            tmpdos = [enew]
            for ll in range(1,lmax+3):
                t = (dos_l_cmplx[ie-1, ll]-dos_l_cmplx[ie+1, ll])*0.5*(0.0+eim*1j)/deltae
                #print ie+1, ll,  dos_l_cmplx[ie, ll], deltae, eim, t, shape(dos_l_cmplx[ie]), lmax
                #tmpdos.append(dos_l_cmplx[ie, ll] + 0.5*(dos_l_cmplx[ie-1, ll]-dos_l_cmplx[ie+1, ll])*(0.+1j*eim)/deltae)
                tmpdos.append(dos_l_cmplx[ie, ll]+t)
            tmpdos = array(tmpdos)
            # build imaginary part (factor -1/2pi is already included)
            tmpdos = array([real(tmpdos[0])]+[imag(ds) for ds in tmpdos[1:]])
            dosnew.append(tmpdos)
        
        # save to big array with all atoms
        dosnew_all_atoms.append(dosnew)
        
        if i1 != npot:
            text = f.readline() # dummy line
            
    dosnew_all_atoms = array(dosnew_all_atoms)
    dos_all_atoms = array(dos_all_atoms)
    
    # close complex.dos file
    f.close()
    
    if return_original:
        return ef, dos_all_atoms, dosnew_all_atoms
    else:
        return ef, dosnew_all_atoms
    
def get_ef_from_potfile(potfile):
    """
    extract fermi energy from potfile
    """
    f = open(potfile)
    txt = f.readlines()
    f.close()
    ef = float(txt[3].split()[1])
    return ef
  
