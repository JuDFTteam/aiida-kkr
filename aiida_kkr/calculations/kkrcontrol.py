#!/usr/bin/python

from __future__ import print_function
import sys

# redefine raw_input for python 3/2.7 compatilbility
if sys.version_info[0] >= 3:
    def raw_input(msg):
        return input(msg)


###############################################################################
### Function definitions

def print_warning(msg):
    """Print warning message, highlighted with '*' around it"""
    n = len(msg)
    print('*'*(n+4))
    print('* '+msg+' *')
    print('*'*(n+4))
    return 0


def abort_with_message(msg):
    """Prints warning message 'msg' and exits"""
    print_warning(msg)
    sys.exit()
    return 0


def cleanup_zpos(basis_vectors):
    """remove offset from z-position and make it positive"""

    #remove offset
    basis_vectors[:, 2] = basis_vectors[:, 2]-basis_vectors[-1, 2]

    #make positive z-component
    for ipos in range(basis_vectors.shape[0]):
        if(basis_vectors[ipos, 2] < 0):
            basis_vectors[ipos, 2] = basis_vectors[ipos, 2] +1

    return basis_vectors


def read_cpainfo(filename='inputcard'):
    """Read CPA information (NATYP, concentrations) from file 'filename'"""
    from scipy import zeros

    #open the template file and read in lines as a list
    tmpl_list = open(filename).read().split('\n')

    #get the number of atoms -- before the modification
    iline_natyp = [('NATYP' in item) for item in tmpl_list].index(True)
    natyp = int(tmpl_list[iline_natyp].split('=')[1].split()[0])

    #get the first line of the atominfo-block
    istart_atomlist = [('<CPA-CONC>' in iline) for iline in tmpl_list].index(True)

    #read in the CPA concentrations
    cpa_conc = zeros(natyp)
    for iat in range(natyp):
        cpa_conc[iat] = float(tmpl_list[istart_atomlist+iat+1].split()[3])

    return cpa_conc


def write_kkr_inputcard_template(bravais, natyp, basis_vectors, nucl_numbers, outfile='inputcard.tmpl'):
    """Write KKR inputcard template (starting point for further scripting)"""

    fout = open(outfile, 'w')

    #Run and test options
    fout.write('RUNOPT\n<RUNOP1><RUNOP2><RUNOP3><RUNOP4><RUNOP5><RUNOP6><RUNOP7><RUNOP8>\n+-------+-------+-------+-------+-------+-------+-------+-------\n')
    fout.write('TESTOPT\n<TSTOP1><TSTOP2><TSTOP3><TSTOP4><TSTOP5><TSTOP6><TSTOP7><TSTOP8>\n<TSTOP9><TSTP10><TSTP11><TSTP12><TSTP13><TSTP14><TSTP15><TSTP16>\n+-------+-------+-------+-------+-------+-------+-------+-------\n')
    fout.write('='*80+'\n')

    #System
    fout.write('NSPIN= [NSPIN]\n')
    fout.write('LMAX=  [LMAX]\n')
    fout.write('='*80+'\n')

    #accuracy
    fout.write('BZDIVIDE= [BZKX]  [BZKY]  [BZKZ]\n')
    fout.write('RCLUSTZ=  [RCLUSTZ]\n')
    fout.write('RCLUSTXY= [RCLUSTXY]\n')
    fout.write('='*80+'\n')

    #energy contour
    fout.write('EMIN= [EMIN]\n')
    fout.write('EMAX= [EMAX]\n')
    fout.write('TEMPR= [TEMPR]\n')
    fout.write('NPOL= [NPOL]\n')
    fout.write('NPT1= [NPT1]\n')
    fout.write('NPT2= [NPT2]\n')
    fout.write('NPT3= [NPT3]\n')
    fout.write('='*80+'\n')

    #convergence
    fout.write('NSTEPS= [NSTEPS]\n')
    fout.write('IMIX= [IMIX]\n')
    fout.write('STRMIX= [STRMIX]\n')
    fout.write('QBOUND= [QBOUND]\n')
    fout.write('='*80+'\n')

    #magnetism
    fout.write('LINIPOL= [LINIPOL]\n')
    #fout.write('XINIPOL= [XINIPOL]\n')
    fout.write('HFIELD= [HFIELD]\n')
    #fout.write('KHFELD= 0     (Apply ext. field in 1st iteration)\n')
    fout.write('='*80+'\n')

    #lattice
    fout.write('  ALATBASIS= [ALATBASIS] lattice constants a (in a.u.)\n')
    fout.write('-'*80+'\n')
    fout.write('   BRAVAIS                   (units of lattice constant)\n')
    fout.write('    %21.15f %21.15f %21.15f\n'%(bravais[0, 0], bravais[0, 1], bravais[0, 2]))
    fout.write('    %21.15f %21.15f %21.15f\n'%(bravais[1, 0], bravais[1, 1], bravais[1, 2]))
    fout.write('    %21.15f %21.15f %21.15f\n'%(bravais[2, 0], bravais[2, 1], bravais[2, 2]))
    fout.write('-'*80+'\n')

    #TODO: create NAEZ!=NATYP keyword and template
    fout.write('NAEZ= [NATYP]\n')
    fout.write('CARTESIAN= [CARTESIAN]\n')
    #TODO: change RBASIS template
    fout.write('<RBASIS>\n')
    for iatom in range(natyp):
        fout.write('  %18.14f  %18.14f  %18.14f\n'%(basis_vectors[iatom, 0], basis_vectors[iatom, 1], basis_vectors[iatom, 2]))
    fout.write('='*80+'\n')

    #slab geometry
    fout.write('INTERFACE= [INTERFACE]\n')
    fout.write('<NLBASIS>= [<NLBASIS>]\n')
    fout.write('<NRBASIS>= [<NRBASIS>]\n')
    fout.write('<RBLEFT>= [<RBLEFT>]\n')
    fout.write('<RBRIGHT>= [<RBRIGHT>]\n')
    fout.write('ZPERIODL= [ZPERIODL]\n')
    fout.write('ZPERIODR= [ZPERIODR]\n')

    #atoms
    fout.write('NATYP= [NATYP]\n')
    #TODO: change ZATOM template
    fout.write('<ZATOM>\n')
    for iatom in range(natyp):
        fout.write('%.1f\n'%(nucl_numbers[iatom]))
    fout.write('='*80+'\n')

    #various
    fout.write('cutoff parameter ewald summation (fcc 7  65)\n')
    fout.write('RMAX= [RMAX]\n')
    fout.write('GMAX= [GMAX]\n')
    fout.write('='*80+'\n')

    #new solver
    fout.write('NPAN_LOG= [NPAN_LOG]\n')
    fout.write('NPAN_EQ= [NPAN_EQ]\n')
    fout.write('NCHEB= [NCHEB]\n')
    fout.write('R_LOG= [R_LOG] \n')
    fout.write('='*80+'\n')

    fout.close()

    return 0


def create_keyword_default_values(**kwargs):
    """Fill default parameters of KKR inputcard to keywords dictionary"""
    default_keywords = dict([('NATYP', [1, '%i']),
                             ('ALATBASIS', [-1, '%f']),
                             ('BZKX', [20, '%i']),
                             ('BZKY', [20, '%i']),
                             ('BZKZ', [20, '%i']),
                             ('RCLUSTZ', [1.0, '%7.2f']),
                             ('RCLUSTXY', [1.0, '%7.2f']),
                             ('EMIN', [-0.4, '%6.3f']),
                             ('EMAX', [1.0, '%6.3f']),
                             ('TEMPR', [800.0, '%7.1f']),
                             ('NPOL', [5, '%4i']),
                             ('NPT1', [3, '%4i']),
                             ('NPT2', [20, '%4i']),
                             ('NPT3', [5, '%4i']),
                             ('NSTEPS', [50, '%4i']),
                             ('IMIX', [0, '%1i']),
                             ('LINIPOL', ['F', '%1s']),
                             #('XINIPOL', [0, '%1i']),
                             ('HFIELD', [0.0, '%5.2f']),
                             ('NPAN_LOG', [15, '%i']),
                             ('NPAN_EQ', [5, '%i']),
                             ('NCHEB', [13, '%i']),
                             ('R_LOG', [1.0, '%5.1f']),
                             ('QBOUND', ['1D-3', '%s']),
                             ('NSPIN', [1, '%i']),
                             ('LMAX', [2, '%i']),
                             ('RMAX', [7, '%i']),
                             ('GMAX', [65, '%i']),
                             ('CARTESIAN', ['F', '%1s']),
                             ('INTERFACE', ['F', '%1s']),
                             ('<NLBASIS>', [1, '%i']),
                             ('<NRBASIS>', [1, '%i']),
                             ('<RBLEFT>', [[0, 0, 0], ['%f'for i in range(3)]]),
                             ('<RBRIGHT>', [[0, 0, 0], ['%f' for i in range(3)]]),
                             ('ZPERIODL', [[0, 0, 0], ['%f' for i in range(3)]]),
                             ('ZPERIODR', [[0, 0, 0], ['%f'for i in range(3)]]),
                             ('NSTEPS', [30, '%i']),
                             ('IMIX', [0, '%i']),
                             ('STRMIX', [0.01, '%f'])])

    for key in kwargs:
        default_keywords[key][0] = kwargs[key]
        
    

    return default_keywords


def read_keywords_from_inputcard(keywords, inputcard='inputcard'):
    """read list of keywords from inputcard and extract values to keywords dict"""
    txt = open(inputcard, 'r').readlines()

    for key in keywords:
        keyvallen = 1
        try:
            if not type(keywords[key][0]) == str:
                keyvallen = len(list(keywords[key][0]))
            else:
                keyvallen = 1
        except:
            keyvallen = 1
        key2 = key.strip()
        if key2[:3] == 'BZK':
            searchkey = 'BZDIVIDE'
            ishift = ' XYZ'.index(key2[3])
        else:
            searchkey = key2
            ishift = 1
        try:
            retrvaltxt = get_value(searchkey, txt, item=ishift, num=keyvallen)
            #print(searchkey, ishift, keyvallen,  retrvaltxt, keywords[key][0], keywords[key][1][-1])
            if keyvallen == 1:
                keywords[key][0] = convert_string2number(retrvaltxt, keywords[key][1][-1])
            else:
                for i in range(len(keywords[key][0])):
                    #print(i,retrvaltxt[i])
                    keywords[key][0][i] = convert_string2number(retrvaltxt[i], keywords[key][1][i][-1])
        except:
            if keyvallen == 1:
                print(('Could not read keyword. Taking default value %s= '+keywords[key][1]+'.')%(key, keywords[key][0]))
            else:
                print(('Could not read keyword. Taking default values %s=')%(key))
                print(keywords[key][0])

    #read in runoptions
    istart = [iline for iline in range(len(txt)) if 'RUNOPT' in txt[iline]][0]
    line = txt[istart+1].strip()
    runopts = [line[i:i+8] for i in range(0, len(line), 8)]

    #read in testioptions (2 lines)
    istart = [iline for iline in range(len(txt)) if 'TESTOPT' in txt[iline]][0]
    line = txt[istart+1].strip()
    testopts = [line[i:i+8] for i in range(0, len(line), 8)]
    line = txt[istart+2].strip()
    testopts = testopts + [line[i:i+8] for i in range(0, len(line), 8)]

    return runopts, testopts

def convert_string2number(str, type):
    """Take input string and convert it to data type 'type'"""
    if type == 'i':
        return int(str)
    elif type == 'f':
        return float(str)
    elif type == 's':
        return str
    else:
        print_warning('Did not recognize type "%s" when converting string "%s"'%(type, str))
        return -1


def get_value(charkey, txt, line=1, item=1, num=1):
    """Search charkey in txt and return value string """
    iline = [ii for ii in range(len(txt)) if charkey in txt[ii]][0]
    txtline = txt[iline]
    chkeq = charkey+'='
    if chkeq in txtline:
        valtxt = txtline.split(chkeq)[1].split()[item-1:item-1+num]
    else:
        nextline = txt[iline+line]
        startpos = txtline.index(charkey)
        valtxt = nextline[startpos:].split()[item-1:item-1+num]
    if num == 1:
        return valtxt[0]
    else:
        return valtxt


def fill_keywords_to_inputcard(keywords, runops=[], testops=[], CPAconc=[], template='inputcard.tmpl', output='inputcard'):
    """Fill new inputcard with keywords/values and rest of template"""
    tmpl = open(template).read()
    tmpl_clean = tmpl
    added = 0
    for key in keywords:
        strkey = '[%s]'%(key)
        #print(strkey, tmpl)
        if strkey in tmpl:
            try:
                arraykey = False
                repltxt = keywords[key][1]%(keywords[key][0])
            except:
                arraykey = True
                repltxt = ''
                #print(strkey,key, keywords[key])
                for i in range(len(keywords[key][1])):
                    repltxt += ' ' + keywords[key][1][i]%(keywords[key][0][i])
            nblanks = len(strkey)-len(repltxt)
            if nblanks < 0 and not arraykey:
                print('WARNING for replacement of keyword %s: too long?'%(key))
            else:
                repltxt = repltxt+' '*nblanks
            tmpl = tmpl.replace(strkey, repltxt)
        else:
            print('WARNING for replacement of keyword %s: not in template!'%(key))
            if not key in tmpl_clean:
                print('... also naked keyword not found. Just appending it.')
                try:
                    repltxt = keywords[key][1]%(keywords[key][0])
                except:
                    repltxt = ''
                    for i in range(len(keywords[key][1])):
                        repltxt += ' ' + keywords[key][1][i]%(keywords[key][0][i])
                tmpl_clean = tmpl_clean+'%s= %s\n'%(key, strkey)
                tmpl = tmpl+'%s= %s\n'%(key, repltxt)
                added = added+1
            else:
                print('but found in template:', tmpl_clean)

    for ic in range(len(CPAconc)):
        strkey = '[cpa_%i]'%(ic+1)
        if strkey in tmpl:
            repltxt = '%7.5f'%(CPAconc[ic])
            tmpl = tmpl.replace(strkey, repltxt)
        elif '<CPA-CONC>' in tmpl:
            print('WARNING for replacement of keyword %s: not in template, but <CPA-CONC> is!'%(strkey))
        else:
            print('Error for replacement of keyword %s: no CPA mode in inputcard detected, but CPAconc-array is given!')
            sys.exit()



    for iop in range(len(runops)):
        strkey = '<RUNOP%i>'%(iop+1)
        repltxt = runops[iop]
        nblanks = len(strkey)-len(repltxt)
        if nblanks < 0:
            print('WARNING for replacement of RUNOPTION %s: too long?'%(repltxt))
        else:
            repltxt = repltxt+' '*nblanks
        tmpl = tmpl.replace(strkey, repltxt)

    for iop in range(len(testops)):
        strkey = '<TSTOP%i>'%(iop+1)
        repltxt = testops[iop]
        nblanks = len(strkey)-len(repltxt)
        if nblanks < 0:
            print('WARNING for replacement of TESTOPTION %s: too long?'%(repltxt))
        else:
            repltxt = repltxt+' '*nblanks
        tmpl = tmpl.replace(strkey, repltxt)

    open(output, 'w').write(tmpl)
    if added > 0:
        open(template, 'w').write(tmpl_clean)


def get_fermi_energy(potfile='potential'):
    """Extract Fermi level from potential file header"""
    return float(open(potfile).readlines()[3].split()[1])


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


def get_valence_min(outfile='out_voronoi'):
    """Construct minimum of energy contour (between valence band bottom and core states)"""
    from scipy import array
    txt = open(outfile).readlines()
    searchstr = 'All other states are above'
    valence_minimum = array([float(line.split(':')[1].split()[0]) for line in txt if searchstr in line])
    return valence_minimum


def check_voronoi_output(potfile, outfile):
    """Read output from voronoi code and create guess of energy contour"""
    from scipy import zeros
    #analyse core levels, minimum of valence band and their difference
    ncore, ecore, lcore = get_corestates_from_potential(potfile=potfile)
    e_val_min = get_valence_min(outfile=outfile)

    #print a table that summarizes the result
    e_core_max = zeros(len(ncore))
    print('pot    Highest core-level     low. val. state    diff')
    for ipot in range(len(ncore)):
        if ncore[ipot] > 0:
            lval, emax, descr = get_highest_core_state(ncore[ipot], ecore[ipot], lcore[ipot])
            e_core_max[ipot] = emax
            print('%3i     %2s %10.6f            %6.2f         %6.2f'%(ipot+1, descr, emax, e_val_min[ipot], e_val_min[ipot]-emax))
        else:
            print('%3i     << no core states >>'%(ipot+1))
            # set to some large negative number for check to not give false positive in case of empty cells
            e_core_max[ipot] = -100

    #get hint for energy integration:
    emin_guess = e_val_min.min()-0.2

    #make a check
    ediff_min = 1.5
    if (emin_guess-e_core_max.max()) < ediff_min:
        print_warning('Highest core level is too close to valence band: check results carefully!')

    return emin_guess


def generate_voronoi_potentials(keywords, dir='voronoi', outfile='out_voronoi', voro_cmd='./voronoi.exe', path_ElementsDatabase='./'):
    """Run voronoi code with parameter from inputcard_template"""
    import os, sys
    from subprocess import call

    #save current working directory
    savecwd = os.getcwd()

    #prepare voronoi-calculation folder
    dest = os.path.curdir+os.path.sep+dir
    if not os.path.exists(dest):
        os.mkdir(dest)

    #prepare the input file
    fill_keywords_to_inputcard(keywords, output=dest+os.path.sep+'inputcard')

    #change cwd
    os.chdir(dest)

    #make symbolic links to voronoi
    for element in ['ElementDataBase', 'voronoi.exe']:
        path = path_ElementsDatabase+os.path.sep+element
        path = os.path.abspath(os.path.expanduser(path))
        dest = os.path.curdir+os.path.sep+element
        if not os.path.exists(dest):
            os.symlink(path, dest)

    #execute voronoi
    voro_cmd = voro_cmd+" &> "+outfile
    try:
        retcode = call(voro_cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal "+str(-retcode), file=sys.stderr)
        else:
            print("Child returned "+str(retcode), file=sys.stderr)
    except OSError as e:
        print("Execution failed: "+str(e), file=sys.stderr)

    #perform some checks on the output
    emin_guess = check_voronoi_output(potfile='output.pot', outfile=outfile)

    os.chdir(savecwd)

    return emin_guess


def run_kkr_step(runcmd='./kkr.x', folderprefix='SYSTEM', runmode=None, outfile='out', save_folder=None, potfile='potential', shapefile='shapefun'):
    """Run KKR calculation step, inputcard, potential and shapefun needs to exist"""
    import os, sys
    from scipy import array
    from subprocess import call
    from shutil import copy2, move

    if save_folder is None:
        n = len(folderprefix)
        folders_exist = array([int(item[n:]) for item in os.listdir('.') if item[:n] == folderprefix])
        if len(folders_exist) == 0:
            nextnum = 1
        else:
            nextnum = folders_exist.max()+1
        save_folder = '%s%i'%(folderprefix, nextnum)

    #prepare kkr input (copy potfile to 'potential' and shapefile to 'shapfun' in workdir)
    if potfile != 'potential':
        copy2(potfile, 'potential')
    if shapefile != 'shapefun' and shapefile != '':
        copy2(shapefile, 'shapefun')

    #execute the kkr program
    try:
        retcode = call(runcmd+' &> %s'%(outfile), shell=True)
        if retcode < 0:
            print("KKR was terminated by signal "+str(-retcode), file=sys.stderr)
        else:
            print("KKR returned "+str(retcode), file=sys.stderr)
    except OSError as e:
        print("Execution of KKR failed: "+str(e), file=sys.stderr)

    #copy potential back
    if runmode == 'SCF':
        if os.path.isfile('out_potential'):
            copy2('out_potential', 'potential')
    elif runmode == 'DOS':
        call('~/sourcecodes/KKRcode/SRC_UTILS/complexdos3.exe', shell=True)
    elif runmode == 'JIJ':
        pass #do nothing
    else:
        print('WARNING: no runmode is given to run_kkr_step')

    #save the results
    files_to_copy = ['inputcard', 'potential', 'nonco_angle.dat', outfile]
    if shapefile != '':
        files_to_copy += ['shapefun']

    files_to_move = ['inputcard_generated.txt', 'out_potential', 'clusters', 'nonco_angle_out.dat'] + [file for file in os.listdir('.') if file[:6] == 'output'] + [file for file in os.listdir('.') if file[:3] == 'Jij'] + ['rimp.dat', 'shells.dat'] + [file for file in os.listdir('.') if 'dos' in file]

    import errno
    try:
        os.makedirs(save_folder)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(save_folder):
            folderprefix = 'SYSTMP'
            n = len(folderprefix)
            folders_exist = array([int(item[n:]) for item in os.listdir('.') if item[:n] == folderprefix])
            if len(folders_exist) == 0:
                nextnum = 1
            else:
                nextnum = folders_exist.max()+1
            save_folder = '%s%i'%(folderprefix, nextnum)
            os.mkdir(save_folder)
        else:
            raise

    for file in files_to_copy:
        if os.path.isfile(file):
            copy2(file, save_folder)
    for file in files_to_move:
        if os.path.isfile(file):
            move(file, save_folder)

    #parse out-file
    #get rms error and other relevant output-data
    results = {}
    outtxt = open(outfile).readlines()
    try:
        results['rms'] = array([float(outtxt[iline].split('=')[1].replace('D', 'E')) for iline in range(len(outtxt)) if 'average rms-error' in outtxt[iline]])
    except:
        results['rms'] = array([-1.0])
    try:
        results['charge_neutrality'] = array([float(outtxt[iline].split('=')[1]) for iline in range(len(outtxt)) if 'charge neutrality in unit cell' in outtxt[iline]])
    except:
        results['charge_neutrality'] = array([-1.0])
    try:
        results['total_magnetic_moment'] = array([float(outtxt[iline].split('=')[1]) for iline in range(len(outtxt)) if 'TOTAL mag. moment in unit cell' in outtxt[iline]])
    except:
        results['total_magnetic_moment'] = array([-1.0])
    try:
        results['EF'] = array([float(outtxt[iline].split('FERMI')[1].split()[0]) for iline in range(len(outtxt)) if 'E FERMI' in outtxt[iline]])
    except:
        results['EF'] = array([-1.0])
    try:
        results['DOS_EF'] = array([float(outtxt[iline].split('=')[1]) for iline in range(len(outtxt)) if 'DOS(E_F)' in outtxt[iline]])
    except:
        results['DOS_EF'] = array([-1.0])
    try:
        results['total_energy'] = array([float(outtxt[iline].split(':')[1]) for iline in range(len(outtxt)) if 'TOTAL ENERGY in ryd.' in outtxt[iline]])
    except:
        results['total_energy'] = array([-1.0])


    return save_folder, results


def scan_valenceband_bottom(keywords, emin, emax=None, runops=[], runcmd='kkr.x', delta_energy=0.05, threshold=1e-2, oldDOSpath=None):
    """Look for bottom of valence band (extracted from DOS run)"""
    import os
    from math import floor
    from scipy import loadtxt

    if oldDOSpath is None:
        if emax is None:
            emax = get_fermi_energy()
        nepts = int((emax-emin)/delta_energy)
        path = calc_DOS(keywords, runops=runops, emin=emin, emax=emax, npts=nepts)
    else:
        path = oldDOSpath

    #analyze new3.dos or valence_band bottom
    file = path+os.path.sep+'new3.dos'
    txt = open(file).readlines()
    nepts = int(txt[2].split('IELAST :')[1])-2

    newdos = loadtxt(file)
    energies = newdos[:nepts, 0]

    totdos = abs(newdos[:, 1].reshape((-1, nepts))).sum(axis=0)
    new_emin = energies[list(totdos > threshold).index(True)-1]
    new_emin_round = floor(new_emin/delta_energy)*delta_energy

    return new_emin_round


def calc_DOS(keywords, runops=[], testops=[], emin=None, emax=None, npts=100, tempr=473.0, runcmd='kkr.x', potfile='potential', shapefile='shapefun'):
    """Set DOS contour in inputcard and run KKR DOS calculation"""
    #save initial settings:
    saved_econt = {}
    for keys_econt in ['EMAX', 'TEMPR', 'NPT1', 'NPT2', 'NPT3', 'NPOL']:
        saved_econt[keys_econt] = keywords[keys_econt][0]
    keywords_tmp = keywords
    if not emin is None:
        keywords_tmp['EMIN'][0] = emin
    else:
        emin = keywords_tmp['EMIN'][0]
    if emax is None:
        emax = (get_fermi_energy(potfile=potfile)-emin)*1.5+emin
    keywords_tmp['EMAX'][0] = emax
    keywords_tmp['TEMPR'][0] = tempr
    keywords_tmp['NPT2'][0] = npts
    keywords_tmp['NPT1'][0] = 0
    keywords_tmp['NPT3'][0] = 0
    keywords_tmp['NPOL'][0] = 0

    fill_keywords_to_inputcard(keywords_tmp, runops=runops, testops=testops)
    path, results = run_kkr_step(runcmd=runcmd, folderprefix='DOS', runmode='DOS', outfile='out', potfile=potfile, shapefile=shapefile)

    #reset to initial settings:
    for keys_econt in ['EMAX', 'TEMPR', 'NPT1', 'NPT2', 'NPT3', 'NPOL']:
        keywords[keys_econt][0] = saved_econt[keys_econt]

    return path


def createNocoAngleFile(natyp, theta, phi, file='nonco_angle.dat'):
    """Write 'nonco_angle.dat' file with input theta and phi"""
    fo = open(file, 'w')
    for iatom in range(natyp):
        if type(theta) == list:
            thetatmp = theta[iatom]
        else:
            thetatmp = theta
        if type(phi) == list:
            phitmp = phi[iatom]
        else:
            phitmp = phi
        fo.write(' %15.11f   %15.11f\n'%(thetatmp, phitmp))
    fo.close()
    return 0

def kkr_dosplot(keywords, dospaths=[], units='Ry', plottitle='', econt=False, noefline=False):
    """Plot the dos of each entry (directory of dos calculation) in dospath list from new3.dos or new3_eV_EF.dos file"""
    from matplotlib.pyplot import plot, xlabel, ylabel, axvline, axhline, twinx, title
    from scipy import loadtxt, pi, array

    for path in dospaths:
        if units == 'Ry':
            d = loadtxt(path+'/new3.dos')
            EF0 = get_fermi_energy(path+'/potential')
        else:
            d = loadtxt(path+'/new3_eV_EF.dos')
            EF0 = 0
        plot(d[:, 0], d[:, 1], label=path)
        if not noefline:
            # vertical line marking EF
            axvline(EF0, color='k')
    # horizontal line showing DOS=0
    axhline(0, color='k')

    if units == 'Ry':
        xlabel('E (Ry)')
        ylabel('DOS (states/Ry)')
    else:
        xlabel('E-EF (eV)')
        ylabel('DOS (states/eV)')
    title(plottitle)

    # plot energy contour on top
    if econt:
        imEmax = pi * 6.333622661443*10**-6 * 2*keywords['NPOL'][0] * keywords['TEMPR'][0]
        econtour = []
        # segment 1: (EMIN, 0) to (EMIN, imEmax) in NPT1 steps
        Nepts = keywords['NPT1'][0]
        for ie in range(Nepts):
            econtour.append([keywords['EMIN'][0], imEmax/Nepts*ie])
        # segment 2: (EMIN, imEmax) to (EF, imEmax) in NPT2 steps
        Nepts = keywords['NPT2'][0]
        for ie in range(Nepts):
            econtour.append([keywords['EMIN'][0]+(EF0-keywords['EMIN'][0])/Nepts*ie, imEmax])
        # segment 3: (EF, imEmax) to (EF, 0) in NPT3 steps
        Nepts = keywords['NPT3'][0]
        for ie in range(Nepts+1):
            econtour.append([EF0, imEmax-imEmax/Nepts*ie])
        econtour = array(econtour)
        # now plotting on top of DOS (separate y-scale)
        ax = twinx()
        plot(econtour[:, 0], econtour[:, 1], 'rx-')
        ylabel('Im(E) (Ry)')
        ylim = ax.get_ylim()
        if keywords['NSPIN'][0] == 2:
            ax.set_ylim(-1.1*ylim[1], 1.1*ylim[1])
        else:
            ax.set_ylim(0, 1.1*ylim[1])


def kkr_plot_scfinfo(resultslist=[], plotrms=False, plotmom=False, plotEtot=False, ylog=False):
    """Plot results from resultslist (list, each entry should be the result of run_kkr_step)"""
    from matplotlib.pyplot import plot, xlabel, ylabel, axvline, gca
    from scipy import array
    
    if plotrms:
        rmslist = list(resultslist[0]['rms'])
        for i in range(len(resultslist)-1):
            rmslist += list(resultslist[i+1]['rms'])
        plot(rmslist, 'o-')
        ylabel('average rms-error')
        xlabel('Iteration')
    elif plotmom:
        momlist = list(resultslist[0]['total_magnetic_moment'])
        for i in range(len(resultslist)-1):
            momlist += list(resultslist[i+1]['total_magnetic_moment'])
        plot(momlist, 'rx-')
        xlabel('Iteration')
        ylabel('total spin moment')
    elif plotEtot:
        Etlist = list(resultslist[0]['total_energy'])
        for i in range(len(resultslist)-1):
            Etlist += list(resultslist[i+1]['total_energy'])
        if ylog:
            Etlist = abs(array(Etlist)-Etlist[-1])
            #Etlist = [0]+list(abs(array(Etlist[1:])-array(Etlist[:-1])))
        plot(Etlist, 'rx-')
        xlabel('Iteration')
        ylabel('total energy')
        if ylog:
            ylabel('change in total energy')
        
    if ylog:
        ax = gca()
        ax.set_yscale('log')
    
    # add vertical lines after each calculation in resultslist
    N0 = 0
    for i in range(len(resultslist)-1):
        N0 += len(resultslist[i]['rms'])
        axvline(N0-0.5)


def reset_emin(keywords, eVBbtm, edist_safe=0.15):
    """Check EMIN value from keywords against eVBbtm and change it to have a distance of edist_save"""
    emin = keywords['EMIN'][0]
    if abs(eVBbtm-emin) < edist_safe:
        print('EMIN too close to eVBbtm, changing it automatically!', emin, eVBbtm, eVBbtm-edist_safe)
        emin = eVBbtm-edist_safe
    else:
        print('EMIN too far away from eVBbtm, changing it automatically!', emin, eVBbtm, eVBbtm-edist_safe)
        emin = eVBbtm-edist_safe
    keywords['EMIN'][0] = emin
    return keywords

###############################################################################
### Eample workflows and workflow-helper functions


def workflow3_from_scratch(keywords, voropath, vorocommand, KKRpath, KKRcommand, runops=[], testops=[], maginit=False, iter1info=(30, 0, 0.01, 10**-3), iter2info=(30, 5, 0.01, 10**-6), iter3info=(100, 5, 0.01, 10**-6)):
    """general workflow taking preconfigured keywords dict and starts kkr calculation from scratch with plotting of DOS, rms, etc. """
    from matplotlib.pyplot import figure, show, subplot, twinx, subplots_adjust
    from time import time
    print('gen_voro')

    # needed for XINIPOL adding
    natyp = keywords['NATYP'][0]
    # write template to inputcard and generate starting potential
    start_time = time()
    emin = generate_voronoi_potentials(keywords, voro_cmd=vorocommand, path_ElementsDatabase=voropath)
    print('Elapsed time:', time()-start_time)
    
    print('starting DOS')

    if maginit:
        # init magnetic state
        print(natyp)
        keywords['LINIPOL'][0] = 'T'
        #TODO: fix XINIPOL template value
        keywords['XINIPOL'] = [[1 for i in range(natyp)], ['%i' for i in range(natyp)]]
        keywords['HFIELD'][0] = 0.02

    # change EMIN to scan larger range
    start_time = time()
    keywords['EMIN'][0] = emin-0.2
    pathDOS0 = calc_DOS(keywords, npts=50, runops=runops, potfile='voronoi/output.pot', shapefile='voronoi/shapefun', runcmd=KKRcommand)
    eVBbtm = scan_valenceband_bottom(keywords, keywords['EMIN'][0], oldDOSpath=pathDOS0)
    keywords = reset_emin(keywords, eVBbtm)
    print('Elapsed time:', time()-start_time)
    

    #save keywords to keywords0 for later reference (needed in plotting of energy contour)
    keywords0 = keywords.copy()

    # plot starting dos
    figure()
    dospaths = [pathDOS0]
    kkr_dosplot(keywords, dospaths=dospaths, units='Ry', plottitle='starting DOS', econt=True)
    show()
    #"""

    #"""
    print('First 30 simple mixing iterations')
    start_time = time()
    keywords['NSTEPS'][0] = iter1info[0]
    keywords['IMIX'][0] = iter1info[1]
    keywords['STRMIX'][0] = iter1info[2]
    keywords['QBOUND'][0] = iter1info[3]
    fill_keywords_to_inputcard(keywords, runops=runops, testops=testops, )
    pathSCF1, resultsSCF1 = run_kkr_step(potfile='voronoi/output.pot', shapefile='voronoi/shapefun', runcmd=KKRcommand, folderprefix='SYSTEM', runmode='SCF', outfile='out')
    print(pathSCF1, resultsSCF1)
    print('Elapsed time:', time()-start_time)

    if maginit:
        #disable init of magnetic moments
        keywords['LINIPOL'][0] = 'F'
        keywords['XINIPOL'][0] = [0 for i in range(natyp)]
        keywords['HFIELD'][0] = 0.00

    print('calculate DOS')
    start_time = time()
    pathDOS1 = calc_DOS(keywords, runops=runops, npts=50, potfile=pathSCF1+'/potential', shapefile='voronoi/shapefun', runcmd=KKRcommand)
    eVBbtm = scan_valenceband_bottom(keywords, keywords['EMIN'][0], oldDOSpath=pathDOS1)
    keywords = reset_emin(keywords, eVBbtm)
    print('Elapsed time:', time()-start_time)

    #save keywords to keywords1 for later reference
    keywords1 = keywords.copy()

    # plot dos, contour and rms error after initial simple mixing
    figure()
    subplot(2, 2, 1)
    kkr_plot_scfinfo(resultslist=[resultsSCF1], plotrms=True, ylog=True)
    if keywords['NSPIN'][0] == 2:
        twinx()
        kkr_plot_scfinfo(resultslist=[resultsSCF1], plotmom=True)
    subplot(2, 2, 2)
    kkr_plot_scfinfo(resultslist=[resultsSCF1], plotEtot=True)
    subplot(2, 2, 3)
    kkr_dosplot(keywords0, dospaths=[pathDOS0], units='Ry', plottitle='initial DOS', econt=True)
    subplot(2, 2, 4)
    kkr_dosplot(keywords1, dospaths=[pathDOS1], units='Ry', plottitle='DOS after simple mixing', econt=True)
    subplots_adjust(hspace=0.5, wspace=0.8)
    show()


    print('second set of iterations')

    start_time = time()
    keywords['NSTEPS'][0] = iter2info[0]
    keywords['IMIX'][0] = iter2info[1]
    keywords['STRMIX'][0] = iter2info[2]
    keywords['QBOUND'][0] = iter2info[3]
    fill_keywords_to_inputcard(keywords, runops=runops, testops=testops)
    pathSCF2, resultsSCF2 = run_kkr_step(runcmd=KKRcommand, folderprefix='SYSTEM', runmode='SCF', outfile='out', potfile=pathSCF1+'/potential')
    print(pathSCF2, resultsSCF2)
    print('Elapsed time:', time()-start_time)

    print('DOS calculation')
    start_time = time()
    pathDOS2 = calc_DOS(keywords, runops=runops, npts=50, potfile=pathSCF2+'/potential', shapefile='voronoi/shapefun', runcmd=KKRcommand)
    eVBbtm = scan_valenceband_bottom(keywords, keywords['EMIN'][0], oldDOSpath=pathDOS2)
    keywords = reset_emin(keywords, eVBbtm)
    print('Elapsed time:', time()-start_time)

    #save keywords to keywords1 for later reference
    keywords2 = keywords.copy()

    # plot dos, contour and rms error after second run
    figure()
    subplot(2, 2, 1)
    kkr_plot_scfinfo(resultslist=[resultsSCF1, resultsSCF2], plotrms=True, ylog=True)
    if keywords['NSPIN'][0] == 2:
        twinx()
        kkr_plot_scfinfo(resultslist=[resultsSCF1, resultsSCF2], plotmom=True)
    subplot(2, 2, 2)
    kkr_plot_scfinfo(resultslist=[resultsSCF1, resultsSCF2], plotEtot=True, ylog=True)
    subplot(2, 2, 3)
    kkr_dosplot(keywords1, dospaths=[pathDOS1], units='Ry', plottitle='previous run', econt=True)
    subplot(2, 2, 4)
    kkr_dosplot(keywords2, dospaths=[pathDOS2], units='Ry', plottitle='last run', econt=True)
    subplots_adjust(hspace=0.5, wspace=0.8)
    show()
    #"""

    #"""
    print('third set of iterations')

    start_time = time()
    keywords['NSTEPS'][0] = iter3info[0]
    keywords['IMIX'][0] = iter3info[1]
    keywords['STRMIX'][0] = iter3info[2]
    keywords['QBOUND'][0] = iter3info[3]

    fill_keywords_to_inputcard(keywords, runops=runops, testops=testops)
    pathSCF3, resultsSCF3 = run_kkr_step(runcmd=KKRcommand, folderprefix='SYSTEM', runmode='SCF', outfile='out', potfile=pathSCF2+'/potential')
    print(pathSCF3, resultsSCF3)
    print('Elapsed time:', time()-start_time)

    print('final DOS calculation with higher accuracy')
    start_time = time()
    keywords['BZKX'][0] = 50
    keywords['BZKY'][0] = 50
    keywords['BZKZ'][0] = 50
    keywords['TEMPR'][0] = 200
    pathDOS3 = calc_DOS(keywords, npts=51, runops=runops, potfile=pathSCF3+'/potential', shapefile='voronoi/shapefun', runcmd=KKRcommand)
    eVBbtm = scan_valenceband_bottom(keywords, keywords['EMIN'][0], oldDOSpath=pathDOS3)
    keywords = reset_emin(keywords, eVBbtm)
    print('Elapsed time:', time()-start_time)

    #save keywords to keywords1 for later reference
    keywords3 = keywords.copy()
    #"""


    # plot some results
    resultslist = [resultsSCF1, resultsSCF2, resultsSCF3]
    dospaths = [pathDOS0, pathDOS1, pathDOS2, pathDOS3]
    allpaths = [pathSCF1, pathSCF2, pathSCF3] + dospaths

    figure()
    subplot(2, 2, 1)
    kkr_plot_scfinfo(resultslist=resultslist, plotrms=True, ylog=True)
    if keywords['NSPIN'][0] == 2:
        twinx()
        kkr_plot_scfinfo(resultslist=resultslist, plotmom=True)
    subplot(2, 2, 2)
    kkr_plot_scfinfo(resultslist=resultslist, plotEtot=True, ylog=True)
    subplot(2, 2, 3)
    kkr_dosplot(keywords3, dospaths=dospaths, units='Ry', noefline=True, econt=True)
    subplot(2, 2, 4)
    kkr_dosplot(keywords, dospaths=dospaths, units='eV_EF')
    subplots_adjust(hspace=0.5, wspace=0.8)
    show()

    return keywords, allpaths, resultslist


def example_bcc110_Fe_2_2_2layer_lmax2_FP_noSOC(KKRpath, KKRcommand, voropath, vorocommand):
    """Example script doing 2-layer bcc-(110) thin film (with 2 vacuum layers on each side) of Fe in lmax 2 and full potential with magnetism"""
    from scipy import array, sqrt

    ### basic setup ###
    # lattice constant in a_Bohr
    alat = 5.416871386 # in a_Bohr
    # number of atom positions in unit cell
    natyp = 6
    # bravais vectors
    bravais = array([[0.7071067812, -0.5, 0.0], [0.7071067812, 0.5, 0.0], [sqrt(2), 0.0, 0.866025404]])

    # atom positions in relative coordinates
    basis_vectors = []
    for iatom in range(natyp):
        tmp = array([0, 0, 0])+iatom*array([0.5, 0.5, bravais[2, 2]])
        tmp[0] = tmp[0]%1
        tmp[1] = tmp[1]%1
        print(iatom, tmp)
        basis_vectors.append(tmp)
    basis_vectors = array(basis_vectors)

    # atomic charges of atoms in unit cell
    nucl_numbers = array([26. for i in range(natyp)])
    nucl_numbers[:2] = 0
    nucl_numbers[-2:] = 0

    print('write inputcard for voro')
    # write template imputcard and fill keywords dict with initial parameters
    write_kkr_inputcard_template(bravais, natyp, basis_vectors, nucl_numbers, outfile='inputcard.tmpl')
    keywords = create_keyword_default_values()
    keywords['NATYP'][0] = natyp
    keywords['ALATBASIS'][0] = alat
    keywords['NSPIN'][0] = 2
    keywords['LMAX'][0] = 2
    # for slab claculation
    keywords['INTERFACE'][0] = 'T'
    keywords['ZPERIODL'][0] = array([-0.5, -0.5, -bravais[2, 2]])
    keywords['ZPERIODR'][0] = array([0.5, 0.5, bravais[2, 2]])
    keywords['<RBLEFT>'][0] = basis_vectors[0]+keywords['ZPERIODL'][0]
    keywords['<RBRIGHT>'][0] = basis_vectors[natyp-1]+keywords['ZPERIODR'][0]
    keywords['BZKZ'][0] = 0
    # choose only coarse energy contour and k-mesh for test purposes
    keywords['NPOL'][0] = 4
    keywords['NPT1'][0] = 3
    keywords['NPT2'][0] = 10
    keywords['NPT3'][0] = 3
    keywords['BZKX'][0] = 10
    keywords['BZKY'][0] = 10
    
    #TODO: shape array for dublicate shapfuns, should be automatized from output of voronoi (atominfo)
    keywords['<SHAPE>'] = [[1 for i in range(natyp)], ['%i' for i in range(natyp)]]
    # run workflow
    keywords, allpaths, results_all = workflow3_from_scratch(keywords, voropath, vorocommand, KKRpath, KKRcommand, maginit=True, iter1info=(30, 0, 0.01, 10**-3), iter2info=(30, 5, 0.01, 10**-6), iter3info=(100, 5, 0.01, 10**-6))

    return keywords, allpaths, results_all


def example_bcc_Fe_lmax2_FP_noSOC(KKRpath, KKRcommand, voropath, vorocommand):
    """Example script doing bulk bcc Fe in lmax 2 and full potential with magnetism"""
    from scipy import array

    ### basic setup ###
    # lattice constant in a_Bohr
    alat = 5.416871386 # in a_Bohr
    # number of atom positions in unit cell
    natyp = 1
    # bravais vectors
    bravais = array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]])

    # atom positions in relative coordinates
    basis_vectors = array([[0, 0, 0]])

    # atomic charges of atoms in unit cell
    nucl_numbers = array([26. for i in range(natyp)])

    print('write inputcard for voro')
    # write template imputcard and fill keywords dict with initial parameters
    write_kkr_inputcard_template(bravais, natyp, basis_vectors, nucl_numbers, outfile='inputcard.tmpl')
    keywords = create_keyword_default_values()
    keywords['NATYP'][0] = natyp
    keywords['ALATBASIS'][0] = alat
    keywords['NSPIN'][0] = 2
    keywords['LMAX'][0] = 2
    # choose only coarse energy contour and k-mesh for test purposes
    keywords['NPOL'][0] = 4
    keywords['NPT1'][0] = 3
    keywords['NPT2'][0] = 10
    keywords['NPT3'][0] = 3
    keywords['BZKX'][0] = 10
    keywords['BZKY'][0] = 10

    # run workflow
    keywords, allpaths, results_all = workflow3_from_scratch(keywords, voropath, vorocommand, KKRpath, KKRcommand, maginit=True, iter1info=(20, 0, 0.01, 10**-3), iter2info=(20, 0, 0.05, 10**-6), iter3info=(100, 5, 0.05, 10**-6))

    return keywords, allpaths, results_all



def example_fcc_Cu_lmax3_FP_noSOC(KKRpath, KKRcommand, voropath, vorocommand):
    """Example script doing bulk Cu in lmax3 and full potential"""
    from scipy import array

    ### basic setup ###
    # lattice constant in a_Bohr
    alat = 6.830000 # in a_Bohr
    # number of atom positions in unit cell
    natyp = 1
    # bravais vectors
    bravais = array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])

    # atom positions in relative coordinates
    basis_vectors = array([[0, 0, 0]])

    # atomic charges of atoms in unit cell
    nucl_numbers = array([29. for i in range(natyp)])

    print('write inputcard for voro')
    # write template imputcard and fill keywords dict with initial parameters
    write_kkr_inputcard_template(bravais, natyp, basis_vectors, nucl_numbers, outfile='inputcard.tmpl')
    keywords = create_keyword_default_values()
    keywords['NATYP'][0] = natyp
    keywords['ALATBASIS'][0] = alat
    keywords['NSPIN'][0] = 1
    keywords['LMAX'][0] = 3
    # choose only coarse energy contour and k-mesh for test purposes
    keywords['NPOL'][0] = 4
    keywords['NPT1'][0] = 3
    keywords['NPT2'][0] = 10
    keywords['NPT3'][0] = 3
    keywords['BZKX'][0] = 10
    keywords['BZKY'][0] = 10

    # run workflow
    keywords, allpaths, results_all = workflow3_from_scratch(keywords, voropath, vorocommand, KKRpath, KKRcommand, iter1info=(30, 0, 0.05, 10**-3), iter2info=(20, 5, 0.05, 10**-6), iter3info=(100, 5, 0.05, 10**-6))

    return keywords, allpaths, results_all
    
    
def example_fcc_Cu_lmax2_ASA_noSOC(KKRpath, KKRcommand, voropath, vorocommand):
    """Exaple script doing a bulk Cu calculation in ASA"""
    from scipy import array

    ### basic setup ###
    # lattice constant in a_Bohr
    alat = 6.830000 # in a_Bohr
    # number of atom positions in unit cell
    natyp = 1
    # bravais vectors
    bravais = array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])

    # atom positions in relative coordinates
    basis_vectors = array([[0, 0, 0]])

    # atomic charges of atoms in unit cell
    nucl_numbers = array([29. for i in range(natyp)])

    print('write inputcard for voro')
    # write template imputcard and fill keywords dict with initial parameters
    write_kkr_inputcard_template(bravais, natyp, basis_vectors, nucl_numbers, outfile='inputcard.tmpl')
    keywords = create_keyword_default_values()
    keywords['NATYP'][0] = natyp
    keywords['ALATBASIS'][0] = alat
    keywords['NSPIN'][0] = 1
    keywords['LMAX'][0] = 2
    # choose only coarse energy contour and k-mesh for test purposes
    keywords['NPOL'][0] = 4
    keywords['NPT1'][0] = 3
    keywords['NPT2'][0] = 10
    keywords['NPT3'][0] = 3
    keywords['BZKX'][0] = 10
    keywords['BZKY'][0] = 10
    # for ASA
    keywords['INS'] = [0, '%i']
    keywords['KSHAPE'] = [0, '%i']

    # run workflow
    keywords, allpaths, results_all = workflow3_from_scratch(keywords, voropath, vorocommand, KKRpath, KKRcommand, iter1info=(30, 0, 0.05, 10**-3), iter2info=(20, 5, 0.05, 10**-6), iter3info=(100, 5, 0.05, 10**-6))

    return keywords, allpaths, results_all
    
    
def example_fcc_Cu_lmax2_ASA_noSOC_continued(keywords, startpath, KKRpath, KKRcommand):
    """Exaple script continuing a bulk Cu calculation in ASA"""

    # double these parameters
    keywords['NPT1'][0] = 6
    keywords['NPT2'][0] = 20
    keywords['NPT3'][0] = 6
    keywords['BZKX'][0] = 30
    keywords['BZKY'][0] = 30
    
    # run scf step
    print('Continued calculation with higher accuracy')
    start_time = time()
    fill_keywords_to_inputcard(keywords)
    pathSCF3, resultsSCF3 = run_kkr_step(runcmd=KKRcommand, folderprefix='SYSTEM', runmode='SCF', outfile='out', potfile=startpath+'/potential', shapefile=startpath+'/shapefun')
    print(pathSCF3, resultsSCF3)
    print('Elapsed time:', time()-start_time)

    print('final DOS calculation with higher accuracy')
    start_time = time()
    keywords['BZKX'][0] = 50
    keywords['BZKY'][0] = 50
    keywords['BZKZ'][0] = 50
    keywords['TEMPR'][0] = 200
    pathDOS3 = calc_DOS(keywords, npts=51, potfile=startpath+'/potential', shapefile=startpath+'/shapefun', runcmd=KKRcommand)
    eVBbtm = scan_valenceband_bottom(keywords, keywords['EMIN'][0], oldDOSpath=pathDOS3)
    keywords = reset_emin(keywords, eVBbtm)
    print('Elapsed time:', time()-start_time)

    #save keywords to keywords1 for later reference
    keywords3 = keywords.copy()
    #"""


    # plot some results
    resultslist = [resultsSCF3]
    dospaths = [pathDOS3]
    allpaths = [pathSCF3] + dospaths

    from matplotlib.pyplot import figure, subplot, twinx, show, subplots_adjust
    figure()
    subplot(2, 2, 1)
    kkr_plot_scfinfo(resultslist=resultslist, plotrms=True, ylog=True)
    subplot(2, 2, 2)
    kkr_plot_scfinfo(resultslist=resultslist, plotEtot=True, ylog=True)
    subplot(2, 2, 3)
    kkr_dosplot(keywords3, dospaths=dospaths, units='Ry', noefline=True, econt=True)
    subplot(2, 2, 4)
    kkr_dosplot(keywords, dospaths=dospaths, units='eV_EF')
    subplots_adjust(hspace=0.5, wspace=0.8)
    show()
    
    return keywords, allpaths, resultslist
###############################################################################

if __name__ == '__main__':
    from shutil import copy2
    from subprocess import call
    from time import time
    
    # code directories etc. (souce commands needed since KKR and voronoi are compiled with this environment)
    KKRpath = './' #'/Users/ruess/sourcecodes/KKRcode/'
    KKRcommand = 'source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64; '+KKRpath+'kkr.x'
    vorocommand = 'source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64; ./voronoi.exe'
    voropath = '/Users/ruess/sourcecodes/voronoi/'
    
    start_time = time()
    # cleanup old DOS and SYSTEM folders
    call('rm -r DOS* SYSTEM*', shell=True)
    """
    print('Test0: fcc Cu bulk lmax3 full potential')    
    call('rm -r inputcard* potential shapefun voronoi/ lmdos.0* *form* gmat gref tmat out* clusters fort.*', shell=True)
    copy2('kkr_lmax3_natyp1_noSOC.x', 'kkr.x')
    keywords0, allpaths0, results_all0 = example_fcc_Cu_lmax3_FP_noSOC(KKRpath, KKRcommand, voropath, vorocommand)
    #"""
    """
    print('Test1: bcc Fe bulk lmax2 full potential magnetic noSOC')    
    call('rm -r inputcard* potential shapefun voronoi/ lmdos.0* *form* gmat gref tmat out* clusters fort.*', shell=True)
    copy2('kkr_lmax2_natyp1_nosoc.x', 'kkr.x')
    keywords1, allpaths1, results_all1 = example_bcc_Fe_lmax2_FP_noSOC(KKRpath, KKRcommand, voropath, vorocommand)
    #"""
    """
    print('Test2: bcc Fe 110 2 layers lmax2 full potential')    
    call('rm -r inputcard* potential shapefun voronoi/ lmdos.0* *form* gmat gref tmat out* clusters fort.*', shell=True)
    copy2('kkr_bcc110_slab.x', 'kkr.x')
    keywords2, allpaths2, results_all2 = example_bcc110_Fe_2_2_2layer_lmax2_FP_noSOC(KKRpath, KKRcommand, voropath, vorocommand)
    #"""
    #"""
    print('Test3: fcc Cu bulk lmax2 ASA noSOC')
    call('rm -r inputcard* potential shapefun voronoi/ lmdos.0* *form* gmat gref tmat out* clusters fort.*', shell=True)
    copy2('kkr_lmax2_natyp1_noSOC_ASA.x', 'kkr.x')
    keywords3, allpaths3, results_all3 = example_fcc_Cu_lmax2_ASA_noSOC(KKRpath, KKRcommand, voropath, vorocommand)
    #"""
    """
    print('Test4: continue Test3 and run with higher accuracy (higher number of kpts, Epts)')
    call('rm -r voronoi/ lmdos.0* *form* gmat gref tmat out* clusters fort.*', shell=True)
    copy2('kkr_lmax2_natyp1_noSOC_ASA.x', 'kkr.x')
    keywords4, allpaths4, results_all4 = example_fcc_Cu_lmax2_ASA_noSOC_continued(keywords3, allpaths3[2], KKRpath, KKRcommand)
    #"""
    
    stop_time = time()
    print('Total runtime:', stop_time-start_time)