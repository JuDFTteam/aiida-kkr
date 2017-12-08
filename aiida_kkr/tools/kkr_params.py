#!/usr/bin/python


#use print('message') instead of print 'message' in python 2.7 as well:
from __future__ import print_function
# redefine raw_input for python 3/2.7 compatilbility
from sys import version_info
if version_info[0] >= 3:
    def raw_input(msg):
        return input(msg)


class kkrparams(object):
    """
    Class for creating and handling the parameter input for a KKR calculation
    Optional keyword arguments are passed to init and stored in values dictionary.

    Example usage:
    params = kkrparams(LMAX=3, BRAVAIS=array([[1,0,0], [0,1,0], [0,0,1]]))

    Alternatively values can be set afterwards either individually with
        params.set_value('LMAX', 3)
    or multiple keys at once with
        params.set_multiple_values(EMIN=-0.5, EMAX=1)

    Other useful functions:
    - print the description of a keyword: params.get_description([key]) where [key] is a string for a keyword in params.values
    - print a list of mandatory keywords: params.get_all_mandatory()
    - print a list of keywords that are set including their value: params.get_set_values()

    Note: KKR-units (e.g. atomic units with energy in Ry, length in a_Bohr) are assumed!
    """


    def __init__(self, **kwargs):
        """
        Initialize class instance with containing the attribute values that also have
        a format, mandatory flags (defaults for KKRcode, changed for example via params_type='voronoi' keyword) and a description.
        """
        if 'params_type' in kwargs:
            self.__params_type = kwargs.pop('params_type')
        else:
            #parameter are set for kkr or voronoi code? (changes mandatory flags)
            self.__params_type = 'kkr' #default value, also possible: 'voronoi', 'kkrimp'
            #TODO: implement update_mandatory for self.__params_type=='kkrimp'

        keyw = self._create_keyword_default_values(**kwargs)

        #values of keywords:
        self.values = {}
        #formatting info
        self.__format = {}
        #mandatory flag
        self._mandatory = {}
        # description of each key
        self.__description = {}

        for key in keyw:
            self.values[key] = keyw[key][0]
            self.__format[key] = keyw[key][1]
            self._mandatory[key] = keyw[key][2]
            self.__description[key] = keyw[key][3]
            

    def get_dict(self, group=None, subgroup=None):
        """
        Returns values dictionary.
        
        Prints values belonging to a certain group only if the 'group' argument
        is one of the following: 'lattice', 'chemistry', 'accuracy', 
                                 'external fields', 'scf cycle', 'other'
          
        Additionally the subgroups argument allows to print only a subset of 
        all keys in a certain group. The following subgroups are available:
        in 'lattice' group:   '2D mode', 'shape functions'
        in 'chemistry' group: 'Atom types', 'Exchange-correlation', 'CPA mode', 
                              '2D mode'
        in 'accuracy' group:  'Valence energy contour', 'Semicore energy contour',
                              'CPA mode', 'Screening clusters', 'Radial solver', 
                              'Ewald summation', 'LLoyd'
        """
        out_dict = self.values
        
        #check for grouping
        group_searchstrings = {'lattice':'Description of lattice', 
                               'chemistry':'Chemistry', 
                               'external fields':'External fields:', 
                               'accuracy':'Accuracy', 
                               'scf cycle':'Self-consistency control:', 
                               'other':['Running and test options', 'Name of potential and shapefun file']}
        subgroups_all = {'lattice':['2D mode', 'shape functions'], 
                         'chemistry':['Atom types', 'Exchange-correlation', 'CPA mode', '2D mode'],
                         'accuracy':['Valence energy contour', 'Semicore energy contour',
                                     'CPA mode', 'Screening clusters', 'Radial solver', 
                                     'Ewald summation', 'LLoyd']}
        if group in ['lattice', 'chemistry', 'accuracy', 'external fields', 'scf cycle', 'other']:
            print('Returning only values belonging to group %s'%group)
            tmp_dict = {}
            for key in out_dict.keys():
                desc = self.__description[key]
                key_in_group = False
                if group_searchstrings[group] != 'other':
                    if group_searchstrings[group] in desc:
                        key_in_group = True
                else:
                    if group_searchstrings[group][0] in desc or group_searchstrings[group][1] in desc:
                        key_in_group = True
                if key_in_group:
                    tmp_dict[key] = self.values[key]
            
            #check for subgrouping and overwrite tmp_dict accordingly
            if group in ['lattice', 'chemistry', 'accuracy']:
                if subgroup in subgroups_all[group]:
                    print('Restrict keys additionally to subgroup %s'%subgroup)
                    tmp_dict2 = {}
                    for key in tmp_dict.keys():
                        desc = self.__description[key]
                        key_in_group = False
                        if subgroup in desc:
                            key_in_group = True
                            if key_in_group:
                                tmp_dict2[key] = self.values[key]
                    tmp_dict = tmp_dict2
                    
            # overwrite out_dict with tmp_dict
            out_dict = tmp_dict
        
        return out_dict

        
    def _get_type_from_string(self, fmtstr):
        """Helper function of get_type"""
        if 'f' in fmtstr:
            keytype = float
        elif 'i' in fmtstr:
            keytype = int
        elif 'l' in fmtstr:
            keytype = bool
        elif 's' in fmtstr:
            keytype = str
        else:
            print('Error: type of keyvalue not found:', fmtstr)
            raise TypeError
        return keytype


    def get_type(self, key):
        """Extract expected type of 'key' from format info"""
        fmtstr = self.__format[key]
        # simple format or complex pattern
        simplefmt = True
        if fmtstr.count('%') > 1:
            simplefmt = False
        if simplefmt:
            keytype = self._get_type_from_string(fmtstr)
        else:
            fmtlist = fmtstr.replace('\n','').replace(' ','').split('%')[1:]
            keytype = []
            for fmtstr in fmtlist:
                keytype.append(self._get_type_from_string(fmtstr))
        return keytype


    def _check_valuetype(self, key):
        """Consistency check if type of value matches expected type from format info"""

        # check if entry is numpy array and change to list automatically:
        try:
            tmpval = self.values[key].flatten().tolist()
        except:
            tmpval = self.values[key]
        tmptype = type(tmpval)

        # get type of value
        if tmptype == list:
            valtype = []
            for val in range(len(tmpval)):
                valtype.append(type(tmpval[val]))
        else:
            valtype = tmptype
        #print(key, valtype, self.get_type(key))

        # check if type matches format info
        cmptypes = self.get_type(key)
        success = True
        #print(key, type(valtype), valtype, cmptypes)
        changed_type_automatically = False
        if valtype == int and cmptypes == float:
            changed_type_automatically = True
            self.values[key] = float(self.values[key])
        elif type(valtype) == list:
            for ival in range(len(valtype)):
                if valtype[ival] == int and cmptypes == float:
                    changed_type_automatically = True
                    self.values[key][ival] = float(self.values[key][ival])
        elif valtype != cmptypes and tmpval is not None:
            success = False
            print('Error: type of value does not match expected type for ', key, self.values[key], cmptypes)
            raise TypeError

        if changed_type_automatically:
            print('Warning: filling value of "%s" with integer but expects float. Converting automatically and continue'%key)

        return success


    def set_value(self, key, value):
        """Sets value of keyword 'key'"""
        if value is None:
            print('Warning setting value None is not permitted!')
            print('Use remove_value funciton instead! Ignore keyword {}'.format(key))
        else:
            self.values[key] = value
            self._check_valuetype(key)
        
        
    def remove_value(self, key):
        """Removes value of keyword 'key', i.e. resets to None"""
        self.values[key] = None


    def set_multiple_values(self, **kwargs):
        """Set multiple values (in example value1 and value2 of keywords 'key1' and 'key2') given as key1=value1, key2=value2"""
        for key in kwargs:
            key2 = key
            if key not in self.values.keys():
                key2 = '<'+key+'>'
            #print('setting', key2, kwargs[key])
            self.set_value(key2, kwargs[key])


    def get_set_values(self):
        """Return a list of all keys/values that are set (i.e. not None)"""
        set_values = []
        added = 0
        for key in self.values.keys():
            if self.values[key] is not None:
                set_values.append([key, self.values[key]])
                added += 1
        if added == 0:
            print('No values set')
        return set_values


    def get_all_mandatory(self):
        """Return a list of mandatory keys"""
        self._update_mandatory()
        mandatory_list = []
        for key in self.values.keys():
            if self.is_mandatory(key):
                mandatory_list.append(key)
        return mandatory_list


    def is_mandatory(self, key):
        """Returns mandatory flag (True/False) for keyword 'key'"""
        return self._mandatory[key]


    def get_description(self, key):
        """Returns description of keyword 'key'"""
        return self.__description[key]


    def _create_keyword_default_values(self, **kwargs):
        """ Creates KKR inputcard keywords dictionary and fills entry if value is given in **kwargs

            entries of keyword dictionary are: 'keyword', [value, format, keyword mandatory, description]
            where
            value can be a single entry or a list of entries
            format contains formatting info
            keyword mandatory is a logical stating if keyword needs to be defined to run a calculation
            description is a string containgin human redable info about the keyword
        """
        default_keywords = dict([# complete list of keywords, detault all that are not mandatory to None
                                # lattice
                                ('ALATBASIS', [None, '%f', True, 'Description of lattice: Length unit in Bohr radii usually conventional lattice parameter']),
                                ('BRAVAIS', [None, '%f %f %f\n%f %f %f\n%f %f %f', True, 'Description of lattice: Bravais vectors in units of [ALATBASIS]']),
                                ('NAEZ', [None, '%i', True, 'Description of lattice: Number of sites in unit cell']),
                                ('<RBASIS>', [None, '%f %f %f', True, 'Description of lattice: Positions of sites in unit cell']),
                                ('CARTESIAN', [None, '%l', False, 'Description of lattice: Interpret the basis vector coordinates as reduced (w. respect to bravais) or as cartesian (in lattice constant units)']),
                                ('INTERFACE', [None, '%l', False, 'Description of lattice, 2D mode: needs to be TRUE for 2D calculation']),
                                ('<NLBASIS>', [None, '%i', False, 'Description of lattice, 2D mode: Number of basis sites forming the half-infinite lattice to the lower (=left) part of the slab.']),
                                ('<RBLEFT>', [None, '%f %f %f', False, 'Description of lattice, 2D mode: Positions of sites forming the basis sites of the half-infinite lattice to the lower (=left) part of the slab.']),
                                ('ZPERIODL', [None, '%f %f %f', False, 'Description of lattice, 2D mode: Lattice vector describing the periodicity perpendicular to the slab-plane for the half-infinite lattice to the lower (=left) part of the slab (plays the role of the 3rd Bravais vector for this half-infinite lattice). The <RBLEFT> vectors are periodically repeated by the ZPERIODL vector.']),
                                ('<NRBASIS>', [None, '%i', False, 'Description of lattice, 2D mode: Number of basis sites forming the half-infinite lattice to the upper (=right) part of the slab.']),
                                ('<RBRIGHT>', [None, '%f %f %f', False, 'Description of lattice, 2D mode: Positions of sites forming the basis sites of the half-infinite lattice to the upper (=right) part of the slab.']),
                                ('ZPERIODR', [None, '%f %f %f', False, 'Description of lattice, 2D mode: Lattice vector describing the periodicity perpendicular to the slab-plane for the half-infinite lattice to the upper (=right) part of the slab (plays the role of the 3rd Bravais vector for this half-infinite lattice). The <RBRIGHT> vectors are periodically repeated by the ZPERIODR vector.']),
                                ('KSHAPE', [None, '%i', False, 'Description of lattice, shape functions: 0 for ASA ([INS]=0), 2 for full potential ([INS]=1)']),
                                ('<SHAPE>', [None, '%i', False, 'Description of lattice, shape functions: Indexes which shape function from the shape-function file to use in which atom. Default is that each atom has its own shape function.']),
                                # chemistry
                                ('<ZATOM>', [None, '%f', True, 'Chemistry, Atom types: Nuclear charge per atom. Negative value signals to use value read in from the potential file.']),
                                ('NSPIN', [None, '%i', True, 'Chemistry, Atom types: Number of spin directions in potential. Values 1 or 2']),
                                ('KVREL', [None, '%i', False, 'Chemistry, Atom types: Relativistic treatment of valence electrons. Takes values 0 (Schroedinger), 1 (Scalar relativistic), 2 (Dirac ; works only in ASA mode)']),
                                ('<SOCSCL>', [None, '%f', False, 'Chemistry, Atom types: Spin-orbit coupling scaling per atom. Takes values between 0. (no spin-orbit) and 1. (full spin-orbit). Works only in combination with the Juelich spin orbit solver (runoption NEWSOSOL)']),
                                ('KEXCOR', [None, '%i', False, 'Chemistry, Exchange-correlation: Type of exchange correlation potential. Takes values 0 (LDA, Moruzzi-Janak-Williams), 1 (LDA, von Barth-Hedin), 2 (LDA, Vosko-Wilk-Nussair), 3 (GGA, Perdew-Wang 91), 4 (GGA, PBE), 5 (GGA, PBEsol)']),
                                ('LAMBDA_XC', [None, '%f', False, 'Chemistry, Exchange-correlation: Scale the magnetic part of the xc-potential and energy. Takes values between 0. (fully suppressed magnetisc potential) and 1. (normal magnetic potential).']),
                                ('NAT_LDAU', [None, '%i', False, 'Chemistry, Exchange-correlation: Numer of atoms where LDA+U will be used']),
                                ('LDAU_PARA', [None, '%i %i %f %f %f', False, 'Chemistry, Exchange-correlation: For each atom where LDA+U should be used, the entries are: [atom type] [angular mom. to apply LDA+U] [Ueff] [Jeff] [Eref] where [atom type] is between 1...[NATYP].']),
                                ('KREADLDAU', [None, '%i', False, "Chemistry, Exchange-correlation: Takes values 0 or 1; if [KREADLDAU]=1 then read previously calculated LDA+U matrix elements from file 'ldaupot'."]),
                                ('NATYP', [None, '%i', False, 'Chemistry, CPA mode: Number of atom types; CPA is triggered by setting [NATYP]>[NAEZ].']),
                                ('<SITE>', [None, '%i', False, 'Chemistry, CPA mode: Takes values 1 < [<SITE>] < [NAEZ] Assigns the position (given by [<RBASIS>]) where the atom-dependent read-in potential is situated. E.g., if the 3rd-in-the-row potential should be positioned at the 2nd <RBASIS> vector, then the 3rd entry of the <SITE> list should have the value 2.']),
                                ('<CPA-CONC>', [None, '%f', False, 'Chemistry, CPA mode: Takes values 0. < [<CPA-CONC>] < 1. Assigns the alloy-concentration corresponding to the atom-dependent read-in potential. Together with the variable <SITE>, <CPA-CONC> assigns the number and concentration of the atom-dependent potentials residing at each site form 1 to [NAEZ]. The sum of concentrations at each site should equal 1.']),
                                ('<KAOEZL>', [None, '%i', False, 'Chemistry, 2D mode: Controls the type of t-matrix at the lower (=left) half-crystal sites in case of embedding as these are given in the left-decimation file (i.e., changes the order compared to the one in the left-decimation file).']),
                                ('<KAOEZR>', [None, '%i', False, 'Chemistry, 2D mode: Controls the type of t-matrix at the upper (=right) half-crystal sites in case of embedding as these are given in the right-decimation file (i.e., changes the order compared to the one in the right-decimation file).']),
                                # external fields
                                ('LINIPOL', [None, '%l', False, 'External fields: If TRUE, triggers an external magn. field per atom in the first iteration.']),
                                ('HFIELD', [None, '%f', False, 'External fields: Value of an external magnetic field in the first iteration. Works only with LINIPOL, XINIPOL']),
                                ('XINIPOL', [None, '%i', False, 'External fields: Integer multiplying the HFIELD per atom']),
                                ('VCONST', [None, '%f', False, 'External fields: Constant potential shift in the first iteration.']),
                                # accuracy
                                ('LMAX', [None, '%i', True, 'Accuracy: Angular momentum cutoff']),
                                ('BZDIVIDE', [None, '%i %i %i', False, 'Accuracy: Maximal Brillouin zone mesh. Should not violate symmetry (e.g cubic symmetry implies i1=i2=i3; terragonal symmetry in xy implies i1=i2; i1=i2=i3 is always safe.)']),
                                ('EMIN', [None, '%f', False, 'Accuracy, Valence energy contour: Lower value (in Ryd) for the energy contour']),
                                ('EMAX', [None, '%f', False, 'Accuracy, Valence energy contour: Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases']),
                                ('TEMPR', [None, '%f', False, 'Accuracy, Valence energy contour: Electronic temperature in K.']),
                                ('NPT1', [None, '%i', False, 'Accuracy, Valence energy contour: Number of energies in the 1st part of the rectangular contour ("going up").']),
                                ('NPT2', [None, '%i', False, 'Accuracy, Valence energy contour: Number of energies in the 2nd part of the rectangular contour ("going right").']),
                                ('NPT3', [None, '%i', False, 'Accuracy, Valence energy contour: Number of energies in the 3rd part of the rectangular contour (Fermi smearing part).']),
                                ('NPOL', [None, '%i', False, 'Accuracy, Valence energy contour: Number of Matsubara poles For DOS calculations, set [NPOL]=0']),
                                ('EBOTSEMI', [None, '%f', False, 'Accuracy, Semicore energy contour: Bottom of semicore contour in Ryd.']),
                                ('EMUSEMI', [None, '%f', False, 'Accuracy, Semicore energy contour: Top of semicore contour in Ryd.']),
                                ('TKSEMI', [None, '%f', False, 'Accuracy, Semicore energy contour: "Temperature" in K controlling height of semicore contour.']),
                                ('NPOLSEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Control of height of semicore contour: Im z = (2 * [NPOLSEMI] * pi * kB * [TKSEMI] ) with kB=0.6333659E-5']),
                                ('N1SEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Number of energies in first part of semicore contour ("going up").']),
                                ('N2SEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Number of energies in second part of semicore contour ("going right").']),
                                ('N3SEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Number of energies in third part of semicore contour ("going down").']),
                                ('FSEMICORE', [None, '%f', False, 'Accuracy, Semicore energy contour: Initial normalization factor for semicore states (approx. 1.).']),
                                ('CPAINFO', [None, '%f %i', False, 'Accuracy, CPA mode: CPA-error max. tolerance and max. number of CPA-cycle iterations.']),
                                ('RCLUSTZ', [None, '%f', False, 'Accuracy, Screening clusters: Radius of screening clusters in units of [ALATBASIS], default is 11 Bohr radii.']),
                                ('RCLUSTXY', [None, '%f', False, 'Accuracy, Screening clusters: If [RCLUSTXY] does not equal [RCLUSTZ] then cylindrical clusters are created with radius [RCLUSTXY] and height [RCLUSTZ].']),
                                ('<RMTREF>', [None, '%f', False, 'Accuracy, Screening clusters: Muffin tin radius in Bohr radii for each site forming screening clusters. Negative value signals automatic calculation by the code.']),
                                ('NLEFTHOS', [None, '%i', False, 'Accuracy, Screening clusters 2D mode: The vectors [<RBLEFT>] are repeated i=1,...,[NLEFTHOS] times, shifted by i*[ZPERIODL], for the later formation of screening clusters.']),
                                ('<RMTREFL>', [None, '%f', False, 'Accuracy, Screening clusters 2D mode: Muffin-tin radius in Bohr radii for each site forming screening clusters in the lower (=left) half-crystal. Negative value signals automatic calculation by the code.']),
                                ('NRIGHTHO', [None, '%i', False, 'Accuracy, Screening clusters 2D mode: The vectors [<RBRIGHT>] are repeated i=1,...,[NRIGHTHO] times, shifted by i*[ZPERIODR], for the later formation of screening clusters.']),
                                ('<RMTREFR>', [None, '%f', False, 'Accuracy, Screening clusters 2D mode: Muffin-tin radius in Bohr radii for each site forming screening clusters in the upper (=right) half-crystal. Negative value signals automatic calculation by the code.']),
                                ('INS', [None, '%i', False, 'Accuracy, Radial solver: Takes values 0 for ASA and 1 for full potential Must be 0 for Munich Dirac solver ([KREL]=2)']),
                                ('ICST', [None, '%i', False, 'Accuracy, Radial solver: Number of iterations in the radial solver']),
                                ('R_LOG', [None, '%f', False, 'Accuracy, Radial solver: Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL']),
                                ('NPAN_LOG', [None, '%i', False, 'Accuracy, Radial solver: Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL']),
                                ('NPAN_EQ', [None, '%i', False, 'Accuracy, Radial solver: Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL']),
                                ('NCHEB', [None, '%i', False, 'Accuracy, Radial solver: Number of Chebyshev polynomials per interval Used in conjunction with runopt NEWSOSOL']),
                                ('<FPRADIUS>', [None, '%f', False, 'Accuracy, Radial solver: Full potential limit per atom (in Bohr radii); at points closer to the nucleus, the potential is assumed spherical. Negative values indicate to use values from potential file. Values larger than the muffin tin indicate to use the muffin tin radius.']),
                                ('RMAX', [None, '%f', True, 'Accuracy, Ewald summation for Madelung potential: Max. radius in [ALATBASIS] for real space Ewald sum']),
                                ('GMAX', [None, '%f', True, 'Accuracy, Ewald summation for Madelung potential: Max. radius in 2*pi/[ALATBASIS] for reciprocal space Ewald sum']),
                                ('<LLOYD>', [None, '%i', False, "Accuracy, LLoyd's formula: Set to 1 in order to use Lloyd's formula"]),
                                ('<DELTAE>', [None, '(%f, %f)', False, "Accuracy, LLoyd's formula: Energy difference for derivative calculation in Lloyd's formula"]),
                                ('<TOLRDIF>', [None, '%f', False, 'Accuracy, Virtual atoms: For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.']),
                                # scf cycle
                                ('NSTEPS', [None, '%i', False, 'Self-consistency control: Max. number of self-consistency iterations. Is reset to 1 in several cases that require only 1 iteration (DOS, Jij, write out GF).']),
                                ('IMIX', [None, '%i', False, "Self-consistency control: Mixing scheme for potential. 0 means straignt (linear) mixing, 3 means Broyden's 1st method, 4 means Broyden's 2nd method, 5 means Anderson's method"]),
                                ('STRMIX', [None, '%f', False, 'Self-consistency control: Linear mixing parameter Set to 0. if [NPOL]=0']),
                                ('ITDBRY', [None, '%i', False, 'Self-consistency control: ow many iterations to keep in the Broyden/Anderson mixing scheme.']),
                                ('FCM', [None, '%f', False, 'Self-consistency control: Factor for increased linear mixing of magnetic part of potential compared to non-magnetic part.']),
                                ('BRYMIX', [None, '%f', False, 'Self-consistency control: Parameter for Broyden mixing.']),
                                ('QBOUND', [None, '%f', False, 'Self-consistency control: Lower limit of rms-error in potential to stop iterations.']),
                                #code options
                                ('RUNOPT', [None, '%s%s%s%s%s%s%s%s', False, 'Running and test options: 8-character keywords in a row without spaces between them']),
                                ('TESTOPT', [None, '%s%s%s%s%s%s%s%s\n%s%s%s%s%s%s%s%s', False, 'Running and test options: optional 8-character keywords in a row without spaces between them plus a secod row of the same.']),
                                #file names
                                ('FILES', [None, '%s', False, 'Name of potential and shapefun file (list of two strings, empty string will set back to default of the one file that is supposed to be changed)'])
                                ])

        for key in kwargs:
            key2 = key
            if key not in default_keywords.keys():
                key2 = '<'+key+'>'
            default_keywords[key2][0] = kwargs[key]

        return default_keywords


    def _update_mandatory(self):
        """Check if mandatory flags need to be updated if certain keywords are set"""
        # initialize all mandatory flags to False and update list afterwards
        for key in self.values.keys():
            self._mandatory[key] = False

        runopts = []
        if self.values['RUNOPT'] is not None:
            for runopt in self.values['RUNOPT']:
                runopts.append(runopt.strip())

        #For a KKR calculation these keywords are always mandatory:
        mandatory_list = ['ALATBASIS', 'BRAVAIS', 'NAEZ', '<RBASIS>', 'NSPIN', 'LMAX', 'RMAX', 'GMAX', '<ZATOM>']

        if self.values['NPOL'] is not None and self.values['NPOL'] != 0:
                mandatory_list += ['EMIN']
        #Mandatory in 2D
        if self.values['INTERFACE']:
            mandatory_list += ['<NLBASIS>', '<RBLEFT>', 'ZPERIODL', '<NRBASIS>', '<RBRIGHT>', 'ZPERIODR']
        #Mandatory in LDA+U
        if 'NAT_LDAU' in self.values.keys() and 'LDAU' in runopts:
            mandatory_list += ['NAT_LDAU', 'LDAU_PARA']
        #Mandatory in CPA
        if self.values['NATYP'] is not None and self.values['NATYP'] > self.values['NAEZ']:
            mandatory_list += ['NATYP', '<SITE>', '<CPA-CONC>']
        #Mandatory in SEMICORE
        if 'EBOTSEMI' in self.values.keys() and 'SEMICORE' in runopts:
            mandatory_list += ['EBOTSEMI', 'EMUSEMI', 'TKSEMI', 'NPOLSEMI', 'N1SEMI', 'N2SEMI', 'N3SEMI', 'FSEMICORE']
        if self.values['INS'] == 1 and 'WRITEALL' not in runopts:
            mandatory_list += ['<SHAPE>']

        for key in mandatory_list:
            self._mandatory[key] = True

        # overwrite if mandatory list needs to be changed (determinded from value of self.__params_type):
        if self.__params_type == 'voronoi':
            self.update_mandatory_voronoi()


    def _check_mandatory(self):
        """Check if all mandatory keywords are set"""
        self._update_mandatory()
        if self.__params_type == 'voronoi':
            self.update_mandatory_voronoi()
        for key in self.values.keys():
            if self._mandatory[key] and self.values[key] is None:
                print('Error not all mandatory keys are set!')
                set_of_mandatory = set(self.get_all_mandatory())
                set_of_keys = set([key[0] for key in self.get_set_values()])
                print(set_of_mandatory-set_of_keys, 'missing')
                raise ValueError


    def _check_array_consistency(self):
        """Check all keys in __listargs if they match their specification (mostly 1D array, except for special cases e.g. <RBASIS>)"""
        from numpy import array, ndarray

        vec3_entries = ['<RBASIS>', '<RBLEFT>', '<RBRIGHT>', 'ZPERIODL', 'ZPERIODR']

        #success = [True]
        for key in self.__listargs.keys():
            if self.values[key] is not None:
                tmpsuccess = True
                #print('checking', key, self.values[key], self.__listargs[key])
                if type(self.values[key]) not in [list, ndarray]:
                    self.values[key] = array([self.values[key]])
                cmpdims = (self.__listargs[key], )
                if key in vec3_entries:
                    cmpdims = (self.__listargs[key], 3)
                    # automatically convert if naez==1 and only 1D array is given
                    if self.__listargs[key] == 1 and len(array(self.values[key]).shape) == 1 and key not in ['ZPERIODL', 'ZPERIODR']:
                        print('Warning: expected 2D array for %s but got 1D array, converting automatically'%key)
                        self.values[key] = array([self.values[key]])
                tmpdims = array(self.values[key]).shape
                if tmpdims[0] != cmpdims[0]:
                    tmpsuccess = False
                if len(tmpdims)==2:
                    if tmpdims[1] != cmpdims[1]:
                        tmpsuccess = False
                #success.append(tmpsuccess)
        
                if not tmpsuccess:
                    print('Error: array input not consistent')
                    print('check consistency:', key, self.values[key], cmpdims, tmpdims, tmpsuccess)
                    raise TypeError


    def _check_input_consistency(self):
        """Check consistency of input, to be done before wrinting to inputcard"""
        from numpy import array

        # first check if all mandatory values are there
        self._check_mandatory()

        # lists of array arguments
        keywords = self.values
        naez = keywords['NAEZ']
        if keywords['NATYP'] is not None:
            natyp = keywords['NATYP']
        else:
            natyp = keywords['NAEZ']
        if keywords['<NLBASIS>'] is not None:
            nlbasis = keywords['<NLBASIS>']
        else:
            nlbasis = 1
        if keywords['<NRBASIS>'] is not None:
            nrbasis = keywords['<NRBASIS>']
        else:
            nrbasis = 1

        listargs = dict([['<RBASIS>', naez], ['<RBLEFT>', nlbasis], ['<RBRIGHT>', nrbasis], ['<SHAPE>', natyp],
                         ['<ZATOM>', natyp], ['<SOCSCL>', natyp], ['<SITE>', natyp], ['<CPA-CONC>', natyp],
                         ['<KAOEZL>', nlbasis], ['<KAOEZR>', nrbasis], ['XINIPOL', natyp], ['<RMTREF>', natyp],
                         ['<RMTREFL>', nlbasis], ['<RMTREFR>', nrbasis], ['<FPRADIUS>', natyp], ['BZDIVIDE', 3],
                         ['<RBLEFT>', nrbasis], ['ZPERIODL', 3], ['<RBRIGHT>', nrbasis], ['ZPERIODR', 3],
                         ['LDAU_PARA', 5], ['CPAINFO', 2], ['<DELTAE>', 2], ['FILES', 2]])
        special_formatting = ['BRAVAIS', 'RUNOPT', 'TESTOPT', 'FILES']

        self.__listargs = listargs
        self.__special_formatting = special_formatting

        # check for consistency of array arguments
        self._check_array_consistency()

        # some special checks
        bulkmode = False
        set_values = [key[0] for key in self.get_set_values()]
        if 'INTERFACE' not in set_values or self.values['INTERFACE']:
            bulkmode = True
            
        bravais = array(self.values['BRAVAIS'])
        if bulkmode and sum(bravais[2]**2)==0:
            print("Error: 'BRAVAIS' matches 2D calculation but 'INTERFACE' is not set to True!")
            raise ValueError


    def fill_keywords_to_inputfile(self, is_voro_calc=False, output='inputcard'):
        """
        Fill new inputcard with keywords/values
        automatically check for input consistency
        if is_voro_calc==True change mandatory list to match voronoi code, default is KKRcode
        """
        from numpy import array
        #rename for easy reference
        keywords = self.values
        keyfmts = self.__format

        sorted_keylist = [#run/testopts
                          'RUNOPT', 'TESTOPT',
                          #lattice:
                          'ALATBASIS', 'BRAVAIS', 'NAEZ', '<RBASIS>', 'CARTESIAN',
                          'INTERFACE', '<NLBASIS>', '<RBLEFT>', 'ZPERIODL', '<NRBASIS>', '<RBRIGHT>', 'ZPERIODR',
                          'KSHAPE', '<SHAPE>',
                          # chemistry
                          '<ZATOM>', 'NSPIN', 'KVREL', '<SOCSCL>',
                          'KEXCOR', 'LAMBDA_XC',
                          'NAT_LDAU', 'LDAU_PARA', 'KREADLDAU',
                          'NATYP', '<SITE>', '<CPA-CONC>',
                          '<KAOEZL>', '<KAOEZR>',
                          # external fields
                          'LINIPOL', 'HFIELD', 'XINIPOL', 'VCONST',
                          # accuracy
                          'LMAX', 'BZDIVIDE', 'EMIN', 'EMAX', 'TEMPR', 'NPT1', 'NPT2', 'NPT3', 'NPOL',
                          'EBOTSEMI', 'EMUSEMI', 'TKSEMI', 'NPOLSEMI', 'N1SEMI', 'N2SEMI', 'N3SEMI', 'FSEMICORE',
                          'CPAINFO',
                          'RCLUSTZ', 'RCLUSTXY',
                          '<RMTREF>', 'NLEFTHOS', '<RMTREFL>', 'NRIGHTHO', '<RMTREFR>',
                          'INS', 'ICST',
                          'R_LOG', 'NPAN_LOG', 'NPAN_EQ', 'NCHEB', '<FPRADIUS>',
                          'RMAX', 'GMAX', '<LLOYD>', '<DELTAE>', '<TOLRDIF>',
                          # scf cycle
                          'NSTEPS', 'IMIX', 'STRMIX', 'ITDBRY', 'FCM', 'BRYMIX', 'QBOUND',
                          #file names
                          'FILES']
        #add everything that was forgotten in sorted_keylist above
        for key in keywords.keys():
            if key not in sorted_keylist:
                sorted_keylist += [key]

        #for key in keyfmts.keys():
        #    keyfmts[key] = keyfmts[key].replace('%l', '%s')

        if is_voro_calc:
            self.__params_type = 'voronoi'

        # check for inconsistencies in input before writing file
        self._check_input_consistency()

        # write all set keys to file
        tmpl = ''
        for key in sorted_keylist:
            if keywords[key] is not None:
                #print(key)
                if (not key in self.__listargs.keys()) and (not key in self.__special_formatting):
                    tmpfmt = (keyfmts[key]).replace('%l', '%s')
                    try:
                        repltxt = tmpfmt%(keywords[key])
                    except:
                        print(key, tmpfmt, keywords[key])
                        repltxt = ''
                        for i in range(len(tmpfmt)):
                            repltxt += ' ' + tmpfmt[i]%(keywords[key][i])
                    tmpl += '%s= %s\n'%(key, repltxt)
                elif key == 'BRAVAIS':
                    self.values[key] = array(self.values[key])
                    tmpl += ('BRAVAIS\n'+self.__format[key]+'\n')%(self.values[key][0, 0], self.values[key][0, 1], self.values[key][0, 2],
                                                                   self.values[key][1, 0], self.values[key][1, 1], self.values[key][1, 2],
                                                                   self.values[key][2, 0], self.values[key][2, 1], self.values[key][2, 2])
                elif key == 'RUNOPT':
                    runops = keywords[key]
                    tmpl += 'RUNOPT\n'
                    for iop in range(len(runops)):
                        repltxt = runops[iop]
                        nblanks = 8 - len(repltxt)
                        if nblanks < 0:
                            print('WARNING for replacement of RUNOPTION %s: too long?'%repltxt)
                            print('RUNOPT %s is ignored and was not set!'%repltxt)
                        else:
                            repltxt = repltxt+' '*nblanks
                        tmpl += repltxt
                    tmpl += '\n'
                elif key == 'TESTOPT':
                    testops = keywords[key]
                    tmpl += 'TESTOPT\n'
                    for iop in range(len(testops)):
                        repltxt = testops[iop]
                        nblanks = 8 - len(repltxt)
                        if nblanks < 0:
                            print('WARNING for replacement of TESTOPTION %s: too long?'%repltxt)
                            print('TESTOPT %s is ignored and was not set!'%repltxt)
                        else:
                            repltxt = repltxt+' '*nblanks
                        tmpl += repltxt
                        if iop==8:
                            tmpl += '\n'
                    tmpl += '\n'
                elif key == 'XINIPOL':
                    tmpl += '%s='%key
                    for ival in range(len(self.values[key])):
                        tmpl += (' %s'%self.__format[key])%self.values[key][ival]
                    tmpl += '\n'
                elif key == 'FILES':
                    files_changed = 0
                    if self.values[key][0]=='':
                        self.values[key][0]='potential'
                    else:
                        files_changed += 1
                    if self.values[key][1]=='':
                        self.values[key][1]='shapefun'
                    else:
                        files_changed += 1
                    if files_changed>0:
                        print('Warning: Changing file name of potential file to "%s" and of shapefunction file to "%s"'%(self.values[key][0], self.values[key][1]))
                        tmpl += 'FILES\n'
                        tmpl += '\n'
                        tmpl += '%s\n'%self.values[key][0]
                        tmpl += '\n'
                        tmpl += '%s\n'%self.values[key][1]
                        tmpl += '\n'
                elif key in self.__listargs.keys():
                    if key in ['<RBASIS>', '<RBLEFT>', '<RBRIGHT>']: # RBASIS needs special formatting since three numbers are filled per line
                        tmpl += '%s\n'%key
                        for ival in range(self.__listargs[key]):
                            tmpl += (self.__format[key]+'\n')%(self.values[key][ival][0], self.values[key][ival][1], self.values[key][ival][2])
                    elif key in ['CPAINFO', '<DELTAE>']:
                        tmpl += '%s= '%key
                        tmpl += (self.__format[key]+'\n')%(self.values[key][0], self.values[key][1])
                    elif key in ['BZDIVIDE', 'ZPERIODL', 'ZPERIODR']:
                        tmpl += '%s= '%key
                        tmpl += (self.__format[key]+'\n')%(self.values[key][0], self.values[key][1], self.values[key][2])
                    elif key in ['LDAU_PARA']:
                        tmpl += '%s= '%key
                        tmpl += (self.__format[key]+'\n')%(self.values[key][0], self.values[key][1], self.values[key][2], self.values[key][3], self.values[key][4])
                    else:
                        #print(key, self.__listargs[key], len(self.values[key]))
                        tmpl += '%s\n'%key
                        for ival in range(self.__listargs[key]):
                            tmpl += (self.__format[key]+'\n')%(self.values[key][ival])
                else:
                    print('Error trying to write keyword %s but writing failed!'%key)
                    raise ValueError

        open(output, 'w').write(tmpl)


    def read_keywords_from_inputcard(self, inputcard='inputcard'):
        """Read list of keywords from inputcard and extract values to keywords dict"""
        """
        txt = open(inputcard, 'r').readlines()
        keywords = self.values
        keyfmts = self.__format

        using _get_value (from below)
        """
        #TODO loop over known keywords and fill with values found in inputcard
        #TODO take different known formatting into account, not only the default that is written here
        print('Not implemented yet!')
        pass

    '''
    def _get_value(self, charkey, txt, line=1, item=1, num=1):
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
    '''

    # redefine _update_mandatory for voronoi code
    def update_mandatory_voronoi(self):
        """Change mandatory flags to match requirements of voronoi code"""
        # initialize all mandatory flags to False and update list afterwards
        for key in self.values.keys():
            self._mandatory[key] = False

        runopts = []
        if self.values['RUNOPT'] is not None:
            for runopt in self.values['RUNOPT']:
                runopts.append(runopt.strip())

        #For a KKR calculation these keywords are always mandatory:
        mandatory_list = ['ALATBASIS', 'BRAVAIS', 'NAEZ', '<RBASIS>', 'NSPIN', 'LMAX', 'RCLUSTZ', '<ZATOM>']

        #Mandatory in 2D
        if self.values['INTERFACE']:
            mandatory_list += ['<NLBASIS>', '<RBLEFT>', 'ZPERIODL', '<NRBASIS>', '<RBRIGHT>', 'ZPERIODR']
        #Mandatory in CPA
        if self.values['NATYP'] is not None and self.values['NATYP'] > self.values['NAEZ']:
            mandatory_list += ['NATYP', '<SITE>', '<CPA-CONC>']

        for key in mandatory_list:
            self._mandatory[key] = True
       
            
    def get_missing_keys(self, use_aiida=False):
        setlist = dict(self.get_set_values()).keys()
        manlist = self.get_all_mandatory()
        missing = []
        for key in manlist:
            if key not in setlist:
                if not use_aiida:
                    missing.append(key)
                else:
                    if key not in ['BRAVAIS', '<RBASIS>', '<ZATOM>', 'ALATBASIS', 'NAEZ']:
                        missing.append(key)
        return missing
