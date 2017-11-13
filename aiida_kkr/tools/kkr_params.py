#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 22:19:00 2017

@author: ruess
"""        
        
class kkrparams(object):
    #TODO: docstring """"""
    
        
    def __init__(self, **kwargs):
        """initialize class instance with attributes values that also have format, mandatory and description information."""
        keyw = self._create_keyword_default_values(**kwargs)
        
        self.values = {}
        self.__format = {}
        self._mandatory = {}
        self.__description = {}
        
        for key in keyw:
            self.values[key] = keyw[key][0]
            self.__format[key] = keyw[key][1]
            self._mandatory[key] = keyw[key][2]
            self.__description[key] = keyw[key][3]
      
    
    def _get_type_from_string(self, fmtstr):
        """helper function of get_type"""
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
        """extract expected type of 'key' from format info"""
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
        """consistency check if type of value matches expected type from format info"""
        
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
        print(key, valtype, self.get_type(key)) 
        
        # check if type matches format info
        cmptypes = self.get_type(key)
        success = True
        if valtype == int and cmptypes == float:
            print('Warning: filling value of "%s" with integer but expects float. Converting automatically and continue'%key)
            self.values[key] = float(self.values[key])
        elif type(valtype) == list:
            for ival in range(len(valtype)):
                if valtype[ival] == int and cmptypes[ival] == float:
                    print('Warning: filling value of "%s" with integer but expects float. Converting automatically and continue'%key)
                    self.values[key][ival] = float(self.values[key][ival])
        elif valtype != cmptypes:
            success = False
            print('Error: type of value does not match expected type')
            raise TypeError
                
        return success
    
        
    def set_value(self, key, value):
        """sets value of keyword 'key'"""
        self.values[key] = value
        type_consistent = self._check_valuetype(key)
        
        
    def set_multiple_values(self, **kwargs):
        """set multiple values (in example value1 and value2 of keywords 'key1' and 'key2') given as key1=value1, key2=value2"""
        for key in kwargs:
            key2 = key
            if key not in self.values.keys():
                key2 = '<'+key+'>'
            print('setting', key2, kwargs[key])
            self.set_value(key2, kwargs[key])
            

    def get_all_values(self):
        set_values = []
        added = 0
        for key in self.values.keys():
            if self.values[key] is not None:
                set_values.append([key, self.values[key]])
                added += 1
        if added == 0:
            print('Not values set')
        return set_values
    
    
    def get_all_mandatory(self):
        self._update_mandatory()
        mandatory_list = []
        for key in self.values.keys():
            if self._mandatory[key]:
                mandatory_list.append(key)
        return mandatory_list
            
            
    def get_mandatory(self, key):
        """returns mandatory flag (True/False) for keyword 'key'"""
        self._update_mandatory()
        return self._mandatory[key]
    
        
    def get_description(self, key):
        """returns description of keyword 'key'"""
        return self.__description[key]

        
    def _create_keyword_default_values(self, **kwargs):
        """ creates KKR inputcard keywords dictionary and fills entry if value is given in **kwargs
        
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
                                ('<RBLEFT>', [None, '%f', False, 'Description of lattice, 2D mode: Positions of sites forming the basis sites of the half-infinite lattice to the lower (=left) part of the slab.']),
                                ('ZPERIODL', [None, '%f', False, 'Description of lattice, 2D mode: Lattice vector describing the periodicity perpendicular to the slab-plane for the half-infinite lattice to the lower (=left) part of the slab (plays the role of the 3rd Bravais vector for this half-infinite lattice). The <RBLEFT> vectors are periodically repeated by the ZPERIODL vector.']),
                                ('<NRBASIS>', [None, '%i', False, 'Description of lattice, 2D mode: Number of basis sites forming the half-infinite lattice to the upper (=right) part of the slab.']),
                                ('<RBRIGHT>', [None, '%f', False, 'Description of lattice, 2D mode: Positions of sites forming the basis sites of the half-infinite lattice to the upper (=right) part of the slab.']),
                                ('ZPERIODR', [None, '%f', False, 'Description of lattice, 2D mode: Lattice vector describing the periodicity perpendicular to the slab-plane for the half-infinite lattice to the upper (=right) part of the slab (plays the role of the 3rd Bravais vector for this half-infinite lattice). The <RBRIGHT> vectors are periodically repeated by the ZPERIODR vector.']),
                                ('KSHAPE', [None, '%i', False, 'Description of lattice, shape functions: 0 for ASA ([INS]=0), 2 for full potential ([INS]=1)']),
                                ('<SHAPE>', [None, '%i', False, 'Description of lattice, shape functions: Indexes which shape function from the shape-function file to use in which atom. Default is that each atom has its own shape function.']),
                                # chemistry
                                ('<ZATOM>', [None, '%f', False, 'Chemistry, Atom types: Nuclear charge per atom. Negative value signals to use value read in from the potential file.']),
                                ('NSPIN', [None, '%i', True, 'Chemistry, Atom types: Number of spin directions in potential. Values 1 or 2']),
                                ('KVREL', [None, '%i', False, 'Chemistry, Atom types: Relativistic treatment of valence electrons. Takes values 0 (Schroedinger), 1 (Scalar relativistic), 2 (Dirac ; works only in ASA mode)']),
                                ('<SOCSCL>', [None, '%f', False, 'Chemistry, Atom types: Spin-orbit coupling scaling per atom. Takes values between 0. (no spin-orbit) and 1. (full spin-orbit). Works only in combination with the Juelich spin orbit solver (runoption NEWSOSOL)']),
                                ('KEXCOR', [None, '%i', False, 'Chemistry, Exchange-correlation: Type of exchange correlation potential. Takes values 0 (LDA, Moruzzi-Janak-Williams), 1 (LDA, von Barth-Hedin), 2 (LDA, Vosko-Wilk-Nussair), 3 (GGA, Perdew-Wang 91), 4 (GGA, PBE), 5 (GGA, PBEsol)']),
                                ('LAMBDA_XC', [None, '%f', False, 'Chemistry, Exchange-correlation: Scale the magnetic part of the xc-potential and energy. Takes values between 0. (fully suppressed magnetisc potential) and 1. (normal magnetic potential).']),
                                ('NAT_LDAU', [None, '%i', False, 'Chemistry, Exchange-correlation: Numer of atoms where LDA+U will be used']),
                                ('LDAU_PARA', [None, '%i %i %f %f %f', False, 'Chemistry, Exchange-correlation: For each atom where LDA+U should be used, the entries are: [atom type] [angular mom. to apply LDA+U] [Ueff] [Jeff] [Eref] where [atom type] is between 1…[NATYP].']),
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
                                ('BZDIVIDE', [None, '%i', False, 'Accuracy: Maximal Brillouin zone mesh. Should not violate symmetry (e.g cubic symmetry implies i1=i2=i3; terragonal symmetry in xy implies i1=i2; i1=i2=i3 is always safe.)']),
                                ('EMIN', [None, '%f', False, 'Accuracy, Valence energy contour: Lower value (in Ryd) for the energy contour']),
                                ('EMAX', [None, '%f', False, 'Accuracy, Valence energy contour: Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases']),
                                ('TEMPR', [None, '%f', False, 'Accuracy, Valence energy contour: Electronic temperature in K.']),
                                ('NPT1', [None, '%i', False, 'Accuracy, Valence energy contour: Number of energies in the 1st part of the rectangular contour (“going up”).']),
                                ('NPT2', [None, '%i', False, 'Accuracy, Valence energy contour: Number of energies in the 2nd part of the rectangular contour (“going right”).']),
                                ('NPT3', [None, '%i', False, 'Accuracy, Valence energy contour: Number of energies in the 3rd part of the rectangular contour (Fermi smearing part).']),
                                ('NPOL', [None, '%i', False, 'Accuracy, Valence energy contour: Number of Matsubara poles For DOS calculations, set [NPOL]=0']),
                                ('EBOTSEMI', [None, '%f', False, 'Accuracy, Semicore energy contour: Bottom of semicore contour in Ryd.']),
                                ('EMUSEMI', [None, '%f', False, 'Accuracy, Semicore energy contour: Top of semicore contour in Ryd.']),
                                ('TKSEMI', [None, '%f', False, 'Accuracy, Semicore energy contour: “Temperature” in K controlling height of semicore contour.']),
                                ('NPOLSEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Control of height of semicore contour: Im z = (2 * [NPOLSEMI] * pi * kB * [TKSEMI] ) with kB=0.6333659E-5']),
                                ('N1SEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Number of energies in first part of semicore contour (“going up”).']),
                                ('N2SEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Number of energies in second part of semicore contour (“going right”).']),
                                ('N3SEMI', [None, '%i', False, 'Accuracy, Semicore energy contour: Number of energies in third part of semicore contour (“going down”).']),
                                ('FSEMICORE', [None, '%f', False, 'Accuracy, Semicore energy contour: Initial normalization factor for semicore states (approx. 1.).']),
                                ('CPAINFO', [None, '%r %i', False, 'Accuracy, CPA mode: CPA-error max. tolerance and max. number of CPA-cycle iterations.']),
                                ('RCLUSTZ', [None, '%f', False, 'Accuracy, Screening clusters: Radius of screening clusters in units of [ALATBASIS], default is 11 Bohr radii.']),
                                ('RCLUSTXY', [None, '%f', False, 'Accuracy, Screening clusters: If [RCLUSTXY] does not equal [RCLUSTZ] then cylindrical clusters are created with radius [RCLUSTXY] and height [RCLUSTZ].']),
                                ('<RMTREF>', [None, '%f', False, 'Accuracy, Screening clusters: Muffin tin radius in Bohr radii for each site forming screening clusters. Negative value signals automatic calculation by the code.']),
                                ('NLEFTHOS', [None, '%i', False, 'Accuracy, Screening clusters 2D mode: The vectors [<RBLEFT>] are repeated i=1,…,[NLEFTHOS] times, shifted by i*[ZPERIODL], for the later formation of screening clusters.']),
                                ('<RMTREFL>', [None, '%f', False, 'Accuracy, Screening clusters 2D mode: Muffin-tin radius in Bohr radii for each site forming screening clusters in the lower (=left) half-crystal. Negative value signals automatic calculation by the code.']),
                                ('NRIGHTHO', [None, '%i', False, 'Accuracy, Screening clusters 2D mode: The vectors [<RBRIGHT>] are repeated i=1,…,[NRIGHTHO] times, shifted by i*[ZPERIODR], for the later formation of screening clusters.']),
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
                                ('NSTEPS', [None, '%i', False, 'Self-consistency controlMax. number of self-consistency iterations. Is reset to 1 in several cases that require only 1 iteration (DOS, Jij, write out GF).']),
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
                                ('FILES', [None, '%s', False, 'Name of potential file'])
                                ])
    
        for key in kwargs:
            default_keywords[key][0] = kwargs[key]        
    
        return default_keywords
        
        
    def _update_mandatory(self):
        runopts = []
        if self.values['RUNOPT'] is not None:
            for runopt in self.values['RUNOPT']:
                runopts.append(runopt.strip())
        
        #For a KKR calculation these keywords are always mandatory:
        mandatory_list = ['ALATBASIS', 'BRAVAIS', 'NAEZ', '<RBASIS>', 'NSPIN', 'LMAX', 'RMAX', 'GMAX']
        
        if 'NPOL' in self.values.keys() and self.values['NPOL'] != 0:
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
        
        
    def _check_mandatory(self):
        self._update_mandatory()
        for key in self.values.keys():
            if self._mandatory[key] and self.values[key] is None:
                print('Error not all mandatory keys are set!')
                raise ValueError
                
                
    def _check_array_consistency(self):
        """Check all keys in __listargs if they match their specification (mostly 1D array, except for special cases e.g. <RBASIS>)"""
        from numpy import array
        
        success = False
        
        for key in self.__listargs.keys():
            if self.values[key] is not None:
                cmpdims = (self.__listargs[key], )
                if key == '<RBASIS>':
                    cmpdims = (self.__listargs[key], 3)
                    # automatically convert if naez==1 and only 1D array is given
                    if self.values['NAEZ'] == 1 and len(array(self.values['<RBASIS>']).shape) == 1:
                        print('Warning: expected 2D array for <RBASIS> but got 1D array, converting automatically')
                        self.values['<RBASIS>'] = array([self.values['<RBASIS>']])
                tmpdims = array(self.values[key]).shape
                if tmpdims[0] == cmpdims[0]:
                    success = True
                if len(tmpdims)==2:
                    if tmpdims[1] != cmpdims[1]:
                        success = False
                        
        if not success:
            print('Error: array input not consistent')
            raise 
    
    
    def _check_input_consistency(self):
        """Check consistency of input, to be done before wrinting to inputcard"""
        
        # check if all mandatory values are there
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
                         ['<RMTREFL>', nlbasis], ['<RMTREFR>', nrbasis], ['<FPRADIUS>', natyp], ['BZDIVIDE', 3]])
        special_formatting = ['BRAVAIS', 'RUNOPT', 'TESTOPT', 'FILES']
        
        self.__listargs = listargs
        self.__special_formatting = special_formatting
        
        # check for consistency of array arguments
        self._check_array_consistency()
        
        
    def fill_keywords_to_inputfile(self, output='inputcard'):
        """Fill new inputcard with keywords/values and rest of template"""
        keywords = self.values
        keyfmts = self.__format
        
        self._check_input_consistency()
        
        tmpl = ''
        for key in keywords:
            if keywords[key] is not None:
                print(key)
                if (not key in self.__listargs.keys()) and (not key in self.__special_formatting):
                    try:
                        repltxt = keyfmts[key]%(keywords[key])
                    except:
                        print(key, keyfmts[key], keywords[key])
                        repltxt = ''
                        for i in range(len(keyfmts[key])):
                            repltxt += ' ' + keyfmts[key][i]%(keywords[key][i])
                    tmpl += '%s= %s\n'%(key, repltxt)
                elif key == 'BRAVAIS':
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
                            print('WARNING for replacement of RUNOPTION %s: too long?'%(repltxt))
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
                            print('WARNING for replacement of TESTOPTION %s: too long?'%(repltxt))
                        else:
                            repltxt = repltxt+' '*nblanks
                        tmpl += repltxt
                        if iop==8:
                            tmpl += '\n'
                    tmpl += '\n'
                elif key == 'FILES':
                    print('Warning: Changing file name of potential file.')
                    tmpl += 'FILES\n'
                    tmpl += '\n'
                    tmpl += '%\n'%self.values[key]
                    tmpl += '\n'
                    tmpl += '\n'
                elif key in self.__listargs.keys():
                    tmpl += '%s\n'%key
                    if key == '<RBASIS>': # RBASIS needs special formatting since three numbers are filled per line
                        print(key, self.__listargs[key], self.__format[key], self.values[key])
                        for ival in range(self.__listargs[key]):
                            tmpl += (self.__format[key]+'\n')%(self.values[key][ival][0], self.values[key][ival][1], self.values[key][ival][2])
                    else:
                        for ival in range(self.__listargs[key]):
                            tmpl += (self.__format[key]+'\n')%(self.values[key][ival])
                else:
                    print('Error trying to write keyword %s but writing failed!'%key)
                    raise ValueError
    
        open(output, 'w').write(tmpl)
        
        
    def read_keywords_from_inputcard(self, keywords, inputcard='inputcard'):
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
                retrvaltxt = self._get_value(searchkey, txt, item=ishift, num=keyvallen)
                #print(searchkey, ishift, keyvallen,  retrvaltxt, keywords[key][0], keywords[key][1][-1])
                if keyvallen == 1:
                    keywords[key][0] = self.convert_string2number(retrvaltxt, keywords[key][1][-1])
                else:
                    for i in range(len(keywords[key][0])):
                        #print(i,retrvaltxt[i])
                        keywords[key][0][i] = self.convert_string2number(retrvaltxt[i], keywords[key][1][i][-1])
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
    
        
    def _convert_string2number(self, str, type):
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
            
            
class voroparams(kkrparams):
    
    # redefine _update_mandatory for voronoi code    
    def _update_mandatory(self):
        # change mandatory flags to 
        for key in self.values.keys():
            self._mandatory[key] = False
        
        runopts = []
        if self.values['RUNOPT'] is not None:
            for runopt in self.values['RUNOPT']:
                runopts.append(runopt.strip())
        
        #For a KKR calculation these keywords are always mandatory:
        mandatory_list = ['ALATBASIS', 'BRAVAIS', 'NAEZ', '<RBASIS>', 'NSPIN', 'LMAX', 'RCLUSTZ']
        
        #Mandatory in 2D
        if self.values['INTERFACE']:
            mandatory_list += ['<NLBASIS>', '<RBLEFT>', 'ZPERIODL', '<NRBASIS>', '<RBRIGHT>', 'ZPERIODR']
        #Mandatory in CPA
        if self.values['NATYP'] is not None and self.values['NATYP'] > self.values['NAEZ']:
            mandatory_list += ['NATYP', '<SITE>', '<CPA-CONC>']
            
        for key in mandatory_list:
            self._mandatory[key] = True

         
if __name__=='__main__':
    p = kkrparams()
    print('mandatory keys:', p.get_all_mandatory())
    from numpy import array
    p.set_multiple_values(NAEZ=1, RBASIS=array([0,0,0]), ALATBASIS=6.83, NSPIN=1, GMAX=65, RMAX=7, LMAX=2, EMIN=-0.4, BRAVAIS=array([[0.5,0.5,0],[0.5,0.0,0.5],[0.0,0.5,0.5]]))
    p.set_value('RCLUSTZ', 1.5)
    p.fill_keywords_to_inputfile()
    #ln -s /Users/ruess/sourcecodes/voronoi/voronoi.exe
    #ln -s /Users/ruess/sourcecodes/voronoi/ElementDataBase/
    #from subprocess import check_output
    #check_output('source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64; ./voronoi.exe | tee out_voro', shell=True)
    #!ln -s /Users/ruess/sourcecodes/KKRcode/kkr.x
    #!cp output.pot potential
    #!source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64; ./kkr.x