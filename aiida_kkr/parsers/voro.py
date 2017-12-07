#-*- coding: utf-8 -*-

from aiida.parsers.parser import Parser
from aiida.orm.data.parameter import ParameterData
from aiida_kkr.tools.voronoi_helper import check_voronoi_output
from aiida_kkr.tools.common_functions import search_string, get_version_info
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.exceptions import InputValidationError


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1"
__contributors__ = ("Jens Broeder", "Philipp Rüßmann")


class VoronoiParser(Parser):
    """
    Parser class for parsing output of voronoi code..
    """

    def __init__(self, calc):
        """
        Initialize the instance of Voronoi_Parser
        """
        # check for valid input
        if not isinstance(calc, VoronoiCalculation):
            raise InputValidationError("Input calc must be a Voronoi Calculation")

        # these files should be present after success of voronoi
        self._default_files = {'outfile': calc._OUTPUT_FILE_NAME, 
                               'atominfo': calc._ATOMINFO, 
                               'radii': calc._RADII}
        
        self._ParserVersion = __version__

        #reuse init of base class
        super(VoronoiParser, self).__init__(calc)

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
        
        #Parse voronoi output, results that are stored in database are in out_dict
        out_dict = {'ParserVersion': self._ParserVersion}
        
        try:
            outfile = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
            code_version, compile_options, serial_number = get_version_info(outfile)
            out_dict['Code_version'] = code_version
            out_dict['Compile_options'] = compile_options
            out_dict['Calculation_serial_number'] = serial_number
        except:
            self.logger.error("Error parsing output of voronoi: Version Info")
            return success, self._get_nodelist(out_dict)
        
        try:
            outfile = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
            potfile = out_folder.get_abs_path(self._calc._OUT_POTENTIAL_voronoi)
            emin = check_voronoi_output(potfile, outfile)
            out_dict['emin'] = emin
            out_dict['units_EMIN'] = 'Ry'
        except:
            self.logger.error("Error parsing output of voronoi: 'EMIN'")
            return success, self._get_nodelist(out_dict)
        
        try:
            Ncls, results = self._get_cls_info(out_folder)
            out_dict['Cluster_number'] = Ncls
            tmpdict_all = []
            for icls in range(Ncls):
                tmpdict = {}
                tmpdict['iatom'] = results[icls][0]
                tmpdict['Refpot'] = results[icls][1]
                tmpdict['RMTref'] = results[icls][2]
                tmpdict['TB-cluster'] = results[icls][3]
                tmpdict['sites'] = results[icls][4]
                tmpdict_all.append(tmpdict)
            out_dict['Cluster_info'] = tmpdict_all
        except:
            self.logger.error("Error parsing output of voronoi: Cluster Info")
            return success, self._get_nodelist(out_dict)
        
        try:
            out_dict['Start_from_jellium_potentials'] = self._startpot_jellium(out_folder)
        except:
            self.logger.error("Error parsing output of voronoi: Jellium startpot")
            return success, self._get_nodelist(out_dict)
        
        try:
            natyp, naez, shapes = self._get_shape_array(out_folder)
            out_dict['shapes'] = shapes
        except:
            self.logger.error("Error parsing output of voronoi: SHAPE Info")
            return success, self._get_nodelist(out_dict)
        
        try:
            Vtot, results = self._get_volumes(out_folder)
            out_dict['Volume_total'] = Vtot
            tmpdict_all = []
            for icls in range(naez):
                tmpdict = {}
                tmpdict['iatom'] = results[icls][0]
                tmpdict['V_atom'] = results[icls][1]
                tmpdict_all.append(tmpdict)
            out_dict['Volume_atoms'] = tmpdict_all
            out_dict['units_Volume'] = 'alat^3'
        except:
            self.logger.error("Error parsing output of voronoi: Volume Info")
            return success, self._get_nodelist(out_dict)
        
        try:
            results = self._get_radii(naez, out_folder)
            tmpdict_all = []
            for icls in range(naez):
                tmpdict = {}
                tmpdict['iatom'] = results[icls][0]
                tmpdict['Rmt0'] = results[icls][1]
                tmpdict['Rout'] = results[icls][2]
                tmpdict['dist_NN'] = results[icls][4]
                tmpdict['Rmt0/Rout'] = results[icls][3]
                tmpdict['Rout/dist_NN'] = results[icls][5]
                tmpdict_all.append(tmpdict)
            out_dict['radii_atoms'] = tmpdict_all
            out_dict['units_radii'] = 'alat'
        except:
            self.logger.error("Error parsing output of voronoi: radii.dat Info")
            return success, self._get_nodelist(out_dict)
        
        
        # some consistency checks comparing lists with natyp/naez numbers
        #TODO implement checks

        success = True

        return success, self._get_nodelist(out_dict)
    
    # here follow the parser functions:
    
    def _get_nodelist(self, out_dict):
        output_data = ParameterData(dict=out_dict)
        link_name = self.get_linkname_outparams()
        node_list = [(link_name, output_data)]
        return node_list
    
    
    def _startpot_jellium(self, out_folder):
        f = open(out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME))
        tmptxt = f.readlines()
        f.close()
        itmp = search_string('JELLSTART POTENTIALS', tmptxt)
        if itmp ==-1:
            return False
        else:
            return True
    
    
    def _get_volumes(self, out_folder):
        f = open(out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME))
        tmptxt = f.readlines()
        f.close()
        
        itmp = search_string('Total volume (alat^3)', tmptxt)
        if itmp>=0:
            Vtot = float(tmptxt.pop(itmp).split()[-1])
        
        itmp = 0
        results = []
        while itmp>=0:
            itmp = search_string(' Volume(alat^3)  :', tmptxt)
            if itmp>=0:
                tmpstr = tmptxt.pop(itmp)
                tmpstr = tmpstr.split()
                tmpstr = [int(tmpstr[2]), float(tmpstr[5])]
                results.append(tmpstr)
        return Vtot, results
    
    
    def _get_cls_info(self, out_folder):
        f = open(out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME))
        tmptxt = f.readlines()
        f.close()
        itmp = 0
        Ncls = 0
        results = []
        while itmp>=0:
            itmp = search_string('CLSGEN_TB: Atom', tmptxt)
            if itmp>=0:
                tmpstr = tmptxt.pop(itmp)
                tmpstr = tmpstr.split()
                tmpstr = [int(tmpstr[2]), int(tmpstr[4]), float(tmpstr[6]), int(tmpstr[8]), int(tmpstr[10])]
                results.append(tmpstr)
                Ncls += 1
        return Ncls, results
    
    
    def _get_shape_array(self, out_folder):
        f = open(out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME))
        txt = f.readlines()
        f.close()
        #naez/natyp number of items either one number (=ishape without cpa or two =[iatom, ishape] with CPA)
        # read in naez and/or natyp and then find ishape array (1..natyp[=naez without CPA])
        itmp = search_string('NAEZ= ', txt)
        if itmp>=0:
            tmp = txt[itmp]
            ipos = tmp.find('NAEZ=')
            naez = int(tmp[ipos+5:].split()[0])
        else:
            naez = -1
        itmp = search_string('NATYP= ', txt)
        if itmp>=0:
            tmp = txt[itmp]
            ipos = tmp.find('NATYP=')
            natyp = int(tmp[ipos+6:].split()[0])
        else:
            natyp = -1
            
        # consistency check
        if naez==-1 and natyp>0:
            naez = natyp
        elif natyp==-1 and naez>0:
            natyp = naez
        elif natyp==-1 and naez==-1:
            raise ValueError('Neither NAEZ nor NATYP found in %s'%self._calc._OUTPUT_FILE_NAME)
        
        # read shape index from atominfo file
        f = open(out_folder.get_abs_path(self._calc._ATOMINFO))
        tmptxt = f.readlines()
        f.close()
        
        itmp = search_string('<SHAPE>', tmptxt) + 1
        ishape = []
        for iatom in range(natyp):
            txt = tmptxt[itmp+iatom]
            if natyp>naez: #CPA option
                ishape.append(int(txt.split()[1]))
            else:
                ishape.append(int(txt.split()[0]))
        
        return natyp, naez, ishape
    
    
    def _get_radii(self, naez, out_folder):
        f = open(out_folder.get_abs_path(self._calc._RADII))
        txt = f.readlines()
        f.close()
        results = []
        for iatom in range(naez):
            # IAT    Rmt0           Rout            Ratio(%)   dist(NN)      Rout/dist(NN) (%)              
            # 1   0.5000001547   0.7071070000       70.71   1.0000003094       70.71
            tmpline = txt[3+iatom].split()
            tmpline = [int(tmpline[0]), float(tmpline[1]), float(tmpline[2]), float(tmpline[3]), float(tmpline[4]), float(tmpline[5])]
            results.append(tmpline)
        return results
    
