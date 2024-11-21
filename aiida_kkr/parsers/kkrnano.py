#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 14:38:59 2021

@author: markus
"""

# -*- coding: utf-8 -*-

from aiida.parsers.parser import Parser
from aiida.orm import Dict, CalcJobNode
from aiida_kkr.calculations.kkrnano import KKRnanoCalculation
from masci_tools.io.common_functions import (search_string, open_general)
import numpy as np
from io import StringIO
from pprint import pprint as pp
import os

__copyright__ = (u'Copyright (c), 2021, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.0.2'
__contributors__ = ('Markus Struckmann', 'Philipp Rüßmann')


class KKRnanoParser(Parser):
    """
    Parser class for parsing output of the KKRnano code
    """

    def __init__(self, calc):
        """
        Initialize
        """
        # these files should be present after successful run of KKRnano
        #         self._default_files = {
        #             #'stdout': KKRnanoCalculation._DEFAULT_OUTPUT_FILE
        #            #'stdout_prep': KKRnanoCalculation._DEFAULT_OUTPUT_PREP_FILE
        #         }

        self._ParserVersion = __version__

        # reuse init of base class
        super(KKRnanoParser, self).__init__(calc)

    # pylint: disable=protected-access

    def _get_lines(self, retrieved_folder, output_file_name):
        """returns list of string lines"""
        with retrieved_folder.open(output_file_name, 'r') as f:
            lines = f.readlines()
        return lines

    def _findSimpleEntries(self, string2find, retrieved_folder, output_file_name, lineindices=[-1], simpleEntry=True):
        """
        read out entries that are simply given at the end of a line preceeded by the string2find.
        returns a list of said entries
        """
        lines = self._get_lines(retrieved_folder, output_file_name)

        returnlist, indexlist = [], []

        #checking if line indices were passed. If so, only the specified lines are searched
        if lineindices[0] == -1:
            lineindices = range(len(lines))

        #looping over indicated or all (default) lines
        for j in lineindices:
            line = lines[j]
            stringposition = line.find(string2find)

            #matching sought string
            if stringposition >= 0:
                indexlist.append(j)
                keyvalue = line[stringposition + len(string2find) + 1:-1]  #.replace(" ", "")
                keyvalue = keyvalue.replace('=', '')
                keyvalue = keyvalue.split(sep='(')[0]

                if simpleEntry:
                    array = np.genfromtxt(StringIO(keyvalue.replace('D', 'e')), delimiter=' ', dtype=None)  #[0]

                try:
                    returnlist.append(array.item()[0])  # pylint: disable=possibly-used-before-assignment
                except TypeError:
                    returnlist.append(array.item())

        return returnlist, indexlist

    def _get_index_list(self, string2find, retrieved_folder, output_file_name='', lines=[]):
        """
        returns list of indicies of lines containting the passed string.
        opens file, if a file name and no string is passed.
        """

        if lines == [] and output_file_name == '':
            print('ERROR in get_index_list: Neither file name nor lines of strings were passed')
        if lines == []:
            lines = self._get_lines(retrieved_folder, output_file_name)

        returnlist, indexlist = [], []

        for j in range(len(lines)):
            line = lines[j]
            stringposition = line.find(string2find)
            if stringposition >= 0:
                indexlist.append(j)
        return indexlist

    def _read_table_block(self, lines, retrieved_folder, output0_file_name, index_multiple_tables=-1):
        """
        reads a table in the output of the --prepare step which is indicated by a borderlines of "---", returns a string array of said table.
        """
        indexlist = self._get_index_list('---------', retrieved_folder, output0_file_name, lines=lines)
        if index_multiple_tables == -1:
            table = lines[indexlist[-2] + 1:indexlist[-1]]
        else:
            table = lines[indexlist[index_multiple_tables * 2] + 1:indexlist[index_multiple_tables * 2 + 1]]
        array = [line.replace('\n', '').split(' ') for line in table]
        for j in range(len(array)):
            for k in range(array[j].count('')):
                array[j].remove('')

        #array=[np.delete(a,np.where(a=="")) for a in array]
        return np.array(array)

    def _find_block(self, lines, string, retrieved_folder, output0_file_name):
        """
        finds a block which contains the indicated string in passed lines.
        Can be used for the output of the --prepare step which is indicated by a borderlines of "===", returns a list of the lineindices in the passed lines.
        """
        pos__output = search_string(string, lines)
        indexlist = self._get_index_list('=======', retrieved_folder, output0_file_name, lines=lines)
        for i in range(len(indexlist) - 1):
            if indexlist[i] < pos__output and indexlist[i + 1] > pos__output:
                return lines[indexlist[i]:indexlist[i + 1]]
        print(f"Warning: Block '{string}' not found!")
        return []

    def _get_total_EnergyeV(self, key, retrieved_folder, output_file_name, stringsInOutputFile):
        return_dict = {}
        _, indexlist = self._findSimpleEntries(stringsInOutputFile[key], retrieved_folder, output_file_name)
        #print(findSimpleEntries("eV  :",output_file_name,retrieved_folder) np.array(indexlist)+1))
        return_dict['total_energy_in_eV'], _ = self._findSimpleEntries(
            'eV  :', retrieved_folder, output_file_name,
            np.array(indexlist) + 1
        )
        return return_dict

    def _get_charge_in_WScell(self, key, retrieved_folder, output_file_name, n_atoms):
        number_of_atoms = 5
        charges, indexlist = self._findSimpleEntries(key, retrieved_folder, output_file_name)
        charge_dict = {}
        charge_dict['atom'] = {}
        num_atoms = n_atoms

        for k in range(num_atoms):
            charge_dict['atom'][k + 1] = [charges[i] for i in np.arange(0 + k, len(charges) + k, num_atoms)]
        return charge_dict

    def _stringFromList(self, stringlist):
        """
        turn a list of strings into a single string
        """
        finalstring = ''
        for string in stringlist:
            finalstring += string
        return finalstring

    def _dict_from_table(self, captions_columns, captions_lines, array):
        """
        returns a dictionary from a table using specified captions for columns and lines
        """
        return_dict = {}
        for caption_index in range(len(captions_columns)):
            return_dict[captions_columns[caption_index]] = dict(zip(captions_lines, array[:, caption_index]))
        return return_dict

    def _extract_l_valence_charges(self, lines, captions):
        """
        extract the l-decomposed valence charges charges from the output file of a KKRnano run:
        reads the used captions and turns the used tables into a dict
        """
        orbitals = [
            's', 'p', 'd', 'f', 'dummy'
        ]  #KKRnano does not write out "orbitals" beyond this, dummy to simplify accounting for varying lengths
        return_dict = {}

        length_top = search_string(
            '----', lines
        )  #this delimiter marks where only total values follow in the output file, varies w.r.t. used LMAX
        length_bottom = len(lines) - length_top

        #after col 21, only floats follow, that are sought to be read in in the following
        blockstring = self._stringFromList([line[21:] for line in lines[:-length_bottom]])

        used_orbitals = orbitals[:length_top]
        used_orbitals[-1] = 'non-spherical'  #last entry is always the non-spherical part

        array = np.genfromtxt(StringIO(blockstring), dtype=float)
        if array.ndim == 1:
            array = np.transpose(
                np.atleast_2d(array)
            )  #make sure that a 2D array is processed in the following, as this would raise an error otherwise
        # Single String for line where the total values are indicated
        totalstring = self._stringFromList(lines[-length_bottom + 1])

        return_dict = self._dict_from_table(captions, used_orbitals, array)

        totalvalues = np.genfromtxt(StringIO(totalstring[21:]), dtype=float)  #again after col 21, only floats follow

        if np.shape(totalvalues) == ():  #accounting for 0D array
            return_dict['total'] = dict(zip(captions, [totalvalues.item()]))
        else:
            return_dict['total'] = dict(zip(captions, np.genfromtxt(StringIO(totalstring[21:]), dtype=float)))

        return return_dict

    def _convert_all_values_2_python_natives(self, sub_dictionary):
        """
        convert all dtypes of numpy to python native data types in a given (potentially nested) dictionary.
        """
        for key, value in sub_dictionary.items():

            if type(value) is dict:
                self._convert_all_values_2_python_natives(value)

            else:
                try:
                    sub_dictionary[key] = value.item()
                except:
                    pass
        return sub_dictionary

    def _identifyBlocks(self, lines):
        """identify the DOS blocks in the respective output files"""
        blockstartlist = [-1]
        index = 0
        for line in lines:

            block_pos = line.find('&')  #find blocks in DOS file
            if block_pos >= 0:
                blockstartlist.append(index)
            index += 1
        blockstartlist.append(len(lines))  #add EOF line index

        #print(blockstartlist)
        return blockstartlist

    def _get_commentlessLineIndices(self, lines):
        nonemptylines = []  #index list
        index = 0
        for line in lines:
            com_pos = line.find('#')  #comment sign position

            if com_pos < 0:
                nonemptylines.append(index)
            index += 1
        return nonemptylines

    def _process_DOS_file(self, retrieved_folder, filename):
        """reading in DOS file output from KKRnano"""
        data = np.array(self._get_lines(retrieved_folder, filename))
        element = data[0][1:3]

        blocklist = self._identifyBlocks(data)
        dict_dos = {}
        multipleSpins = False
        if len(blocklist) > 1:
            multipleSpins = True
            dict_dos['spin_directions'] = 2
            first_spin = data[0].split(sep='SPIN')[1].split(sep=' ')[1]
            if first_spin == 'DOWN':
                spin_list = ['spin_down', 'spin_up']
            elif first_spin == 'UP':
                spin_list = ['spin_up', 'spin_down']
            else:
                dict_dos['spin_directions'] = 1
        #DOS_blocks=[]

        DOS_captions = ['energy_in_ryd', 's', 'p', 'd', 'f', 'non-spherical', 'total_DOS']
        for i in range(len(blocklist) - 1):
            datablockindices=np.array(self._get_commentlessLineIndices( \
                    data[blocklist[i]+1:blocklist[i+1]-1]))+blocklist[i]+1
            datablockcontent = self._stringFromList(data[datablockindices])
            #print(datablockcontent)
            DOS_block = np.genfromtxt(StringIO(datablockcontent))

            block_dict = {}
            final_DOS_captions = DOS_captions[:np.shape(DOS_block)[1] - 2]
            final_DOS_captions.append(DOS_captions[-2])
            final_DOS_captions.append(DOS_captions[-1])

            for p in range(np.shape(DOS_block)[1]):
                block_dict[final_DOS_captions[p]] = DOS_block[:, p]

            if multipleSpins:
                dict_dos[spin_list[i]] = block_dict  # pylint: disable=used-before-assignment
            else:
                dict_dos = block_dict
            #DOS_blocks.append(DOS_block)
        return element, dict_dos  #DOS_blocks

    def parse(self, debug=False, **kwargs):
        """
        Parse output data folder, store results in database.

        :param retrieved: a dictionary of retrieved nodes, where
          the key is the link name
        :returns: nothing if everything is fine or an exit code defined in the voronoi calculation class
        """

        success = True
        node_list = ()

        # Get retrieved folders
        try:
            retrieved_folder = self.retrieved
            print(retrieved_folder)
        except:
            print('OUT FOLDER NOT FOUND')
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = retrieved_folder.list_object_names()

        calc_node = retrieved_folder.get_incoming(node_class=CalcJobNode).first().node

        #Check if a StrucWithPotData object was used as input and take then corresponding structure or
        #if a parent_structure has to be found
        #         if hasattr(calc_node.inputs, 'strucwithpot') and not hasattr(calc_node.inputs, 'parent_folder'):
        #             struc = calc_node.inputs.strucwithpot.structure
        #         else:
        #             struc = find_parent_structure(calc_node)

        #Also for the convert step, these are the files that are supposed to be parsed
        output0_file_name = KKRnanoCalculation._DEFAULT_OUTPUT_PREP_FILE
        output_file_name = KKRnanoCalculation._DEFAULT_OUTPUT_FILE

        # initialize out_dict and parse output files
        out_dict_final = {'parser_version': self._ParserVersion}

        #from --prepare output file
        out0_dict = {}
        #"""
        lines0 = self._get_lines(retrieved_folder, output0_file_name)

        # number of atoms
        try:
            out0_dict['num_atoms'] = int(lines0[search_string('atoms in rbasis.xyz', lines0)].split(sep=' ')[1])

        except:
            print('Number of atoms was not read!')

        # lattice constant

        try:
            alat_string = lines0[search_string('Lattice constants :  ALAT',
                                               lines0)].split(sep=' = ')[1].split(sep='   ')[0].replace(' ', '')
            out0_dict['alat_internal'] = float(alat_string)
            pass
        except:
            print('ALAT was not read!')
        out0_dict['alat_internal_unit'] = 'a_Bohr'

        # Reading k-mesh details

        kmesh_dict = {}
        #kmesh_dict['number_different_kmeshes']=int(lines0[search_string("number of different k-meshes",lines0)].split(sep=" : ")[1])

        try:
            kmesh_dict['number_different_kmeshes'] = int(
                lines0[search_string('number of different k-meshes', lines0)].split(sep=' : ')[1]
            )

        except:
            print('Number of diff. k-meshes was not read!')

        kmesh_caption_KKRnano = 'k-mesh NofKs N kx N ky N kz vol BZ'
        kmesh_caption_aiida_KKRhost = ['number_of_kpts', 'n_kx', 'n_ky', 'n_kz']

        table = self._read_table_block(
            self._find_block(lines0, kmesh_caption_KKRnano, retrieved_folder, output0_file_name), retrieved_folder,
            output0_file_name
        )

        kmesh_dict['number_kpoints_per_kmesh'] = {}

        #filling dictionary with retrieved data
        try:
            for keyindex in range(len(kmesh_caption_aiida_KKRhost)):
                key = kmesh_caption_aiida_KKRhost[keyindex]
                kmesh_dict['number_kpoints_per_kmesh'][key] = np.array(table[:, 1 + keyindex], dtype=int)

            out0_dict['kmesh_group'] = kmesh_dict
        except:
            print('Number of different k-meshes was not read due to unusal format.')

        # reciprocal Bravais matrix
        bravais_caption = 'Reciprocal lattice cell vectors'
        table = self._read_table_block(
            self._find_block([line[:45] for line in lines0], bravais_caption, retrieved_folder, output0_file_name),
            retrieved_folder, output0_file_name
        )

        table = np.array(table[:, 1:], dtype=float)
        out0_dict['reciprocal_bravais_matrix'] = table[:, :3]
        out0_dict['reciprocal_bravais_matrix_unit'] = '2*pi / alat'

        # Bravais matrix
        bravais_caption = 'Direct lattice cell vectors'
        table = self._read_table_block(
            self._find_block([line[:45] for line in lines0], bravais_caption, retrieved_folder, output0_file_name),
            retrieved_folder,
            output0_file_name,
            index_multiple_tables=0
        )

        table = np.array(table[:, 1:], dtype=float)
        out0_dict['direct_bravais_matrix'] = table[:, :3]
        out0_dict['direct_bravais_matrix_unit'] = 'alat'
        #"""

        # read entries from the main output

        stringsInOutputFile = {
            'total_energy_in_ryd': 'TOTAL ENERGY in ryd. :',
            'rms_all_iterations': 'v+ + v-',
            'rms_minus_all_iterations': 'v+ - v-',
            'fermi_energy_in_ryd': 'Fermi energy =',
            'charge_neutrality_in_e': 'charge neutrality in unit cell =',
            'total_magn_moment_in_unit_cell': 'TOTAL mag. moment in unit cell ='
        }

        out_dict = {}
        for key in stringsInOutputFile:
            out_dict[key], _ = self._findSimpleEntries(stringsInOutputFile[key], retrieved_folder, output_file_name)

        out_dict = {
            **out_dict,
            **self._get_total_EnergyeV('total_energy_in_ryd', retrieved_folder, output_file_name, stringsInOutputFile)
        }

        # Get charges in WS cell
        dict_WScell_keys = {
            'charge_in_e': 'charge in Wigner Seitz cell =',
            'spin_moment': 'spin moment in Wigner Seitz cell =',
            'nuclear_charge_in_e': 'nuclear charge',
            'core_charge_in_e': 'core charge'
        }
        WScell_dict = {}
        for key in dict_WScell_keys:
            WScell_dict[key] = self._get_charge_in_WScell(
                dict_WScell_keys[key], retrieved_folder, output_file_name, out0_dict['num_atoms']
            )

        # Extract the l-decomposed valence charges information for all "orbitals"
        # and all iterations and add them to the dict
        lines = self._get_lines(retrieved_folder, output_file_name)

        #find block where the valence charges are indicated
        string2find = 'l-decomposed valence charges'
        indexlist = np.array(self._get_index_list(string2find, retrieved_folder, output_file_name, lines)) + 1

        #using a subdictionary to store the information
        dict_orbitals = {}
        for m in range(len(indexlist)):
            #find iteration block to process
            index = indexlist[m]
            nextindex = indexlist[(m + 1) % len(indexlist)]
            if nextindex <= index:
                nextindex = -1  #use EOF as nextindex, if necessary

            #read in the captions of the table and convert them to a format that is easier to process
            captions_table = np.genfromtxt(StringIO(lines[index + 2]), dtype=str, delimiter='  ')
            captions_table = np.delete(captions_table, np.where(captions_table == ''))
            captions_table = [item.replace(' ', '_') for item in captions_table][1:]

            # removing leading "_" (occurs for some spins)
            for s in range(len(captions_table)):
                if captions_table[s].find('_') == 0:
                    captions_table[s] = captions_table[s][1:]

            #find string list with the lines to process
            blockend = index + 2 + search_string('#########', lines[index + 2:nextindex])
            blocklines = lines[index + 3:blockend]

            #             print("blocklines",blocklines)
            #             print("indexlist", indexlist)
            #             print("index", index)
            #             print("lines", lines[:5])
            #             print("blockend", blockend)
            #             print("nextindex",nextindex)

            atomblocks = np.array(
                self._get_index_list('===', retrieved_folder, lines=blocklines)
            ) + index + 3  #retrieved folder is actually not needed, but is passed for keeping it simple

            if len(atomblocks) > 1:
                atomblocklength = atomblocks[1] - atomblocks[0]
            else:
                atomblocklength = blockend - atomblocks[0]

            #loop over the atomblocks to read in the information for each atom

            dict_atoms = {}
            dict_atoms['atom'] = {}
            for j in range(len(atomblocks)):
                dict_atoms['atom'][j + 1] = self._extract_l_valence_charges(
                    lines[atomblocks[0] + j * atomblocklength + 1:atomblocks[0] + (j + 1) * atomblocklength + 1],
                    captions_table
                )  #all blocks should have the same length

            #add the atom dic-tionary to the one for the iterations
            dict_orbitals[m + 1] = dict_atoms

        # identify DOS-files
        list_of_DOS_files = []
        for file in list_of_files:
            if file.find('DOS') > -1:
                list_of_DOS_files.append(file)

        #process DOS files if necessary
        out_dict_dos = {}
        if len(list_of_DOS_files) > 0:
            for file_index in range(len(list_of_DOS_files)):
                file = list_of_DOS_files[file_index]

                element, DOSblocks = self._process_DOS_file(retrieved_folder, file)
                atomname = f'atom {file_index + 1}'
                out_dict_dos[atomname] = DOSblocks
                out_dict_dos[atomname]['element'] = element

            out_dict_final['DOS'] = out_dict_dos

        dict_orbitals = {'iterations': dict_orbitals}

        out_dict_final['prepare'] = out0_dict
        out_dict_final['WS_charges'] = WScell_dict
        out_dict_final['l_decomposed_charges'] = dict_orbitals

        out_dict_final = {**out_dict, **out_dict_final}

        # create output node and link
        self.out('output_parameters', Dict(dict=out_dict_final))

        if not success:
            return self.exit_codes.ERROR_PARSING_FAILED
