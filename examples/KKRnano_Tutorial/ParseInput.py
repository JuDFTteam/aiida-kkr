#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from io import StringIO
import numpy as np


def f2p(string):
    string = str.replace(string, 'D', 'e')
    string = str.replace(string, 'd', 'e')
    return string


p2f = lambda string: str.replace(string, 'e', 'D')


def check_for_bool_value(string):
    if type(string) is str:
        if string == 't':
            return True
        elif string == 'f':
            return False
    else:
        return string


def unpack_from_arrays(array):
    """ change into python natives """
    if array.shape == ():
        item = array.item()
        if type(item) is bytes:
            item = item.decode('utf-8')
        return item
    else:
        return array


def remove_leading_Whitespace(linepart):
    while linepart.find(' ') == 0:
        linepart = linepart[1:]
    return linepart


def get_inputdict_nano(filename):
    """
    returns an inputdict from a KKRnano input.conf file

    filename: path for file to read in from
    """

    with open(filename, 'r') as reader:
        data = reader.readlines()
        reader.close()
    inputdict = {}
    for line in data:
        eq_pos = line.find('=')  #equal sign position
        hash_pos = line.find('#')  # Hash position
        if (eq_pos > 0 and hash_pos < 0) or (eq_pos > 0 and hash_pos > eq_pos):
            lineparts = line.split('=')
            parameter = lineparts[0].replace(' ', '')
            lineparts_nocomment = f2p(lineparts[1].split('#')[0].replace('\n', ''))

            lineparts_nocomment = remove_leading_Whitespace(lineparts_nocomment)
            value = np.genfromtxt(StringIO(lineparts_nocomment), dtype=None)  #read in all entries as strings
            value = unpack_from_arrays(value)
            value = check_for_bool_value(value)
            # print("{},{}".format(parameter,value))
            # print(type(value))
            #        if np.any(np.isnan(value)):
            #            value=lineparts_nocomment
            #value=lineparts_nocomment
            inputdict.update({parameter: value})

    protected_keys = ['alat', 'basisscale', 'bravais_a', 'bravais_b', 'bravais_c', 'cartesian']
    for key in protected_keys:
        try:
            del inputdict[key]
        except KeyError:
            pass

    inputdict_nano = {}
    for key in inputdict:
        inputdict_nano[key] = {'value': inputdict[key]}
    return inputdict_nano
