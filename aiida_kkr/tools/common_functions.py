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