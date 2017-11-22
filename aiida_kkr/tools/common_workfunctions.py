#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 10:06:23 2017

@author: ruess
"""


from aiida.common.exceptions import InputValidationError
from aiida.work import workfunction as wf
from aiida.orm import DataFactory
from aiida_kkr.tools.kkr_params import kkrparams

#define aiida structures from DataFactory of aiida
ParameterData = DataFactory('parameter')

@wf
def update_params_wf(parameternode, updatenode):
    """
    Work function to update a KKR input parameter node.
    Stores new node in database and creates a link from old parameter node to new node
    Returns updated parameter node using update_params function

    :note: Input nodes need to be valid aiida ParameterData objects.
    
    :param parameternode: Input aiida ParameterData node cotaining KKR specific parameters
    :param updatenode: Input aiida ParameterData node containing a dictionary with the parameters that are supposed to be changed.
    
    :note: If 'nodename' is contained in dict of updatenode the string corresponding to this key will be used as nodename for the new node. Otherwise a default name is used
    :note: Similar for 'nodedesc' which gives new node a description
        
    :example: updated_params = ParameterData(dict={'nodename': 'my_changed_name', 'nodedesc': 'My description text', 'EMIN': -1, 'RMAX': 10.})
              new_params_node = update_params_wf(input_node, updated_params)
    """
    updatenode_dict = updatenode.get_dict()
    if 'nodename' in updatenode_dict.keys():
        # take nodename out of dict (should only contain valid KKR parameter)
        nodename = updatenode_dict.pop('nodename')
    else:
        nodename = None
    if 'nodedesc' in updatenode_dict.keys():
        # take nodename out of dict (should only contain valid KKR parameter)
        nodedesc = updatenode_dict.pop('nodedesc')
    else:
        nodedesc = None
    
    # do nothing if updatenode is empty
    if len(updatenode_dict.keys())==0:
        print('Input node is empty, do nothing!')
        raise InputValidationError('Nothing to store in input') 
    # 
    new_parameternode = update_params(parameternode, nodename=nodename, **updatenode_dict)
    return new_parameternode
    
        
def update_params(node, nodename=None, **kwargs):
    """
    Update parameter node given with the values given as kwargs.
    Returns new node.
    
    :param node: Input parameter node (needs to be valid KKR input parameter node).
    :param **kwargs: Input keys with values as in kkrparams.
    :param linkname: Input linkname string. Give link from old to new node a name . 
                     If no linkname is given linkname defaults to 'updated parameters'
                     
    :return: parameter node
    
    :example usage: OutputNode = KkrCalculation.update_params(InputNode, EMIN=-1, NSTEPS=30)
    
    :note: Keys are set as in kkrparams class. Check documentation of kkrparams for further information.
    """    
    # check if node is a valid KKR parameters node
    if not isinstance(node, ParameterData):
        print('Input node is not a valid ParameterData node')
        raise InputValidationError('update_params needs valid parameter node as input')
    
    print('input kwargs', kwargs)
    
    #initialize temporary kkrparams instance containing all possible KKR parameters
    params = kkrparams()
    
    # extract input dict from node
    inp_params = node.get_dict()
    
    # check if input dict contains only values for KKR parameters
    for key in inp_params:
        if key not in params.values.keys():
            print('Input node contains unvalid key "{}"'.format(key))
            raise InputValidationError('unvalid key "{}" in input parameter node'.format(key))
    
    # copy values from input node
    for key in inp_params:
        value = inp_params[key]
        if value is not None:
            params.set_value(key, value)
            
    # check if values are given as **kwargs (otherwise return input node)
    if len(kwargs)==0:
        print('No additional input keys given, return input node')
        return node
    else:
        for key in kwargs:
            params.set_value(key, kwargs[key])
            
    # set linkname with input or default value
    if nodename is None or type(nodename) is not str:
        nodename = 'updated KKR parameters'
        
    # create new node
    ParaNode = ParameterData(dict=params.values)
    ParaNode.label = nodename
    
    return ParaNode