"""
Here we implement a Structure Data type that is specific to the KKR code
It can store empty spheres, i.e can store 'dummy atom' with charge 0
"""
#from aiida.orm.data.structure import *
from __future__ import absolute_import
from aiida.orm import Data
from aiida.plugins import DataFactory

StructureData = DataFactory('structure')

class KkrstructureData(StructureData):
    """
    Extention of AiiDA StructureData type that can store empty spheres/dummy atoms.
    it also provides some methods to convert to and from an aiida StructreData node.
    """

    def __init__(self, *args, **kwargs):
        super(KkrstructureData, self).__init__(*args, **kwargs)

    pass
