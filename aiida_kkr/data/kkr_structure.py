"""
Here we implement a Structure Data type that is specific to the KKR code
It can store empty spheres, i.e can store 'dummy atom' with charge 0
"""
from aiida.orm.data import StructureData

class KkrstructureData(StructureData):
    """
    Extention of AiiDA StructureData type that can store empty spheres/dummy atoms.
    it also provides some methods to convert to and from an aiida StructreData node.
    """
    pass


