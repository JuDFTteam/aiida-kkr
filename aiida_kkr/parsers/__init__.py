"""
Parsers provided by the aiida-kkr plugin

Register parsers via the "aiida.parsers" entry point in setup.json.
"""
# import all calculations here to expose them in `aiida_kkr.parsers` directly
from .voro import VoronoiParser
from .kkr import KkrParser
from .kkrimp import KkrimpParser
# broken at the moment
#from .kkrimporter import KkrImporterParser
