"""
Calculations provided by plugin

Register calculations via the "aiida.calculations" entry point in setup.json.
"""

# import all calculations here to expose them in `aiida_kkr.calculations` directly
from .voro import VoronoiCalculation
from .kkr import KkrCalculation
from .kkrimp import KkrimpCalculation
# broken at the moment:
#from .kkrimporter import KkrimporterCalculation
