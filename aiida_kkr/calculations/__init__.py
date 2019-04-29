"""
Calculations provided by plugin

Register calculations via the "aiida.calculations" entry point in setup.json.
"""

# import all calculations here to expose them in `aiida_kkr.calculations` directly
from .voro import VoronoiCalculation
from .kkr import KkrCalculation
# broken at the moment:
#from .kkrimp import KkrimpCalculation
#from .kkrimporter import KkrimporterCalculation
