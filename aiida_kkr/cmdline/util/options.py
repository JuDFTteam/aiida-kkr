# -*- coding: utf-8 -*-
"""
Options commonly used throughout the aiida-kkr command line interface.
To standardize at least of the options, they are kept as close as possible to
aiida-fleur. aiida-core and aiida-quantumespresso
"""
import click
from aiida.cmdline.params import types
from aiida.cmdline.params.options import OverridableOption
from .types import StructureNodeOrFileParamType
from .defaults import get_voro, get_kkr, get_kkrimp, get_kkrpara_defaults

STRUCTURE_OR_FILE = OverridableOption(
    '-s',
    '--structure',
    type=StructureNodeOrFileParamType(),
    help='StructureData node, given by pk or uuid or file in any for mat which will be converted.'
)

STRUCTURE = OverridableOption(
    '-s',
    '--structure',
    type=types.DataParamType(sub_classes=('aiida.data:core.structure',)),
    help='StructureData node, given by pk or uuid.'
)

VORO = OverridableOption(
    '-v',
    '--voro',
    type=types.CodeParamType(entry_point='kkr.voro'),
    default=get_voro,
    show_default=True,
    help='A code node or label for an voronoi executable.'
)

KKR = OverridableOption(
    '-k',
    '--kkr',
    type=types.CodeParamType(entry_point='kkr.kkr'),
    default=get_kkr,
    show_default=True,
    help='A code node or label for a KKR executable.'
)

KKRIMP = OverridableOption(
    '-ki',
    '--kkr-imp',
    type=types.CodeParamType(entry_point='kkr.kkrimp'),
    default=get_kkrimp,
    show_default=True,
    help='A code node or label for a KKRimp executable.'
)

PARAMETERS = OverridableOption(
    '-p',
    '--parameters',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    default=get_kkrpara_defaults,
    show_default=True,
    help='Dict with calculation parameters which will be given to the calculation or workchain.'
)

KPOINTS = OverridableOption(
    '-kpt',
    '--kpoints',
    type=types.DataParamType(sub_classes=('aiida.data:core.array.kpoints',)),
    help='An aiida KpointsData node.'
)

POTENTIAL_OVERWRITE = OverridableOption(
    '--potential-overwrite',
    type=types.DataParamType(sub_classes=('aiida.data:core.singlefile',)),
    help='Use a node that specifies the potential which is used instead of the voronoi output potential'
)

IMPURITY_INFO = OverridableOption(
    '-i',
    '--impurity-info',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help=
    'Dict containing parameters that specify properties for a following impurity calculation (e.g. setting of impurity cluster in scoef file that is automatically created'
)

WF_PARAMETERS = OverridableOption(
    '-wf',
    '--wf-parameters',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help='Dict containing parameters given to the workchain.'
)

VORO_PARAMETERS = OverridableOption(
    '-vp',
    '--voro-parameters',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help='Dict containing parameters for the auxiliary voronoi starting potential workflow.'
)

OPTION_NODE = OverridableOption(
    '-opt',
    '--option-node',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help=
    'Dict, an option node for the workchain containing for instance: {"withmpi": False, "max_wallclock_seconds": 6000, "resources": {"num_machines": 1, "num_mpiprocs_per_machine": 1}}'
)

MAX_NUM_MACHINES = OverridableOption(
    '-N',
    '--max-num-machines',
    type=click.INT,
    default=1,
    show_default=True,
    help='The maximum number of machines (nodes) to use for the calculations.'
)

MAX_WALLCLOCK_SECONDS = OverridableOption(
    '-W',
    '--max-wallclock-seconds',
    type=click.INT,
    default=1800,
    show_default=True,
    help='The maximum wallclock time in seconds to set for the calculations.'
)

NUM_MPIPROCS_PER_MACHINE = OverridableOption(
    '-M',
    '--num-mpiprocs-per-machine',
    type=click.INT,
    default=12,
    show_default=True,
    help='Run the simulation with so many num-mpi-procs-per-machine.'
)

WITH_MPI = OverridableOption(
    '-I', '--with-mpi', is_flag=True, default=False, show_default=True, help='Run the calculations with MPI enabled.'
)

QUEUE_NAME = OverridableOption(
    '-Q',
    '--queue-name',
    default='',
    show_default=True,
    help='The queue_name to be used in the submission script of the job.'
)

PARENT_FOLDER = OverridableOption(
    '-P',
    '--parent-folder',
    'parent_folder',
    type=types.DataParamType(sub_classes=('aiida.data:core.remote',)),
    show_default=True,
    required=False,
    help='The PK of a parent remote folder (for restarts).'
)

DAEMON = OverridableOption(
    '-d',
    '--daemon',
    is_flag=True,
    default=False,
    show_default=True,
    help='Submit the process to the daemon instead of running it locally. -d flag does not need any argument'
)

WF_PARAMETERS = OverridableOption(
    '-wf',
    '--wf-parameters',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help='Dict containing parameters given to the workchain.'
)

VORO_PARAMETERS = OverridableOption(
    '-vp',
    '--voro-parameters',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help='Dict containing parameters for the auxiliary voronoi starting potential workflow.'
)

OPTION_NODE = OverridableOption(
    '-opt',
    '--option-node',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help='Dict, an option node for the workchain.'
)

MAX_NUM_MACHINES = OverridableOption(
    '-N',
    '--max-num-machines',
    type=click.INT,
    default=1,
    show_default=True,
    help='The maximum number of machines (nodes) to use for the calculations.'
)

MAX_WALLCLOCK_SECONDS = OverridableOption(
    '-W',
    '--max-wallclock-seconds',
    type=click.INT,
    default=1800,
    show_default=True,
    help='The maximum wallclock time in seconds to set for the calculations.'
)

NUM_MPIPROCS_PER_MACHINE = OverridableOption(
    '-M',
    '--num-mpiprocs-per-machine',
    type=click.INT,
    default=12,
    show_default=True,
    help='Run the simulation with so many num-mpi-procs-per-machine.'
)

WITH_MPI = OverridableOption(
    '-I', '--with-mpi', is_flag=True, default=False, show_default=True, help='Run the calculations with MPI enabled.'
)

PARENT_FOLDER = OverridableOption(
    '-P',
    '--parent-folder',
    'parent_folder',
    type=types.DataParamType(sub_classes=('aiida.data:core.remote',)),
    show_default=True,
    required=False,
    help='The PK of a parent remote folder (for restarts).'
)

DAEMON = OverridableOption(
    '-d',
    '--daemon',
    is_flag=True,
    default=False,
    show_default=True,
    help='Submit the process to the daemon instead of running it locally.'
)

NOCO_ANGLES = OverridableOption(
    '--noco-angles',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    help='Dict containing the non-collinear angles for the magnetic moments. See KkrCalculation for details.'
)
