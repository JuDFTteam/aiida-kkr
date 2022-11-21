# -*- coding: utf-8 -*-
'''
Module with CLI commands to launch for calcjob and workflows of aiida-kkr.
'''
import click

from aiida.cmdline.params import types
from aiida.cmdline.params.options import OverridableOption
from aiida.orm import Code, load_node, Dict
from aiida.plugins import WorkflowFactory
from aiida.plugins import CalculationFactory
from masci_tools.io.kkr_params import kkrparams

from aiida_kkr.tools.dict_util import clean_nones
from ..util import options
from ..util.utils import launch_process
from ..util import defaults

# TODO: command for kkrimporter
# not sure if this should be a launch command
# TODO: commands for workchains: voro_start, eos, gf_writeout, kkr_imp, kkr_imp_dos, kkr_imp_sub
# Check_para_convergence, check_magnetic_state. base_restart_calc?


@click.command('voro')
@options.STRUCTURE_OR_FILE(default=defaults.get_cu_bulk_structure, show_default=True)
@options.VORO()
@options.PARAMETERS()
@options.PARENT_FOLDER()
@options.POTENTIAL_OVERWRITE()
@options.MAX_WALLCLOCK_SECONDS()
@options.QUEUE_NAME()
@options.DAEMON()
def launch_voro(
    structure,
    voro,
    parameters,
    parent_folder,
    potential_overwrite,
    max_wallclock_seconds,
    queue_name,
    daemon,
):
    """
    Launch an Voronoi calcjob on given input
    """
    # TODO?: maybe allow for additional metadata to be given.
    process_class = CalculationFactory('kkr.voro')

    inputs = {
        'structure': structure,
        'code': voro,
        'parameters': parameters,
        'parent_kkr': parent_folder,
        'potential_overwrite': potential_overwrite,
        'metadata': {
            'options': {
                'withmpi': False,
                'queue_name': queue_name,
                'max_wallclock_seconds': max_wallclock_seconds,
                'resources': {
                    'num_machines': 1,
                    'num_mpiprocs_per_machine': 1,
                }
            }
        }
    }
    inputs = clean_nones(inputs)
    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)


@click.command('kkr')
@options.KKR()
@options.PARAMETERS()
@options.PARENT_FOLDER(required=True)
@options.IMPURITY_INFO()
@options.KPOINTS()
@options.DAEMON()
@options.WITH_MPI()
@options.NUM_MPIPROCS_PER_MACHINE()
@options.MAX_WALLCLOCK_SECONDS()
@options.MAX_NUM_MACHINES()
@options.QUEUE_NAME()
def launch_kkr(
    kkr,
    parameters,
    parent_folder,
    impurity_info,
    kpoints,
    daemon,
    with_mpi,
    num_mpiprocs_per_machine,
    max_wallclock_seconds,
    max_num_machines,
    queue_name,
):
    """
    Launch an KKRhost calcjob on given input
    """
    # TODO?: maybe allow for additional metadata to be given.
    process_class = CalculationFactory('kkr.kkr')

    inputs = {
        'code': kkr,
        'parameters': parameters,
        'parent_folder': parent_folder,
        'impurity_info': impurity_info,
        'metadata': {
            'options': {
                'withmpi': with_mpi,
                'max_wallclock_seconds': max_wallclock_seconds,
                'queue_name': queue_name,
                'resources': {
                    'num_machines': max_num_machines,
                    'num_mpiprocs_per_machine': num_mpiprocs_per_machine,
                },
            }
        }
    }
    inputs = clean_nones(inputs)
    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)


@click.command('kkrimp')
@options.KKRIMP()
@options.PARAMETERS()
@options.PARENT_FOLDER(required=True)
@options.IMPURITY_INFO()
@options.DAEMON()
@options.WITH_MPI()
@options.NUM_MPIPROCS_PER_MACHINE()
@options.MAX_WALLCLOCK_SECONDS()
@options.MAX_NUM_MACHINES()
def launch_kkrimp(
    kkrimp,
    parameters,
    parent_folder,
    impurity_info,
    daemon,
    with_mpi,
    num_mpiprocs_per_machine,
    max_wallclock_seconds,
    max_num_machines,
):
    """
    Launch an KKRimp calcjob on given input
    """
    # TODO?: maybe allow for additional metadata to be given.
    process_class = CalculationFactory('kkr.kkrimp')

    inputs = {
        'code': kkrimp,
        'parameters': parameters,
        'parent_folder': parent_folder,
        'impurity_info': impurity_info,
        'metadata': {
            'options': {
                'withmpi': with_mpi,
                'max_wallclock_seconds': max_wallclock_seconds,
                'resources': {
                    'num_machines': max_num_machines,
                    'num_mpiprocs_per_machine': num_mpiprocs_per_machine,
                }
            }
        }
    }
    inputs = clean_nones(inputs)
    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)


@click.command('dos')
@options.KKR()
@options.WF_PARAMETERS()
@options.OPTION_NODE()
@options.PARENT_FOLDER()
@options.DAEMON()
def launch_dos(kkr, wf_parameters, option_node, parent_folder, daemon):
    """
    Launch an KKRhost density of states workflow
    """
    # TODO?: maybe allow for additional metadata to be given.
    process_class = WorkflowFactory('kkr.dos')

    inputs = {
        'kkr': kkr,
        'remote_data': parent_folder,
        'options': option_node,
        'wf_parameters': wf_parameters,
    }
    inputs = clean_nones(inputs)
    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)


@click.command('bs')
@options.KKR()
@options.WF_PARAMETERS()
@options.KPOINTS()
@options.OPTION_NODE()
@options.PARENT_FOLDER()
@options.DAEMON()
def launch_bs(kkr, wf_parameters, kpoints, option_node, parent_folder, daemon):
    """
    Launch an KKRhost bandstructure workflow with required inputs
    (kkr code, remote_data, options, wf-parameters, daemon)
    """
    process_class = WorkflowFactory('kkr.bs')

    inputs = {
        'kkr': kkr,
        'remote_data': parent_folder,
        'kpoints': kpoints,
        'options': option_node,
        'wf_parameters': wf_parameters,
    }
    inputs = clean_nones(inputs)
    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)


@click.command('scf')
@options.KKR()
@options.VORO()
@options.STRUCTURE_OR_FILE(default=defaults.get_cu_bulk_structure, show_default=True)
@options.PARAMETERS()
@options.PARENT_FOLDER()
@options.DAEMON()
@options.OPTION_NODE()
@options.WF_PARAMETERS()
@options.POTENTIAL_OVERWRITE()
@options.NOCO_ANGLES()
def launch_scf(
    kkr,
    voro,
    structure,
    parameters,
    parent_folder,
    daemon,
    option_node,
    wf_parameters,
    potential_overwrite,
    noco_angles,
):
    """
    Launch an KKRhost self-consistency workflow
    """
    # TODO?: maybe allow for additional metadata to be given.
    process_class = WorkflowFactory('kkr.scf')

    inputs = {
        'kkr': kkr,
        'voronoi': voro,
        'structure': structure,
        'calc_parameters': parameters,
        'remote_data': parent_folder,
        'options': option_node,
        'wf_parameters': wf_parameters,
        'startpot_overwrite': potential_overwrite,
        'initial_noco_angles': noco_angles,
    }
    inputs = clean_nones(inputs)
    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)


PARENT_FOR_IMP = OverridableOption(
    '-P',
    '--parent-folder',
    'parent_folder',
    type=types.DataParamType(sub_classes=('aiida.data:core.remote',)),
    show_default=True,
    required=True,
    help=
    'The remote folder of a parent calculation (either the converged calculation or a GF writeout calculation which prevents recalculating the host Greens function).'
)

PARAMS_HOST_GF = OverridableOption(
    '-p',
    '--parameters-hostGF',
    'parameters_hostgf',
    type=types.DataParamType(sub_classes=('aiida.data:core.dict',)),
    show_default=True,
    help='Set of parameters that are overwritten in the host GF writeout step.'
)

IMP_STARTPOT = OverridableOption(
    '--impurity-startpot',
    'impurity_startpot',
    type=types.DataParamType(sub_classes=('aiida.data:core.singlefile',)),
    help=
    'Use this as the starting potential for the impurity calculation. Needs to match the settings in the impurity info node.'
)


def _check_parent_calc_type(parent_folder):
    """
    check if parent folder was a host GF writeout calculation
    """
    from aiida_kkr.calculations import KkrCalculation
    calc = parent_folder.get_incoming(node_class=KkrCalculation).first()
    calc = calc.node
    ret = calc.outputs.retrieved
    return 'kkrflex_green' in ret.list_object_names()


@click.command('kkrimpscf')
@options.KKR()
@options.VORO()
@options.KKRIMP()
@options.IMPURITY_INFO(required=True, show_default=True)
@PARENT_FOR_IMP()
@options.OPTION_NODE()
@options.DAEMON()
@PARAMS_HOST_GF()
@options.WF_PARAMETERS()
@IMP_STARTPOT()
def launch_kkrimp_scf(
    kkr,
    voro,
    kkr_imp,
    impurity_info,
    parent_folder,
    option_node,
    daemon,
    parameters_hostgf,
    wf_parameters,
    impurity_startpot,
):
    """
    Launch an kkr calcjob on given input
    """
    # TODO?: maybe allow for additional metadata to be given.
    process_class = WorkflowFactory('kkr.imp')

    inputs = {
        'kkr': kkr,
        'kkrimp': kkr_imp,
        'voronoi': voro,
        'options': option_node,
        'wf_parameters': wf_parameters,
        'impurity_info': impurity_info,
        'params_kkr_overwrite': parameters_hostgf,
        'startpot': impurity_startpot,
    }
    inputs = clean_nones(inputs)

    if _check_parent_calc_type(parent_folder):
        # this means we can reuse the GF writeout from the parent calculation
        inputs['remote_data_gf'] = parent_folder
        click.echo('The parent_folder already has the written out GF of the host, we can therefore skip that step')
    else:
        # this means the calculation was not a GF writeout KKR calculation
        inputs['remote_data_host'] = parent_folder
        click.echo('The parent_folder does not contain the host GF, we calculate it.')

    builder = process_class.get_builder()
    builder.update(inputs)
    launch_process(builder, daemon)
