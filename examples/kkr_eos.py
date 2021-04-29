"""
Example calculation of fcc Cu making use of the SCF workchain for the kkrhost
code.
"""
import ase.spacegroup
from aiida import orm
from aiida.engine import submit
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.workflows.eos import kkr_eos_wc
from aiida_kkr.workflows.voro_start import kkr_startpot_wc
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc


def create_fcc_cu(a_lat=lambda: orm.Float(3.6142)) -> orm.StructureData:
    """
    Function to generate a fcc Cu structure

    The function also makes sure that the formula and the space group number
    are stored in the extras so that they can be more easily queried via the
    QueryBuilder.

    :param a_lat: lattice parameter for fcc Cu, defaults to lambda:orm.Float(3.6142)
    :type a_lat: [type], optional
    """
    # Get the lattice distance for the fcc lattice
    alat = 0.5 * a_lat.value
    structure = orm.StructureData(cell=[[alat, alat, 0.0], [alat, 0.0, alat], [0.0, alat, alat]])
    # Add the positions of the atoms in the unitcell
    structure.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')

    # Get the spacegroup information of the unitcell
    sg = ase.spacegroup.get_spacegroup(structure.get_ase())

    # Storing the formula as an extra so that it can easily be found
    structure.set_extra('formula', structure.get_formula(mode='count'))
    # Storing the spacegroup number so that it can be easily found
    structure.set_extra('space_group_number', sg.no)

    return structure


def generate_kkr_parameters(input_parameters):
    """
    Function to generate a dictionary with the parameters needed for a
    kkrhost run. It overwrites the default parameters with the provided
    parameters
    """
    kkr_calculation_parameters = kkrparams(**input_parameters)

    return kkr_calculation_parameters


def main():
    """
    Main function to run an SCF run using the kkrhost code
    """
    # Define the KKR code via the code string
    kkr_code = orm.load_code('jukkr-host-intel@moggie')
    # Define the Voronoi code via the code string
    voronoi_code = orm.load_code('jukkr-voronoi-intel@moggie')

    # Generate the fcc Cu structure that will be used for the calculation
    fcc_structure = create_fcc_cu(a_lat=orm.Float(3.6142))

    # Generate a dictionary with the necessary parameters for the kkrhost
    # calculation
    kkr_calculation_parameters = generate_kkr_parameters(
        input_parameters={
            'LMAX': 2,
            'RMAX': 7,
            'GMAX': 65,
            'NSPIN': 1,
            'RCLUSTZ': 1.9
        }
    )

    # For the SGE scheduler
    if kkr_code.computer.scheduler_type in ['sge']:
        # Setting parameters for the submission of the job
        options = {
            'resources': {
                # Total number of mpi processes
                'tot_num_mpiprocs': 16,
                # Name of the parallel environment
                'parallel_env': 'mpi',
            },
            # Maximum allowed execution time in seconds
            'max_wallclock_seconds': 18000,
            # Whether to run in parallel
            'withmpi': True,
            'custom_scheduler_commands': '#$ -l slot_type=hpc'
        }

    # For slurm and pbs scheduler types
    if kkr_code.computer.scheduler_type in ['slurm', 'pbspro']:
        # Setting parameters for the submission of the job
        options = {
            'resources': {
                # Total number of mpi processes
                'num_machines': 1,
            },
            # Maximum allowed execution time in seconds
            'max_wallclock_seconds': 18000,
            # Whether to run in parallel
            'withmpi': True,
        }

    wf_parameters = {
        # Minimum and maximum volume that will be studied
        'scale_range': [0.90, 1.10],
        # Number of volume steps in the eos determination
        'nsteps': 11,
        # create and return a structure which has the ground state volume
        # determined by the fit used
        'ground_state_structure': True,
        # use seekpath to get primitive structure after scaling to reduce
        # computational time
        'use_primitive_structure': True,
        # fitfunction used to determine ground state volume
        # (see ase.eos.EquationOfState class for details)
        'fitfunction': 'birchmurnaghan',
        # settings for kkr_startpot behavior
        'settings_kkr_startpot': kkr_startpot_wc.get_wf_defaults(silent=True),
        # settings for kkr_scf behavior
        'settings_kkr_scf': kkr_scf_wc.get_wf_defaults(silent=True)[0]
    }

    # Setting a label for the calculation
    label = 'JM-KKR_host_fcc_Cu'
    # Provide a more lengthy description of the actual calculation
    description = 'Test FCC Cu SCF calculation with the JM-KKR host code'

    # Populate the inputs of the kkrhost code
    builder = kkr_eos_wc.get_builder()
    builder.kkr = kkr_code
    builder.voronoi = voronoi_code
    builder.options = orm.Dict(dict=options)
    builder.structure = fcc_structure
    builder.wf_parameters = orm.Dict(dict=wf_parameters)
    builder.calc_parameters = orm.Dict(dict=kkr_calculation_parameters)
    builder.metadata.label = label
    builder.metadata.description = description

    # Submit the calculation to the daemon
    submission = submit(builder)
    print('EOS workchain calculation submitted with pk {}'.format(submission.pk))


if __name__ == '__main__':
    main()
