# Workflow for STM swf_parameterson around a magnetic impurity

from aiida.engine import WorkChain, ToContext, if_, calcfunction
from aiida.orm import Dict, RemoteData, Code, CalcJobNode, WorkChainNode, Float, Bool, XyData, SinglefileData, List
from aiida.orm import Group, load_group
from aiida.common import LinkType
from aiida_kkr.workflows.kkr_imp_dos import kkr_imp_dos_wc
from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc
from aiida_kkr.tools.find_parent import get_calc_from_remote
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode

__copyright__ = (u'Copyright (c), 2024, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.1'
__contributors__ = (u'Raffaele Aliberti', u'David Antognini Silva', u'Philipp Rüßmann')
_VERBOSE_ = True


class kkr_STM_wc(WorkChain):
    """
    Workchain for the Simulation of a (superconducting) STM tip around a magnetic impurity.
    The workchain uses the (converged) impurity calculation of a host system and combines them
    with vacuum sites in positions

    inputs::
        :param options: (Dict), computer options
        :param tip_position: (Dict), specify the position of the STM tip
        :param impurity_info: (Dict), information about the impurity cluster
        :param host_calc: (RemoteData), information about the host structure of the sample
        :param wf_parameters: (Dict), parameters that are used to run the workflow
        :param kkr: (Code), KKR host code for the writing out kkrflex files
        :param kkrimp: (Code), KKR impurity code for the normal state impurity scf and BdG impurity DOS calculation
        :param gf_writeout.params_kkr_overwrite (Dict), overwrite parameters for the GF calculation
        :param kkr_imp_sub.params_kkr_overwrite (Dict), overwrite parameters for the impurity calculation

     returns::

        :return workflow_info: (Dict), Information of workflow results
                            like success, last result node, list with convergence behavior
        :return STM_dos_data: (XYData), Returns the plot of the lmDOS of the calculation
        :retrun STM_lmdos_data: (XYData), Returns the interpolated lmDOS of the calculation"""

    # TO DO: Add the initialize step.
    # TO DO: Add a better creation of the impurity cluster.
    # TO DO: Add check that between the ilayer and the actual number of layers in the structure.
    # TO DO: Add to the outputs the calculated imp_info and imp_potential.
    # TO DO: Add BdG options for the builder run

    _wf_version = __version__
    _wf_label = 'STM_wc'
    _wf_description = 'Workflow for simulating an STM measurement'

    _options_default = {
        'queue_name': '',  # Queue name to submit jobs too
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 48,
            'num_cores_per_mpiproc': 1
        },  # resources to allocate for the job
        'max_wallclock_seconds': 3600 * 2,  # walltime after which the job gets killed (gets parsed to KKR)}
        'custom_scheduler_commands': '',  # some additional scheduler commands
        'withmpi': True
    }  # execute KKR with mpi or without

    _wf_default = {
        'jij_run': False,  # calculate Jij's energy resolved
        'lmdos':
        True,  # calculate also (l,m) or only l-resolved DOS, for the STM wf this is alwyas set on True as default
        'retrieve_kkrflex': False,  # retrieve kkrflex files to repository or leave on remote computer only
    }
    # add defaults of dos_params since they are passed onto that workflow
    _wf_default['dos_params'] = kkr_imp_dos_wc.get_wf_defaults()['dos_params']

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(kkr_STM_wc, cls).define(spec)

        spec.input('kkr', valid_type=Code, required=False, help='KKRhost code, needed if gf_dos_remote is not given.')
        spec.input('kkrimp', valid_type=Code, required=True, help='KKRimp code, always needed.')
        spec.input(
            'options',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help='Computer options (resources, quene name, etc.).'
        )
        spec.input(
            'wf_parameters',
            valid_type=Dict,
            required=False,
            default=lambda: Dict(dict=cls._wf_default),
            help='Workflow parameter (see `kkr_dos_wc.get_wf_defaults()`).'
        )
        spec.input(
            'tip_position',
            valid_type=Dict,
            required=False,
            # Find a way to create an area of this position.
            default=lambda: Dict({
                'ilayer': 0,
                'nx': 0,
                'ny': 0
            }),
            # In the previous line we have set to study the first layer
            # nx is the number of (symmetric) steps that we take in the x-direction starting from the impurity
            # ny is the number of (symmetric) steps that we take in the y-direction starting from the impurity
            # (0,0) correspond to calculate the DOS only on the impurity site
            help=
            'How many sites will be scanned in the da and db direction (Bravais Lattice). And the layer that is being scanned.'
        )
        spec.input(
            'imp_info',
            valid_type=Dict,
            required=True,
            help='Information of the impurity like position in the unit cell, screening cluster, atom type.'
        )
        spec.input(
            'host_calc',
            valid_type=RemoteData,
            required=False,
            help='The information about the clean host structure is required in order to continue the cluster'
            'Inside a bigger host structure with empty sites.'
        )
        spec.input(
            'host_remote',
            valid_type=RemoteData,
            required=True,
            help='Remote Data containing the remote folder from the outputs of the host calculation',
        )
        spec.input(
            'imp_potential_node',
            valid_type=SinglefileData,
            required=True,
            help='Impurity potential node',
        )
        spec.input(
            'remote_data',
            valid_type=RemoteData,
            required=False,
            help='Remote data from a converged kkr calculation, required for the gf writeout step',
        )
        spec.input(
            'kkrflex_files',
            valid_type=RemoteData,
            required=False,
            help='with this input we can directly load the gf_dos files without calculating them'
        )

        # Here we expose the inputs for the GF calculations step.
        # One parameter which is crucial is the NSHELD, which determines the impurity cluster radius.
        spec.expose_inputs(kkr_flex_wc, namespace='gf_writeout', include=('params_kkr_overwrite', 'options'))

        # Here we expose the BdG calculations from the kkr_imp_dos_wc
        spec.expose_inputs(kkr_imp_sub_wc, namespace='BdG', include=('params_overwrite'))
        spec.expose_inputs(kkr_imp_sub_wc, include=('initial_noco_angles', 'rimpshift'))

        # Specify the possible outputs
        spec.output('tip_position', valid_type=Dict)
        spec.output('STM_dos_data', valid_type=XyData, required=True)
        spec.output('STM_dos_data_lmdos', valid_type=XyData, required=True)
        #spec.output('workflow_info', valid_type=Dict)
        spec.output('kkrflexfiles', valid_type=RemoteData)
        spec.output('combined_imp_info', valid_type=Dict)
        spec.output('combined_imp_potential', valid_type=SinglefileData)

        # Define all possible error messages

        spec.exit_code(100, 'ERROR_STM_POSITION_NOT_VALID', 'The position provided for the STM probe are incorrect')
        spec.exit_code(101, 'ERROR_IMP_INFO_NOT_CORRECT', 'The node provided for the impurity info is not valid')
        spec.exit_code(102, 'ERROR_NO_IMP_POT_SFD', 'No impurity node has been given in the intput')
        spec.exit_code(103, 'ERROR_NO_IMPURITY_INFO', 'No impurity info has been given in the input')
        spec.exit_code(
            104, 'ERROR_NO_DATA_FOR_THE_GF_STEP', """Neither the kkrflex files nor the KKR builder have been given.
Please provide already converged kkrflex files, or the kkr builder to evaluate them"""
        )
        spec.exit_code(201, 'ERROR_IMP_SUB_WORKFLOW_FAILURE', 'A step in the kkr_imp_dos workflow has failed')

        spec.outline(
            # For initializing workflow
            cls.start,
            # We first aggregate all the impurity data
            # The gf is then used to evaluate the STM lmdos
            cls.STM_lmdos_run,
            cls.results
        )

    def combine_potentials(self, host_structure, impurity_to_combine, da, db):
        from aiida_kkr.tools.tools_STM_scan import get_imp_info_add_position
        import numpy as np  # TO DO: optimize this call, only need append from numpy
        """
        Here we want to combine the impurity information and the host information
        """
        tip_position = {}

        tip_position['ilayer'] = self.inputs.tip_position['ilayer']
        tip_position['da'] = da
        tip_position['db'] = db
        imp_info = self.inputs.imp_info  #(impurity to combine)

        combined_imp_info = get_imp_info_add_position(Dict(tip_position), host_structure, imp_info)
        # Since the objects in AiiDA are immutable we have to create a new dictionary and then convert
        # it to the right AiiDA type

        #new_combined_imp_info = {}
        # Add check to see if imp_cls is there
        if 'imp_cls' in impurity_to_combine:

            for key in impurity_to_combine.keys():
                if key == 'Zimp':
                    impurity_to_combine[key].append(combined_imp_info[key][-1])
                else:
                    impurity_to_combine[key] = np.append(impurity_to_combine[key], [combined_imp_info[key][-1]], axis=0)

            # Convert to an AiiDA Dictionary
            new_combined_imp_info = impurity_to_combine

        else:

            new_combined_imp_info = combined_imp_info

        return new_combined_imp_info

    def combine_nodes(self, host_calc, node_to_combine, da, db):
        from aiida_kkr.tools.tools_STM_scan import create_combined_potential_node
        """
        Here we create a combined potential node from the host potential (no impurity)
        and from the impurity potential
        """

        # Since the objects in AiiDA are immutable we have to create a new dictionary and then convert
        # it to the right AiiDA type
        tip_position = {}
        # for now we require that the z position remains the same.
        tip_position['ilayer'] = self.inputs.tip_position['ilayer']
        tip_position['da'] = da
        tip_position['db'] = db

        combined_node = create_combined_potential_node(tip_position, host_calc, node_to_combine)
        return combined_node

    def start(self):
        """
        Initialise context and some parameters
        """

        self.report(f'INFO: started STM workflow version {self._wf_version}')
        if _VERBOSE_:
            self.report(f'inputs: {self.inputs}')

        # Input both wf and options parameters
        # We check if the inputs for the wc where given
        # Otherwise we assign the default values
        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
            # empty dictionary evaluate as false to Python
            if options_dict == {}:
                options_dict = self._options_default
                self.report('INFO: using default options parameters')

        if 'wf_parameters' in self.inputs:
            wf_param_dict = self.inputs.wf_parameters.get_dict()
            # empty dictionary evaluate as false to Python
            if wf_param_dict == {}:
                options_dict = self._wf_default
                self.report('INFO: usign defalut wf parameters')

        # In this section we assign the computational resources to the builder
        self.ctx.withmpi = options_dict.get('withmpi', self._options_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get(
            'max_wallclock_seconds', self._options_default['max_wallclock_seconds']
        )
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get(
            'custom_scheduler_commands', self._options_default['custom_scheduler_commands']
        )
        self.ctx.options_params_dict = Dict({
            'withmpi': self.ctx.withmpi,
            'resources': self.ctx.resources,
            'max_wallclock_seconds': self.ctx.max_wallclock_seconds,
            'queue_name': self.ctx.queue,
            'custom_scheduler_commands': self.ctx.custom_scheduler_commands
        })

        # Set workflow parameters for the KKR imputrity calculations
        # This part is really important, this should always be set to True for an STM calculation
        self.ctx.lmdos = wf_param_dict.get('lmdos', self._wf_default['lmdos'])

        self.ctx.retrieve_kkrflex = wf_param_dict.get('retrieve_kkrflex', self._wf_default['retrieve_kkrflex'])

        self.ctx.dos_params_dict = wf_param_dict.get('dos_params', self._wf_default['dos_params'])

        # fill missing key, value pairs with defaults
        for k, v in self._wf_default['dos_params'].items():
            if k not in self.ctx.dos_params_dict.keys():
                self.ctx.dos_params_dict[k] = v

        # set workflow parameters for the KKR impurity calculation
        self.ctx.jij_run = wf_param_dict.get('jij_run', self._wf_default['jij_run'])

        # set workflow label and description
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)

        if _VERBOSE_:
            message = f"""
INFO: use the following parameter:
withmpi: {self.ctx.withmpi}
Resources: {self.ctx.resources}
Walltime (s): {self.ctx.max_wallclock_seconds}
queue name: {self.ctx.queue}
scheduler command: {self.ctx.custom_scheduler_commands}
description: {self.ctx.description_wf}
label: {self.ctx.label_wf}
                      """
            self.report(message)

    def validate_input(self):
        """Check if inputs are valid"""
        inputs = self.inputs

        if not 'imp_potential_node' in inputs:
            return self.exit_codes.ERROR_NO_IMP_POT_SFD  # pylint: disable=no-member

        if not 'imp_info' in inputs:
            return self.exit_codes.ERROR_NO_IMP_INFO  # pylint: disable=no-member

        if not 'kkrflex_files' and 'kkr' in inputs:
            return self.exit_codes.ERROR_NO_DATA_FOR_THE_GF_STEP  # pylint: disable=no-member

    def impurity_cluster_evaluation(self):
        """
        Create the combined impurity cluster and potential for the impurity region
        used in self-consistency + the additional scanning sites.
        """
        from aiida_kkr.tools import find_parent_structure
        if _VERBOSE_:
            from time import time

            # measure time at start
            t_start = time()

        # Here we create an impurity cluster that has inside all the positions on which the STM will scan

        impurity_info = self.inputs.imp_info  # for the first step we combine the impurity info from the input
        imp_potential_node = self.inputs.imp_potential_node  # for the first step we combine the impurity node from the input

        host_remote = self.inputs.host_remote
        host_calc = host_remote.get_incoming(node_class=CalcJobNode).first().node
        host_structure = find_parent_structure(host_remote)

        # now find all the positions we need to scan
        coeff = self.get_scanning_positions(host_remote)

        if _VERBOSE_:
            # timing counters
            t_imp_info, t_pot = 0., 0.

        # construct impurity potential and imp_info for the impurity cluster + scanning area
        for element in coeff:
            if _VERBOSE_:
                t0 = time()
            tmp_imp_info = self.combine_potentials(host_structure, impurity_info, element[0], element[1])
            impurity_info = tmp_imp_info

            if _VERBOSE_:
                t_imp_info += time() - t0
                t0 = time()

            # Aggregation the impurity nodes
            tmp_imp_pot = self.combine_nodes(host_calc, imp_potential_node, element[0], element[1])
            imp_potential_node = tmp_imp_pot

            if _VERBOSE_:
                t_pot += time() - t0

        if _VERBOSE_:
            # report elapsed time for cluster generation
            self.report(f'time for cluster generation (s): {time()-t_start}, imp_info={t_imp_info}, pot={t_pot}')

        return impurity_info, imp_potential_node

    def get_scanning_positions(self, host_remote):
        """
        Extract scanning positions either from input 'scan_positions' or from 'nx', 'ny' + symmetry analysis

        If  'scan_positions' is found in 'tip_position' input dict we use these positions which should
        be 2D array of integers with the positions in units of the structure's in-plane Bravais matrix.

        Otherwise we use the 'nx', 'ny' input to define a scanning region where an automated symmetry
        analysis is done to reduce the scanning area to the irreducible part.
        """
        from aiida_kkr.tools import tools_STM_scan

        generate_scan_positions = True
        if 'scan_positions' in self.inputs.tip_position:
            coeff = self.inputs.tip_position['scan_positions']
            if coeff is not None:
                # check if coefficients exists and are valid
                # TODO: improve the validity check
                generate_scan_positions = False

        if generate_scan_positions:

            # Information of the host structure
            struc_info, symm_matrices = tools_STM_scan.STM_pathfinder(host_remote)

            # We now want to iterate over several in-plane positions.
            # These are the number of vectors in which we want to move the STM tip.
            x = self.inputs.tip_position['nx']
            y = self.inputs.tip_position['ny']

            # Path creation step. (The the identity operator is present, but will be excluded)
            unused_pos, used_pos = tools_STM_scan.lattice_generation(x, y, symm_matrices, struc_info['plane_vectors'])

            # Since the combine tools use the element already in the units of da and db, we use a helper function
            # to have the indices of the linear combination of the used position vectors in the base of the Bravais lattice.
            coeff = tools_STM_scan.find_linear_combination_coefficients(struc_info['plane_vectors'], used_pos)

        return coeff

    def STM_lmdos_run(self):
        """In this part of the worflow we want to simulate the lmdos which a STM is able to measure """

        # First we would like to distinguish between an impurity dos and a normal state calculation
        builder = kkr_imp_dos_wc.get_builder()

        # Code loading
        builder.kkrimp = self.inputs.kkrimp  # needed for the kkr_imp_dos_wc

        # Builder options
        builder.options = self.ctx.options_params_dict

        # Check if the kkrflex files are already given in the outputs
        if 'kkrflex_files' in self.inputs:
            builder.gf_dos_remote = self.inputs.kkrflex_files
            message = f'Remote host function is given in the outputs from the node: {self.inputs.kkrflex_files}'
            self.report(message)
        else:
            builder.kkr = self.inputs.kkr  # needed to evaluate the kkr_flex files in the DOS step

        # NSHELD is the parameter that controls the radius in the impurity cluster.
        # The bigger the scanning position, the greater it must be set.
        if 'gf_writeout' in self.inputs:
            if 'params_kkr_overwrite' in self.inputs.gf_writeout:
                builder.gf_writeout.params_kkr_overwrite = self.inputs.gf_writeout.params_kkr_overwrite  # pylint: disable=no-member
            if 'options' in self.inputs.gf_writeout:
                builder.gf_writeout.options = self.inputs.gf_writeout.options  # pylint: disable=no-member
        else:
            # This is a big value of NSHELD to make sure that most calculations work
            builder.gf_writeout.params_kkr_overwrite = Dict(dict={'NSHELD': 1500})  # pylint: disable=no-member

        # Update the BdG parameters if they are inserted in the workflow
        if 'BdG' in self.inputs:
            if 'params_kkr_overwrite' in self.inputs.BdG:
                builder.BdG.params_overwrite = self.inputs.BdG.params_kkr_overwrite  # pylint: disable=no-member

        self.ctx.kkrimp_params_dict = Dict(
            dict={
                'nsteps': 1,  # redundant because this is already set inside the kkr_imp_dos workchain?!
                'kkr_runmax': 1,  # redundant because this is already set inside the kkr_imp_dos workchain?!
                'dos_run': True,  # redundant because this is already set inside the kkr_imp_dos workchain?!
                'retrieve_kkrflex': self.ctx.retrieve_kkrflex,
                'lmdos': self.ctx.lmdos,
                'jij_run': self.ctx.jij_run,
                'dos_params': self.ctx.dos_params_dict
            }
        )

        # We want to set the energy to the Fermi level
        if 'emin' not in self.ctx.dos_params_dict:
            self.ctx.kkrimp_params_dict['dos_params']['emin'] = 0 - 0.005
        if 'emax' not in self.ctx.dos_params_dict:
            self.ctx.kkrimp_params_dict['dos_params']['emax'] = 0 + 0.005

        # Finally we overwrite the number of energy points to 1
        # This is because we want many epoints around the impurity position
        if 'nepts' not in self.ctx.dos_params_dict:
            self.ctx.kkrimp_params_dict['dos_params'][
                'nepts'] = 7  # Here 7 because of the interpolated files that aren't generated

        builder.wf_parameters = self.ctx.kkrimp_params_dict
        # Host remote files that will be used for the actual plot step.
        builder.host_remote = self.inputs.host_remote

        # Here we create the impurity cluster for the STM scanning tool
        impurity_info, imp_pot_sfd = self.impurity_cluster_evaluation()

        # impurity info for the workflow
        builder.impurity_info = impurity_info
        builder.imp_pot_sfd = imp_pot_sfd

        # submit calculation
        calc = self.submit(builder)

        # print report
        message = f"""INFO: running DOS step for an STM measurement (pk: {calc.pk}) at position (ilayer: {self.inputs.tip_position['ilayer']})"""
        if 'params_kkr_overwrite' in self.inputs.BdG:
            if self.inputs.BdG.params_kkr_overwrite:
                message += f'\nINFO: runnig DOS step (pk: {calc.pk}) BdG is present'
        self.report(message)

        # Save the calculated impurity cluster and impurity info in the context
        self.ctx.impurity_info = impurity_info
        self.ctx.imp_pot_sfd = imp_pot_sfd

        return ToContext(STM_data=calc)

    def results(self):
        """Collect results and return output nodes"""

        if not self.ctx.STM_data.is_finished_ok:
            self.report('ERROR: sub workflow for STM calculation failed')
            return self.exit_codes.ERROR_IMP_SUB_WORKFLOW_FAILURE  # pylint: disable=no-member
        else:
            # Declaring the output
            self.out('STM_dos_data', self.ctx.STM_data.outputs.dos_data)
            self.out('STM_dos_data_lmdos', self.ctx.STM_data.outputs.dos_data_lm)
            self.out('tip_position', self.inputs.tip_position)
            if 'gf_dos_remote' in self.ctx.STM_data.outputs:
                self.out('kkrflexfiles', self.ctx.STM_data.outputs.gf_dos_remote)
            self.out('combined_imp_info', self.ctx.impurity_info)
            self.out('combined_imp_potential', self.ctx.imp_pot_sfd)

        self.report('INFO: created output nodes for KKR STM workflow.')

        self.report(
            '\n'
            '|------------------------------------------------------------------------------------------------------------------|\n'
            '|-----------------------------------------| Done with the STM workflow! |------------------------------------------|\n'
            '|------------------------------------------------------------------------------------------------------------------|'
        )
