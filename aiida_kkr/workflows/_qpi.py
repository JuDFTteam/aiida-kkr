# Workflow for real-space quasiparticle interference calculation around an impurity

from aiida import orm, engine
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode

__copyright__ = (u'Copyright (c), 2023, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = (u'Philipp Rüßmann')

class qpi_wc(engine.WorkChain):
    """
    AiiDA-KKR workchain for calculation of real-space quasiparticle interference (QPI) images
    (like Fig. 2a Ref. [1]).
    
    The idea is to combine an impurity calculation with its converged potential in the impurity cluster
    with additional "scanning positions" (in analogy to an STM measurement). Typically this is an impurity 
    at a surface and the scanning are is a region above the surface where an STM tip would scan.

    The real-space QPI image is calculated as the change in the local density of states on the 
    scanning positions due to the presence of the nearby impurity. According to the Tersof-Hamann model
    this change is expressed as the local density of states at the scanning sites. We calculate this in
    the usual way of computing a DOS for a single energy point with an impurity cluster that consists of
    the converged impurity region plus additional scanning positions where the host potential is taken.
    The change in the density at the scanning positions is then not a change in the potential but the
    standing wave that forms due to scattering of electrons at the defect
    (i.e. these are Friedel oscillations around the defect).


    inputs:

    *   workflow settings:
        *   energy-point (eV units relative to EF)
        *   accuracy parameters (via kkr-params overwrite?): k-point density, should somehow be inverse-proportional to scanning range
        *   (additional) scanning positions, batch size?
    *   kkrhost code: `host.code`, `host.options`
    *   kkrimp code: `imp.code`, `imp.options`
    *   impurity scf: `imp.potential`, hostparent taken from impurity scf or overwirte from `host.parent_remote`
    *   general options (overwritten by imp, host namespaces)
    *   expose optional inputs for `gf_writeout` step `host.kkrparams_overwrite`

    outputs:

    *   QPI (DOS at scanning positions)
    *   expose kkrimp calculations, host gf

    outline:

    *   check inputs
    *   create `impurity_info` for unique positions, do this for all batches (reuse stuff from combine_imps workchain)
    *   create impurity potential combining impurity cluster + host potential from additional scanning positions
    *   run GF writeout step
    *   run KKRimp step
    *   collect data from batches

    References:

    [1] P. Rüßmann, P. Mavropoulos, and S. Blügel, *Ab Initio Theory of Fourier-Transformed Quasiparticle Interference Maps and Application to the Topological Insulator Bi2Te3*. Phys. Status Solidi B **258**, 2000031 (2021). [https://doi.org/10.1002/pssb.202000031](https://doi.org/10.1002/pssb.202000031)
    """

    _wf_version = __version__
    #TODO change Dict input to own class for this input that inherits from Dict?
    _wf_default = {
        'energy_point': 0.0, # energy point where QPI is calculated, relative to Fermi energy in eV units
        'temperature': 100., # temperature smearing for the QPI calculation
        'batch_size': -1, # size of the batches in which the QPI positions are grouped. -1 means only a single batch
        'RCLUSTZ': None, # size of screening cluster, if `None` take value from host calculation.
        'k-mesh': None, # size of k-point mesh, if `None` take value from host calculation.
    }

    @classmethod
    def get_wf_defaults(self, silent=False):
        """
        Return the default values of the workflow parameters (wf_parameters input node)
        """
        if not silent:
            print(f'Version of the kkrimp_BdG_wc workflow: {self._wf_version}')
        return self._wf_default.copy()

    @classmethod
    def define(cls, spec):
        """
        Layout of the workflow, defines the input nodes and the outline of the workchain
        """
        super(qpi_wc, cls).define(spec)

        # define input ports
        #TODO define input ports
        spec.input(
            'node',
            valid_type=,
            required=False,
            default=lambda: Dict(dict=cls._options_default),
            help=''
        )

        # define outputs
        #TODO define output ports
        spec.output('workflow_info', valid_type=Dict, required=False)

        # Outline of the workchain
        #TODO define outline
        spec.outline(
            # initialization and input validation
            cls.start,
            cls.create_imp_info,
            cls.create_imp_pot,
            cls.write_host_gf,
            cls.run_kkrimp,
            cls.results
        )

        # Define all possible error messages
        #TODO define exit codes
        spec.exit_code(
            100, 'ERROR_KKRCODE_NOT_CORRECT', 'The code you provided for kkr does not use the plugin kkr.kkr'
        )
        spec.exit_code(
            101, 'ERROR_KKRIMPCODE_NOT_CORRECT', 'The code you provided for kkrimp does not use the plugin kkr.kkrimp'
        )

    ###########################################################
    # Now here we define all the functions from the outline

    def start(self):
        """
        Set up context of the workflow and test input nodes for consistency
        """
        self.report(f'INFO: started KKR BdG impurity version {self._wf_version}')
    
        #TODO check inputs with custom validators instead of in here

        # validate for kkr code
        if 'kkr' in self.inputs:
            try:
                test_and_get_codenode(self.inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT  # pylint: disable=no-member

        # validate for kkrimp code, always done because this is a required input
        try:
            test_and_get_codenode(self.inputs.kkrimp, 'kkr.kkrimp', use_exceptions=True)
        except ValueError:
            return self.exit_codes.ERROR_KKRIMPCODE_NOT_CORRECT  # pylint: disable=no-member

    def create_imp_info(self):
        """
        create `impurity_info` for unique positions, do this for all batches (reuse stuff from combine_imps workchain)
        """
        #TODO implement
        raise NotImplementedError("'create_imp_info' step not impemented yet in qpi workchain")

    def create_imp_pot(self):
        """
        create impurity potential combining impurity cluster + host potential from additional scanning positions
        """
        #TODO implement
        raise NotImplementedError("'write_host_gf' step not impemented yet in qpi workchain")

    def write_host_gf(self):
        """
        run GF writeout step using KKRhost for all batches
        """
        #TODO implement
        raise NotImplementedError("'run_kkrimp' step not impemented yet in qpi workchain")

    def run_kkrimp(self):
        """
        run KKRimp step to calculate impurity DOS for all batches
        """
        #TODO implement
        raise NotImplementedError("'run_kkrimp' step not impemented yet in qpi workchain")

    def results(self):
        """
        Collect results of the calculations and return results node
        """
        #TODO implement
        raise NotImplementedError("'results' step not impemented yet in qpi workchain")

    ###########################################################
    # Now here we define all the helper functions used in the outline steps

    #TODO add if necessary