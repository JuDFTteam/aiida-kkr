#!/usr/bin/env python
# coding: utf-8

# In[18]:

# from __future__ import print_function, division, absolute_import
from six.moves import range
from aiida.orm import Code,  Dict, RemoteData, StructureData, Float
from aiida.orm import XyData, WorkChainNode, load_node, CalcJobNode, ArrayData
from aiida.engine import WorkChain, if_, ToContext, submit
from aiida.engine import CalcJob, calcfunction
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import test_and_get_codenode, get_parent_paranode
from aiida_kkr.tools.common_workfunctions import update_params_wf, get_inputs_kkr
from aiida_kkr.calculations.kkr import KkrCalculation
from aiida_kkr.calculations.voro import VoronoiCalculation
from aiida.common.exceptions import InputValidationError, ConfigurationError

from aiida_kkr.tools.save_output_nodes import create_out_dict_node
from aiida.tools.data.array.kpoints import get_explicit_kpoints_path
__copyright__ = (u"This is the Copyright")
__license__ = "Here is the license"
__version__ = ("Here is the version")
__contributors__ = "using technique from Philpp RÜssmann's DOS wf"


# In[19]:



class kkr_BS_wf(WorkChain):
    """
    Workchain for BandStructure calculation, starting from the previous converged KKR calculation remote folder data

    ####### 001 with kkr starting from the remote data node
    ####### of a previous calculation (either Voronoi or KKR)

    inputs:
    :param wf_parameters: (Dict); Workchain Specifications, contains npt2, tempr, emin(ev), emax(ev), kmesh
                          keys. The Energy emin and emax are the energy difference from the fermi energy.

    :param options: (Dict); Computer Specifications
    :param kpoints: (KpointsData); Kpoints data type from thestructure,
                                   but not mendatory as it can be extrruct from structure internaly from the remote data
    :param remote_data: (RemoteData); From the previous kkr-converged calcualtion.
    :param kkr: (Code); KKR code specifiaction

    returns:
    :out result_wf: python dict for the work_chain_specifications
    """

    _wf_version = __version__
    _wf_label = 'kkr_BandStructure_wc'
    _wf_description = """Workflow for a bandstructure calculation starting eithe from a structure with automatic voronoi'
                        calculation or a valid RemoteData of a previous calculation."""
    _wf_params_default = {
                        'EMIN':10,  # Energy below the fermi surface for energy contour in ev unit
                        'EMAX': 5,     # Energy over the fermi surface for energy contour in ev unit
                        "NPT2": 200,   # Energy points in the energy contour
                        "kmesh": [30, 30, 30],
                        'TEMPR': 50.   # temperature
                          }
    _wf_option_default = {'max_wallclock_seconds': 36000,
                          'resources': {'tot_num_mpiprocs': 48, 'num_machines': 1},
                          'custom_scheduler_commands':
                          '#SBATCH --account=jara0191\n\nulimit -s unlimited; export OMP_STACKSIZE=2g',
                          'withmpi': True,
                          'queue_name': ''
                         }
    @classmethod
    def get_wf_params_defaults(self, silent=False):
        if not silent: print('Version of the wf pparameters  {}'.format(self._wf_version))
        return self._wf_params_default

    @classmethod
    def define(cls, spec):
        super(kkr_BS_wf, cls).define(spec)
        # here inputs are defined
        spec.input("wf_parameters", valid_type = Dict, required = False,default=lambda: Dict(self._wf_params_default))
        spec.input("kpoints", valid_type=KpointsData, required = False)
        spec.input("options", valid_type=Dict, required=False,default=lambda: Dict(dict=cls._wf_option_default))
        spec.input("remote_data", valid_type=RemoteData, required=True)
        spec.input("kkr", valid_type=Code, required=True)
        spec.input("label",  required=False)
        spec.input("description", required=False)
        spec.input("fromula", required=False)
        
        # Here outputs are defined
        spec.output("results_wf", valid_type=Dict, required=True)        
        spec.output('BS_data', valid_type=ArrayData, required=True)

        # Here outlines are being specified
        spec.outline(
            # For initialiging workflow
            cls.start,
            if_(cls.validate_input)(
                cls.set_params_BS,
                cls.get_BS),
            cls.return_results
        )
        # definition of exit code in case something goes wrong in this workflow
        spec.exit_code(161, 'ERROR_NO_INPUT_REMOTE_DATA', 'No remote_data was provided as Input')
        spec.exit_code(162, 'ERROR_KKRCODE_NOT_CORRECT', 'The code you provided for kkr does not use the plugin kkr.kkr')
        spec.exit_code(163, 'ERROR_CALC_PARAMETERS_INVALID', 'calc_parameters given are not consistent! Hint: did you give an unknown keyword?')
        spec.exit_code(164, 'ERROR_CALC_PARAMETERS_INCOMPLETE', 'calc_parameters not complete')
        spec.exit_code(165, 'ERROR_BS_CALC_FAILED', 'KKR Band Structure calculation failed')
        spec.exit_code(166, 'ERROR_NO_KPOINTS_EXTRACTED', 'No K-POINTS was extracted from the structure data')

    def start(self):
        self.report('INFO: started KKR Band Structure workflow version {}'
                    ''.format(self._wf_version))
        wf_dict = self.inputs.wf_parameters.get_dict()
        options_dict = self.inputs.options.get_dict()
        if wf_dict == {}:
            self.report('INFO: Using default wf parameter')
            wf_dict = self._wf_params_default
        if options_dict == {}:
            self.report('INFO: Using default wf Options')
            options_dict = self._wf_option_default
        self.ctx.withmpi = options_dict.get('withmpi', self._wf_option_default['withmpi'])
        self.ctx.resources = options_dict.get('resources', self._wf_option_default['resources'])
        self.ctx.max_wallclock_seconds = options_dict.get('max_wallclock_seconds', self._wf_option_default['max_wallclock_seconds'])
        self.ctx.queue = options_dict.get('queue_name', self._wf_option_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._wf_option_default['custom_scheduler_commands'])
        self.ctx.BS_params_dict = wf_dict
        self.ctx.BS_kkrparams = None # is set in set_params_BS
        self.ctx.BS_kpoints = None
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)
        self.report('INFO: use the following parameter:\n'
                    'withmpi: {}\n'
                    'Resources: {}\n'
                    'Walltime (s): {}\n'
                    'queue name: {}\n'
                    'scheduler command: {}\n'
                    'description: {}\n'
                    'label: {}\n'

                    'BS_params: {}\n'.format(self.ctx.withmpi, self.ctx.resources, self.ctx.max_wallclock_seconds,
                                              self.ctx.queue, self.ctx.custom_scheduler_commands,
                                              self.ctx.description_wf, self.ctx.label_wf,
                                              self.ctx.BS_params_dict))
        self.ctx.successful = True
        self.ctx.errors = []

    def validate_input(self):
        """
        # validate input and find out which path (1, or 2) to take
        # return True means run voronoi if false run kkr directly
        """
        inputs = self.inputs
        if 'remote_data' in inputs:
            input_ok = True
        else:
            input_ok = False
            return self.exit_codes.ERROR_NO_INPUT_REMOTE_DATA
        input_remote = self.inputs.remote_data
        parents = input_remote.get_incoming(node_class=CalcJobNode)
        nparents = len(parents.all_link_labels())

        if nparents!=1:
            # extract parent workflow and get uuid of last calc from output node
            parent_workflow = input_remote.inputs.last_RemoteData
            if not isinstance(parent_workflow, WorkChainNode):
                raise InputValidationError("Input remote_data node neither output of a KKR/voronoi calculation nor of kkr_scf_wc workflow")

            parent_workflow_out = parent_workflow.outputs.output_kkr_scf_wc_ParameterResults
            uuid_last_calc = parent_workflow_out.get_dict().get('last_calc_nodeinfo').get('uuid')
            last_calc = load_node(uuid_last_calc)

            if not isinstance(last_calc, KkrCalculation) and not isinstance(last_calc, VoronoiCalculation):
                raise InputValidationError("Extracted last_calc node not of type KkrCalculation: check remote_data input node")

            # overwrite remote_data node with extracted remote folder
            output_remote = last_calc.outputs.remote_folder

            self.inputs.remote_data = output_remote

        struc_kkr, remote_voro = VoronoiCalculation.find_parent_structure(self.inputs.remote_data)
        kpts = get_explicit_kpoints_path(struc_kkr).get('explicit_kpoints')
        self.ctx.BS_kpoints = kpts
        if 'kkr' in inputs:
            try:
                test_and_get_codenode(inputs.kkr, 'kkr.kkr', use_exceptions=True)
            except ValueError:
                input_ok = False
                return self.exit_codes.ERROR_KKRCODE_NOT_CORRECT

        # set self.ctx.input_params_KKR
        self.ctx.input_params_KKR = get_parent_paranode(self.inputs.remote_data)
        return input_ok

    def set_params_BS(self):

        params = self.ctx.input_params_KKR
        input_dict = params.get_dict()
        para_check = kkrparams()
        try:
            for key, val in input_dict.items():
                para_check.set_value(key, val, silent=True)
        except:
            return self.exit_codes.ERROR_CALC_PARAMETERS_INVALID
        label = ''
        descr = '(pk - {}, and uuid - {})'.format(self.inputs.remote_data.pk, self.inputs.remote_data.uuid)

        missing_list = para_check.get_missing_keys(use_aiida=True)

        if missing_list != []:
            kkrdefaults = kkrparams.get_KKRcalc_parameter_defaults()[0]
            kkrdefaults_updated = []
            for key_default, val_default in list(kkrdefaults.items()):
                if key_default in missing_list:
                    para_check.set_value(key_default, val_default)#kkrdefaults.get(key_default), silent=True)
                    kkrdefaults_updated.append(key_default)
                    missing_list.remove(key_default)
            if len(missing_list)>0:
                self.report('ERROR: calc_parameters misses keys: {}'.format(missing_list))
                return self.exit_codes.ERROR_CALC_PARAMETERS_INCOMPLETE
            else:
                self.report('updated KKR parameter node with default values: {}'.format(kkrdefaults_updated))
                label = 'add_defaults_'
                descr = 'added missing default keys, '
##+++ Starts to add the NTP2, EMAX and EMIN from the DOS
        econt_new = self.ctx.BS_params_dict
        kkr_calc = self.inputs.remote_data.get_incoming().first().node
        ef = kkr_calc.outputs.output_parameters.get_dict()['fermi_energy']
        self.ctx.fermi_energy = ef ## in ev unit
        try:
            for key, val in econt_new.items():
                if key=='kmesh':
                    key = 'BZDIVIDE'
                    para_check.set_value(key, val, silent=True)
                elif key=='nepts':
                    key = 'NPT2'
                    # add IEMXD which has to be big enough
                    print('setting IEMXD', val)
                    para_check.set_value('IEMXD', val, silent=True)
                elif key=='emin':
                    key = 'EMIN'
                    new_val = ef - ev_To_Ry(val)
                    para_check.set_value(key, new_val, silent=True)

                    para_check.set_value(key, new_val, silent=True)
                elif key=='emax':
                    key = 'EMAX'
                    new_val = ef + ev_To_Ry(val)
                    para_check.set_value(key, new_val, silent=True)
                elif key=='tempr':
                    key = 'TEMPR'
                    para_check.set_value(key, val, silent=True)

        # set params
        except:
            return self.exit_codes.ERROR_DOS_PARAMS_INVALID
        para_check.set_multiple_values(NPT1=0, NPT3=0, NPOL=0,)
        updatenode = Dict(dict=para_check.get_dict())
        updatenode.label = label+'KKRparam_BS'
        updatenode.description = 'KKR parameter node extracted from remote_folder' + descr + ' as well as wf_parameter input node.'

        paranode_BS = update_params_wf(self.ctx.input_params_KKR, updatenode)
        self.ctx.BS_kkrparams = paranode_BS

    def get_BS(self):
        label = 'KKR BS calc.'
        BS_dict = self.ctx.BS_params_dict
        key_list = list(BS_dict)
        description = 'User defined BandStructure parameters '
        for index in range(len(key_list)):
            description += '{}= {} ,'.format(key_list[index], BS_dict[key_list[index]])

        code = self.inputs.kkr
        remote = self.inputs.remote_data
        params = self.ctx.BS_kkrparams
        kpoints = self.ctx.BS_kpoints
        options = {"max_wallclock_seconds": self.ctx.max_wallclock_seconds,
                   "resources": self.ctx.resources,
                   "queue_name" : self.ctx.queue}#,
        if self.ctx.custom_scheduler_commands:
            options["custom_scheduler_commands"] = self.ctx.custom_scheduler_commands

        create_builder = KkrCalculation.get_builder()
        create_builder.parent_folder = remote
        create_builder.code = code
        create_builder.parameters = params
        create_builder.kpoints = kpoints
        create_builder.metadata.options = options

        BS_run = self.submit(KkrCalculation, **create_builder)
        self.ctx.last_calc = BS_run


        return ToContext(BS_run=BS_run)

    def return_results(self):
        """
        Collect results, parse DOS output and link output nodes to workflow node
        """

        # check wether or not calculation was taked from cached node
        caching_info = "INFO: cache_source of dos calc node: {}".format(self.ctx.BS_run.get_cache_source)
        print(caching_info)
        self.report(caching_info)

        if not self.ctx.BS_run.is_finished_ok:
            self.ctx.successful = False
            error = ('ERROR: DOS calculation failed somehow it is '
                    'in state {}'.format(self.ctx.BS_run.process_state))
            print(error)
            from pprint import pprint
            pprint(self.ctx.BS_run.attributes)
            print('stdout', self.ctx.BS_run.get_scheduler_stdout())
            print('stderr', self.ctx.BS_run.get_scheduler_stderr())
            self.report(error)
            self.ctx.errors.append(error)
            return self.exit_codes.ERROR_BS_CALC_FAILED

        # create dict to store results of workflow output
        outputnode_dict = {}
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._wf_version
        outputnode_dict['withmpi'] = self.ctx.withmpi
        outputnode_dict['resources'] = self.ctx.resources
        outputnode_dict['max_wallclock_seconds'] = self.ctx.max_wallclock_seconds
        outputnode_dict['queue_name'] = self.ctx.queue
        outputnode_dict['custom_scheduler_commands'] = self.ctx.custom_scheduler_commands
        outputnode_dict['BS_params'] = self.ctx.BS_params_dict
        try:
            outputnode_dict['nspin'] = self.ctx.BS_run.res.nspin
        except:
            error = "ERROR: nspin not extracted"
            self.report(error)
            self.ctx.successful = False
            self.ctx.errors.append(error)
        outputnode_dict['successful'] = self.ctx.successful
        outputnode_dict['list_of_errors'] = self.ctx.errors

        # create output node with data-provenance
        outputnode = Dict(dict=outputnode_dict)
        outputnode.label = 'kkr_scf_wc_results'
        outputnode.description = ''

        self.report("INFO: create Banstructure results nodes")
        try:
            self.report("INFO: create Bandstructure results nodes. BS calc retrieved node={}".
                        format(self.ctx.BS_run.outputs.retrieved))
            has_BS_run = True
        except AttributeError as e:
            self.report("ERROR: no Bandstructure calc retrieved node found")
            self.report("Caught AttributeError {}".format(e))
            has_BS_run = False
        
        if has_BS_run:
            BS_retrieved = self.ctx.BS_run.outputs.retrieved
            
        ef = self.ctx.fermi_energy

        outdict = {}
        if has_BS_run:
            ArraData = parse_BS_data(BS_retrieved, Float(ef))
            outdict['BS_data'] = ArraData['BS_data']
        
        link_nodes = outdict.copy() # link also dos output nodes
        if has_dosrun: link_nodes['BS_calc_remote'] =self.ctx.BS_run.outputs.remote_folder
        outdict['results_wf'] = create_out_dict_node(outputnode, **link_nodes)

    # create links to output nodes
        for link_name, node in outdict.items():
            self.out(link_name, node)

        self.report("INFO: done with BS_workflow!\n")

    def ev_To_Ry(value):
        return float(value/13.6056980659)
    

@calcfunction
def parse_BS_data(retrieved_folder, ef):
    from masci_tools.io.common_functions import get_Ry2eV
    from aiida.plugins import DataFactory

    eVscale = get_Ry2eV   
    
    retrieved_list = retrieved_folder.list_object_names()
    qdos_file_list = [i for i in retrieved_list if 'qdos.' in i]
    q_vec_file = 'qvec.dat'
    if q_vec_file in retrieved_list:
        file_open = retrieved_folder.open(q_vec_file)
        q_vec = file_open.readlines()
        file_open.close()
    else:
        print('{} file is not present in the retrived folder'.format(q_vec_file))

    num_dos_files = len(qdos_file_list)
    total_dos = np.loadtxt(qdos_file_list[0])

    no_q_vec = int(q_vec[0])

    for i in qdos_file_list[1:]:
        i_atom = i.split('.')[1]
        load_file = np.loadtxt(i)
        total_dos[:,5:] += load_file[:,5:]    
    
    ef = ef.value
    total_dos[:,0] = total_dos[:,0]*eVscale - ef*np.ones(np.shape(total_dos[:,0]))
    eng_points = set(total_dos[:,0])
    no_eng_points = len(eng_points)

    dos_intensity = np.ndarray(shape = (no_eng_points, no_q_vec))
    for ne in range(np.shape(dos_intensity)[0]):
        nk = np.shape(dos_intensity)[1]

        dos_intensity[ne,:] = total_dos[ ne*nk : (ne+1)*nk, 5 ]/eVscale 

    dos_intensity = dos_intensity.T # setting eng-kpts corresponding x-y asix
    q_vec = np.asarray(q_vec) # converting q_vec into array
    eng_points = np.asarray(eng_points) # converting eng_popints into array

    # To save into the ArrayData
    ArrayData = DataFactory('array')
    array = ArrayData()
    array.set_array('BS_Matrix', dos_intensity)
    array.set_array('Kpts', q_vec)
    array.set_array('Eng_levels',list(eng_points))

    return {'BS_Data' : array}

                                         


# In[20]:


##++ Starts submision of kkr_BS_wf
from aiida import load_profile, get_profile
load_profile()
print(get_profile().name)
from aiida.engine import run
from aiida.orm import Code, load_node, load_code
from aiida.plugins import DataFactory

Dict = DataFactory('dict')
metadata_option_1 = {'max_wallclock_seconds': 36000,'resources':
              {'tot_num_mpiprocs': 48, 'num_machines': 1},
              'custom_scheduler_commands':
              '#SBATCH --account=jara0191\n\nulimit -s unlimited; export OMP_STACKSIZE=2g',
              'withmpi': True}
Dict_metadata_option_1 = Dict(dict=metadata_option_1)

calc_dos = load_node(203)
calc_kkr_converged = load_node(180)
remote_folder_kkr = calc_kkr_converged.outputs.remote_folder
# remote_folder_dos = calc_dos.outputs.remote_folder

# Parameters_from_DOS_calc_dict = calc_dos.inputs.parameters
# Parameters_from_DOS_calc = remote_folder_dos.get_incoming().first().node.inputs.parameters
wf_parameters_dict = {'NPT2' : 200,'EMIN' : 10, 'EMAX' : 5, 'TEMPR' : 50, 'kmesh' : [30,30,30]}
wf_parameters = Dict(dict=wf_parameters_dict)#calc_kkr_converged.inputs.parameters#.get_dict()
kkr = Code.get_from_string('kkrhost@claix18')
inputs = { 'wf_parameters':wf_parameters, 'options':Dict_metadata_option_1, 'remote_data':remote_folder_kkr,
           'kkr':kkr }
submition_kkr_BS_wf = run(kkr_BS_wf, **inputs)
##-- Ends submision of the kkr_BS_wf
                                     


# In[ ]:




