#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module you find the base workflow for a impurity DOS calculation and
some helper methods to do so with AiiDA
"""
from __future__ import print_function, absolute_import
from aiida.orm import Code, load_node, CalcJobNode
from aiida.plugins import DataFactory
from aiida.engine import if_, ToContext, WorkChain
from aiida.common import LinkType
from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc

__copyright__ = (u"Copyright (c), 2019, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.4"
__contributors__ = (u"Fabian Bertoldo", u"Philipp Ruessmann")

#TODO: improve workflow output node structure
#TODO: generalise search for imp_info and conv_host from startpot


Dict = DataFactory('dict')
RemoteData = DataFactory('remote')
SinglefileData = DataFactory('singlefile')


class kkr_imp_dos_wc(WorkChain):
    """
    Workchain of a DOS calculation for an impurity system starting from a
    converged impurity calculation or workflow

    :param options: (Dict), computer options
    :param wf_parameters: (Dict), specifications for the DOS
    :param kkr: (Code), mandatory: KKR code for gf_writeout step
    :param kkrimp: (Code), mandatory: KKRimp code for DOS calculation
    :param imp_host_pot: (SinglefileData), mandatory: impurity startpotential

    :return workflow_info: (Dict), Information on workflow results
    :return last_calc_output_parameters: (Dict), output parameters of
                                         the last called calculation
    :return last_calc_info: (Dict), information of the last called calculation
    """

    _workflowversion = __version__
    _wf_label = 'kkr_imp_dos_wc'
    _wf_description = 'Workflow for a KKR impurity DOS calculation'


    _options_default = {'queue_name' : '',                        # Queue name to submit jobs too
                        'resources': {"num_machines": 1},         # resources to allocate for the job
                        'max_wallclock_seconds' : 60*60,          # walltime after which the job gets killed (gets parsed to KKR)}
                        'custom_scheduler_commands' : '',         # some additional scheduler commands
                        'use_mpi' : True}                         # execute KKR with mpi or without

    _wf_default = {'ef_shift': 0. ,                                  # set costum absolute E_F (in eV)
                   'dos_params': {'nepts': 61,                       # DOS params: number of points in contour
                                  'tempr': 200, # K                  # DOS params: temperature
                                  'emin': -1, # Ry                   # DOS params: start of energy contour
                                  'emax': 1,  # Ry                   # DOS params: end of energy contour
                                  'kmesh': [30, 30, 30]},            # DOS params: kmesh for DOS calculation (typically higher than in scf contour)
                   'non_spherical': 1,                      # use non-spherical parts of the potential (0 if you don't want that)
                   'spinorbit' : False,                     # SOC calculation (True/False)
                   'newsol' : False }                       # new SOC solver is applied



    @classmethod
    def get_wf_defaults(self):
        """
        Print and return _wf_defaults dictionary. Can be used to easily create set of wf_parameters.
        returns _wf_defaults
        """

        print('Version of workflow: {}'.format(self._workflowversion))
        return self._wf_default


    @classmethod
    def define(cls, spec):
        """
        Defines the outline of the workflow
        """

        # Take input of the workflow or use defaults defined above
        super(kkr_imp_dos_wc, cls).define(spec)

        spec.input("kkr", valid_type=Code, required=True)
        spec.input("kkrimp", valid_type=Code, required=True)
        spec.input("host_imp_pot", valid_type=SinglefileData, required=True)
        spec.input("impurity_info", valid_type=Dict, required=False)
        spec.input("host_remote", valid_type=RemoteData, required=False)
        spec.input("options", valid_type=Dict, required=False,
                       default=Dict(dict=cls._options_default))
        spec.input("wf_parameters", valid_type=Dict, required=False,
                       default=Dict(dict=cls._wf_default))

        # Here the structure of the workflow is defined
        spec.outline(
            cls.start,                  # start and initialise workflow
            if_(cls.validate_input)(    # validate the given input
                cls.run_gfstep,         # run GF step with DOS energy contour
                cls.run_imp_dos),       # run DOS for the impurity problem
            cls.return_results          # terminate workflow and return results
            )

        # Define possible exit codes for the workflow
        spec.exit_code(220, "ERROR_UNKNOWN_PROBLEM", 
            message="Unknown problem detected.")
        spec.exit_code(221, "ERROR_NO_PARENT_FOUND", 
            message="Unable to find the parent remote_data node that led to "
                    "the input impurity calculation. You need to specify "
                    "`host_remote` and `impurity_info` nodes.")
        spec.exit_code(222, "ERROR_GF_WRITEOUT_UNSUCCESFUL", 
            message="The gf_writeout workflow was not succesful, cannot continue.")

        # specify the outputs
        spec.output('workflow_info', valid_type=Dict)
        spec.output('last_calc_output_parameters', valid_type=Dict)
        spec.output('last_calc_info', valid_type=Dict)


    def start(self):
        """
        Initialise context and some parameters
        """

        self.report('INFO: started KKR impurity DOS workflow version {}'
                    ''.format(self._workflowversion))

        # input both wf and options parameters
        if 'wf_parameters' in self.inputs:
            wf_dict = self.inputs.wf_parameters.get_dict()
            if wf_dict == {}:
                wf_dict = self._wf_default
                self.report('INFO: using default wf parameters')

        if 'options' in self.inputs:
            options_dict = self.inputs.options.get_dict()
            if options_dict == {}:
                options_dict = self._options_default
                self.report('INFO: using default options parameters')


        # set values, or defaults
        self.ctx.use_mpi = options_dict.get('use_mpi', self._options_default['use_mpi'])
        self.ctx.resources = options_dict.get('resources', self._options_default['resources'])
        self.ctx.walltime_sec = options_dict.get('max_wallclock_seconds', self._options_default['max_wallclock_seconds'])
        self.ctx.queue = options_dict.get('queue_name', self._options_default['queue_name'])
        self.ctx.custom_scheduler_commands = options_dict.get('custom_scheduler_commands', self._options_default['custom_scheduler_commands'])
        self.ctx.options_params_dict = Dict(dict={'use_mpi': self.ctx.use_mpi, 'resources': self.ctx.resources,
                                                           'max_wallclock_seconds': self.ctx.walltime_sec, 'queue_name': self.ctx.queue,
                                                           'custom_scheduler_commands': self.ctx.custom_scheduler_commands})

        # set workflow parameters for the KKR imputrity calculations
        self.ctx.ef_shift = wf_dict.get('ef_shift', self._wf_default['ef_shift'])
        self.ctx.dos_params_dict = wf_dict.get('dos_params', self._wf_default['dos_params'])

        # set workflow parameters for the KKR impurity calculation
        self.ctx.nsteps = 1 # always only one step for DOS calculation
        self.ctx.kkr_runmax = 1 # no restarts for DOS calculation
        self.ctx.non_spherical = wf_dict.get('non_spherical', self._wf_default['non_spherical'])
        self.ctx.spinorbit = wf_dict.get('spinorbit', self._wf_default['spinorbit'])
        self.ctx.newsol = wf_dict.get('newsol', self._wf_default['newsol'])

        # set workflow label and description
        self.ctx.description_wf = self.inputs.get('description', self._wf_description)
        self.ctx.label_wf = self.inputs.get('label', self._wf_label)


        self.report('INFO: use the following parameter:\n'
                    'use_mpi: {}\n'
                    'Resources: {}\n'
                    'Walltime (s): {}\n'
                    'queue name: {}\n'
                    'scheduler command: {}\n'
                    'description: {}\n'
                    'label: {}\n'.format(self.ctx.use_mpi, self.ctx.resources, self.ctx.walltime_sec,
                                         self.ctx.queue, self.ctx.custom_scheduler_commands,
                                         self.ctx.description_wf, self.ctx.label_wf))

        # return para/vars
        self.ctx.successful = True
        self.ctx.errors = []
        self.ctx.formula = ''


    def validate_input(self):
        """
        Validate input and catch possible errors
        """

        inputs = self.inputs
        inputs_ok = True

        if 'host_imp_pot' in inputs:
            # check if input potential has incoming return link
            if len(inputs.host_imp_pot.get_incoming(link_type=LinkType.RETURN).all()) < 1:
                self.report("input potential not from kkrimp workflow: take remote_data folder of host system from input")
                if 'impurity_info' in inputs and 'host_remote' in inputs:
                    self.ctx.imp_info = inputs.impurity_info
                    self.ctx.conv_host_remote = inputs.host_remote
                else:
                    self.report('WARNING: startpot has no parent and can not find '
                                'a converged host RemoteData node')
                    if 'impurity_info' not in inputs:
                        self.report("`impurity_info` optional input node not given but needed in this case.")
                    if 'host_remote' not in inputs:
                        self.report("`host_remote` optional input node not given but needed in this case.")
                    inputs_ok = False
                    self.ctx.errors.append(1)
            else:
                # if return ink is found get input nodes automatically
                self.report('INFO: get converged host RemoteData node and '
                            'impurity_info node from database')
                self.ctx.kkr_imp_wf = inputs.host_imp_pot.created_by.called_by
                self.report('INFO: found underlying kkr impurity workflow '
                            '(pk: {})'.format(self.ctx.kkr_imp_wf.pk))
                self.ctx.imp_info = self.ctx.kkr_imp_wf.inputs.impurity_info
                self.report('INFO: found impurity_info node (pk: {})'.format(
                            self.ctx.imp_info.pk))
                try:
                    self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inputs.remote_data_gf.inputs.remote_folder.inputs.parent_calc_folder.inputs.remote_folder.outputs.remote_folder
                    self.report('INFO: imported converged_host_remote (pk: {}) and '
                                'impurity_info from database'.format(self.ctx.conv_host_remote.pk))
                except:
                    self.ctx.conv_host_remote = self.ctx.kkr_imp_wf.inputs.gf_remote.inputs.remote_folder.inputs.parent_calc_folder.inputs.remote_folder.outputs.remote_folder
                    self.report('INFO: imported converged_host_remote (pk: {}) and '
                                'impurity_info from database'.format(self.ctx.conv_host_remote.pk))

        self.report('INFO: validated input successfully: {}'.format(inputs_ok))

        return inputs_ok


    def run_gfstep(self):
        """
        Start GF writeout step with DOS energy contour
        """

        options = self.ctx.options_params_dict
        kkrcode = self.inputs.kkr
        converged_host_remote = self.ctx.conv_host_remote
        imp_info = self.ctx.imp_info

        wf_params_gf = Dict(dict={'ef_shift':self.ctx.ef_shift, 'dos_run':True,
                                           'dos_params':self.ctx.dos_params_dict})
        label_gf = 'GF writeout for imp DOS'
        description_gf = 'GF writeout step with energy contour for impurity DOS'

        builder = kkr_flex_wc.get_builder()
        builder.metadata.label = label_gf
        builder.metadata.description = description_gf
        builder.kkr = kkrcode
        builder.options = options
        builder.wf_parameters = wf_params_gf
        builder.remote_data = converged_host_remote
        builder.impurity_info = imp_info

        future = self.submit(builder)

        self.report('INFO: running GF writeout (pid: {})'.format(future.pk))

        return ToContext(gf_writeout=future)


    def run_imp_dos(self):
        """
        Use previous GF step to calculate DOS for the impurity problem
        """

        if not self.ctx.gf_writeout.is_finished_ok:
            return self.exit_codes.ERROR_GF_WRITEOUT_UNSUCCESFUL
            

        options = self.ctx.options_params_dict
        kkrimpcode = self.inputs.kkrimp
        gf_writeout_wf = self.ctx.gf_writeout
        gf_writeout_calc = load_node(self.ctx.gf_writeout.outputs.workflow_info.get_dict().get('pk_flexcalc'))
        gf_writeout_remote = gf_writeout_wf.outputs.GF_host_remote
        impurity_pot = self.inputs.host_imp_pot
        imps = self.ctx.imp_info

        nspin = gf_writeout_calc.outputs.output_parameters.get_dict().get('nspin')
        self.ctx.nspin = nspin
        self.report('nspin: {}'.format(nspin))
        self.ctx.kkrimp_params_dict = Dict(dict={'nspin': nspin,
                                                          'nsteps': self.ctx.nsteps,
                                                          'kkr_runmax': self.ctx.kkr_runmax, 'non_spherical': self.ctx.non_spherical,
                                                          'spinorbit': self.ctx.spinorbit, 'newsol': self.ctx.newsol,
                                                          'dos_run': True})
        kkrimp_params = self.ctx.kkrimp_params_dict

        label_imp = 'KKRimp DOS (GF: {}, imp_pot: {}, Zimp: {}, ilayer_cent: {})'.format(
                    gf_writeout_wf.pk, impurity_pot.pk, imps.get_dict().get('Zimp'), imps.get_dict().get('ilayer_center'))
        description_imp = 'KKRimp DOS run (GF: {}, imp_pot: {}, Zimp: {}, ilayer_cent: {}, R_cut: {})'.format(
                    gf_writeout_wf.pk, impurity_pot.pk, imps.get_dict().get('Zimp'), imps.get_dict().get('ilayer_center'),
                    imps.get_dict().get('Rcut'))

        builder = kkr_imp_sub_wc.get_builder()
        builder.metadata.label = label_imp
        builder.metadata.description = description_imp
        builder.kkrimp = kkrimpcode
        builder.options = options
        builder.wf_parameters = kkrimp_params
        builder.remote_data = gf_writeout_remote
        builder.host_imp_startpot = impurity_pot

        future = self.submit(builder)

        self.report('INFO: running DOS step for impurity system (pid: {})'.format(future.pk))

        return ToContext(kkrimp_dos=future)


    def return_results(self):
        """
        Return the results and create all of the output nodes
        """

        if self.ctx.errors!=[]:
            if self.ctx.errors[0]==1:
                return self.exit_codes.ERROR_NO_PARENT_FOUND
            else:
                return self.exit_codes.ERROR_UNKNOWN_PROBLEM

        self.report('INFO: creating output nodes for the KKR imp DOS workflow ...')

        last_calc_pk = self.ctx.kkrimp_dos.outputs.workflow_info.get_dict().get('last_calc_nodeinfo')['pk']
        last_calc = load_node(last_calc_pk)
        last_calc_output_params = last_calc.outputs.output_parameters
        last_calc_info = self.ctx.kkrimp_dos.outputs.workflow_info

        outputnode_dict = {}
        outputnode_dict['impurity_info'] = self.ctx.imp_info.get_dict()
        outputnode_dict['workflow_name'] = self.__class__.__name__
        outputnode_dict['workflow_version'] = self._workflowversion
        outputnode_dict['used_subworkflows'] = {'gf_writeout': self.ctx.gf_writeout.pk,
                                                'impurity_dos': self.ctx.kkrimp_dos.pk}
        outputnode_t = Dict(dict=outputnode_dict)
        outputnode_t.label = 'kkr_imp_dos_wc_inform'
        outputnode_t.description = 'Contains information for workflow'

        # interpol dos file and store to XyData nodes
        dos_retrieved = last_calc.outputs.retrieved
        if 'out_ldos.interpol.atom=01_spin1.dat' in dos_retrieved.list_object_names():
            kkrflex_writeout = load_node(self.ctx.gf_writeout.outputs.workflow_info.get_dict().get('pk_flexcalc'))
            parent_calc_kkr_converged = kkrflex_writeout.inputs.parent_folder.get_incoming(node_class=CalcJobNode).first().node
            ef = parent_calc_kkr_converged.outputs.output_parameters.get_dict().get('fermi_energy')
            natom = last_calc_output_params.get_dict().get('number_of_atoms_in_unit_cell')
            dosXyDatas = parse_impdosfiles(dos_retrieved, natom, self.ctx.nspin, ef)
            dos_extracted = True
        else:
            dos_extracted = False
        if dos_extracted:
            self.out('dos_data', dosXyDatas[0])
            self.out('dos_data_interpol', dosXyDatas[1])
 

        self.report('INFO: workflow_info node: {}'.format(outputnode_t.uuid))

        self.out('workflow_info', outputnode_t)
        self.out('last_calc_output_parameters', last_calc_output_params)
        self.out('last_calc_info', last_calc_info)

        self.report('INFO: created output nodes for KKR imp DOS workflow.')
        self.report('\n'
                    '|------------------------------------------------------------------------------------------------------------------|\n'
                    '|-------------------------------------| Done with the KKR imp DOS workflow! |--------------------------------------|\n'
                    '|------------------------------------------------------------------------------------------------------------------|')

def parse_impdosfiles(dos_retrieved, natom, nspin, ef):
    from masci_tools.io.common_functions import get_Ry2eV, get_ef_from_potfile
    from numpy import loadtxt, array
    from aiida.plugins import DataFactory
    XyData = DataFactory('array.xy')

    # read dos files
    dos, dos_int = [], []
    for iatom in range(1, natom+1):
        for ispin in range(1, nspin+1):
            with dos_retrieved.open('out_ldos.atom=%0.2i_spin%i.dat'%(iatom, ispin)) as dosfile:
                tmp = loadtxt(dosfile)
                dos.append(tmp)
            with dos_retrieved.open('out_ldos.interpol.atom=%0.2i_spin%i.dat'%(iatom, ispin)) as dosfile:
                tmp = loadtxt(dosfile)
                dos_int.append(tmp)
    dos, dos_int = array(dos), array(dos_int)

    # convert to eV units
    eVscale = get_Ry2eV()
    dos[:,:,0] = (dos[:,:,0]-ef)*eVscale
    dos[:,:,1:] = dos[:,:,1:]/eVscale
    dos_int[:,:,0] = (dos_int[:,:,0]-ef)*eVscale
    dos_int[:,:,1:] = dos_int[:,:,1:]/eVscale

    # create output nodes
    dosnode = XyData()
    dosnode.label = 'dos_data'
    dosnode.description = 'Array data containing uniterpolated DOS (i.e. dos at finite imaginary part of energy). 3D array with (atoms, energy point, l-channel) dimensions.'
    dosnode.set_x(dos[:,:,0], 'E-EF', 'eV')

    name = ['tot', 's', 'p', 'd', 'f', 'g']
    name = name[:len(dos[0,0,1:])-1]+['ns']

    ylists = [[],[],[]]
    for l in range(len(name)):
        ylists[0].append(dos[:,:,1+l])
        ylists[1].append('dos '+name[l])
        ylists[2].append('states/eV')
    dosnode.set_y(ylists[0], ylists[1], ylists[2])

    # node for interpolated DOS
    dosnode2 = XyData()
    dosnode2.label = 'dos_interpol_data'
    dosnode2.description = 'Array data containing iterpolated DOS (i.e. dos at finite imaginary part of energy). 3D array with (atoms, energy point, l-channel) dimensions.'
    dosnode2.set_x(dos_int[:,:,0], 'E-EF', 'eV')
    ylists = [[],[],[]]
    for l in range(len(name)):
        ylists[0].append(dos_int[:,:,1+l])
        ylists[1].append('interpolated dos '+name[l])
        ylists[2].append('states/eV')
    dosnode2.set_y(ylists[0], ylists[1], ylists[2])

    return dosnode, dosnode2
