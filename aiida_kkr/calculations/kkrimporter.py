# -*- coding: utf-8 -*-
"""
Plug-in to import a KKR calculation. This is based on the PwImmigrantCalculation of the aiida-quantumespresso plugin.
"""

from __future__ import absolute_import
import os
from aiida.plugins import DataFactory
from aiida.engine.calculation.job import _input_subfolder
from aiida.common.utils import classproperty
from aiida.common.folders import SandboxFolder
from aiida.common.exceptions import InputValidationError, InvalidOperation
from aiida.common.datastructures import calc_states
from aiida.common.links import LinkType
from .kkr import KkrCalculation
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.tools.common_workfunctions import structure_from_params
from six.moves import range


#define aiida structures from DataFactory of aiida
RemoteData = DataFactory('remote')
Dict = DataFactory('dict')
StructureData = DataFactory('structure')


__copyright__ = (u"Copyright (c), 2017, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.2"
__contributors__ = ("Philipp Rüßmann")



class KkrImporterCalculation(KkrCalculation):
    """
    Importer dummy calculation for a previous KKR run

    :param remote_workdir: Absolute path to the directory where the job was run.
        The transport of the computer you link ask input to the calculation is
        the transport that will be used to retrieve the calculation's files.
        Therefore, ``remote_workdir`` should be the absolute path to the job's
        directory on that computer.
    :type remote_workdir: str
    :param input_file_names: The file names of the job's input file.
    :type input_file_name: dict with str entries
    :param output_file_name: The file names of the job's output file (i.e. the
        file containing the stdout of kkr.x).
    :type output_file_name: dict with str entries
    """

    def _init_internal_params(self):
        """
        Init internal parameters at class load time
        """
        # reuse base class function
        super(KkrImporterCalculation, self)._init_internal_params()

        # calculation plugin version
        self._CALCULATION_PLUGIN_VERSION = __version__

        # parser
        self._default_parser = 'kkr.kkrimporterparser'

    @classproperty
    def _use_methods(cls):
        """
        Add use_structure method for KKRimporter calculations.
        """
        use_dict = KkrCalculation._use_methods
        use_dict.update({
            "structure": {
                'valid_types': StructureData,
                'additional_parameter': None,
                'linkname': 'structure',
                'docstring':
                ("Use a node that specifies the input crystal structure ")
                }
            })
        return use_dict

    def create_input_nodes(self, open_transport, input_file_names=None,
                           output_file_names=None, remote_workdir=None):
        """
        Create calculation input nodes based on the job's files.

        :param open_transport: An open instance of the transport class of the calculation's computer. See the tutorial for more information.
        :type open_transport: aiida.transport.plugins.local.LocalTransport or aiida.transport.plugins.ssh.SshTransport

        This method parses the files in the job's remote working directory to
        create the input nodes that would exist if the calculation were
        submitted using AiiDa. These nodes are:
        * a ``'parameters'`` Dict node, based on the namelists and
        their variable-value pairs;
        * ...;
        and can be retrieved as a dictionary using the ``get_inputs_dict()``
        method. *These input links are cached-links; nothing is stored by this
        method (including the calculation node itself).*

        **Keyword arguments**
        .. note:: These keyword arguments can also be set when instantiating the
        class or using the ``set_`` methods (e.g. ``set_remote_workdir``).
        Offering to set them here simply offers the user an additional
        place to set their values. *Only the values that have not yet been
        set need to be specified.*

        :param input_file_names: The file name of the job's input files (inputcard, shapefun).
        :type input_file_names: dict of str values
        :param output_file_names: The file name of the job's output files (i.e.
            the files containing the stdout of KKR, the output.000.txt, output.0.txt, output.2.txt, out_timing.txt, output potential, output nonco_angle file).
        :type output_file_names: dict of str values
        :param remote_workdir: Absolute path to the directory where the job
            was run. The transport of the computer you link ask input to the
            calculation is the transport that will be used to retrieve the
            calculation's files. Therefore, ``remote_workdir`` should be the
            absolute path to the job's directory on that computer.
        :type remote_workdir: str
        :raises aiida.common.exceptions.InputValidationError: if
            ``open_transport`` is a different type of transport than the
            computer's.
        :raises aiida.common.exceptions.InvalidOperation: if
            ``open_transport`` is not open.
        :raises aiida.common.exceptions.InputValidationError: if
            ``remote_workdir``, ``input_file_names``, and/or ``output_file_names``
            are not set prior to or during the call of this method.
        :raises aiida.common.exceptions.FeatureNotAvailable: if the input file
            uses anything which is not currently implimented in aiida-kkr.
        :raises aiida.common.exceptions.ParsingError: if there are issues
            parsing the input file.
        :raises IOError: if there are issues reading the input file.
        """
        # Make sure the remote workdir and input + output file names were
        # provided either before or during the call to this method. If they
        # were just provided during this method call, store the values.
        if remote_workdir is not None:
            self.set_remote_workdir(remote_workdir)
        elif self.get_attr('remote_workdir', None) is None:
            raise InputValidationError(
                'The remote working directory has not been specified.\n'
                'Please specify it using one of the following...\n '
                '(a) pass as a keyword argument to create_input_nodes\n'
                '    [create_input_nodes(remote_workdir=your_remote_workdir)]\n'
                '(b) pass as a keyword argument when instantiating\n '
                '    [calc = KkrImporterCalculation(remote_workdir=your_remote_workdir)]\n'
                '(c) use the set_remote_workdir method\n'
                '    [calc.set_remote_workdir(your_remote_workdir)]'
            )

        # check if necessary input file names are given and set filenames from dict
        self.set_input_file_names(input_file_names)

        # set output file names instead of using defaults
        self.set_output_file_names(output_file_names)

        # Check that open_transport is the correct transport type.
        if type(open_transport) is not self.get_computer().get_transport_class():
            raise InputValidationError(
                "The transport passed as the `open_transport` parameter is "
                "not the same transport type linked to the computer. Please "
                "obtain the correct transport class using the "
                "`get_transport_class` method of the calculation's computer. "
                "See the tutorial for more information."
            )

        # Check that open_transport is actually open.
        if not open_transport._is_open:
            raise InvalidOperation(
                "The transport passed as the `open_transport` parameter is "
                "not open. Please execute the open the transport using it's "
                "`open` method, or execute the call to this method within a "
                "`with` statement context guard. See the tutorial for more "
                "information."
            )

        # Copy the input file and psuedo files to a temp folder for parsing.
        with SandboxFolder() as folder:

            # Copy the input file to the temp folder.
            remote_path = os.path.join(self._get_remote_workdir(),
                                       self._INPUT_FILE_NAME)
            open_transport.get(remote_path, folder.abspath)

            # Parse the input file.
            local_path = os.path.join(folder.abspath, self._INPUT_FILE_NAME)
            parameters = kkrparams(params_type='kkr')
            parameters.read_keywords_from_inputcard(inputcard=local_path)

            # Create Dict node based on the namelist and link as input.
            aiida_parameters = Dict(dict=parameters.get_dict())
            self.use_parameters(aiida_parameters)

            # Create a StructureData node from inputs
            success, structuredata = structure_from_params(parameters)
            if not success:
                raise InputValidationError('Something went wrong when creating the an aiida structure from your input file. Check you input')
            self.use_structure(structuredata)

        self._set_attr('input_nodes_created', True)

    def _prepare_for_retrieval(self, open_transport):
        """
        Prepare the calculation for retrieval by daemon.

        :param open_transport: An open instance of the transport class of the
            calculation's computer.
        :type open_transport: aiida.transport.plugins.local.LocalTransport or
            aiida.transport.plugins.ssh.SshTransport

        Here, we
        * manually set the files to retrieve
        * store the calculation and all it's input nodes
        * copy the input file to the calculation's raw_input_folder in the
        * store the remote_workdir as a RemoteData output node

        """

        # Manually set the files that will be copied to the repository and that
        # the parser will extract the results from. This would normally be
        # performed in self._prepare_for_submission prior to submission.

        natom = 1000 # maximal number of atom-resolved files that are retrieved
        # TODO take actual natom value (maybe extract from number of files that are there)
        self._set_attr('retrieve_list',
                       [self._INPUT_FILE_NAME, # inputcard needed for parsing
                        self._DEFAULT_OUTPUT_FILE, # out_kkr, std shell output
                        self._NONCO_ANGLES_OUT, # nonco angles files
                        self._OUTPUT_0_INIT, self._OUTPUT_000, self._OUTPUT_2, # output files in new(er) style
                        'output.0','output.1a','output.1b','output.1c','output.2', # try to import old style output
                        self._OUT_TIMING_000, # timing file
                        self._OUT_POTENTIAL, self._SHAPEFUN  # make sure to retrieve potential and shapefun as well
                        # add Jij files etc for other run options
                        ] + [self._SHELLS_DAT] + [self._Jij_ATOM%iatom for iatom in range(1,natom+1)]
                        )
        self._set_attr('retrieve_singlefile_list', [])

        # Make sure the calculation and input links are stored.
        self.store_all()

        # Store the original input file in the calculation's repository folder.
        remote_path = os.path.join(self._get_remote_workdir(),
                                   self._INPUT_FILE_NAME)
        raw_input_folder = self.folder.get_subfolder(_input_subfolder,
                                                     create=True)
        open_transport.get(remote_path, raw_input_folder.abspath)

        # Manually add the remote working directory as a RemoteData output
        # node.
        self._set_state(calc_states.SUBMITTING)
        remotedata = RemoteData(computer=self.get_computer(),
                                remote_path=self._get_remote_workdir())
        remotedata.add_link_from(self, label='remote_folder',
                                 link_type=LinkType.CREATE)
        remotedata.store()

    def prepare_for_retrieval_and_parsing(self, open_transport):
        """
        Tell the daemon that the calculation is computed and ready to be parsed.

        :param open_transport: An open instance of the transport class of the
            calculation's computer. See the tutorial for more information.
        :type open_transport: aiida.transport.plugins.local.LocalTransport
            or aiida.transport.plugins.ssh.SshTransport

        The next time the daemon updates the status of calculations, it will
        see this job is in the 'COMPUTED' state and will retrieve its output
        files and parse the results.
        If the daemon is not currently running, nothing will happen until it is
        started again.
        This method also stores the calculation and all input nodes. It also
        copies the original input file to the calculation's repository folder.

        :raises aiida.common.exceptions.InputValidationError: if
            ``open_transport`` is a different type of transport
            than the computer's.
        :raises aiida.common.exceptions.InvalidOperation: if
            ``open_transport`` is not open.
        """

        # Check that the create_input_nodes method has run successfully.
        if not self.get_attr('input_nodes_created', False):
            raise InvalidOperation(
                "You must run the create_input_nodes method before calling "
                "prepare_for_retrieval_and_parsing!"
            )

        # Check that open_transport is the correct transport type.
        if type(open_transport) is not self.get_computer().get_transport_class():
            raise InputValidationError(
                "The transport passed as the `open_transport` parameter is "
                "not the same transport type linked to the computer. Please "
                "obtain the correct transport class using the "
                "`get_transport_class` method of the calculation's computer. "
                "See the tutorial for more information."
            )

        # Check that open_transport is actually open.
        if not open_transport._is_open:
            raise InvalidOperation(
                "The transport passed as the `open_transport` parameter is "
                "not open. Please execute the open the transport using it's "
                "`open` method, or execute the call to this method within a "
                "`with` statement context guard. See the tutorial for more "
                "information."
            )

        # Prepare the calculation for retrieval
        self._prepare_for_retrieval(open_transport)

        # Manually set the state of the calculation to "COMPUTED", so that it
        # will be retrieved and parsed the next time the daemon updates the
        # status of calculations.
        self._set_state(calc_states.COMPUTED)

    def set_remote_workdir(self, remote_workdir):
        """
        Set the job's remote working directory.

        :param remote_workdir: Absolute path of the job's remote working
            directory.
        :type remote_workdir: str
        """
        # This is the functionality as self._set_remote_workir, but it bypasses
        # the need to have the calculation state set as SUBMITTING.
        self._set_attr('remote_workdir', remote_workdir)

    def set_input_file_names(self, input_file_names):
        """
        Set the file names of the job's input file (e.g. ``'inputcard'`` etc.).
        :param input_file_names: The file names of the job's input files (inputcard, potential, etc.).
        :type input_file_names: dict

        Keys of input_file_names dict should be one or more of:

        * ``'input_file'``
        * ``'potential_file'``
        * ``'shapefun_file'``

        :note: The keys 'input_file' and 'potential_file' are mandatory!
        """
        if self._INPUT_FILE_NAMES is None:
            self.logger.info('setting input_file_names: {}'.format(input_file_names))
            self._set_attr('input_file_names', input_file_names)
            self._set_filenames_in(input_file_names)

    def set_output_file_names(self, output_file_names):
        """
        Set the file names of the job's output files (e.g. ``'output.000.txt'`` etc.).

        :param output_file_names: The dictionary of file names of file containing the job's
            outputs.
        :type output_file_names: dict

        Keys of output_file_names dict should be one or more of:

        * ``'out_file'``
        * ``'out_potential_file'``
        * ``'out_nonco_angles_file'``
        * ``'output_0_file'``
        * ``'output_000_file'``
        * ``'output_2_file'``
        * ``'timing_file'``

        :note: In a filename is not given, the default value of aiida_kkr.calculations.kkr are used
        """
        if self._OUTPUT_FILE_NAMES is None:
            if output_file_names is None: #set default values if nothing is chosen
                self._OUT_POTENTIAL = self._POTENTIAL
                self._DEFAULT_OUTPUT_FILE = KkrCalculation()._DEFAULT_OUTPUT_FILE
                self._NONCO_ANGLES_OUT = KkrCalculation()._NONCO_ANGLES_OUT
                self._OUTPUT_0_INIT = KkrCalculation()._OUTPUT_0_INIT
                self._OUTPUT_000 = KkrCalculation()._OUTPUT_000
                self._OUTPUT_2 = KkrCalculation()._OUTPUT_2
                self._OUT_TIMING_000 = KkrCalculation()._OUT_TIMING_000
                output_file_names = {'out_potential_file': self._OUT_POTENTIAL,
                                     'out_file': self._DEFAULT_OUTPUT_FILE,
                                     'out_nonco_angles_file': self._NONCO_ANGLES_OUT,
                                     'output_0_file': self._OUTPUT_0_INIT,
                                     'output_000_file': self._OUT_TIMING_000,
                                     'output_2_file': self._OUTPUT_2,
                                     'timing_file': self._OUT_TIMING_000}
            else:
                self.logger.info('setting output_file_names: {}'.format(output_file_names))
                self._set_filenames_out(output_file_names)
             # finally set ourput_file_names as class attribute
            self._set_attr('output_file_names', output_file_names)

    def _set_filenames_in(self, input_file_names):
        """
        helper function to set filenames from input_file_names dict
        also check if at least inputcard and potential are given
        """
        has_input_files = False
        if input_file_names is not None:
            has_input_files = True
            if input_file_names.get('input_file') is not None:
                self._INPUT_FILE_NAME = input_file_names.get('input_file')
            else:
                self.logger.info('KKRimporter: input_file not given!')
                has_input_files = False
            if input_file_names.get('potential_file') is not None:
                self._POTENTIAL = input_file_names.get('potential_file')
            else:
                self.logger.info('KKRimporter: potential_file not given!')
                has_input_files = False
            if input_file_names.get('shapefun_file') is not None:
                self._SHAPEFUN = input_file_names.get('shapefun_file')
            else:
                self._SHAPEFUN = KkrCalculation()._SHAPEFUN
        if not has_input_files:
            raise InputValidationError(
                'The input file_names has not been specified correctly.\n'
                'Please specify it using one of the following...\n '
                '(a) pass as a keyword argument to create_input_nodes\n'
                '    [create_input_nodes(input_file_names=your_file_names_dictionary)]\n'
                '(b) pass as a keyword argument when instantiating\n '
                '    [calc = KkrImporterCalculation(input_file_names={"input_file":"inputcard", "potential_file":"potential", ...}]\n'
                '(c) use the set_input_file_names method\n'
                '    [calc.set_input_file_name(your_file_names)]'
            )

    def _set_filenames_out(self, output_file_names):
        """
        helper function that sets the output file names from the output_file_names dictionary
        """
        if output_file_names is not None and output_file_names.get('out_potential_file') is not None:
            self._OUT_POTENTIAL = output_file_names.get('out_potential_file')
        else:
            self._OUT_POTENTIAL = self._POTENTIAL
        #raise InputValidationError('test: {}'.format(self._OUT_POTENTIAL))
        if output_file_names is not None:
            if output_file_names.get('out_file') is not None:
                self._DEFAULT_OUTPUT_FILE = output_file_names.get('out_file')
            else:
                self._DEFAULT_OUTPUT_FILE = KkrCalculation()._DEFAULT_OUTPUT_FILE
            if output_file_names.get('out_nonco_angles_file') is not None:
                self._NONCO_ANGLES_OUT = output_file_names.get('out_nonco_angles_file')
            else:
                self._NONCO_ANGLES_OUT = KkrCalculation()._NONCO_ANGLES_OUT
            if output_file_names.get('output_0_file') is not None:
                self._OUTPUT_0_INIT = output_file_names.get('output_0_file')
            else:
                self._OUTPUT_0_INIT = KkrCalculation()._OUTPUT_0_INIT
            if output_file_names.get('output_000_file') is not None:
                self._OUTPUT_000 = output_file_names.get('output_000_file')
            else:
                self._OUTPUT_000 = KkrCalculation()._OUTPUT_000
            if output_file_names.get('output_2_file') is not None:
                self._OUTPUT_2 = output_file_names.get('output_2_file')
            else:
                self._OUTPUT_2 = KkrCalculation()._OUTPUT_2
            if output_file_names.get('timing_file') is not None:
                self._OUT_TIMING_000 = output_file_names.get('timing_file')
            else:
                self._OUT_TIMING_000 = KkrCalculation()._OUT_TIMING_000


    # These value are set as class attributes in the parent class,
    # BasePwInputGenerator, but they will be different for a job that wasn't
    # run using aiida, and they will likely vary from job to job. Therefore,
    # we override the parent class's attributes using properties, whose
    # setter methods store the values as db attributes, and whose getter
    # methods retrieve the stored values from the db.

    @property
    def _INPUT_FILE_NAMES(self):
        return self.get_attr('input_file_names', None)

    @_INPUT_FILE_NAMES.setter
    def _INPUT_FILE_NAMES(self, value):
        self._set_attr('input_file_names', value)

    @property
    def _OUTPUT_FILE_NAMES(self):
        return self.get_attr('output_file_names', None)

    @_OUTPUT_FILE_NAMES.setter
    def _OUTPUT_FILE_NAMES(self, value):
        self._set_attr('output_file_names', value)
