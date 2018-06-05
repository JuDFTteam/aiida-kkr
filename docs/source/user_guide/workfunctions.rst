=============
Workfunctions
=============

Here the workfunctions provided by the aiida-kkr plugin are presented. The workfunctions are 
small tools useful for small tasks performed on aiida nodes that keep the provenance in the 
database. 


update_params_wf
++++++++++++++++

The workfunktion ``aiida_kkr.tools.common_workfunctions.update_params_wf`` takes as an input a
*ParameterData* node (``parameternode``) containing a KKR parameter set (i.e. created using the ``kkrparams`` class)
and updates the parameter node with new values given in the dictionary of the second 
*ParameterData* input node (``updatenode``).

Input nodes:
    * ``parameternode`` (ParameterData): aiida node of a KKR parameter set
    * ``updatenode`` (ParameterData): aiida node containing parameter names with new values

Output node:
    * ``updated_parameter_node`` (ParameterData): new parameter node with updated values
    
.. note:: If the ``updatenode`` contains the keys ``nodename`` and/or ``nodedesc`` then the 
    label and/or description of the output node will be set accordingly.

Example Usage::

    # initial KKR parameter node
    input_node = ParameterData(dict=kkrparams(LMAX=3, EMIN=0))
    input_node.store()
    # update some values (e.g. change EMIN)
    updated_params = ParameterData(dict={'nodename': 'my_changed_name', 'nodedesc': 'My description text', 'EMIN': -1, 'RMAX': 10.})
    new_params_node = update_params_wf(input_node, updated_params)

    
neworder_potential_wf
+++++++++++++++++++++

The workfunction ``aiida_kkr.tools.common_workfunctions.neworder_potential_wf`` creates a 
*SingleFileData* node that contains the new potential based in a potential file in the 
*RemoteFolder* input node (``settings_node``) which is braught to a new order according to 
the workfunction settings in the *ParameterData* input node (``parent_calc_folder``).

Input nodes:
    * ``settings_node`` (ParameterData): Settings like filenames and neworder-list
    * ``parent_calc_folder`` (RemoteData): folder where initial potential file is found
    * ``parent_calc_folder2`` (RemoteData, optional): folder where second potential is found

Output node:
    * ``potential_file`` (SingleFileData): output potential in new order
    
.. note:: 
    The settings_dict should contain the following keys:
        * ``pot1``, mandatory: *<filename_input_potential>*
        * ``out_pot``, mandatory: *<filename_output_potential>*
        * ``neworder``, mandatory: *[list of intended order in output potential]* 
        * ``pot2``, mandatory if ``parent_calc_folder2`` is given as input node: *<filename_second_input_file>*
        * ``replace_newpos``, mandatory if ``parent_calc_folder2`` is given as input node: *[[position in neworder list which is replace with potential from pot2, position in pot2 that is chosen for replacement]]*
        * ``label``, optional: *label_for_output_node*
        * ``description``, optional: *longer_description_for_output_node*


prepare_VCA_structure_wf
++++++++++++++++++++++++

.. warning:: Not implemented yet!



prepare_2Dcalc_wf
+++++++++++++++++

.. warning:: Not implemented yet!


