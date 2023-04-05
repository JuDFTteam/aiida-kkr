"""
these are calcfunctions used to store nodes which ensures the data provenance
"""

from aiida.engine import calcfunction

__copyright__ = (u'Copyright (c), 2020, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = (u'Philipp Ruessmann')


@calcfunction
def create_out_dict_node(out_node, **input_nodes):
    """
    This calcfunction makes sure the data provenance is kept for the output nodes of workflows.
    Stores the output node and links to inputs nodes if given.

    :param out_node: output Dict node containing results of a workflow
    :type out_node: :class: aiida.orm.node

    :param input_nodes: dict of link_label, input_node pairs: `input_nodes = {'link_label', Node}`
    :type node_dict: dict
    ..note::
        The nodes in input_nodes cannot be stored already

    :return output_node: output Node (stored via calcfunction for data-provenance)
    :type output_node: :class: aiida.orm.node
    ..note::
        The output node is a clone of the input `out_node` and will be doubled in the database
    """

    # need to return a clone of the node to not have a cyclic graph
    return out_node.clone()
