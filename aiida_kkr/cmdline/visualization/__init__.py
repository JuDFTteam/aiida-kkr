# -*- coding: utf-8 -*-
"""
Module with CLI commands for various visualizations of data types.
"""
import click
from aiida.cmdline.utils import decorators
from aiida.cmdline.params import arguments


@click.command(name='plot',)
@click.option('-f', 'filename', type=click.File('r'), default=None)
@click.option(
    '-o',
    '--option',
    type=str,
    default=None,
    multiple=True,
    help=
    "Option in that is passed as key=value pair to plot_kkr function. Multiple options are allowed. Should be given like this: -o marker=x -o lw=2. If TEXT does not contain '=', the option is ignored."
)
@arguments.NODES('nodes')
#@click.pass_context
def cmd_plot(nodes, option, filename):
    """
    Invoke the plot_kkr command on given nodes
    """

    # TODO add kwargs input or add some options
    # TODO: make iatom etc selectors
    from aiida_kkr.tools.plot_kkr import plot_kkr

    nodes = list(nodes)
    nodesf = []
    if filename is not None:
        nodesf = filename.read().split()
        filename.close()
    nodes = nodes + nodesf

    kwargs = {}
    for item in option:
        if '=' not in item:
            # if we cannot split the string to extract the key, value pair this input is ignored
            print(f"Warning: could not split option '{item}' because it does not contain a '='.")
            pass
        else:
            key, val = item.split('=')
            # check if we need to convert the datatype of some argument
            if key in _mappings_dtype:
                if _mappings_dtype[key] == 'bool':
                    val = _get_bool(val)
                elif _mappings_dtype[key] == 'float':
                    val = float(val)
                elif _mappings_dtype[key] == 'list:int':
                    val = [int(i) for i in val.split()]
            kwargs[key] = val
    print(kwargs)
    plot_kkr(nodes, **kwargs)


# helper functions

_mappings_dtype = {
    'nolegend': 'bool',
    'silent': 'bool',
    'noshow': 'bool',
    'show_empty_atoms': 'bool',
    'filled': 'bool',
    'nofig': 'bool',
    'strucplot': 'bool',
    'logscale': 'bool',
    'newfig': 'bool',
    'interpol': 'bool',
    'all_atoms': 'bool',
    'l_channels': 'bool',
    'sum_spins': 'bool',
    'switch_xy': 'bool',
    'switch_sign_spin2': 'bool',
    'dos_only': 'bool',
    'xscale': 'float',
    'yscale': 'float',
    'xshift': 'float',
    'label': 'str',
    'only': 'str',
    'ptitle': 'str',
    'iatom': 'list:int',
    'subplot': 'list:int',
}


def _get_bool(value):
    """
    convert str to bool, i.e. we check if value is a string and then see if we can recognize True in there, otherwise we set it to False
    """
    # do nothing if the input is a string
    if not type(value) == str:
        return value
    # make string lowercase and search for 't' in there (this way we allow 'True', '.true.', and 'T')
    return 't' in value.lower()
