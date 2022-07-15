# -*- coding: utf-8 -*-
"""
Helper functions used in KKR and voronoi parser. Needed to properly open file with the context manager.
"""


def open_context_to_stack(stack, out_folder, outfile_name, rw_mode=u'r'):
    """Open a file in a context manager which are collected
    in a context manager stack. The file is only opened if
    the filename is not None.

    :param stack:  a contextlib.ExitStack to which the opened file context is added
    :param out_folder: a retrieved folder for the calculation where the files are found
    :param outfile_name: name of the file that is supposed to be opened from the out_folder

    :returns: the file handle or None if the filename was None
    """
    if outfile_name is not None:
        outfile = stack.enter_context(out_folder.open(outfile_name, rw_mode))
    else:
        outfile = None
    return outfile


def open_files_in_context(stack, out_folder, *filenames):
    """Open files and add to context manager stack.

    :param stack:  a contextlib.ExitStack to which the opened file context is added
    :param out_folder: a retrieved folder for the calculation where the files are found
    :param filnames: list of filenames (should be strings) that are openend in out_folder

    :returns: a list of file handles
    """
    fhandles = []
    for fname in filenames:
        fhandles.append(open_context_to_stack(stack, out_folder, fname))

    return fhandles
