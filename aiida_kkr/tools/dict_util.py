# -*- coding: utf-8 -*-
"""
This contains code snippets and utility useful for dealing with parameter data nodes
commonly used by the plugin and workflows
"""

import math
import numbers
import re

import numpy as np


def clean_nones(dict_to_clean):
    """Recursively remove all keys which values are None from a nested dictionary
    return the cleaned dictionary

    :param dict_to_clean: (dict): python dictionary to remove keys with None as value
    :return: dict, cleaned dictionary
    """
    new_dict = {}
    for key, val in dict_to_clean.items():
        if isinstance(val, dict):
            new_val = clean_nones(val)
        else:
            new_val = val
        if new_val is not None:  # currently we keep empty dicts
            new_dict[key] = new_val

    return new_dict


def sanitize_nonfinite(out_dict):
    """Replace every non-finite float (NaN, inf, -inf) in a nested parser output by None, in place.

    AiiDA refuses to store non-finite floats in a Dict node, so a parser output that contains them
    makes the calculation except at store time and lose its exit code. None is used instead of a
    number on purpose: e.g. 0.0 would read as a perfectly converged rms or a real total energy.
    Lists keep their full length, so the position of a replaced value stays recoverable.

    Dicts and lists are modified in place, tuples and numpy arrays are rebuilt (arrays as lists).

    :param out_dict: (dict) parser output dictionary, nested to any depth
    :return: list of str, paths of the replaced values, e.g. 'convergence_group.rms_all_iterations[27]'
    """
    paths = []
    _sanitize_value(out_dict, '', paths)
    return paths


def _sanitize_value(val, path, paths):
    """Recursive worker of sanitize_nonfinite, returns the (possibly replaced) value."""
    if isinstance(val, dict):
        for key in val:
            val[key] = _sanitize_value(val[key], f'{path}.{key}' if path else str(key), paths)
        return val
    if isinstance(val, np.ndarray):
        val = val.tolist()
    if isinstance(val, (list, tuple)):
        new = [_sanitize_value(item, f'{path}[{i}]', paths) for i, item in enumerate(val)]
        if isinstance(val, tuple):
            return tuple(new)
        val[:] = new
        return val
    if isinstance(val, numbers.Real) and not isinstance(val, numbers.Integral) and not math.isfinite(val):
        paths.append(path)
        return None
    return val


def record_nonfinite(out_dict, paths):
    """Record in the parser output which values sanitize_nonfinite replaced.

    Does nothing if `paths` is empty. Otherwise adds
      * `nonfinite_values`: the list of replaced paths,
      * `convergence_group.first_nonfinite_iteration_index`: 0-based index of the first replaced
        entry of `convergence_group.rms_all_iterations` (only if that list has one),
      * a warning in `parser_warnings`.

    :param out_dict: (dict) parser output dictionary, already sanitized
    :param paths: list of str, return value of sanitize_nonfinite
    """
    if not paths:
        return
    out_dict['nonfinite_values'] = list(paths)
    rms_indices = [
        int(match.group(1))
        for match in (re.fullmatch(r'convergence_group\.rms_all_iterations\[(\d+)\]', p) for p in paths)
        if match
    ]
    if rms_indices:
        out_dict['convergence_group']['first_nonfinite_iteration_index'] = min(rms_indices)
    out_dict.setdefault('parser_warnings', []).append(
        f'Warning! {len(paths)} non-finite values (NaN/inf) in parser output replaced by None, '
        "see 'nonfinite_values'."
    )
