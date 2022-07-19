# -*- coding: utf-8 -*-
"""
This contains code snippets and utility useful for dealing with parameter data nodes
commonly used by the plugin and workflows
"""


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
