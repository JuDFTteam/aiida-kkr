#!/usr/bin/env python
# coding: utf-8
"""Test extracted data.
In this current test, the test will try to populate the MultiImpuritiesData
obj with the data from combine_imps_wc nodes in all possible ways.
The final obj will construct the dict maping {'data_headers':[data list]}.
In this test typicaly three pandas data frame have been created

1. Instanciating the MultiImpuritiesData obj with validating data node
    '22d72273-676a-4243-9d41-f20273d7941b'.

2. Instantiating class method ExtractDictMultipleNode(node_list), where
   the list contains the validating node. Later append further nodes using
   the AppendDataMultipleNode method.

3. Instantiating the class method ExtractDictMultipleNode(node_list).
    After populate the obj with the node list containg the validating node.

Later, the data represented in 'pandas data frame' will be validated
w.r.t. the combine_imps_wc node '22d72273-676a-4243-9d41-f20273d7941b'
which is the first node in data_list_1. So, technically this node is
one equipment of the test.
"""

import pytest
import os
import pandas as pd
from aiida.orm import load_node
from aiida.manage.tests.pytest_fixtures import clear_database_after_test, clear_database_before_test
from aiida_kkr.tools.multi_imps_data_extract import MultiImpuritiesData
from ..conftest import import_with_migration


def test_multi_imp_to_data_frame(clear_database_before_test, clear_database_after_test):
    """ Test Pandas DataFrame.

    Data is being tested under the validating node
    '22d72273-676a-4243-9d41-f20273d7941b'
    """

    TEST_DIR = os.path.dirname(os.path.realpath(__file__))
    import_with_migration(os.path.join(TEST_DIR, '../files/db_dump_imps_data_extract_tools.tar.gz'))

    data_list_1 = [
        '22d72273-676a-4243-9d41-f20273d7941b', '5f8bf3a8-fc7d-45d5-8140-9a1e7f566207',
        '7ce9e9dc-3891-4922-88ae-0a63e5893026'
    ]
    data_list_1 = [load_node(i) for i in data_list_1]

    data_list_2 = [
        '967ecb19-75b9-4f72-9ef8-04da8b7376f0', 'f3669dc2-a94d-4eee-962e-2720dcd444cd',
        'cc2de321-52ce-482d-962d-76bbf0ed7d84'
    ]
    data_list_2 = [load_node(i) for i in data_list_2]

    node_1 = data_list_1[0]

    data_aggr_1 = MultiImpuritiesData(node_1)
    data_aggr_1.AppendDataMultipleNode(data_list_1[1:])
    data_frame_1 = pd.DataFrame(data_aggr_1.GetDataDict())

    data_aggr_2 = MultiImpuritiesData.ExtractDictMultipleNode(data_list_1[::-1])
    data_aggr_2.AppendDataMultipleNode(data_list_2)
    data_frame_2 = pd.DataFrame(data_aggr_2.GetDataDict())

    data_aggr_3 = MultiImpuritiesData.ExtractDictMultipleNode(data_list_2)
    data_aggr_3.AppendDataMultipleNode(data_list_1)
    data_frame_3 = pd.DataFrame(data_aggr_3.GetDataDict())

    # data validation from combine_imps_wc node
    # '22d72273-676a-4243-9d41-f20273d7941b'
    calc_imps = 2
    J_data_uuid = '5c5cff9f-f037-4ab2-ad48-f4ce3743c3df'

    imp0 = 21
    offset0 = 0
    ilayer0 = 3

    imp1 = 22
    offset1 = 1
    ilayer1 = 3

    i = 0
    j = 1
    Z_i = 21
    Z_j = 22

    J = 1.24270853e-02
    D = 2.95299862e-03
    Dx = -3.34314299e-04
    Dy = -2.02302285e-03
    Dz = -2.12504431e-03

    rx = -0.41851489357247
    ry = 0.0
    rz = 0.0
    r = abs(rx)

    mom1 = 1.8567015891424e-06
    mom2 = 0.81527563023296
    tot_mom = mom1 + mom2
    for panda_df in [data_frame_1, data_frame_2, data_frame_3]:
        probing_data = panda_df.loc[panda_df['J_data_uuid'] == '5c5cff9f-f037-4ab2-ad48-f4ce3743c3df'].iloc[0]

        assert probing_data['calc_imps'] == calc_imps
        assert probing_data['imp0'] == imp0
        assert probing_data['offset0'] == offset0
        assert probing_data['ilayer0'] == ilayer0
        assert probing_data['imp1'] == imp1
        assert probing_data['offset1'] == offset1
        assert probing_data['ilayer1'] == ilayer1
        assert probing_data['i'] == i
        assert probing_data['j'] == j
        assert probing_data['Z_i'] == Z_i
        assert probing_data['Z_j'] == Z_j
        assert (probing_data['J'] - J) <= 1e-7
        assert (probing_data['D'] - D) <= 1e-7
        assert (probing_data['Dx'] - Dx) <= 1e-7
        assert (probing_data['Dy'] - Dy) <= 1e-7
        assert (probing_data['Dz'] - Dz) <= 1e-7
        assert (probing_data['rx'] - rx) <= 1e-7
        assert (probing_data['ry'] - ry) <= 1e7
        assert (probing_data['rz'] - rz) <= 1e-7
        assert (probing_data['r'] - r) <= 1e-7
        assert (probing_data['mom1'] - mom1) <= 1e-7
        assert (probing_data['mom2'] - mom2) <= 1e-7
        assert (probing_data['tot_mom'] - tot_mom) <= 1e-7
