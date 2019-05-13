#!/usr/bin/env python

from __future__ import print_function
from masci_tools.io.common_functions import search_string
with open('out_rmq_check') as f:
    txt = f.readlines()
    i = search_string('Sent', txt)
    t_sent = txt.pop(i)
    t_sent = t_sent.split('Hello World! ')[1]
    j = search_string(t_sent, txt)
    print('rmq check passed?', t_sent in txt[j])
    assert t_sent in txt[j]

