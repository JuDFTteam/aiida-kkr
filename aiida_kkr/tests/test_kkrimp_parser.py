#!/usr/bin/env python

import pytest

# some global settings

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_parser():
    """
    Tests for the kkrimp calculation
    """
    
    def test_parse_kkrimp_calc(self):
        """
        simple Cu noSOC, FP, lmax2
        """
        from aiida.orm import load_node
        from aiida_kkr.parser.kkrimp import KkrimpParser
        from aiida_kkr.parsers.kkrimp import KkrimpParser
        kkrimp_calc = load_node(179)
        parser = KkrimpParser(kkrimp_calc)
        success, outnodes = parser.parse_from_calc()
        assert success


if __name__=='__main__':
    t = Test_kkrimp_parser()
    t.test_parse_kkrimp_calc()
