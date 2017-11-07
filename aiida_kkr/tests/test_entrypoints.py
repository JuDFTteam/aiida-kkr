#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest

@pytest.mark.usefixtures("aiida_env")
class TestAiida_fleur_entrypoints:
    """
    tests all the entry points of the Kkr plugin. Therefore if the plugin is reconized by AiiDA 
    and installed right. 
    """
    
    # Calculation
    def test_kkrimpcalculation_entry_point(aiida_env):
        from aiida.orm import CalculationFactory
        from aiida_kkr.calculations.kkrimp import KkrimpCalculation

        kkrimp_calculation = CalculationFactory('kkr.kkrimp')
        assert kkrimp_calculation == KkrimpCalculation

    
    # Data
    def test_kkrstructuredata_entry_point(aiida_env):
        from aiida.orm import DataFactory, Data
        from aiida_kkr.data.kkrstructure import KkrstructureData
        
        StructureData = DataFactory('structure')
        kkrstruc = DataFactory('kkr.kkrstructure')
        assert kkrstruc == KkrstructureData
        assert isinstance(kkrstruc(), Data)
        assert isinstance(kkrstruc(), StructureData)

    # Parsers

    def test_kkrimp_parser_entry_point(aiida_env):
        from aiida.parsers import ParserFactory
        from aiida_kkr.parsers.kkrimp import KkrimpParser

        parser = ParserFactory('kkr.kkrimpparser')
        assert parser == KkrimpParser

    # Workchains


