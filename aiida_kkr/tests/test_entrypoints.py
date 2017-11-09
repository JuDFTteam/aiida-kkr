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
    def test_kkrcalculation_entry_point(aiida_env):
        from aiida.orm import CalculationFactory
        from aiida_kkr.calculations.kkr import KkrCalculation

        kkr_calculation = CalculationFactory('kkr.kkr')
        assert kkr_calculation == KkrCalculation

    
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

    def test_kkr_parser_entry_point(aiida_env):
        from aiida.parsers import ParserFactory
        from aiida_kkr.parsers.kkr import KkrParser

        parser = ParserFactory('kkr.kkrparser')
        assert parser == KkrParser

    # Workchains


