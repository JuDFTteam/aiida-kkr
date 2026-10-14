#!/usr/bin/env python
"""kkr_imp_wc takes ALATBASIS from the host GF write-out when this daemon's Ang->Bohr constant differs.

KKRimp stops if a potential's ALAT differs from the GF write-out's by more than 1e-11. masci-tools
picks ANG_BOHR_KKR per process (MASCI_TOOLS_USE_OLD_CONSTANTS: old / interim / new, about 4e-10
apart), so an auxiliary Voronoi run that recomputes alat from the structure under another constant
than the stored GF write-out produced an impurity potential KKRimp refused.
"""

import re
from types import SimpleNamespace

import pytest

from aiida_kkr.workflows.kkr_imp import get_alat_override

CONST = {'old': 1.8897261254578281, 'interim': 1.8897261249935897, 'new': 1.8897261246257702}
ALAT_GF_OLD = 4.82381975000676  # ALATBASIS of a GF write-out run under the 'old' constant (fcc Cu)


@pytest.fixture
def cu_fcc(aiida_profile):
    from aiida.orm import StructureData
    structure = StructureData(cell=[[1.805, 1.805, 0.0], [1.805, 0.0, 1.805], [0.0, 1.805, 1.805]])
    structure.append_atom(position=[0, 0, 0], symbols='Cu')
    return structure


def _use_constant(monkeypatch, name):
    import masci_tools.util.constants
    monkeypatch.setattr(masci_tools.util.constants, 'ANG_BOHR_KKR', CONST[name])


def _gf_calc(inputcard):
    """stand-in for a KkrCalculation node; inputcard=None means the file is not in the repository"""

    def get_object_content(name):
        if inputcard is None:
            raise FileNotFoundError(name)
        return inputcard

    return SimpleNamespace(base=SimpleNamespace(repository=SimpleNamespace(get_object_content=get_object_content)))


def _card_values(path):
    txt = path.read_text()
    values = {key: re.search(rf'^{key}=\s*(\S+)', txt, re.MULTILINE).group(1) for key in ('ALATBASIS', 'RMAX', 'GMAX')}
    values['BRAVAIS'] = txt.split('BRAVAIS\n')[1].split('\n')[0].split()
    return values


@pytest.mark.parametrize('constant', ['interim', 'new'])
def test_get_alat_override_when_constants_differ(monkeypatch, cu_fcc, constant):
    _use_constant(monkeypatch, constant)
    assert get_alat_override(_gf_calc(f'ALATBASIS=      {ALAT_GF_OLD}\n'), cu_fcc) == ALAT_GF_OLD


def test_get_alat_override_none_when_constants_agree(monkeypatch, cu_fcc):
    """the write-out's 12-decimal rounding alone must not trigger an override (keeps inputs and caching unchanged)"""
    _use_constant(monkeypatch, 'old')
    assert get_alat_override(_gf_calc('ALATBASIS=        4.823819750007\n'), cu_fcc) is None


def test_get_alat_override_none_without_inputcard(cu_fcc):
    assert get_alat_override(_gf_calc(None), cu_fcc) is None


@pytest.mark.parametrize('constant', ['old', 'interim', 'new'])
def test_alat_override_writes_alatbasis_only(monkeypatch, cu_fcc, tmp_path, constant):
    """ALATBASIS takes the override; BRAVAIS, RMAX and GMAX are the same as without it under any constant"""
    from aiida.orm import Dict
    from aiida_kkr.tools.common_workfunctions import generate_inputcard_from_structure

    _use_constant(monkeypatch, constant)
    params = {'LMAX': 2, 'NSPIN': 1, 'RMAX': 7.0, 'GMAX': 65.0, 'RCLUSTZ': 2.3}
    generate_inputcard_from_structure(Dict(params), cu_fcc, str(tmp_path / 'plain'), isvoronoi=True)
    generate_inputcard_from_structure(
        Dict({
            **params, 'alat_override': ALAT_GF_OLD
        }), cu_fcc, str(tmp_path / 'override'), isvoronoi=True
    )
    plain, override = _card_values(tmp_path / 'plain'), _card_values(tmp_path / 'override')

    assert override['ALATBASIS'] == '4.823819750007'
    assert {
        k: v for k, v in override.items() if k != 'ALATBASIS'
    } == {
        k: v for k, v in plain.items() if k != 'ALATBASIS'
    }
    assert override['BRAVAIS'] == ['0.707106781187', '0.707106781187', '0.000000000000']


def test_update_params_keeps_alat_override(aiida_profile):
    """kkr_startpot_wc rebuilds its parameters with update_params, which drops keys kkrparams does not know"""
    from aiida.orm import Dict
    from aiida_kkr.tools.common_workfunctions import update_params

    node = update_params(Dict({'LMAX': 2, 'alat_override': ALAT_GF_OLD}), RMAX=8.0)
    assert node['alat_override'] == ALAT_GF_OLD
    assert update_params(Dict({'LMAX': 2}), alat_override=ALAT_GF_OLD)['alat_override'] == ALAT_GF_OLD
    assert 'alat_override' not in update_params(Dict({'LMAX': 2}), RMAX=8.0).get_dict()
