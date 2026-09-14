#!/usr/bin/env python

from builtins import object
import pytest
from aiida import orm
from ..conftest import import_with_migration, test_dir

# tests


def test_parse_kkrimp_calc(aiida_profile):
    """
    simple Cu noSOC, FP, lmax2
    """
    from aiida_kkr.parsers.kkrimp import KkrimpParser
    group_pk = import_with_migration('data_dir/kkrimp_full_wc.aiida')
    kkrimp_calc = [
        i for i in orm.load_group(group_pk).nodes if i.label == 'KKRimp calculation step 2 (IMIX=0, Zimp: [30.0])'
    ][0]
    print(kkrimp_calc, kkrimp_calc.label)

    print(kkrimp_calc.outputs.retrieved.list_object_names())
    parser = KkrimpParser(kkrimp_calc)
    out = parser.parse(debug=False)
    print(out)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    assert out_dict['parser_errors'] == []


def test_parse_kkrimp_calc_complex(aiida_profile):
    """
    complex magnetic impurity with SOC
    """
    from aiida_kkr.parsers.kkrimp import KkrimpParser
    import_with_migration('files/export_kkrimp_calc.aiida')
    kkrimp_calc = orm.load_node('b9a2b29a-e250-4992-ae6a-579b733ad1f8')
    parser = KkrimpParser(kkrimp_calc)
    out = parser.parse(debug=False, doscalc=False)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    assert out_dict['parser_errors'] == []
    assert 'nonfinite_values' not in out_dict


def test_parse_kkrimp_nonfinite(aiida_localhost, tmp_path):
    """
    diverged KKRimp calculation whose output contains NaN values

    Retrieved files of a real run (Ba impurity in Cu, calc uuid 90a0e26b-71c3-430b-a160-473f5951cac8,
    created 2026-09-12, aiida-kkr fc98f34, masci-tools 487413bb whose KKR parser files are identical to
    3392746f). The run blows up in its 28th and last iteration. The calculation was killed before
    retrieval, so the six output files are shipped whole as a tarball: cutting out_log.000.txt removes
    the diverged iteration. masci-tools parses these files with success and no errors, so the non-finite
    values have to be caught independently of the parser's verdict.
    """
    import tarfile
    from aiida.common.links import LinkType
    from aiida_kkr.parsers.kkrimp import KkrimpParser

    with tarfile.open(test_dir / 'files/kkrimp_parser/nonfinite_ba_cu_rung0.tar.xz') as tar:
        tar.extractall(tmp_path)

    kkrimp_calc = orm.CalcJobNode(computer=aiida_localhost, process_type='aiida.calculations:kkr.kkrimp')
    kkrimp_calc.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
    kkrimp_calc.store()
    retrieved = orm.FolderData(tree=tmp_path)
    retrieved.base.links.add_incoming(kkrimp_calc, link_type=LinkType.CREATE, link_label='retrieved')
    retrieved.store()

    parser = KkrimpParser(kkrimp_calc)
    exit_code = parser.parse(debug=False, doscalc=False)
    assert exit_code.status == 303

    # storing is where the unsanitized output raised 'nan and inf/-inf can not be serialized to the database'
    output_parameters = parser.outputs.output_parameters.store()
    out_dict = output_parameters.get_dict()

    assert out_dict['parser_version'] == KkrimpParser(kkrimp_calc)._ParserVersion
    assert out_dict['parser_errors'] == []
    assert len(out_dict['nonfinite_values']) == 54
    assert 'convergence_group.rms' in out_dict['nonfinite_values']
    assert 'energy' in out_dict['nonfinite_values']
    assert 'convergence_group.total_spin_moment_all_iterations[1][513][2]' in out_dict['nonfinite_values']
    convergence = out_dict['convergence_group']
    assert convergence['first_nonfinite_iteration_index'] == 27
    assert len(convergence['rms_all_iterations']) == 28
    assert convergence['rms_all_iterations'][27] is None
    assert all(isinstance(rms, float) for rms in convergence['rms_all_iterations'][:27])
    assert convergence['rms'] is None and out_dict['energy'] is None
    assert any('nonfinite_values' in warning for warning in out_dict['parser_warnings'])
