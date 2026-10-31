#!/usr/bin/env python

if __name__ != '__main__':
    import pytest
    from ..conftest import kkrimp_local_code, kkrhost_local_code, test_dir, data_dir
import numpy as np
from aiida.orm import load_node, Dict, load_group
from aiida.engine import run_get_node
from aiida_kkr.workflows import combine_imps_wc, kkr_scf_wc
from ..conftest import import_with_migration

# activate debug mode?
_debug = True


def write_graph(node, label=''):
    #if create_graph_file:
    from aiida.tools.visualization import Graph
    builder = node.get_builder_restart()
    name = str(builder.process_class).split('.')[-1].strip("'>")
    name += label
    graph = Graph(engine='dot', node_id_type='uuid')
    graph.recurse_ancestors(node, depth=None, link_types=(), annotate_links='both', include_process_outputs=True)
    graph.recurse_descendants(node, depth=None, link_types=(), annotate_links='both', include_process_inputs=True)
    output_file_name = graph.graphviz.render(str(data_dir) + '/' + name, format='pdf', view=False, cleanup=True)
    print('wrote graph to', output_file_name)


def get_single_imp_inputs():
    # import single imp calculations
    group_pk = import_with_migration(test_dir / 'data_dir/kkrimp_sub_wc.aiida')
    if _debug:
        print(group_pk, [i.label for i in load_group(group_pk).nodes])
    for node in load_group(group_pk).nodes:
        if 'KKRimp calculation step 4' in node.label:
            imp1 = node
    if _debug:
        print(imp1, list(imp1.outputs), list(imp1.inputs))
    imp1_out = imp1.outputs.output_parameters
    imp2_out = imp1_out  # use the same impurity and create a dimer

    return imp1_out, imp2_out


def get_builder_basic(label, kkrhost_local_code, kkrimp_local_code):

    imp1_out, imp2_out = get_single_imp_inputs()

    #set up combine_imps_wc workflow
    builder = combine_imps_wc.get_builder()
    builder.metadata.label = label
    builder.impurity1_output_node = imp1_out
    builder.impurity2_output_node = imp2_out
    builder.offset_imp2 = Dict({'index': 1})
    # set code
    builder.host_gf.kkr = kkrhost_local_code  # should not be required if gf_host_remote is given, seems to be a problem of aiida-testing's run_with_cache
    builder.scf.kkrimp = kkrimp_local_code
    # set computer options
    options = {
        'queue_name': '',
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    builder.scf.options = Dict(options)

    scf_settings = kkr_scf_wc.get_wf_defaults(silent=True)[0]
    scf_settings.kkr_runmax = 1
    scf_settings.nsteps = 2
    scf_settings.do_final_cleanup = False  # this is needed to allow for caching
    builder.scf.wf_parameters = Dict(scf_settings)
    builder.scf.params_overwrite = Dict({'TOL_ALAT_CHECK': 1e-8})

    builder.host_gf.wf_parameters = Dict({'retrieve_kkrflex': True})  # this is needed to allow for caching

    builder.host_gf.options = builder.scf.options
    builder.wf_parameters_overwrite = Dict({'global': {'allow_unconverged_inputs': True}})

    return builder


def test_combine_imps(
    aiida_profile_clean,
    kkrhost_local_code,
    kkrimp_local_code,
    enable_archive_cache,
    ndarrays_regression,
    nopytest=False
):
    """
    test for combine_imps_wc (place two imps next to each other)
    """

    builder = get_builder_basic('test_combine_imps', kkrhost_local_code, kkrimp_local_code)

    # now submit
    print(builder, type(builder))
    if not nopytest:
        with enable_archive_cache(data_dir / 'combine_imps.aiida'):
            out, node = run_get_node(builder)
    else:
        out, node = run_get_node(builder)

    print((out, node))
    write_graph(node)

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert 'remote_data_gf' in out  # make sure GF writeout step was done
    check_dict = {'rms_per_atom': node.outputs.last_calc_output_parameters['convergence_group']['rms_per_atom']}
    ndarrays_regression.check(check_dict)


def test_combine_imps_params_kkr_overwrite(
    aiida_profile_clean,
    kkrhost_local_code,
    kkrimp_local_code,
    enable_archive_cache,
    ndarrays_regression,
    nopytest=False
):
    """
    test for combine_imps_wc overwriting the k-mesh with hte params_kkr_overwrite input to the gf writeout step
    """

    builder = get_builder_basic('test_combine_imps_params_kkr_overwrite', kkrhost_local_code, kkrimp_local_code)
    # increase k-mesh for GF writeout step with params_kkr_overwrite
    builder.host_gf.params_kkr_overwrite = Dict({'BZDIVIDE': [20, 20, 20]})

    # now submit
    print(builder, type(builder))
    if not nopytest:
        with enable_archive_cache(data_dir / 'combine_imps_params_kkr_overwrite.aiida'):
            out, node = run_get_node(builder)
    else:
        out, node = run_get_node(builder)
    print((out, node))
    write_graph(node, '_params_kkr_overwrite')

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert 'remote_data_gf' in out  # make sure GF writeout step was done
    check_dict = {'rms_per_atom': node.outputs.last_calc_output_parameters['convergence_group']['rms_per_atom']}
    ndarrays_regression.check(check_dict)


def test_combine_imps_reuse_gf(
    aiida_profile_clean,
    kkrhost_local_code,
    kkrimp_local_code,
    enable_archive_cache,
    ndarrays_regression,
    nopytest=False
):
    """
    test for combine_imps_wc reusing the host gf from a previous calculation
    """

    # import previous combine_imps workflow and reuse the host GF
    group_pk = import_with_migration(test_dir / 'data_dir/combine_imps.aiida')
    for node in load_group(group_pk).nodes:
        if node.label == 'test_combine_imps':
            imp_combine_imported = node
    host_gf_remote = imp_combine_imported.outputs.remote_data_gf

    builder = get_builder_basic('test_combine_imps_reuse_gf', kkrhost_local_code, kkrimp_local_code)
    # set pre-calculated host_gf
    builder.gf_host_remote = host_gf_remote

    # now submit
    print(builder, type(builder))
    if not nopytest:
        with enable_archive_cache(data_dir / 'combine_imps_reuse_gf.aiida'):
            out, node = run_get_node(builder)
    else:
        out, node = run_get_node(builder)
    print((out, node))
    write_graph(node, '_reuse_gf')

    # check outcome
    results = out['workflow_info'].get_dict()
    print(results)
    assert 'remote_data_gf' not in out  # make sure GF writeout step was skipped
    check_dict = {'rms_per_atom': node.outputs.last_calc_output_parameters['convergence_group']['rms_per_atom']}
    ndarrays_regression.check(check_dict)


def _fake_kkr_imp_wc(settings_LDAU=None):
    """stand-in for a kkr_imp_wc node: kkr_imp_wc -> kkr_imp_sub_wc -> KkrimpCalculation"""
    from types import SimpleNamespace
    from aiida.common import AttributeDict
    from aiida_kkr.workflows import kkr_imp_wc, kkr_imp_sub_wc

    def outgoing(node):
        get_outgoing = lambda node_class: SimpleNamespace(first=lambda: SimpleNamespace(node=node))
        return SimpleNamespace(links=SimpleNamespace(get_outgoing=get_outgoing))

    calc = SimpleNamespace(outputs=SimpleNamespace(retrieved='retrieved_of_kkrimp_calc'))
    inputs = AttributeDict({'settings_LDAU': settings_LDAU} if settings_LDAU is not None else {})
    sub = SimpleNamespace(process_class=kkr_imp_sub_wc, inputs=inputs, base=outgoing(calc))
    return SimpleNamespace(process_class=kkr_imp_wc, base=outgoing(sub))


@pytest.mark.parametrize('imp1_ldau', [None, 'settings_LDAU_of_imp1'])
def test_get_ldau_combined_kkr_imp_wc_inputs(monkeypatch, imp1_ldau):
    """
    get_ldau_combined must handle kkr_imp_wc impurity inputs (crashed with UnboundLocalError before)
    """
    from types import SimpleNamespace
    import aiida_kkr.workflows._combine_imps as combine_imps_module

    calls = []
    combined = SimpleNamespace(uuid='uuid_combined', get_dict=dict)
    monkeypatch.setattr(combine_imps_module, 'combine_settings_ldau', lambda **kwargs: calls.append(kwargs) or combined)
    fake_self = SimpleNamespace(
        ctx=SimpleNamespace(imp1=_fake_kkr_imp_wc(imp1_ldau), imp2=_fake_kkr_imp_wc(), kickout_info='kickout_info'),
        report=print,
    )

    add_ldau, settings = combine_imps_wc.get_ldau_combined(fake_self)

    if imp1_ldau is None:
        assert (add_ldau, settings, calls) == (False, {}, [])
    else:
        assert (add_ldau, settings) == (True, combined)
        assert calls == [{
            'settings_LDAU1': imp1_ldau,
            'retrieved1': 'retrieved_of_kkrimp_calc',
            'kickout_info': 'kickout_info'
        }]


def _fake_combine_imps(**inputs):
    """stand-in for a combine_imps_wc instance that has run `start` and `update_params` with these inputs"""
    from types import SimpleNamespace
    from plumpy.utils import AttributesFrozendict
    from aiida.common import AttributeDict

    inputs.setdefault('scf', {})
    inputs = AttributesFrozendict({
        k: AttributesFrozendict(v) if k in ('scf', 'host_gf') else v for k, v in inputs.items()
    })
    fake_self = SimpleNamespace(
        inputs=inputs, ctx=AttributeDict(), report=print, _wf_default=combine_imps_wc._wf_default
    )
    if 'wf_parameters_overwrite' in inputs:
        fake_self.ctx.wf_parameters_overwrite = inputs.wf_parameters_overwrite.get_dict()
    fake_self.ctx.run_options = dict(combine_imps_wc._wf_default)
    combine_imps_wc.update_params(fake_self)
    return fake_self


def test_update_params_nested_overwrite_reaches_sub_workflows(aiida_profile, monkeypatch):
    """nested wf_parameters_overwrite reaches the wf_parameters given to kkr_imp_sub_wc and kkr_flex_wc (#212)"""
    from types import SimpleNamespace
    import aiida_kkr.workflows._combine_imps as combine_imps_module

    fake_self = _fake_combine_imps(
        scf={
            'wf_parameters': Dict({
                'nsteps': 2,
                'strmix': 0.1
            }),
            'kkrimp': 'kkrimp_code'
        },
        host_gf={'wf_parameters': Dict({'retrieve_kkrflex': True})},
        gf_host_remote='gf_remote',
        wf_parameters_overwrite=Dict({
            'global': {
                'jij_run': True
            },
            'scf': {
                'nsteps': 50
            },
            'host_gf': {
                'retrieve_kkrflex': False
            },
        }),
    )
    assert fake_self.ctx.scf_wf_parameters == {'nsteps': 50, 'strmix': 0.1}
    assert fake_self.ctx.wf_parameters_flex == {'retrieve_kkrflex': False}
    assert fake_self.ctx.run_options['jij_run'] is True
    assert combine_imps_wc._wf_default['jij_run'] is False  # class defaults are not changed

    # run_kkrimp_scf passes the merged parameters on
    builder = SimpleNamespace(metadata=SimpleNamespace())
    monkeypatch.setattr(combine_imps_module.kkr_imp_sub_wc, 'get_builder', lambda: builder)
    fake_self.ctx.update(
        imp1=SimpleNamespace(label='imp1'),
        imp2=SimpleNamespace(label='imp2'),
        imp_info_combined='imp_info',
        combined_potentials='pot'
    )
    fake_self.get_ldau_combined = lambda: (False, {})
    fake_self.submit = lambda b: SimpleNamespace(uuid='uuid')
    combine_imps_wc.run_kkrimp_scf(fake_self)
    assert builder.wf_parameters.get_dict() == {'nsteps': 50, 'strmix': 0.1}


def test_update_params_flat_overwrite(aiida_profile):
    """flat keys keep working: scf keys go to the scf parameters, run options to the run options"""
    fake_self = _fake_combine_imps(
        scf={'wf_parameters': Dict({'nsteps': 2})},
        host_gf={'wf_parameters': Dict({'retrieve_kkrflex': True})},
        wf_parameters_overwrite=Dict({
            'nsteps': 50,
            'jij_run': True,
            'retrieve_kkrflex': False
        }),
    )
    assert fake_self.ctx.scf_wf_parameters == {'nsteps': 50}
    assert fake_self.ctx.wf_parameters_flex == {'retrieve_kkrflex': False}
    assert fake_self.ctx.run_options['jij_run'] is True


def test_update_params_optional_inputs_omitted(aiida_profile):
    """neither scf.wf_parameters nor the host_gf namespace are needed (#212)"""
    from aiida_kkr.workflows._combine_imps import _wf_parameters_node

    fake_self = _fake_combine_imps(gf_host_remote='gf_remote')
    assert fake_self.ctx.scf_wf_parameters == {}
    assert fake_self.ctx.wf_parameters_flex == {}
    # nothing to pass on: the sub-workflows use their defaults
    assert _wf_parameters_node(fake_self.inputs.scf, fake_self.ctx.scf_wf_parameters) is None


def test_wf_parameters_node_reuses_unchanged_input(aiida_profile):
    """unchanged parameters pass the input node itself on (keeps provenance and caching), changed ones a new Dict"""
    from aiida_kkr.workflows._combine_imps import _wf_parameters_node

    given = Dict({'nsteps': 2})
    assert _wf_parameters_node({'wf_parameters': given}, {'nsteps': 2}) is given
    new = _wf_parameters_node({'wf_parameters': given}, {'nsteps': 3})
    assert new is not given and new.get_dict() == {'nsteps': 3}


def test_import_tools_combine_imps_alone():
    """aiida_kkr.tools.combine_imps can be imported in a fresh interpreter (no circular import with workflows)"""
    import subprocess
    import sys
    subprocess.run([sys.executable, '-c', 'from aiida_kkr.tools.combine_imps import get_host_structure'], check=True)


def _cluster_with_ties():
    """single-impurity cluster (alat units) with two pairs of equidistant atoms"""
    pos = np.array([[0., 0., 0.], [0.5, 0., 0.], [-0.5, 0., 0.], [0., 0.6, 0.], [0., -0.6, 0.], [0.3, 0.3, 0.6]])
    dist = np.linalg.norm(pos, axis=1)
    return np.column_stack([pos, np.zeros(len(pos)), np.full(len(pos), 29.), dist])


def _gf_retrieved_with_scoef(scoef, tmp_path):
    from aiida.orm import FolderData
    from aiida_kkr.tools.tools_kkrimp import write_scoef
    write_scoef(scoef, tmp_path / 'scoef')
    retrieved = FolderData()
    retrieved.base.repository.put_object_from_file(tmp_path / 'scoef', 'scoef')
    return retrieved


def test_order_like_gf_scoef(aiida_profile, tmp_path):
    """the rebuilt cluster is put in the atom order of the host GF's scoef, also if that is in other alat units"""
    from aiida_kkr.tools.imp_cluster_tools import order_like_gf_scoef

    clust = _cluster_with_ties()
    # scoef of the GF writeout: mirror pairs swapped and in units of a 0.41851 times smaller lattice constant
    order_gf = [0, 2, 1, 4, 3, 5]
    scoef = clust[order_gf].copy()
    scoef[:, :3] /= 0.41851
    scoef[:, -1] /= 0.41851
    reordered = order_like_gf_scoef(clust, _gf_retrieved_with_scoef(scoef, tmp_path))
    assert np.array_equal(reordered, clust[order_gf])

    # a different cluster cannot be mapped
    with pytest.raises(ValueError):
        order_like_gf_scoef(clust[:-1], _gf_retrieved_with_scoef(scoef, tmp_path))
    shifted = scoef.copy()
    shifted[3, 0] += 0.1
    with pytest.raises(ValueError):
        order_like_gf_scoef(clust, _gf_retrieved_with_scoef(shifted, tmp_path))


def _jij_retrieved(natom):
    """out_Jijmatrix of a Jij run with natom sites, Jraw_kl(i, j) = -(i + j) for all components"""
    from io import StringIO
    from aiida.orm import FolderData
    jij = np.zeros((3, natom, natom, 3, 3))
    for i in range(natom):
        for j in range(natom):
            jij[:, i, j] = -(i + j)
    retrieved = FolderData()
    txt = StringIO()
    np.savetxt(txt, jij.reshape(-1, 3))
    retrieved.base.repository.put_object_from_bytes(txt.getvalue().encode(), 'out_Jijmatrix')
    return retrieved


@pytest.mark.parametrize(
    'zimp, rimp_rel, expected_pairs',
    [
        # three impurities of the same species (Mn-Mn + Mn), all three must be in JijData
        ([25, 25, 25], [[0, 0, 0], [0.5, 0, 0], [0, 0.6, 0]], [(0, 1), (0, 3), (1, 3)]),
        # two impurities in a host that contains the impurity species: no host sites in JijData
        ([29, 29], [[0, 0, 0], [0, 0.6, 0]], [(0, 3)]),
    ],
)
def test_parse_jij_finds_impurities_by_position(aiida_profile, zimp, rimp_rel, expected_pairs):
    from aiida_kkr.workflows._combine_imps import parse_Jij

    clust = _cluster_with_ties()
    for r in rimp_rel:
        clust[np.argmin(np.linalg.norm(clust[:, :3] - r, axis=1)), 4] = zimp[0]
    impurity_info = Dict({'imp_cls': clust, 'Zimp': zimp, 'Rimp_rel': rimp_rel})
    # call the plain function, the calcfunction would need a running broker
    out = parse_Jij.__wrapped__(_jij_retrieved(len(clust)), impurity_info)
    jijdata = out['Jijdata'].get_array('JijData')
    assert len(jijdata) == len(expected_pairs)
    for row, (i, j) in zip(jijdata, expected_pairs):
        assert np.allclose(row[:3], clust[j, :3] - clust[i, :3])
    assert 'between impurities' in out['info']['text']


def test_r_offset_without_index(aiida_profile, monkeypatch):
    """offset_imp2 with r_offset only: no KeyError and no false 'neighbor too small' error"""
    from types import SimpleNamespace
    from aiida.common import AttributeDict
    from aiida.orm import Bool
    import aiida_kkr.workflows._combine_imps as combine_imps_module

    calls = {}

    def fake_create_combined_imp_info_cf(*args, **kwargs):
        calls['args'], calls['kwargs'] = args, kwargs
        return {'imp_info_combined': 'imp_info', 'kickout_info': 'kickout'}

    monkeypatch.setattr(combine_imps_module, 'create_combined_imp_info_cf', fake_create_combined_imp_info_cf)
    impinfo = Dict({'Zimp': 25, 'ilayer_center': 0, 'Rcut': 1.0})
    fake_self = SimpleNamespace(
        ctx=AttributeDict({
            'imp1': SimpleNamespace(inputs={}, name='imp1'),
            'imp2': SimpleNamespace(inputs={}, name='imp2'),
            'single_vs_single': True,
            'host_structure': 'host'
        }),
        inputs=AttributeDict({'offset_imp2': Dict({'r_offset': [0.5, 0., 0.]})}),
        report=print,
        exit_codes=combine_imps_wc.exit_codes,
        get_impinfo_from_hostGF=lambda imp: impinfo,
        extract_imps_info_exact_cluster=lambda: (0, {
            'Zimps': [25, 25],
            'ilayers': [0, 0],
            'offset_imps': [0, None]
        }),
        get_and_check_zimp_list=lambda impinfo: ([25], True),
        get_host_gf_retrieved=lambda imp: 'gf_retrieved_' + imp.name,
    )
    assert combine_imps_wc.create_big_cluster(fake_self) is None
    assert calls['kwargs'] == {'gf_retrieved1': 'gf_retrieved_imp1', 'gf_retrieved2': 'gf_retrieved_imp2'}

    # r_offset on top of the first impurity is still refused
    fake_self.inputs = AttributeDict({'offset_imp2': Dict({'r_offset': [0., 0., 0.]})})
    assert combine_imps_wc.create_big_cluster(fake_self) == combine_imps_wc.exit_codes.ERROR_INPLANE_NEIGHBOR_TOO_SMALL


# run manual:
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    from aiida.orm import Code
    import pathlib

    test_dir = pathlib.Path('.')
    #define codes
    try:
        # on mac
        kkrhost_local_code = Code.get_from_string('kkrhost_intel19@localhost')
        kkrimp_local_code = Code.get_from_string('kkrimp_intel19@localhost')
    except:
        # on iff desktop
        kkrhost_local_code = Code.get_from_string('kkrhost@localhost')
        kkrimp_local_code = Code.get_from_string('kkrimp@localhost')

    # run test
    #test_combine_imps(None, kkrhost_local_code, kkrimp_local_code, None, nopytest=True)
    test_combine_imps_reuse_gf(None, kkrhost_local_code, kkrimp_local_code, None, nopytest=True)
