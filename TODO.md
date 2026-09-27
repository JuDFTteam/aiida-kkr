# TODO

Open work on `develop`: pull requests awaiting a decision, and known failures that are
deferred. Keep this file short; an item that is done is deleted here, and its record lives in
the merged PR or closed issue it links to. Last updated 2026-09-27.

## Open pull requests

- [ ] **[#115](https://github.com/JuDFTteam/aiida-kkr/pull/115) — Implement base restart functionality** (open since 2022-12).
  Not reviewed as part of the current work.

## Release and branch layout

- [ ] **Rename `master` to `main`, bring it up to date, and cut a release that reaches PyPI.**
  The default branch was switched to `develop` on 2026-09-20; `master` still exists, holds zero
  commits that `develop` does not, and its tip is the 2024-12-05 release. Steps, in order: rename
  `master` to `main` on GitHub (which keeps redirects for old links and retargets any open pull
  request), merge `develop` into `main` — a fast-forward, since nothing diverged — bump the version
  in `aiida_kkr/__init__.py`, `pyproject.toml` and `.bumpversion.cfg`, then tag `vX.Y.Z`. The tag is
  what publishes: `cd.yml` fires on `v*` and uploads to PyPI only if the tests pass. Then **check
  the upload actually arrived**, rather than assuming a green workflow means a release.
  Two edits travel with it: [README.md](README.md) line 81 names `master` in the release steps, and
  `2.3.1` has been the version since before the last release, so the next number is a decision
  (2.4.0 is the honest one — `find_cluster_radius` returns different radii now, see the entry
  below). **Why it matters beyond tidiness:** the [#182](https://github.com/JuDFTteam/aiida-kkr/issues/182)
  fix is on `develop` only, so anyone installing from PyPI still gets the version that crashes on
  open lattices and mis-sizes clusters on off-origin structures.

## Known bugs and limitations, not yet fixed

- [ ] **Cluster radii computed by aiida-kkr >= 1.1.12 may be too small where a structure has no
  site at the origin.** [#182](https://github.com/JuDFTteam/aiida-kkr/issues/182), fixed by
  [#187](https://github.com/JuDFTteam/aiida-kkr/pull/187) going forward; already-stored radii are
  not corrected by the fix. The trigger is `natom_in_cls_min > 0` alone — `voro_start` then always
  computes the radius, and an `RCLUSTZ` in `calc_parameters` is no protection: it is read only in
  the `natom_in_cls_min <= 0` branch, and where both exist the larger wins, which is usually the
  computed one. One production database was counted on 2026-09-20 and has none, because every
  off-origin run there predates the rewrite — census in the issue comment. Any other database on
  aiida-kkr >= 1.1.12 needs its own count.
- [ ] **[#195](https://github.com/JuDFTteam/aiida-kkr/issues/195) — KKRimp leaves non-spherical channels below `QBOUND` out of `out_potential`, so the
  simple-mixing step hands a truncated potential to the Anderson (`IMIX 5`) step.** It runs at
  `QBOUND = threshold_aggressive_mixing`, and KKRimp also zeroes channels below `POT_NS_CUTOFF`
  (default 0.1 × `QBOUND`) in memory. From a converged start potential this costs 16–23 Anderson
  iterations instead of 1–2 (19 of 20 Cu embeddings, vimp-prediction pilot). Stop-gap on branch
  `feature/kkrimp-pass-input-pot` (PR to follow): after a simple-mixing step that converged in one
  iteration, hand on its input potential. Switch `pass_input_pot_after_quick_simple_mixing`,
  default on. Root fix in progress elsewhere: a KKRimp/KKRhost key `POT_NS_WRITE_CUTOFF` (jukkr,
  local), and masci-tools `kkrparams` support for it and `POT_NS_CUTOFF` (merged in
  [masci-tools#255](https://github.com/JuDFTteam/masci-tools/pull/255)). aiida-kkr can then
  set both keys on the simple-mixing step. `kkr_scf_wc` has the same structure; not measured.
- [ ] **`kkr_imp_sub_wc` falls back to simple mixing on every parameter update.**
  `update_kkrimp_params` reads the previous `IMIX` from a freshly created `kkrparams`, which is
  always `None`, so any update that is not an explicit switch (for example a reduced mixing factor
  after an unconverged Anderson step) restarts with `IMIX 0`. This is the restart path the false
  convergence of #192 (fixed in #199) rode on. Found while tracing it; not fixed.
- [ ] **`kkr_imp_sub_wc`: a calculation that reaches `QBOUND` in one iteration is never a restart
  anchor.** After an unconverged step that is not on track, `update_kkrimp_params` restarts from the
  latest earlier calculation whose last rms is below its first; with one iteration they are equal.
  With a high `threshold_aggressive_mixing` the first calculation always takes one iteration, so an
  Anderson (`IMIX 5`) step whose rms rises ends the workchain with exit 127 instead of retrying.
  Seen in workchain 881846 (Fe:Cu, threshold 1e3) while testing #192. Not fixed.
- [ ] **In-flight `kkr_imp_wc` cannot be resumed across the
  [#184](https://github.com/JuDFTteam/aiida-kkr/pull/184) upgrade**, because plumpy persists the
  outline position as a bare index. Processes must be allowed to finish or be restarted.

## Known bugs found while fixing #177, not yet addressed

Details in [#177](https://github.com/JuDFTteam/aiida-kkr/issues/177), section "Found along the way".

- [ ] **`ignore_nan` default mismatch.** `KkrimpParser.parse` defaults to `True`,
  `KkrimpParserFunctions.parse_kkrimp_outputfile` in masci-tools to `False`.
- [ ] **`Parser.parse_from_node` cannot drive any aiida-kkr parser.** All five parsers declare
  `def __init__(self, calc)`, a pre-AiiDA-2.x signature. aiida-core instantiates parsers two ways:
  the daemon uses the positional form `parser_class(self.node)`, which works, while
  `Parser.parse_from_node` uses the keyword form `cls(node=node)`, which raises
  `TypeError: __init__() got an unexpected keyword argument 'node'`. So production parsing has
  never been affected — only the documented re-parse and testing entry point, which is why this
  survived the 2.x migration unnoticed. (Statements quoted rather than line numbers, which drift
  between aiida-core versions.) Same fault class as
  [#180](https://github.com/JuDFTteam/aiida-kkr/issues/180): a pre-2.x signature left behind,
  invisible because the common path avoids it. Fixing it would let a parser change be validated by
  replaying real stored calculations instead of archived fixtures — which is also why the #178
  tests import an archive and hand-build a `CalcJobNode` rather than replaying a node.
  Found while validating #178 against real stored failures. Not filed upstream.
- [ ] **masci-tools requirement is too loose.** The PyPI release masci-tools 0.15.0 rejects the
  `ignore_nan` and `doscalc` arguments the KKRimp parser passes; only masci-tools `develop` works,
  while `pyproject.toml` requires `masci-tools >= 0.4.8.dev5`.

## Known CI failures, deferred

- [ ] **[#198](https://github.com/JuDFTteam/aiida-kkr/issues/198) — CI pins masci-tools to `c90e5815`, the `develop` commit before masci-tools#255, until the
  test archives are re-exported.** #255 added three `kkrparams` keys. `get_dict()` returns every key,
  unset ones as `None`, so every parameter `Dict` the tests build changes its hash and misses the
  archive cache. That fell through to the fake `kkr.x` in 5 tests. The re-export below lifts the pin too.
- [ ] **[#198](https://github.com/JuDFTteam/aiida-kkr/issues/198) — aiida-core cannot move past 2.5.x until the test archives are re-exported.** Measured
  2026-09-20 across four CI runs: on aiida-core 2.9.2 with `aiida-test-cache@main`, 16 of the 20
  whitelisted workflow tests fail; on aiida-core 2.5.2 with the PyPI `aiida-test-cache` 0.0.1, 19
  of 20 pass. Identical on Python 3.10, 3.11 and 3.12. The failures are archive-cache misses —
  `kkr.voro` calculations hit the cache 18/18 while `kkr.kkr` calculations miss 6/6 — which fall
  through to the zero-byte fake `kkr.x` and surface as `KeyError` on a workflow output. The
  archives in `tests/data_dir/` are already at the current `export_version` (`main_0001`) and
  `tests/migrate_exports.py` migrates them every run, so this is **not** a schema-migration
  problem; their `metadata.json` records `aiida_version 2.5.2`, and they need re-**exporting**
  under the new aiida-core, which is the follow-up PR. Two constraints travel with this: PyPI
  `aiida-test-cache` 0.0.1 pins `aiida-core < 2.6` and needs `setuptools < 81` for
  `pkg_resources`, while `aiida-test-cache@main` needs aiida-core >= 2.6 despite declaring
  `>= 2.1` (it calls `NodeCaching.compute_hash`, added in 2.6).
- [ ] **[#185](https://github.com/JuDFTteam/aiida-kkr/issues/185) — `workflows/test_stm.py` is
  `xfail`ed.** `IndexError: list index out of range` in
  `tools_STM_scan.lattice_generation`, from `kkr_STM.get_scanning_positions`. The failing
  expression is the bounds test `p[0] < xmax and … p[1] < ymax`, where
  `p = [i * x + j * y for x, y in zip(vec[0], vec[1])]` is shorter than 2 when the plane vectors
  are not what the code assumes. `feature/kkr-bdg-workflow` carries the identical line, so no
  unmerged branch fixes it. Same failure class as [#182](https://github.com/JuDFTteam/aiida-kkr/issues/182)
  (`find_cluster_radius`): unchecked geometry indexing in `aiida_kkr/tools`. STM is maintained by
  its contributors. Marked `strict=False` so an upstream fix does not turn the suite red.
- [ ] **pre-commit.ci** reports the local `pylint` hook as skipped even though that hook is listed
  under `ci: skip`; cosmetic.
