# TODO

Open work on `develop`: pull requests awaiting a decision, and known failures that are
deferred. Keep this file short; an item that is done is deleted here, and its record lives in
the merged PR or closed issue it links to. Last updated 2026-09-20.

## Open pull requests

- [ ] **[#186](https://github.com/JuDFTteam/aiida-kkr/pull/186) — Make the `tests` legs green.**
  Stacked on #179. Decouples the two pytest passes in `run_all.sh` (the unit pass used to abort the
  script, so the workflow pass never ran), pins `matplotlib < 3.11` instead of regenerating the
  plot baselines, `xfail`s `test_stm.py` for [#185](https://github.com/JuDFTteam/aiida-kkr/issues/185),
  fixes `np.float_` in `kkrnano.py` (removed in NumPy 2.0, so KKRnano input files could not be
  written at all), and imports `get_format` from `aiida.tools.archive.abstract`. Expected result:
  92 passed in the unit pass, 19 passed and 1 xfailed in the workflow pass, per leg. **The
  aiida-core 2.9 bump is deliberately not in it** — see the deferred item below.
- [ ] **[#179](https://github.com/JuDFTteam/aiida-kkr/pull/179) — Fix CI setup: `reentry`, unparsable
  `tools_STM_scan.py`, and the flynt crash on pre-commit.ci.** Rebuilt on `develop` (`c47af80`) on
  2026-09-20, so the earlier rebase conflict with #178's formatting commit is gone. Removes
  `reentry` / `reentry scan` from the test job in **both** `ci.yml` and `cd.yml` (it is aiida-core 1.x
  machinery and now dies on `ModuleNotFoundError: No module named 'pkg_resources'`), pins
  `setuptools < 81` for tests because `aiida-test-cache` 0.0.1 imports `pkg_resources`, fixes the
  over-indented `_version_ = 0.1` that makes `tools_STM_scan.py` untokenizable plus the undefined
  names pylint then exposes, and bumps the flynt hook from 1.0.1 to 1.0.6. The flynt bump is what
  fixes pre-commit.ci: its runner image is on Python 3.14, where `ast.Str` was removed, and flynt
  1.0.1 touches `ast.Str` at import. The in-repo `pre-commit` job never saw this because it runs
  Python 3.12. With #179, `pre-commit`, `pre-commit.ci` and `docs` go green and the `tests` legs
  reach the suite; the `tests` legs stay red on the known failures below, which the follow-up PR
  handles together with the aiida-core bump.
- [ ] **[#175](https://github.com/JuDFTteam/aiida-kkr/pull/175) — pre-commit.ci autoupdate** (bot, since 2025-08).
  Its flynt 1.0.1 → 1.0.6 bump is now carried by #179 instead, so #175 is no longer on the critical
  path; it still bumps every other hook and should be rebased or closed once #179 lands.
- [ ] **[#115](https://github.com/JuDFTteam/aiida-kkr/pull/115) — Implement base restart functionality** (open since 2022-12).
  Not reviewed as part of the current work.

Merged recently: [#183](https://github.com/JuDFTteam/aiida-kkr/pull/183) (issue
[#180](https://github.com/JuDFTteam/aiida-kkr/issues/180), closed), `kkr_flex_wc` now expands
`{username}` in the Green's function upload path — verified live against a templated computer.
[#184](https://github.com/JuDFTteam/aiida-kkr/pull/184) (issue
[#181](https://github.com/JuDFTteam/aiida-kkr/issues/181), closed), `kkr_imp_wc` returns exit codes
145/146 instead of excepting on a failed sub-workflow; no end-to-end live check, and **in-flight
`kkr_imp_wc` cannot be resumed across that upgrade** — plumpy persists the outline position as a
bare index.

Earlier: [#178](https://github.com/JuDFTteam/aiida-kkr/pull/178), non-finite parser output
gets exit codes 303/304 and `kkr_imp_sub_wc` exit code 134 (issue
[#177](https://github.com/JuDFTteam/aiida-kkr/issues/177), closed). #177 also lists smaller existing
problems that #178 did not address.

- [ ] **#178 is only half verified against real data.** Exit **304** (KKRhost parser) is confirmed
  on a live profile: two stored calculations that diverged in 2024 were replayed, parsed to 304 with
  `nonfinite_values` populated (300 and 139 entries), `parser_version` 0.9.0 read back off the output
  node, and the output `Dict` **stored successfully** — the store being the operation that used to
  raise. Exit **303** (KKRimp parser) and `kkr_imp_sub_wc` exit **134** versus the generic 130 have
  never been run against real data and remain covered only by unit tests.

## Known bugs, filed but not yet fixed

- [ ] **[#182](https://github.com/JuDFTteam/aiida-kkr/issues/182) — `find_cluster_radius`: unchecked
  index, and neighbour distances wrong off-origin.** Two defects two lines apart. It searches a
  fixed 10 Å radius and then indexes `dist_all[nclsmin - 2]` without checking enough neighbours were
  found, so it raises `IndexError` on any lattice open enough that the requested cluster does not
  fit — this blocked all 60 embeddings of one host in a production batch. Separately, it measures
  `np.linalg.norm(n.coords)`, the neighbour's absolute position rather than its displacement from
  the central site, which is correct only for a site at the origin; off-origin it counts the central
  site's own image at distance 0 and returns a radius that is **too small, silently**. Stored radii
  for multi-site cells may already be affected; how many is not known and no count is scheduled.

## Known bugs found while fixing #177, not yet addressed

Details in [#177](https://github.com/JuDFTteam/aiida-kkr/issues/177), section "Found along the way".

- [ ] **`kkr_imp_sub_wc.condition()` returns an `ExitCode` from a `while_` predicate.** plumpy
  only warns and treats it as `True`, so that branch never aborts the workchain.
- [ ] **`kkr_imp_sub_wc`: `ctx.kkr_step_success` is never set to `False` after the first step.**
  `update_kkrimp_params` checks it to decide on reducing the mixing factor, but `inspect_kkrimp`
  sets `ctx.kkrimp_step_success` instead.
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

- [ ] **aiida-core cannot move past 2.5.x until the test archives are re-exported.** Measured
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

Resolved 2026-09-20: the stale plot baselines. They were not stale — `tests/files/baseline_images/`
still matches matplotlib 3.10.9 exactly for five of eight comparisons and within tolerance for the
rest. Only matplotlib 3.11 disagrees, by RMS 6.5-18 against tolerances of 2-8, and matplotlib 3.11
requires Python >= 3.11 so the 3.10 leg cannot use it anyway. Pinned `matplotlib < 3.11` in the
`testing` extra instead of regenerating.
