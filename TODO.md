# TODO

Open work on `develop`: pull requests awaiting a decision, and known failures that are
deferred. Keep this file short; an item that is done is deleted here, and its record lives in
the merged PR or closed issue it links to. Last updated 2026-09-20.

## Open pull requests

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

Seen in the CI run of #179, the first run in a while that got past setup.

- [ ] **Plot baseline images are stale.** On Python 3.11 and 3.12 (matplotlib 3.11.2), the seven
  image comparisons in `tests/test_plot_kkr.py` fail against `tests/files/baseline_images/`
  (last updated 2022-11). On 3.10 (matplotlib 3.10.9) they pass. Decide whether to regenerate the
  baselines or pin matplotlib for tests. This used to also hide the workflow pass on 3.11 and
  3.12, because `run_all.sh` let the unit pass abort the whole script; the two passes are now
  decoupled, so the workflow results show up alongside these failures rather than behind them.
- [ ] **`workflows/test_stm.py` fails.** `IndexError: list index out of range` in
  `tools_STM_scan.lattice_generation` (line 623), called from `kkr_STM.get_scanning_positions`.
  This test could not run before #179, because the module did not import. STM is maintained by its
  contributors: first check whether an unmerged STM or BdG branch (including forks) already fixes
  this, together with the syntax error #179 repairs.
- [ ] **pre-commit.ci** reports the local `pylint` hook as skipped even though that hook is listed
  under `ci: skip`; cosmetic. Its flynt 1.0.1 crash is fixed by #179.
