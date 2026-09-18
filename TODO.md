# TODO

Open work on `develop`: pull requests awaiting a decision, and known failures that are
deferred. Keep this file short; an item that is done is deleted here, and its record lives in
the merged PR or closed issue it links to. Last updated 2026-09-18.

## Open pull requests

- [ ] **[#179](https://github.com/JuDFTteam/aiida-kkr/pull/179) — Fix CI: test setup and unparsable `tools_STM_scan.py`.**
  Removes `reentry` from the CI test job, pins `setuptools < 81` for tests (for `aiida-test-cache`),
  and fixes the syntax error and undefined names in `aiida_kkr/tools/tools_STM_scan.py` plus
  pre-commit formatting of the STM and BdG files. With it, CI reaches the tests again, but is still red
  because of the known failures below. Needs rebasing onto `develop` after #178: it touches the same
  files as #178's pre-commit.ci formatting commit.
- [ ] **[#175](https://github.com/JuDFTteam/aiida-kkr/pull/175) — pre-commit.ci autoupdate** (bot, since 2025-08).
  Among other bumps it moves flynt from 1.0.1 to 1.0.6, which fixes the pre-commit.ci crash
  `module 'ast' has no attribute 'Str'`.
- [ ] **[#115](https://github.com/JuDFTteam/aiida-kkr/pull/115) — Implement base restart functionality** (open since 2022-12).
  Not reviewed as part of the current work.

Merged recently: [#178](https://github.com/JuDFTteam/aiida-kkr/pull/178), non-finite parser output
gets exit codes 303/304 and `kkr_imp_sub_wc` exit code 134 (issue
[#177](https://github.com/JuDFTteam/aiida-kkr/issues/177), closed). #177 also lists smaller existing
problems that #178 did not address.

## Known bugs, fixed on a branch and awaiting verification

Both are committed and carry a stub regression test. Neither reproduces without a real AiiDA
profile on a computer with a templated work directory, so both wait on a live check before a
pull request is opened. Each workflow stamps `workflow_version` into its output node, so a
verification result records which code actually ran.

- [ ] **[#180](https://github.com/JuDFTteam/aiida-kkr/issues/180) — `kkr_flex_wc` does not
  substitute `{username}` into the work directory.** On a computer whose `workdir` is a template,
  the Green's function is written to a directory named literally `{username}` — silently where the
  parent is writable, or with a misleading `EACCES` where an earlier run already created it. The
  inner `KkrCalculation` finishes with exit 0 first, so each attempt pays for a full write-out and
  loses it. Fixed on `fix/gf-writeout-username` by expanding the template from the transport the
  step already opens; `kkr_flex_wc` 0.5.6 → 0.6.0.
- [ ] **[#181](https://github.com/JuDFTteam/aiida-kkr/issues/181) — `kkr_imp_wc` excepts instead
  of returning an exit code when a sub-workflow fails.**
  `construct_startpot` recorded a failure in `ctx.exit_code` and carried on into code assuming
  success, so a failed `kkr_startpot_wc` ended as `excepted` with `ValueError: max() arg is an
  empty sequence`, and a failed `kkr_flex_wc` with `NotExistentAttributeError` on a missing
  `workflow_info`. Same fault as #177/#178: an exception escaping where an exit code was already
  available. Fixed on `fix/kkr-imp-startpot-exit-code`; adds exit code 146, deletes an inverted
  predicate, `kkr_imp_wc` 0.9.3 → 0.10.0. **Deployment note:** the outline gained a step, so
  in-flight `kkr_imp_wc` instances cannot be resumed across this upgrade — plumpy persists the
  outline position as a bare index.

## Known bugs found while fixing #177, not yet addressed

Details in [#177](https://github.com/JuDFTteam/aiida-kkr/issues/177), section "Found along the way".

- [ ] **`kkr_imp_sub_wc.condition()` returns an `ExitCode` from a `while_` predicate.** plumpy
  only warns and treats it as `True`, so that branch never aborts the workchain.
- [ ] **`kkr_imp_sub_wc`: `ctx.kkr_step_success` is never set to `False` after the first step.**
  `update_kkrimp_params` checks it to decide on reducing the mixing factor, but `inspect_kkrimp`
  sets `ctx.kkrimp_step_success` instead.
- [ ] **`ignore_nan` default mismatch.** `KkrimpParser.parse` defaults to `True`,
  `KkrimpParserFunctions.parse_kkrimp_outputfile` in masci-tools to `False`.
- [ ] **masci-tools requirement is too loose.** The PyPI release masci-tools 0.15.0 rejects the
  `ignore_nan` and `doscalc` arguments the KKRimp parser passes; only masci-tools `develop` works,
  while `pyproject.toml` requires `masci-tools >= 0.4.8.dev5`.

## Known CI failures, deferred

Seen in the CI run of #179, the first run in a while that got past setup.

- [ ] **Plot baseline images are stale.** On Python 3.11 and 3.12 (matplotlib 3.11.2), the seven
  image comparisons in `tests/test_plot_kkr.py` fail against `tests/files/baseline_images/`
  (last updated 2022-11). On 3.10 (matplotlib 3.10.9) they pass. Decide whether to regenerate the
  baselines or pin matplotlib for tests. Because `run_all.sh` stops at the first failure, the
  workflow tests do not run on 3.11 and 3.12.
- [ ] **`workflows/test_stm.py` fails.** `IndexError: list index out of range` in
  `tools_STM_scan.lattice_generation` (line 623), called from `kkr_STM.get_scanning_positions`.
  This test could not run before #179, because the module did not import. STM is maintained by its
  contributors: first check whether an unmerged STM or BdG branch (including forks) already fixes
  this, together with the syntax error #179 repairs.
- [ ] **pre-commit.ci** still fails on flynt 1.0.1 (fixed by #175), and it still reports the
  local `pylint` hook although that hook is listed under `ci: skip`.
