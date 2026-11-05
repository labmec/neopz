---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: insufficient-evidence
severity: major
evidence-commit: 6ffd38b12
tags:
  - neopz
  - testing
  - working-tree
---

# 5 of 40 unit-test suites crash in the local working-tree build (pin itself is CI-green)

## Runtime evidence [run, 2026-07-02, ctest on `NeoPZ_divfree/build` (Release, BUILD_UNITTESTING=ON)]
`ctest -j3`: 35/40 pass in 16.6 s; **failures**: `TestCondensedHangingNodes` (SIGTRAP — reaches "Testing Hanging Nodes… test 0 EH1", then traps ⇒ consistent with a `DebugStop`), `TestReduced`, `TestErrorAnalysis`, `TestHangingNode`, `TestSBFem` (all Bus error). Full log: job tmp `ctest-run.log`.

## Attribution analysis
- Stale-ABI hypothesis **eliminated**: all failing binaries timestamp Jun 30 16:29, same session as `libpz.dylib` (16:28) [repo ls].
- Build identity: stamped `PZ_REVISION "852a5116c"`; built at 16:28, i.e. **before** HEAD commit `4de234fae` (16:47) — the build = 852a5116c + then-uncommitted working-tree edits (which became the MultiplyByScalar-virtual commit).
- **Upstream GitHub Actions "Run Unit Tests" is `success` on develop at both `6ffd38b` (the analysis pin) and `852a511`** [web: api.github.com runs list]. So the pinned develop is green in CI; the crashes belong to the *local working-tree state and/or this machine's configuration*.
- Suggestive pattern: 4 of 5 failing suites exercise condensation / hanging-node constraints / reduced spaces — the exact area refactored by the recent commits (`GroupElements`/`CondenseElements` split, virtual changes) — [inference, not proof].

## Classification
**Insufficient evidence / requires a rebuild to attribute** (options: (a) uncommitted-edit breakage now embodied in `4de234fae`, (b) machine/config-specific (Release+macOS+Accelerate vs CI images), (c) flaky memory bug surfacing locally). *Not* counted against the pinned develop, which is CI-green. Severity major because a current-branch user sees 12.5% of suites crashing.

## What would resolve
Rebuild the build tree at `4de234fae` (or current HEAD) and rerun the 5 suites; if still crashing, bisect the 3-commit delta; run one suite under lldb to get the trap site. (Fresh builds were out of this engagement's authorized scope.)

## Related
[[refinement-hanging-nodes]] · [[static-condensation]] · [[error-estimation-convergence]] · [[finding-hybridelasticity2d-missing-rhs-at-pin]] (same delta window)
