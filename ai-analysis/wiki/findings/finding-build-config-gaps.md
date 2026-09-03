---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: major
evidence-commit: 6ffd38b12
tags:
  - neopz
  - build
  - cpp
---

# Build-configuration gaps: default build gets neither debug checks nor release paths; 14 dead DEBUG blocks; warnings off

## Repository evidence (verified first-hand)
1. `CMakeLists.txt:82-87`: `PZNODEBUG` only under `$<CONFIG:RELEASE>`, `PZDEBUG` only under `$<CONFIG:DEBUG>` (plus RELEASE-only `ZERO_INTERNAL_RESIDU`/`MAKEINTERNAL`). The **default build type is RelWithDebInfo** (`cmake/StandardPZSettings.cmake:1-14` [agent, consistent with observed builds]) — it matches neither generator expression ⇒ the recommended default compiles with **no `PZDEBUG` assertions and no `PZNODEBUG`/release fast-paths**; 142 headers gate on PZDEBUG [agent count]. (Observed: the user's own `NeoPZ_divfree_build` is RelWithDebInfo.)
2. **14 dead `#ifdef DEBUG` blocks** (macro never defined by the build — project macro is `PZDEBUG`): verified sample `Mesh/TPZCompElHDivDuplConnects.cpp:62-66` (a bounds check + DebugStop that never compiles); count verified via grep across Mesh/Shape/Material. Includes checks in `pzelchdiv.cpp`, `pzelchdivbound2.cpp`, `TPZMixedElasticityND.cpp`, `TPZShapeHDivOptimized.cpp` [agent].
3. **Warnings effectively off**: zero `-Wall/-Wextra/-Werror` in the build (verified grep); only `-Wsuggest-override`, `-Wno-narrowing`, `-Wno-alloc-size-larger-than`; and on Apple, `XCODE_ATTRIBUTE_WARNING_CFLAGS ""` explicitly silences Xcode warnings (`CMakeLists.txt:88-90`, verified).

## Why it matters
These three interact: the checks that would catch defects (PZDEBUG guards, dead DEBUG guards, compiler warnings) are all disabled in the configuration users actually build. Several defects found in this assessment (dead guards around exactly the duplicated-connect code; narrowing suppressions) would likely have surfaced earlier with `-Wall` + working assert config.

## Classification
**Confirmed defects (accidental)** — none is a domain tradeoff; all are hours-scale fixes: add a RelWithDebInfo case (or key on `NDEBUG`), global `DEBUG`→`PZDEBUG` rename, `-Wall -Wextra` at least in CI.

## Related
[[finding-debugstop-throws-release]] · [[finding-approx-creator-hygiene]] · CPP_TECHNICAL_REVIEW §2-H5 · TESTING_AND_VALIDATION_REVIEW (CI angle)
