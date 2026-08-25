---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - application
---

# App repo: voronoi_mixed_elas dereferences null gAnalytic on its active path

## Repository evidence
`divfreebubbles/targets/voronoi_mixed_elas.cpp` [repo]: `:84` `TPZAnalyticSolution *gAnalytic = 0;`; the only assignments (`gAnalytic = elas;`) are commented out (`:93`, `:103`); `:143` executes `an.SetExact(gAnalytic->ExactSolution());` → null dereference when the target runs. The existing binary (built Mar 27) predates or was built from a different source state.

## Assessment
- **Classification: confirmed bug — application repo, not NeoPZ.** Severity minor (research driver; crashes at startup of the error-setup stage, no silent wrong numbers).
- Value for the NeoPZ assessment: (a) evidence that app drivers drift out of compilable/runnable state without CI (`BUILD_TESTS=OFF`, no app CI) — same theme as the `EMHMSparse` drift (OQ6); (b) API-design observation: `TPZAnalysis::SetExact` taking a callable forces the *caller* to null-check the provider — a fluent guard (accepting a provider object) would fail softer → Phase 5 note.

## What would resolve
Compile+run the target (out of scope: fresh builds not authorized) or user confirmation. Fix is app-side: uncomment the intended `gAnalytic = elas` (or guard the SetExact call).

## Related
[[divfree-support-lib]] · [[TPZAnalysis]] · [[flow-mhm-hdivconstant]] (drift sibling)
