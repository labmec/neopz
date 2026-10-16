---
type: concept
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
---

# H1-conforming spaces

**Idea.** Piecewise-polynomial spaces continuous across element boundaries; the natural home of primal formulations (Poisson, elasticity displacement). Conformity requirement: function traces match on shared faces/edges/vertices.

**In NeoPZ (Session-2 trace, [[element-families]]).** `TPZCompElH1<TSHAPE> : TPZIntelGen : TPZInterpolatedElement`; continuity by *shared connects*: one connect per topological side including corners; corner connects clamped to order 1, side/internal orders per connect (hierarchical `TPZShapeH1` basis) → native p-adaptivity. `EffectiveSideOrder` = max over contained sub-sides (face ≥ edge order rule). Hanging nodes fully supported through the generic scalar `RestrainSide` L2 projection ([[TPZConnect]]).
**H1Family resolved**: `{EH1Standard, EH1WidePrism}` — the family only matters for prisms and only at element creation (template argument `TPZShapeWidePrism` vs `TPZShapePrism`, `pzcreateapproxspace.cpp:111-121`); the stored `fh1fam` is never branched on at runtime.
Creators: `SetAllCreateFunctionsContinuous` (+`WithMem`), problem-level `TPZH1ApproxCreator` (hybrid variants, [[hybridization]]). Broken-H1 (disconnected build) doubles as the p>0 L² realization ([[discontinuous-l2-dg]]).

**Downstream evidence.** H1 is a primary research surface, not a baseline: GFEM builds enriched fracture spaces by subclassing `TPZCompElH1` ([[app-gfem]]); ErrorEstimation reconstructs conforming potentials in H1 ([[app-error-estimation]]); wann uses H1 companion meshes as error references ([[app-wann]]).

**Validated in-tree.** `TestH1ApproxSpaceCreator` (constant/linear exact representation), De Rham pairs H1→HCurl ([[de-rham-complex]]), hanging-node suites. Gap (Phase 7): no convergence-*rate* tests.

**Reference anchors.** Devloo–Bravo–Rylo 2009 (systematic shape construction); Szabó–Babuška (p-version); Ern–Guermond.

Related: [[element-families]] · [[shape-functions]] · [[hybridization]] · [[de-rham-complex]] · [[TPZConnect]]
