---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - adaptivity
---

# h-refinement, refinement patterns & hanging nodes

**Idea.** Subdivide elements (h-refinement) possibly non-uniformly → "hanging" nodes on interfaces between refinement levels; conformity restored by constraining hanging DOFs to coarse-side DOFs (dependency/constraint matrices), or by pattern-conforming closures.

**In NeoPZ (distinctive design).** *Runtime-defined refinement patterns*: `Refine/TPZRefPattern.h` + `TPZRefPatternDataBase` + 71 `.rpt` data files (`Refine/RefPatterns/`) describe father→sons subdivisions as little meshes; `TPZGeoElRefPattern<TGeo>` applies them ([[TPZGeoMesh]]). Uniform refinement via per-topology `TPZRef*` classes. Hanging-node constraints live on connects: `TPZConnect` dependency matrices ([[TPZCompMesh]]), built by `TPZInterpolatedElement` restraint logic; validated by `TestHangingNode`, `TestCondensedSpace` ("Constrained Space"), and a refinement suite [agent]. README claims hp-adaptivity + hanging-node support as headline features [repo:README.md:16-19].

**Traced (Session 2, [[geometry-refinement-maps]] §3 + [[TPZConnect]]).** The full path is now line-cited: geometric `Divide` → `CreateMidSideConnect` consults `LowerLevelElementList` (ancestor walk) → `RestrainSide` builds the L2 projection of the coarse trace through `SideTransform3` → `TPZConnect::AddDependency`; constraints resolve per element in `ApplyConstraints` (complex-correct congruence transform), topologically ordered. H(div)/H(curl) have their **own** `RestrainSide` overrides (HCurl DebugStops when small-side order < large order). Dependency closure enforced by `BuildDependencyOrder` fixpoint. Pattern compatibility fails loudly (midnode mismatch DebugStop). Downstream: directional refinement to wells/crack tips ([[app-wann]], [[app-gfem]]); custom patterns as macro-element constructors ([[app-mixed-elasticity]]); hand-built `AddDependency` for cross-dimensional coupling ([[app-wann]]).
**Remaining Phase-7 gap.** Constrained H(div)/H(curl) under nonuniform refinement is still the thin *test* area (the code paths exist).

**Reference anchors.** Devloo's early adaptivity papers (Devloo–Oden 1987-89 line); Demkowicz hp book (constrained approximation); Šolín et al. as contrast.

Related: [[TPZGeoMesh]] · [[TPZCompMesh]] · [[shape-functions]] · [[hp-adaptivity]]
