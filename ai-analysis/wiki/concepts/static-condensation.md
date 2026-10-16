---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
---

# Static condensation

**Idea.** Eliminate element-interior DOFs at element level (Schur complement of the interior block) so the global solve only sees interface/skeleton unknowns; recover interior afterwards.

**In NeoPZ.** `Mesh/pzcondensedcompel.h` (`TPZCondensedCompEl`) is a decorator over `fReferenceCompEl` exposing only `fActiveConnectIndexes` (NConnects/ConnectIndex overridden; condensed connects listed separately; `Unwrap()` reverses; internal Schur container = library `TPZMatRed` via `pzmatred.h` include; `fKeepMatrix` flag controls memory retention) — verified [repo pzcondensedcompel.h:20-110]. It wraps an element/group; `Mesh/pzelementgroup.h` groups elements pre-condensation; connect "Lagrange levels" order which DOFs are condensable (e.g. keep one pressure per element for zero-mean constraints — rigid-body spaces `fIsRBSpaces` in [[approx-space-creators]] exist exactly to make internal blocks invertible [repo:TPZApproxCreator.h:67-68]). `TPZSubCompMesh` provides the coarser-grained variant (whole submesh condensed) used by [[mhm]]. Matrix-side: `Matrix/pzmatred.h` `TPZMatRed` (K11/K01 reduction container) [agent]; app-side `TPZMatRedSolver` drives reductions iteratively ([[divfree-support-lib]]).
`TPZCompMesh::NEquations()` (condensed) vs `Solution().Rows()` (full) — observed in iter_elast [repo:257-272].

**Invariants to check (Phase 4).** Invertibility of condensed blocks (pivoting? symmetric LDLt assumptions); consistency of recovery step (`LoadSolution` path through condensed wrappers); interaction with equation filters and renumbering; correctness under `SetShouldCondense(false)` + later `GroupAndCondenseElements` (iter_elast does exactly this sequence [repo:224-233]).

**Session 2:** full code-level trace (Resequence partition rules, K11Reduced/UGlobal round trip, SetKeepMatrix memory mode, submesh Schur exposure, ordering constraints) now in [[condensation-groups-submeshes]]; connect-level mechanics in [[TPZConnect]]. Downstream: manual group+condense in [[app-iterative-saddle-point]] and [[app-mixed-elasticity]].

Related: [[hybridization]] · [[mixed-methods]] · [[structural-matrices]] · [[mhm]] · [[flow-iter-elast]] · [[condensation-groups-submeshes]]
