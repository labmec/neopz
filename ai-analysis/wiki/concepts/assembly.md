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

# Assembly (element → global system)

**Idea.** Loop elements; compute element stiffness `ek` and load `ef` by quadrature over material `Contribute`; scatter into the global matrix/rhs via connect → equation numbering; apply constraints (hanging nodes, condensation) on the way.

**In NeoPZ.** Pipeline (to be traced precisely in Phase 2/4):
`TPZAnalysis::Assemble` → [[structural-matrices|TPZStructMatrix]]`::Assemble` → per-element `TPZCompEl::CalcStiff(ek,ef)` → `TPZInterpolatedElement` gathers [[shape-functions]] + [[geometric-mappings]] into `TPZMaterialData` → [[material-system|material]]`::Contribute` at each [[quadrature]] point → `ek/ef` constrained (connect dependencies, condensation wrappers) → scatter through `TPZConnect` sequence numbers + `TPZBlock` offsets into the chosen matrix storage ([[matrix-and-solvers]]).
Parallel variants (OR/OT/TBB-flow) partition the element loop; equality with serial assembly is unit-tested (`TestMultithreading` [agent]).

**Invariants to check.** Constraint application order (dependencies before scatter); symmetric-storage assembly writes only one triangle (sym sparse `TPZSYSMPMatrix` — delta file caution); block offsets vs connect sequence renumbering ([[TPZAnalysis]] RenumType); thread-safety of shared `TPZMaterialData`/materials (materials are shared across threads — `Contribute` must be const/reentrant? → Phase 5).

Related: [[structural-matrices]] · [[material-system]] · [[static-condensation]] · [[quadrature]] · [[TPZCompMesh]]
