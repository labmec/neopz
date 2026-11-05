---
type: app-survey
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: "app @ 2025-12-18 HEAD; embedded neopz Australia25 @ 182a80985 (2025-12-16)"
tags:
  - neopz
  - downstream
  - gfem
  - fracture
  - sbfem
---

# App survey: GFEM (~/GitHub/GFemResearch)

> Downstream-usage evidence (Session 2). Claims [agent], load-bearing ones spot-verified [✓]. Citations refer to the app repo.

**What it is.** **Generalized FEM for fracture mechanics**: enriched H1 approximations carrying displacement-jump discontinuities across a crack and singular fields at the crack tip, 2D and 3D (Darcy + elasticity physics). Targets: `GFem2D`, `GFem3D` (penny-shaped crack), `SBFem2D` (derives the singular enrichment modes via a scaled-boundary eigen-analysis), `TestGFemEnrichment`, 2 Catch2 test targets, all over a shared `GFEM_library` of NeoPZ extensions.

**Library surface exercised:**
- **Element-level extension — the headline**: `TPZGFemCompElH1<TSHAPE> : TPZCompElH1<TSHAPE>` [✓ `NeoPZ_extensions/TPZGFemCompElH1.h:9`] overrides `ComputeShape`/`ComputeRequiredData` to multiply standard H1 shapes by an enrichment/jump function per integration point (product rule for gradients) — partition-of-unity enrichment injected through the documented `ComputeShape` seam, instantiated for linear/triangle/quad/tetra shapes. Partner subclass `TPZGFemCompMesh : TPZCompMesh` carries a `std::map<int64_t, GFemShapeFunctionType>` associating enrichment functions to **connects** [✓ `TPZGFemCompMesh.h:25,32`].
- **Multiphysics as space composition**: background H1 + enriched GFem (+ optional H1 "mirror") atomic meshes combined via `TPZMultiphysicsCompMesh::BuildMultiphysicsSpace`; the custom combined-spaces materials (`TPZGFemDarcyFlow`, `TPZGFemElasticity2D/3D` — each inheriting *both* the library single-space physics class and `TPZMatCombinedSpacesT` + `TPZMatErrorCombinedSpaces`) sum sub-mesh contributions in `Solution()`. No HDiv/HCurl/L2 anywhere — a pure H1/multiphysics stress test.
- **SBFem as a production tool**: `TPZBuildSBFem` + `TPZSBFemElementGroup::EigenValues()/LoadEigenVector()` [✓ `2D/SBFem2D.cpp:494,592`] extract crack-tip singular modes from a JSON crack definition — the library's scaled-boundary eigen machinery used to *generate basis functions* for another method ([[sbfem]]).
- **Custom StrMatrix/Matrix pair**: `TPZSSpMatRedStructMatrix : TPZStructMatrixT<TVar>, TPar` [✓ `NeoPZ_extensions/TPZSSpMatRedStructMatrix.h:15`] + app-side `TPZSparseMatRed : TPZMatrix` produce a K00(background)/K11(enrichment) block reduction solved by CG preconditioned with direct LDLt of K11 (`3D/GFem3D.cpp:978-1008`) — the same MatRed pattern as divfreebubbles, independently reimplemented for a conditioning problem specific to enrichment methods.
- **Conditioning research on the connect layer**: `TPZGFemOrthogonal` orthogonalizes enrichment DOFs against the background space patch-by-patch using local dense `SolveEigenProblem` + `TPZMatRed` — near-linear-dependence control that only enrichment methods need.
- **Geometry**: gmsh input; `gRefDBase.InitializeRefPatterns` + **`TPZRefPatternTools::RefineDirectional` toward the crack tip** [✓ `2D/GFem2D.cpp:254`]; `TPZArc3D` + `TPZGeoBlend` curved crack-tip sectors; programmatic node-by-node geometry for SBFem domains.
- Embedded neopz on branch `Australia25` with app-motivated fixes (eigenvector loading resize; linear-element subdivision) — same downstream-drives-library pattern as wann.

**What it teaches about the library.** (1) The `TPZCompElH1`/`ComputeShape` pipeline is open enough for arbitrary per-point basis enrichment without touching the library — the cleanest evidence for element-layer extensibility (contrast the ~500-line new-family cost noted in `CPP_TECHNICAL_REVIEW.md` §6: *modifying* a family is much cheaper than *adding* one). (2) `TPZMultiphysicsCompMesh` composes same-physics spaces (background+enrichment), not just different fields. (3) SBFem is live downstream, not a dormant breadth item. (4) MatRed-style block reduction is a recurring *user pattern* across independent apps — a strong argument for first-class library support.

Related: [[element-families]] · [[multiphysics-composition]] · [[sbfem]] · [[refinement-hanging-nodes]] · [[matrix-and-solvers]] · [[apps-overview]]
