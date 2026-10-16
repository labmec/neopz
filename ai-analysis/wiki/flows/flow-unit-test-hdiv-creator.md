---
type: flow
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - testing
  - hdiv
  - de-rham
---

# Flow: NeoPZ unit-test slice — TestHDivApproxSpaceCreator + TestDeRham

In-library validation flows (what NeoPZ itself proves about its spaces). Files: `UnitTest_PZ/TestHDivApproxSpaceCreator/TestHDivApproxSpaceCreator.cpp` (1081 lines) and `UnitTest_PZ/TestDeRham/TestDeRham.cpp` (721 lines) [repo].

## TestHDivApproxSpaceCreator — "linear solution representation" grid
- One `TEST_CASE("HDiv Approx Space Creator","[hdiv_linear_solution_representation]")` (:152) fanning a nested `GENERATE` grid [repo :152-216]:
  `HDivFamily {EHDivConstant, EHDivStandard, EHDivOptimized}` × `ProblemType {EDarcy, EElastic}` × `MMeshType {Quad, Tri, Tetra, Hexa}` × `pOrder {1,2}` × `extraporder {0,1}` × `HybridizationType {ENone, EStandard, ESemi}` × `isRBSpaces {0,1}` × `isCondensed {0,1}` × `isRef {0,1}` (+ an MHM variant instantiating `TPZMHMHDivApproxCreator` — **the library has its own MHM creator**, cf. same-named app-side class, :519-524).
- Per combination `TestHdivApproxSpaceCreator(...)` (:451): builds gmesh → creator → multiphysics cmesh → solve → checks:
  `CheckIntegralOverDomain` (flux/pressure integrals vs analytic values), `CheckError` (errors ≈ expected for a linear exact solution), `TestKnownSol` (constant-solution reproduction), `CheckNEqCondensedProb` (equation counts after condensation), `PostProcessVTK` (:56-72 decls).
- **What it proves**: the creator pipeline produces spaces that exactly represent linear/constant solutions across families × hybridizations × condensation × refinement, with correct condensed equation counts. **What it does not prove**: convergence *rates*, curved geometry, 3D prism/pyramid, stability constants → [[error-estimation-convergence]], Phase 7.

## TestDeRham — basis-level exact-sequence checks
- `TEMPLATE_TEST_CASE("Dimension Compatibility","[derham_tests]", dim=2|3)` + `"Inclusion"`; LAPACK-gated (SVD) (:75-120).
- `CheckRankKerDim<left,right,dim>(k)`: builds mass-like matrices of `op(φ_i)` and verifies `rank(M_left) = ker(M_right)` for pairs H1→HCurl, HCurl→L2 (2D) / HCurl→HDiv (3D), HDiv→L2, HDivConst→L2; `CheckInclusion` verifies range(op(left)) ⊆ span(right) via block-matrix rank identity `rank([A B; C D]) = rank(D)` (:49-73 header comments quoted).
- Dedicated pair materials (`TPZMatDeRhamH1HCurl` etc., 7 helper materials in dir listing [repo]) assemble the mixed Gram matrices.
- **What it proves**: dimension-level exactness/inclusion of the discrete sequence on the tested meshes/orders (k=1..3). **What it does not prove**: commuting-diagram/interpolation properties, exactness on curved or distorted elements → [[de-rham-complex]].

## Flow anatomy (both)
gmesh (`TPZGeoMeshTools`) → atomic cmeshes → `TPZMultiphysicsCompMesh` / basis matrices → assembly (`TPZFStructMatrix`/creators) → LAPACK SVD or solve → Catch2 `REQUIRE(... Approx ...)`. Central classes: `TPZHDivApproxCreator`, `TPZMHMHDivApproxCreator` (lib), `TPZLinearAnalysis`, `TPZFMatrix::SVD`.

Related: [[de-rham-complex]] · [[hdiv-space]] · [[approx-space-creators]] · [[mixed-methods]] · [[flow-dupl-connects]]
