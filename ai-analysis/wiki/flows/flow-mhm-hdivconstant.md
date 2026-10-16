---
type: flow
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - mhm
  - hdiv
---

# Flow: MHM_HDivConstant — multiscale hybrid-mixed Darcy on polygonal partition

Driver: `divfreebubbles/targets/main_MHM_HDivConstant.cpp` (read :60-240 [repo]). Binary from **Mar 27** — older than current sources (see drift note below).

## Shallow trace
1. **Problem**: Laplace/Darcy (`TLaplaceExample1::EX` exact) on a domain partitioned into polygonal subdomains (UNSW quadtree file), solved with [[mhm]] using `EHDivConstant` local spaces and **rigid-body spaces on** (needed for total condensation of subdomain interiors).
2. **Start**: `RunMHM<tshape>(xdiv,pOrder)` (:67).
3-4. **Mesh**: `ReadUNSWQuadtreeMesh("polygon00.txt", elpartition, scalingcenterindices)` (app-side `Common.cpp`) → `TPZAutoPointer<TPZGeoMesh>`; polygons triangulated via `mhm_gcreator.CreateTriangleElements(gmesh, matmap, partition, scalingcenters)` — triangles fanned around scaling centers, element-partition vector tracks coarse cell ids (:92-123) → [[TPZGeoMesh]].
5-6. **Spaces**: `TPZMHMHDivApproxCreator mhm_ccreator(mhm_gcreator, gmesh)` (app-side, [[divfree-support-lib]]): `HdivFamily()=EHDivConstant`, `ProbType()=EDarcy`, `IsRigidBodySpaces()=true`, `SetShouldCondense(true)`, `HybridType()=ENone`, `SetPOrderSkeleton(pOrder)` → `BuildMultiphysicsCMesh()`; then **`PutinSubstructures(*multiCmesh)` + `CondenseElements(*multiCmesh)`** — coarse cells become `TPZSubCompMesh` substructures, interiors condensed (:135-164) → [[mhm]], [[static-condensation]].
7. **Materials**: inserted via `mhm_ccreator.InsertMaterialObjects(LaplaceExact)` (analytic-driven) (:146).
8-9. **Assembly/solve**: `TPZLinearAnalysis` + `TPZSSpStructMatrix<STATE,TPZStructMatrixOR>` (10 threads) + `TPZStepSolver ELDLt`; the `ESemi` branch would use `TPZMatRedSolver(...,EMHMSparse)` but is dead code here (:166-201).
10-11. **Outputs**: equation/element counts to stdout; **always-on debug block** (`if(1)`) writes `mphysics2.txt` + `cmesh_multi.vtk` (geo-mesh-of-cmesh via `TPZVTKGeoMesh::PrintCMeshVTK`) into CWD (:149-162) → [[vtk-output]]. No error computation active.
12. **Central classes**: `TPZMHMGeoMeshCreator`, `TPZMHMHDivApproxCreator`, `TPZSubCompMesh`, `TPZVTKGeoMesh`, `TPZSSpStructMatrix`.
13. **Research before judging**: MHM skeleton/subdomain construction correctness; rigid-body-space role in subdomain invertibility; relation of app-side creators to the in-library `Pre/TPZMHMHDivApproxCreator` (same name, two repos? verify which is used — include path decides) → Phase 4.

## Drift finding — CONFIRMED (Phase 4)
Resolved via enum git history (`git log -L`): `EDefault/ESparse/EMHMSparse` were all removed at HEAD `fbe9696` (the commit that introduced `ProblemOrigin` + the elasticity mode); `main_MHM_HDivConstant.cpp:184` additionally uses a 3-arg constructor that no longer exists (`TPZMatRedSolver.h:17-19` has only default + 2-arg). Target still registered in `targets/CMakeLists.txt:7-8` ⇒ **does not compile at HEAD**; binary is a Mar 27 relic. → [[finding-mhm-target-uncompilable]].
Also noted: `new TPZLinearAnalysis(...)` never deleted (:168) [app-side hygiene].

Related: [[mhm]] · [[static-condensation]] · [[hdiv-space]] · [[divfree-support-lib]] · [[vtk-output]]
