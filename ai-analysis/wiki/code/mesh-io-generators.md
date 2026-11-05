---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - neopz
  - pre-processing
  - mesh-io
---

# Pre/ — mesh input, generation & analytic solutions

## Responsibility
Everything that happens *before* the computational mesh exists: import meshes, generate structured grids, provide analytic benchmark solutions; plus the [[approx-space-creators]] (own page) and [[mhm]] controllers.

## Key files [repo paths]
- `Pre/TPZGmshReader.h` — **in-tree parser** of Gmsh `.msh` (v3/v4?) files with physical-group → material-id mapping; gmsh is *not* a linked dependency. Used by divfreebubbles (`main_1element.cpp`, `main_2fractures.cpp`).
- `Pre/TPZGenGrid2D.h`, `TPZGenGrid3D.h`, `TPZAcademicGeoMesh.h`, `Mesh/TPZGeoMeshTools.h` (`CreateGeoMeshOnGrid`, `CreateGeoMeshSingleEl` — used by iter_elast [repo]) — structured generators.
- `Pre/TPZReadGIDGrid.h`, `TPZGMSHReadMesh.h` — older importers [agent].
- `Pre/TPZAnalyticSolution.h` — family of manufactured solutions: `TElasticity2DAnalytic` (used by iter_elast [repo:66]), `TElasticity3DAnalytic`, `TLaplaceExample1`, `TStokesAnalytic`… each provides exact solution + forcing consistent with a chosen problem type → backbone of convergence validation ([[error-estimation-convergence]]).
- `Pre/pzbuildmultiphysicsmesh.h` — utilities to combine/transfer atomic meshes ↔ multiphysics mesh (`TPZBuildMultiphysicsMesh::TransferFromMultiPhysics` etc.).
- Hybridization utilities + MHM controllers → [[approx-space-creators]], [[mhm]].

## Related
[[TPZGeoMesh]] · [[approx-space-creators]] · [[error-estimation-convergence]] · [[flow-dfreebubbles-1el]] · [[flow-mhm-hdivconstant]]

## Open questions
- Which `.msh` format versions the reader supports (v2? v4?) and how robust it is — relevant to app-repo mesh assets (Phase 4).
- `TPZAnalyticSolution` uses FAD (auto-diff) to derive gradients/forcings? (README mentions forward AD [repo]; confirm mechanism.)
