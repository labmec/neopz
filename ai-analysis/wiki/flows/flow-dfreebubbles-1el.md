---
type: flow
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - hdiv
  - vtk
---

# Flow: dFreeBubbles1el — minimal manual mixed H(div) with error + VTK

Driver: `divfreebubbles/targets/main_1element.cpp` (read :100-405 [repo]; 980 lines incl. many commented variants). Binary rebuilt Jul 2 — actively used.

## Shallow trace
1. **Problem**: Darcy/Laplace on a single-element 2D domain (well-type physical groups), mixed H(div)×L² formulation; the closest surviving driver to the repo's original "div-free bubbles" intent (H1 and DFB comparison paths present but commented, :194-275).
2. **Start**: straight-line `main` (:115).
3. **Mesh**: `TPZGmshReader` with explicit physical-name→matid map (`stringtoint[dim]["Surface"]=1` etc.), reads `MESHDIR + "1element.msh"`; immediately dumps `gmesh.vtk` via `TPZVTKGeoMesh::PrintGMeshVTK` (:124-145) → [[mesh-io-generators]], [[TPZGeoMesh]].
4. **Geometric elements**: one 2D quad element + 1D boundary lines + points (from the .msh physical groups).
5. **Computational meshes — the manual pattern** (pre-ApproxCreator style):
   - `FluxCMesh`: `TPZNullMaterial` per matid (space placeholder — physics lives only in the multiphysics mesh), `cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim)`, `SetDefaultOrder(pOrder=4)`, `AutoBuild()` (:307-332) — **layer-1 factory used directly** ([[approx-space-creators]]).
   - `PressureCMesh` (analogous, L²/discontinuous [not yet read — Phase 4]).
   - `MultiphysicCMesh`: combines {flux, pressure} into `TPZMultiphysicsCompMesh` with real material (`TPZMixedDarcyFlow`-family; body TBC) (:167-172).
   - The commented `FluxCMeshDFB` shows the *fully manual* kernel-H(div) construction: per-element `new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel)` + neighbor walking (`TPZGeoElSide::Neighbour()`) to instantiate wrap/point elements (:335-405) — direct evidence of the element-level API and of [[TPZCompElHDiv]] kernel variants.
6. **Space selection**: explicit `SetAllCreateFunctionsHDiv` (+ commented HDivConstant/HDivKernel alternates) — family switching at factory level.
7. **Materials**: `TPZNullMaterial` in atomic meshes; physical material in multiphysics mesh (body Phase 4).
8-9. **Assembly/solve**: `TPZLinearAnalysis an(cmesh); SolveProblemDirect(an,cmesh)` (helper wraps struct-matrix + `ELDLt` [agent]) (:182-183).
10. **Errors**: `an.SetExact(exactSolError2,2); util.ComputeError(an, "postprocessHdiv.txt")` (:189-191) → [[error-estimation-convergence]] — this slice exercises the error leg.
11. **VTK**: `PrintResultsMultiphysic(dim, meshvector, an, cmesh)` → multiphysics VTK output (helper wraps `TPZVTKGenerator` [agent; verify Phase 4]) (:186) → [[vtk-output]] — exercises the VTK leg.
12. **Central classes**: `TPZGmshReader`, `TPZNullMaterial`, `TPZCreateApproximationSpace` (via `ApproxSpace()`), `TPZMultiphysicsCompMesh`, `TPZCompElKernelHDiv` (commented path), `TPZKernelHdivUtils`.
13. **Research before judging**: how `TPZBuildMultiphysicsMesh` wires atomic→multiphysics connects; what `exactSolError2` actually is; wrap/point element roles in kernel-HDiv constructions.

## Quirks noted (app-side)
- `exactSol` computes a 4-well log-potential then **overwrites with `u=x, ∇u=(1,0)`** (:288-293) — dead code above live code inside one lambda; the effective exact solution is linear (space reproduces it exactly → expected ~machine-zero errors).
- Physics-free `TPZNullMaterial` carries `SetBigNumber(1e10)` (:316) — big-number penalty convention surfacing even in placeholder materials; understand `BigNumber`'s role in BC imposition (Phase 4, [[material-system]]).

Related: [[mesh-io-generators]] · [[approx-space-creators]] · [[TPZCompElHDiv]] · [[error-estimation-convergence]] · [[vtk-output]]
