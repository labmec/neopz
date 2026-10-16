---
type: flow
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - hdiv
  - darcy
---

# Flow: dupl_connects2 — semi-hybrid mixed H(div) Darcy benchmark

Driver: `divfreebubbles/targets/dupl_connects.cpp` (main read :140-340 [repo]). Binary `build/targets/dupl_connects2` (Jun 30).

## Shallow trace
1. **Problem**: mixed Darcy (flux×pressure) on unit square/cube, `HDivFamily::EHDivConstant`, **semi-hybridization** (`HybridizationType::ESemi`), benchmark of `TPZMatRedSolver(EDarcyHDiv)` vs direct; DIM from argv (default 3D, idivs {2,8,12,16,32}) (:155-214).
2-4. **Mesh**: same `CreateGeoMesh<tshape>` → `TPZGeoMeshTools::CreateGeoMeshOnGrid` pattern as iter_elast; quads (2D) / hexes (3D).
5-6. **Spaces**: **library-side** `TPZHDivApproxCreator hdivCreator(gmesh)`: `HdivFamily()=EHDivConstant`, `SetProbType(EDarcy)`, `IsRigidBodySpaces()=false`, `SetExtraInternalOrder(0)`, `SetShouldCondense(true)`, `SetHybridType(ESemi)` → `CreateApproximationSpace()` → `TPZMultiphysicsCompMesh` (:206-237). Contrast with [[flow-iter-elast]]: creator + condensation handled fully in-library here → cleanest Phase 4 window into `Pre/TPZHDivApproxCreator.cpp` and the duplicated-connects machinery behind ESemi ([[hybridization]]).
7. **Material**: `TPZMixedDarcyFlow(EDomain, DIM)`, constant permeability 1, exact + forcing lambdas (order 4); Dirichlet BC type 0 from exact (:223-234) → [[material-system]], [[mixed-methods]].
8-9. **Assembly/solve**: identical benchmark switch as iter_elast: `TPZMatRedSolver(an, EDarcyHDiv).Solve()` vs `TPZSSpStructMatrix`/`Mumps` + ELDLt (:266-323).
10-11. **Outputs**: equation counts to `results_Harmonic2D.txt`, time/memory to `results_memory_time.txt`; error & VTK legs commented out (same pattern as iter_elast).
12. **Central classes**: `TPZHDivApproxCreator`, `TPZMixedDarcyFlow`, `TPZMultiphysicsCompMesh`, `TPZMatRedSolver`, `TPZSSpStructMatrix(Mumps)`.
13. **RESOLVED (Phase 4 trace)**: ESemi duplicates each interior facet connect into even=constant-flux + odd=higher-order pairs (`TPZCompElHDivDuplConnects`, ratio verified: even connect gets 1 shape function, odd gets nshape−1, `TPZCompElHDivDuplConnects.cpp:58-138`); `SemiHybridizeDuplConnects` rebinds only the even connect on the `sideOrient==−1` side to the wrap element (`TPZHDivApproxCreator.cpp:1239-1299`); multiplier submesh order 0; glue = shared connect + `TPZMultiphysicsInterfaceElement`/`TPZLagrangeMultiplierCS`. `EHDivConstant` (order-0 pressure) + per-facet constant-flux multipliers is what makes total condensation well-posed (cf. the singular-K00 guard `TPZHDivApproxCreator.cpp:85-89`). Full details: [[hybridization]].

## Quirks noted (app-side)
- The analytic pair is inconsistent as a manufactured solution: `exactSol` = 3D harmonic (`sin·sin·sinh`, Δu=0, and ≡0 on z=0 for 2D runs) yet `forcefunction` is nonzero (and uses `cosh(√2πx)` — x, not z) (:144-153). Harmless for the timing benchmark (error computation disabled), but any re-enabled error check would be meaningless. [hypothesis: copy-paste vestige]

Related: [[hdiv-space]] · [[hybridization]] · [[static-condensation]] · [[mixed-methods]] · [[flow-iter-elast]]
