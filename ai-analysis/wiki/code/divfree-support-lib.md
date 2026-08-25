---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - application
  - hdiv
---

# divfree/ — application support library (repo: ../divfreebubbles)

> Application-repo code (branch `3DKernelHdiv`), **not** part of NeoPZ. It is the vehicle for the execution-flow analysis and shows how downstream research code extends the library.

## Contents [repo paths under `divfreebubbles/divfree/`; purposes agent-cited, key ones verified via iter_elast]
- **Materials**: `TPZMatDivFreeBubbles.h` (namesake Laplace-with-divfree-bubbles), `TPZMixedDarcyH1.h`, `TPZMixedDarcyFlowHybrid.h`, `TPZMixedDarcyFlowOrtotropic.h` (+`TPZOrtotropicPermeability`), `TPZMatCurlDotCurl.h`.
- **Creators**: `TPZH1HybridApproxCreator.h` — derives NeoPZ `Pre/TPZH1ApproxCreator` (a develop-delta file); adds `ComputeOrthogonalizingRestraints`, `HybridizeLowOrderFluxes`, `GroupAndCondenseElements` (called in `targets/iter_elast.cpp:229-233` [repo]). `TPZMHMGeoMeshCreator.h`, `TPZMHMHDivApproxCreator.h`, `TPZMixedElasticityCMeshCreator.h`.
- **Custom elements** (excluded from current build [agent]): `TPZCompElHDivDuplConnects{,Bound}.h` (duplicated connects for semi-hybridization), `TPZCompElConstFluxHybrid.h`, `TPZAlgebraicInterface.h`.
- **Solvers**: `TPZMatRedSolver.h` — Schur-complement/matrix-reduction driver, modes `EDefault/EDarcyHDiv/EDarcyH1Hybrid/EMHMSparse`; `TPZSparseMatRed.h`, `TPZDoubleMatRed.h`. Relationship to NeoPZ's own `TPZMatRed` → open question in [[matrix-and-solvers]].
- **Utils**: `TPZKernelHdivUtils.h` (print/solve/error helpers used by most targets), `TPZKernelHdivHybridizer.h` (excluded from build), `Common.{h,cpp}` (enums, `LaplaceExact`, UNSW quadtree reader), vendored `JSON.hpp`, generated `divfree_config.h` (absolute `MESHDIR`).
- `deprecated/TPZHDivApproxSpaceCreator.h` — superseded 50KB space creator [agent].

## Why it matters for the NeoPZ assessment
- Shows the *extension surface* actually used by researchers: derive approx creators, add materials, wrap solvers → extensibility evidence for Phase 5.
- The semi-hybrid / duplicated-connects / kernel-H(div) line here mirrors in-library counterparts (`Mesh/TPZCompElHDivDuplConnects*`, `TPZCompElKernelHDiv*`) — migration path app→lib visible in git history [inference].

## Related
[[approx-space-creators]] · [[TPZCompElHDiv]] · [[hybridization]] · [[flow-iter-elast]] · [[flow-dupl-connects]] · [[flow-mhm-hdivconstant]] · [[apps-overview]] (Session 2: five more downstream apps surveyed — this repo is one data point of six)
