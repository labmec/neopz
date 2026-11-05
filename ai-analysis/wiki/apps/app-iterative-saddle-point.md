---
type: app-survey
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: "app main @ 317b0c6 (2026-03-27); embedded neopz develop @ d366830a5 (2026-03-20)"
tags:
  - neopz
  - downstream
  - saddle-point
  - solvers
---

# App survey: Iterative-Saddle_Point (~/GitHub/IterativeResearch)

> Downstream-usage evidence (Session 2). Claims are [agent] from a read-only survey with the load-bearing ones spot-verified first-hand [✓]. Citations refer to the app repo, **not** the NeoPZ pin.

**What it is.** Research app for **iteratively solving saddle-point (mixed/incompressible) systems via a compressibility perturbation** (augmented-Lagrangian/Uzawa-style): add `−α` to the pressure diagonal so the system becomes SPD and Cholesky-factorizable once, then outer-iterate to recover the incompressible solution [✓ `sources/TPZMixedCompressibleDarcyFlow.cpp:233` — `ek(phrq+ip,phrq+ip) += -fAlpha*weight`]. Targets: `iterative-no-condense-darcy`, `iterative-condensed-darcy`, `iterative-condensed-stokes` (`targets/CMakeLists.txt:5-11`); shell scripts sweep mesh size × α ∈ [1e+1…1e−6].

**Library surface exercised** (beyond the divfreebubbles axis):
- **Spaces**: HDiv (`EHDivConstant`/`EHDivStandard`) flux × discontinuous L2 pressure × continuous H1 traction, combined in `TPZMultiphysicsCompMesh`. *Both* construction idioms: `TPZHDivApproxCreator` (Darcy, `iterative-condensed-darcy.cpp:448-472`) and fully **manual atomic-mesh assembly** for Stokes (4 atomic-mesh builders; `BuildMultiphysicsSpace`, `iterative-condensed-stokes.cpp:664-757`).
- **Manual hybridization of Stokes**: tangential-velocity/traction Lagrange spaces (`SetLagrangeMultiplier`), geometric-element surgery via `TPZGeoElBC` + `BuildConnectivity`, explicit `TPZMultiphysicsInterfaceElement` creation (`iterative-condensed-stokes.cpp:203-413`), then **manual `TPZElementGroup` + `TPZCondensedCompElT<STATE>`** condensation (`:784-901`) — the do-it-yourself counterpart of what `TPZApproxCreator` automates ([[hybridization]], [[condensation-groups-submeshes]]).
- **Low-level Matrix/Solvers usage** [✓]: hand-assembled divergence operator Bᵀ as CSR `TPZFYsmpMatrix::SetData(iBT,jBT,valBT)` (`iterative-condensed-darcy.cpp:318`, `...-stokes.cpp:1066`); matrix-free `MultAdd` residual updates; `Clone()` of the global `TPZSYsmpMatrixPardiso` + direct diagonal mutation; `SetDefPositive(true)`; one factorization reused across all outer iterations. The "iterative" method is a bespoke outer loop around a direct inner solve — no `TPZMatRed`, no library Krylov.
- **Materials as the extension layer**: 4 custom materials, all on `TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>, …>` [✓ `TPZMixedCompressibleDarcyFlow.h:25-26`; also `TPZHybridStokes`, `TPZHybridCompressibleStokes`, `TPZInterfaceStokes` (+`TPZLagrangeMultiplierBase`)]. No custom elements/struct-matrices/analyses. One material calls **`cblas_dgemm` directly inside `Contribute`** [✓ `TPZMixedCompressibleDarcyFlow.cpp:122,143,163`].
- **Geometry/config**: `TPZGeoMeshTools::CreateGeoMeshOnGrid` structured grids (hex/tet/prism); `TPZGmshReader` Poiseuille meshes; JSON-driven problem config (`nlohmann::json`). Renumbering `EMetis` (direct) vs `ENone` (iterative). No refinement, no curved maps.

**What it teaches about the library.** (1) The Matrix/Solvers layer is a *user-facing research surface*, not just plumbing — researchers hand-assemble coupling operators in CSR, clone and mutate Pardiso matrices, and exploit factorization reuse. (2) The manual Stokes path documents exactly what the creator layer abstracts away (and that vector-valued dim−1 multiplier spaces work). (3) Material-layer-only extension suffices for a whole solver-methodology study.

Related: [[matrix-and-solvers]] · [[hybridization]] · [[condensation-groups-submeshes]] · [[mixed-methods]] · [[apps-overview]]
