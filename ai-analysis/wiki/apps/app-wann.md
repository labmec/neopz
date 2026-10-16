---
type: app-survey
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: "app @ 2025-09-23 HEAD; embedded neopz develop @ f3b4000be (2025-09-18)"
tags:
  - neopz
  - downstream
  - darcy
  - multiscale-coupling
---

# App survey: wann (~/GitHub/WannResearch)

> Downstream-usage evidence (Session 2). Claims [agent], load-bearing ones spot-verified [✓]. Note: agent citations said `sources/`; actual dir is `src/` (corrected below where verified).

**What it is.** "Wellbore flow ANalysis using Neural Networks": a coupled **3D-reservoir + 1D-wellbore mixed-Darcy simulator** whose post-processed output (position → productivity index) trains a downstream PyTorch model. Targets: `wann3d` (JSON-configured coupled solve + VTK + ANN-training export), `wann3dRef` (adds an H1 companion mesh + H1-vs-mixed error estimator driving adaptive h-refinement), `oldwann3d`, `hanging-nodes-test`, `test-divide`.

**Library surface exercised:**
- **Dimensionally heterogeneous multiphysics**: 3D HDiv reservoir + 2D H1 "pressure skin" on the well cylinder + 1D HDiv wellbore, glued in one `TPZMultiphysicsCompMesh` with `TPZLagrangeMultiplierCS` + `TPZMultiphysicsInterfaceElement` — hybridization/condensation machinery explicitly *disabled* (`SetShouldCondense(false)` [✓ `targets/old-wann3d.cpp:152`]).
- **Direct TPZConnect surgery** [✓ `src/TPZWannApproxTools.cpp:242-244,546`]: `AllocateNewConnect`, `SetConnectIndex`, `SetLagrangeMultiplier`, and hand-built **`AddDependency` restraints** to weld spaces of different dimension along the well — downstream code programs the connect/dependency layer directly ([[TPZConnect]]).
- **Curved geometry a-posteriori**: mesh nodes projected onto the well cylinder, then `TPZChangeEl::ChangeToCylinder` (exact `TPZCylinderMap`) with neighbor `ChangeToGeoBlend` transition elements [✓ `src/TPZWannGeometryTools.cpp:157,187`] — the SpecialMaps/blend workflow on an imported mesh ([[geometry-refinement-maps]]).
- **Refinement breadth**: uniform (`TPZCheckGeom::UniformRefine`), **directional refinement** toward well heel/toe (`gRefDBase.InitializeRefPatterns` + `TPZRefPatternTools::RefineDirectional` [✓ `src/TPZWannGeometryTools.cpp:34-36`]), and estimator-driven adaptive `gel->Divide` loops — refinement patterns exercised as a live research tool, not a legacy feature.
- **Custom material**: `TPZNonlinearWell : TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>>` [✓ `src/TPZNonlinearWell.h:23`] — a nonlinear (friction/Forchheimer-type, |Q|^{3/4}) wellbore law hand-assembling a Newton tangent inside `Contribute`; compiled but not yet wired to a target (WIP).
- **Estimation**: hand-rolled H1-vs-mixed energy-norm comparison (two discretizations of the same problem as mutual error estimators), `std::thread`-parallel — a downstream pattern the [[app-error-estimation]] repo industrializes.
- Solvers: direct only (`TPZSSpStructMatrix` + `ELDLt`/`ECholesky`). Everything real-valued (`STATE`); no eigenproblems, no HCurl, no complex scalars [agent, exhaustive grep].

**What it teaches about the library.** (1) NeoPZ supports genuinely multi-dimensional coupled problems (3D/2D/1D in one mesh) — but at the cost of manual connect/dependency programming; there is no creator-level support for this pattern. (2) The connect-restraint (`AddDependency`) machinery doubles as a general-purpose coupling tool beyond hanging nodes. (3) Refinement patterns + cylinder/blend maps are actively used downstream. (4) The embedded neopz carries wann-motivated fixes on develop (boundary-element orientation, side-orient-aware dependency verification, 2025-09) — downstream needs drive library evolution.

Related: [[TPZConnect]] · [[geometry-refinement-maps]] · [[refinement-hanging-nodes]] · [[multiphysics-composition]] · [[apps-overview]]
