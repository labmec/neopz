---
type: app-survey
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: "app development @ d6c0496 (2025-11-05); embedded neopz Australia25 @ 85cafdd8c (2025-11-23)"
tags:
  - neopz
  - downstream
  - error-estimation
  - adaptivity
  - sbfem
---

# App survey: ErrorEstimation (~/GitHub/ErrorEstimationResearch)

> Downstream-usage evidence (Session 2). Claims [agent], load-bearing ones spot-verified [✓]. Citations refer to the app repo.

**What it is.** An **a-posteriori error-estimation and adaptivity research suite**: solve (mostly Darcy/Poisson, newly elasticity), **reconstruct** an improved conforming solution in an auxiliary space, use the difference as an element-wise estimator, and drive h/hp-adaptive loops. Four estimator families: HDiv-mixed potential reconstruction, Hybrid-H1 reconstruction (H1 and HDiv flavors), MHM estimation, and patch-based partition-of-unity flux reconstruction. ~15 targets over an `ErrorEstimationLib`; active line = `ErrorNaca` (+2DSqrt/Lshaped/contrasting-permeability benchmarks), the HybridH1 reconstructions, and a new SBFem/elasticity thread (Nov 2025).

**Library surface exercised:**
- **The reconstruction layer** (the core novelty): `TPZHDivErrorEstimator` owns a `TPZMultiphysicsCompMesh fPostProcMesh` and a `TPZHybridizeHDiv` [✓ `ErrorEstimation/TPZHDivErrorEstimator.h:30,72`]; `PotentialReconstruction()` clones the solved pressure mesh, builds skeleton geo+comp elements, computes edge and nodal pressure averages, and produces a conforming potential — post-processable in **H(div)++ or H1** (`fPostProcesswithHDiv`). Sibling drivers `TPZHybridH1CreateH1Reconstruction`/`...HDivReconstruction` do the same from hybrid-H1 solutions. `TPZPostProcessError` [✓ `.h:85`] runs colored **patch solves keyed on partition-of-unity connects**. This whole layer is built from library primitives: mesh cloning, `TPZHybridizeHDiv`, `TPZElementGroup`+`TPZCondensedCompElT`, `TPZSubCompMesh` traversal, order-increase utilities (`IncreaseSideOrders` — HDiv++ as an estimation device).
- **Closed adaptive loops**: `ErrorNaca.cpp` runs a 13-step loop [✓ `Hrefinement`/`HPrefinement` :487-489] — estimator → per-element refinement indicator → `gel->Divide` / p-order maps → re-solve. `Tools::hAdaptivity`, `RandomRefinement`, MHM skeleton division (`DivideSkeletonElements`) round out the refinement API usage. This is the missing in-library "adaptive driver" of [[hp-adaptivity]] — it lives here, downstream.
- **Singularity-resolution geometry**: mesh styles `{ETraditional, ECollapsed, EQuarterPoint, ESBFem}` [✓ `ErrorNaca.cpp:241`] — **quarter-point elements** [✓ `CreateQuarterPointElements` :149,395], collapsed elements, and SBFem element groups as alternatives around the NACA trailing edge; custom exact map `TPZNacaProfile : TPZBlendNACA` [✓ `tpznacaprofile.h:31`] (the vendored NACA SpecialMap subclassed downstream). SBFem extensions: `TPZBuildSBFemHybrid : TPZBuildSBFem`, `TPZSBFemElementGroupPostProcess : TPZElementGroup`.
- **Materials as estimation logic**: ~12 custom materials subclass `TPZMixedPoisson`/`TPZDarcyFlow`/`TPZHybridDarcyFlow`/`TPZElasticity2D`, mixing in `TPZMatCombinedSpacesT`+`TPZMatErrorCombinedSpaces` so that `Contribute` computes estimator contributions — including a template mixin `TPZMixedErrorEstimate<MixedMat>` parameterized over the wrapped material. The material interface carries estimation logic, not just weak forms.
- **Cross-repo lineage flag**: an app-side `TPZHybridElasticity2D : TPZElasticity2D, TPZMatCombinedSpacesT, TPZMatErrorCombinedSpaces` [✓ `ErrorEstimation/Material/TPZHybridElasticity2D.h:25`, edited 2025-11-05] shares its name with the library's `Material/Elasticity/TPZHybridElasticity2D` — one of the 5 working-tree delta files of this engagement's pin. The SemiHybridElasticity work spans both repos; same-name-class migration risk (cf. `CPP_TECHNICAL_REVIEW.md` §6) applies here too.
- Spaces: H1, HDiv (incl. ++enrichment), discontinuous L2, multiphysics — manual `SetAllCreateFunctions*` composition plus a custom space factory `TPZCreateHybridH1Space` (pre-dates/parallels the library `TPZApproxCreator` layer). Solvers: direct only (`TPZSSpStructMatrix`, skyline, `TPZParFrontStructMatrix` frontal).

**What it teaches about the library.** (1) Error estimation and adaptivity are a *platform capability*: everything the estimators need (cloning, hybridization, condensation, submesh traversal, order manipulation) already exists as library primitives — but the drivers live downstream, confirming the [[hp-adaptivity]] page's hypothesis. (2) Mesh cloning + material swapping is a supported (if undocumented) workflow. (3) The geometry layer's exotic corners (quarter-point, collapsed, blend-NACA, SBFem) are exercised together in one adaptive driver. (4) A third independent app running its own MatRed-free reconstruction solves per patch shows small-dense local solves are a common downstream idiom.

Related: [[error-estimation-convergence]] · [[hp-adaptivity]] · [[hybridization]] · [[sbfem]] · [[condensation-groups-submeshes]] · [[geometry-refinement-maps]] · [[apps-overview]]
