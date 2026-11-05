---
type: app-survey
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: "surveys of five app repos at their 2025-2026 HEADs (see each page)"
tags:
  - neopz
  - downstream
  - usage-breadth
---

# Downstream usage survey — the five most recent NeoPZ applications

**Purpose (Session 2).** The Session-1 assessment leaned on one application (divfreebubbles: mixed/hybrid H(div) Darcy + hybrid H1 elasticity). To judge the *library* rather than one usage of it, this survey reads the five most recently active application repos under `~/GitHub`, each embedding its own NeoPZ copy (all near-develop; branches noted per page). Method: one read-only explorer per repo, load-bearing claims spot-verified first-hand [✓]. Evidence class: app-repo evidence — **not** statements about the analysis pin.

| App (last active)                        | Problem domain                                                       | Spaces used                                                    | Space-construction idiom                | Distinctive library surface                                                                                                                          |
| ---------------------------------------- | -------------------------------------------------------------------- | -------------------------------------------------------------- | --------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| [[app-iterative-saddle-point]] (2026-03) | Uzawa/augmented-Lagrangian iteration for mixed Darcy & Stokes        | HDiv(Constant/Standard) × L2 × H1 traction, multiphysics       | creator (Darcy) + fully manual (Stokes) | hand-built CSR operators, Pardiso matrix cloning/mutation, factorization reuse; manual Stokes hybridization + manual condensation                    |
| [[app-gfem]] (2025-12)                   | GFEM fracture mechanics (jump + singular enrichment)                 | H1 only + multiphysics composition                             | manual + custom `TPZCompElH1` subclass  | element-level shape enrichment via `ComputeShape` override; connect→function maps; SBFem eigenmodes as enrichment generator; custom StrMatrix/MatRed |
| [[app-error-estimation]] (2025-11)       | a-posteriori estimation + h/hp-adaptivity                            | H1, HDiv(++), L2, multiphysics                                 | manual + own space factory              | solution reconstruction (clone/swap/average), patch solves, closed adaptive loops, quarter-point/collapsed/SBFem/NACA-blend geometry                 |
| [[app-wann]] (2025-09)                   | 3D reservoir + 1D wellbore coupled Darcy → ANN training data         | HDiv × H1 skin × L2, multiphysics                              | creator + manual connect surgery        | 3D/2D/1D dimensionally heterogeneous coupling via `AddDependency`; cylinder/blend curved maps; directional refinement                                |
| [[app-mixed-elasticity]] (2025-05)       | stress-based (Hellinger–Reissner) elasticity, weak & strong symmetry | tensor-valued HDiv × L2 × rotation × RB constants (3–7 fields) | manual + creator + `TPZHybridizeHDiv`   | weak-symmetry multipliers; Johnson–Mercier macro-elements from custom refpatterns; MHM controllers with `TPZSubCompMesh`; frontal solver             |

## Cross-cutting observations (feed the Session-2 deliverable revisions)

1. **Every space family except H(curl) is heavily used downstream.** H1 is not a "lesser" space here — two of five apps are H1-centric (GFEM, much of ErrorEstimation). Discontinuous/L2 spaces appear in every mixed app (pressure, rotation, constants), including a downstream subclass of `TPZCompElDisc` (`TPZCompElDiscScaled`). H(curl)/electromagnetics and complex scalars did **not** appear in these five — that usage lives elsewhere (e.g. the older WGMAResearch line; noted as a boundary, not absence of capability).
2. **Manual space construction is alive and load-bearing.** All five apps use `SetAllCreateFunctions*` + atomic-mesh composition somewhere; three also use `TPZ*ApproxCreator`; two wrap their own space factories. The creator layer is a convenience roof, not the foundation — reviews must treat the low-level path as first-class API.
3. **The connect/dependency layer is a public research surface.** wann programs `AllocateNewConnect`/`AddDependency` directly for cross-dimensional coupling; GFEM keys enrichment functions to connect indices; MixedElasticity hand-stitches connect continuity for macro-elements. TPZConnect is not an internal detail ([[TPZConnect]]).
4. **MatRed-style block reduction is a recurring independent idiom** (divfreebubbles' `TPZSparseMatRed`, GFEM's `TPZSSpMatRedStructMatrix`+`TPZSparseMatRed`, Iterative's hand-built Schur loop) — three separate reimplementations of the same pattern argue for a first-class library citizen.
5. **Refinement patterns are used three ways**: adaptivity (ErrorEstimation), directional/geometric grading (wann, GFEM crack tips), and **space construction** (MixedElasticity's Johnson–Mercier splits). The runtime-pattern design earns its keep.
6. **The material mixin layer absorbs most extension needs**: ~20 downstream material classes across the five apps, typically `TPZMatBase<STATE, TPZMatCombinedSpacesT, TPZMatErrorCombinedSpaces, …>` or subclasses of concrete physics materials; several carry estimation or nonlinear-Newton logic. Only GFEM needed a computational-element subclass; nobody needed to touch Topology/Shape/Geom internals.
7. **App→lib migration is the growth mechanism — and a standing risk.** Observed same-name classes across app/lib: `TPZMixedElasticityND`, `TPZHybridElasticity2D` (the ErrorEstimation copy edited 2025-11 vs the library delta file in this very working tree), `TPZMHMHDivApproxCreator`, `TPZSparseMatRed`. Include paths silently arbitrate which is compiled (cf. `CPP_TECHNICAL_REVIEW.md` §6).
8. **Embedded NeoPZ copies confirm the user's "few changes" description**: all five sit on develop or short-lived branches (`Australia25` ×2) whose recent commits are app-motivated library fixes (HDiv side-orient dependency checks for wann; eigenvector-loading and element-subdivision fixes for GFEM) — downstream needs drive library evolution in small increments.
9. **Solver usage is direct-dominated** (Pardiso/LDLt/Cholesky, one frontal user); the only Krylov uses are bespoke (GFEM's preconditioned CG on a Schur complement; Iterative's outer Uzawa loop). Eigen-analysis appears via SBFem groups and matrix-level `SolveEigenProblem`, not `TPZEigenAnalysis`.

Related: [[divfree-support-lib]] (the sixth data point) · [[element-families]] · [[TPZConnect]] · [[multiphysics-composition]] · [[condensation-groups-submeshes]]
