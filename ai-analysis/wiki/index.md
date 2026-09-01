---
type: index
status: reviewed
updated: 2026-07-06
tags:
  - neopz
  - index
---

# Wiki index

Analyzed commit: `develop @ 6ffd38b12`. Evidence rules and chronology: [[log]] (`log.md`).
Top-level deliverables live one directory up (`ai-analysis/*.md`); this wiki is the working knowledge graph feeding them.
**Session 2 (2026-07-06)** rebalanced the knowledge base toward the library itself: new deep-dive pages for the non-HDiv element families, TPZConnect, multiphysics composition, condensation/groups/submeshes, and the geometry/refinement/maps layer, plus a survey of the five most recently active downstream applications (`apps/`).

## Concepts (`concepts/`)
- [[h1-space]] — continuous conforming spaces; where NeoPZ realizes them. → shape-functions, hybridization
- [[hdiv-space]] — normal-trace-conforming vector spaces; NeoPZ families (Standard/Constant/Kernel). → TPZCompElHDiv, mixed-methods
- [[hcurl-space]] — tangential-trace-conforming spaces; "NoGrads" variant. → de-rham-complex
- [[de-rham-complex]] — exact-sequence property; NeoPZ's own SVD-based exactness tests. → TestDeRham
- [[mixed-methods]] — saddle-point two-field formulations; multiphysics machinery. → approx-space-creators
- [[hybridization]] — Lagrange-multiplier continuity; ENone/EStandard/EStandardSquared/ESemi taxonomy. → flow-iter-elast
- [[static-condensation]] — interior-DOF elimination; condensed elements/groups/submeshes; rigid-body spaces.
- [[mhm]] — multiscale hybrid-mixed method; controllers + creators; substructures. → flow-mhm-hdivconstant
- [[sbfem]] — scaled-boundary FEM sub-framework (breadth item).
- [[refinement-hanging-nodes]] — runtime refinement patterns (.rpt), connect dependency constraints.
- [[hp-adaptivity]] — per-connect orders + h-refinement; in-tree mechanisms vs downstream drivers.
- [[geometric-mappings]] — master→physical maps, blend + exact curved maps; axes convention.
- [[piola-transformations]] — vector-basis mapping requirement; NeoPZ's realization = key open question.
- [[quadrature]] — per-topology rules, long-double internals, order sufficiency questions.
- [[assembly]] — element→global pipeline across StrMatrix/materials/connects.
- [[error-estimation-convergence]] — SetExact/PostProcessError path; reproduction- vs rate-testing; downstream estimator suite.
- [[vtk-output]] — legacy .vtk writer model, subdivision of high-order fields.
- [[discontinuous-l2-dg]] — *(Session 2)* L² spaces (TPZCompElDisc vs broken-H1), interface elements, the compositional DG path.

## Code (`code/`)
- [[TPZGeoMesh]] — geometric mesh container, element/side/neighbor model, refinement genealogy.
- [[TPZCompMesh]] — computational mesh: connects, materials, solution block; Mesh↔Pre coupling.
- [[material-system]] — TPZMaterial/TPZMatBase variadic mixin design; physics dirs; needrefactor legacy layer.
- [[approx-space-creators]] — 2-layer space creation (factory + problem-level creators); hybridization data. (Correction C1 lives here.)
- [[TPZCompElHDiv]] — H(div)/H(curl)/kernel element families and conformity mechanics.
- [[shape-functions]] — static per-topology shape engine; families; hierarchical bases.
- [[topology-module]] — master-element sides, transforms, permutations.
- [[structural-matrices]] — assembly strategies, parallel schemes, equation filter.
- [[matrix-and-solvers]] — storage zoo + direct/iterative/eigen solvers; decomposition state.
- [[TPZAnalysis]] — solve orchestrator; renumbering; preconditioner factory.
- [[post-processing-vtk]] — graph-mesh legacy + TPZVTKGenerator (NGSolve-derived).
- [[persistence]] — TPZSavable/ClassId/chunk translators; thin test coverage.
- [[mesh-io-generators]] — gmsh reader (in-tree), grid generators, analytic solutions.
- [[TPZAutoPointer]] — house ref-counted pointer; dual raw/smart ownership conventions.
- [[divfree-support-lib]] — ../divfreebubbles application-side extensions (creators, Schur solvers, custom elements).

Session-2 deep-dive pages:
- [[element-families]] — H1 / H(curl) / discontinuous / interface elements; family enums resolved; dispatch table inventory.
- [[TPZConnect]] — DOF anatomy, dependency (restraint) machinery, Lagrange levels, SaddlePermute, renumbering strata.
- [[multiphysics-composition]] — TPZMultiphysicsCompMesh mechanics: AddElements/AddConnects, datavec stacking, solution transfer.
- [[condensation-groups-submeshes]] — TPZElementGroup / TPZCondensedCompEl / TPZSubCompMesh; ordering constraints; rigid-body modes.
- [[geometry-refinement-maps]] — TPZGeoEl hierarchy, uniform vs pattern refinement, .rpt format, genealogy→restraints bridge, special/blend maps, TPZGeoElMapped.

## Apps (`apps/`) — Session-2 downstream usage survey
- [[apps-overview]] — method, comparison table, nine cross-cutting observations.
- [[app-iterative-saddle-point]] — Uzawa/augmented-Lagrangian saddle-point iteration; low-level Matrix/Pardiso usage.
- [[app-gfem]] — GFEM fracture enrichment via TPZCompElH1 subclass; SBFem-derived singular modes.
- [[app-error-estimation]] — reconstruction-based estimators + closed hp-adaptive loops; quarter-point/SBFem/NACA geometry.
- [[app-wann]] — 3D/2D/1D coupled reservoir-wellbore Darcy; connect surgery; cylinder maps; directional refinement.
- [[app-mixed-elasticity]] — tensor-valued mixed elasticity, weak/strong symmetry, 3–7-field multiphysics, MHM controllers.

## Flows (`flows/`)
- [[flow-iter-elast]] — mandated slice: 2D hybrid-squared H1 elasticity benchmark; MatRedSolver vs direct; exercises the develop-delta files.
- [[flow-dupl-connects]] — semi-hybrid mixed H(div) Darcy benchmark (`ESemi`, `EHDivConstant`), lib-side creator.
- [[flow-mhm-hdivconstant]] — MHM on polygonal partition; substructures + condensation; app-source drift finding (OQ6).
- [[flow-dfreebubbles-1el]] — minimal manual mixed H(div); the active error + VTK legs.
- [[flow-unit-test-hdiv-creator]] — in-library validation: creator grid test + De Rham rank/kernel tests.

## Findings (`findings/`)
NeoPZ (library):
- [[finding-hybridelasticity2d-missing-rhs-at-pin]] — **major/confirmed**: body-force RHS absent from TPZHybridElasticity2D::Contribute at the pin (fixed upstream 2 commits later); no test would catch it.
- [[finding-hdivconstant-fad-index]] — minor/confirmed: FAD branch facet-count inconsistency in TPZShapeHDivConstant (latent: curved × HDivConstant × variable order).
- [[finding-approxspace-copy-drops-families]] — minor/confirmed: copy ops drop HDiv/H1/HCurl family flags.
- [[finding-approx-creator-hygiene]] — minor/confirmed cluster: tautological guard, mis-scoped if, stubs, dead interface machinery in the creator layer.

NeoPZ C++/architecture (Phase 5):
- [[finding-debugstop-throws-release]] — **major/confirmed**: assertion macro throws messageless bad_exception in all configs (3k call sites, 9 catches).
- [[finding-global-state-cluster]] — major/confirmed: per-TU gTolerance (SetTolerance is a silent no-op), gRefDBase mutability, legacy statics.
- [[finding-mesh-lifetime-ownership]] — major/confirmed pattern: gmesh-first destruction dangles cmesh; raw-owner copy hazards (TPZAnalysis::fSolver).
- [[finding-thread-shared-materials]] — major/math-risk: parallel assembly relies on unenforced material statelessness (non-const Contribute).
- [[finding-build-config-gaps]] — major/confirmed: RelWithDebInfo gets neither PZDEBUG nor PZNODEBUG; 14 dead #ifdef DEBUG blocks; warnings off.
- [[finding-local-test-crashes-workingtree]] — major/insufficient-evidence: 5/40 suites crash in local working-tree build (pin is CI-green).

divfreebubbles (application):
- [[finding-matred-solver-mode-mislabel]] — minor/confirmed, measurement-relevant: elasticity benchmarked with Darcy-shaped preconditioner.
- [[finding-rusage-memory-units]] — minor/confirmed: ru_maxrss bytes-vs-KB ⇒ 1024× memory overstatement on macOS benchmark tables.
- [[finding-mhm-target-uncompilable]] — minor/confirmed: registered target uses removed enum + ctor.
- [[finding-voronoi-null-ganalytic]] — minor/confirmed: null gAnalytic dereference on active path.

## Sources (`sources/`)
- [[devloo-1997-pz-environment]] — founding architecture paper (design intent baseline).
- [[devloo-group-shape-construction]] — 2009/2013 shape-construction papers = published spec of Shape/ H(div)/H(curl) bases.
- [[avancini-2025-double-hybrid-elasticity]] — CMAME 2025 primal double-hybrid elasticity (↔ `EStandardSquared`, iter_elast) + Taraschi–Correa 2026 analysis.
- [[carvalho-2024-semi-hybrid-stokes]] — IJNME 2024 semi-hybrid-mixed (↔ `ESemi`, duplicated connects).
- [[araya-2013-mhm]] — MHM origin (SINUM 2013) + scalable-implementation companion.
- [[devloo-mhm-elasticity-polygonal]] — MHM elasticity on polygonal meshes (Devloo group).
- [[boffi-brezzi-fortin-2013]] — canonical mixed-FEM text (traces, inf-sup, Piola, exactness).
- [[cockburn-2009-unified-hybridization]] — canonical hybridization frame.
- [[devloo-hdiv-variants-accuracy]] — H(div) flavors on curved/hp meshes + divergence-accuracy remark.

## Deliverables (../) — all phases complete 2026-07-02
- `CODEBASE_ATLAS.md` — Phase 1 ✔
- `EXECUTION_FLOWS.md` — Phases 2+4 ✔
- `DOMAIN_PRIMER.md` — Phase 3 ✔
- `ALGORITHM_NOTES.md` — Phase 4 ✔
- `CPP_TECHNICAL_REVIEW.md` — Phase 5 ✔
- `FINDINGS_AND_ROADMAP.md` — Phase 6 ✔ (+roadmap)
- `TESTING_AND_VALIDATION_REVIEW.md` — Phase 7 ✔
- `NEOPZ_TECHNICAL_ASSESSMENT.md` — Phase 9 final report ✔
