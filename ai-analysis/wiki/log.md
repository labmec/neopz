---
type: log
status: reviewed
updated: 2026-07-02
tags:
  - neopz
  - log
---

# Analysis log

Chronological record of major steps, discoveries, corrections, contradictions and open questions.
Newest entries at the bottom. Every entry lists evidence class: [repo] = repository evidence, [ref] = reference evidence, [agent] = subagent report (not yet independently re-verified), [run] = runtime observation.

## 2026-07-02 — Session 1

### Phase 0 — setup and pinning
- Engagement started. Full plan approved by user (autonomous execution, phases 1→9).
- [repo] Repo = clone of `git@github.com:labmec/neopz.git`; repo default branch is `main` (`refs/remotes/origin/HEAD`), review canon is `develop` per user instruction.
- **Pinned analysis commit: local `develop` @ `6ffd38b12`** (2026-06-12, "Making some methods virtual").
- [repo] Working tree = branch `SemiHybridElasticity` @ `4de234fae` = develop + 3 commits touching exactly 5 files: `Material/Elasticity/TPZHybridElasticity2D.cpp`, `Matrix/TPZSYSMPMatrix.h`, `Matrix/pzmatrix.h`, `Pre/TPZH1ApproxCreator.{h,cpp}`. Rule: cross-check these via `git show develop:<path>` before citing.
- [repo] Local develop is 2 commits behind origin/develop @ `852a5116c` (last fetch 2026-06-19); those 2 commits belong to the same SemiHybridElasticity line.
- [repo] Runtime artifacts: divfreebubbles links `find_package(NeoPZ)` → `../neopz_install` prefix; `libpz.dylib` rebuilt 2026-06-30 16:29, `iter_elast` binary 16:43 same day. Runtime traces = evidence about that installed build (≈ develop + 5-file delta), labeled [run], never conflated with the develop pin.
- User decisions: (1) pin as above, (2) autonomous cadence with per-phase checkpoints, (3) run existing artifacts only, cwd-in-tmp; never overwrite `divfreebubbles/build/targets/results_*.txt` (live data, last written 2026-07-02 22:13).

### Phase 1 — cartography (this session)
- Three read-only Explore agents swept: (a) NeoPZ module structure, (b) divfreebubbles app repo, (c) NeoPZ tests/CI/docs. Reports archived in conversation; load-bearing claims re-verified in main thread before entering wiki.
- Verified directly [repo]: CMake ≥3.14, `project(PZ)`, C++17, single `add_library(pz ...)` target (CMakeLists.txt:3,8,13,52); 24 `option(...)` flags; module `add_subdirectory` order (CMakeLists.txt:320-343); 5 GitHub workflow files; 4 `.h.h` template-body files (`Geom/pznoderep.h.h`, `Mesh/TPZGeoElement.h.h`, `Mesh/pzgeoelrefless.h.h`, `Mesh/tpzgeoelrefpattern.h.h`); `Publications/` = 3 H(div)-paper companion sources; `Material/needrefactor/` = 19 entries + `REAL/` subdir with 108 files.
- Verified directly [repo]: core headers read — `Material/TPZMatBase.h` (variadic mixin `TPZMatBase<TVar, Interfaces...>`), `Pre/TPZApproxCreator.h` (`HybridizationType {ENone,EStandard,EStandardSquared,ESemi}`, `ProblemType {ENone,EElastic,EDarcy,EStokes}`), `Mesh/pzcmesh.h`, `Mesh/pzgmesh.h`, `StrMatrix/TPZStructMatrix.h`, `Analysis/TPZAnalysis.h`, `Post/TPZVTKGenerator.h` (adapted from NGSolve, attributed), `Util/tpzautopointer.h` (atomic ref-count, Devloo), `README.md`, `Pre/pzcreateapproxspace.h` (function-pointer factory, Devloo 2009).
- **CORRECTION**: explorer agent (a) claimed `pzcreateapproxspace.h` lives in `Mesh/` and flagged a cross-module split as a "surprise". Wrong — it is `Pre/pzcreateapproxspace.h` [repo, verified by find]. The low-level factory and high-level creators are both in `Pre/`. Lesson: agent-only claims stay marked [agent] until re-verified.
- [repo] First-hand read of `divfreebubbles/targets/iter_elast.cpp` (full file): 2D hybrid elasticity benchmark, `TPZH1HybridApproxCreator` (app-side, `divfree/`), `HybridizationType::EStandardSquared`, orthogonalizing restraints + `HybridizeLowOrderFluxes` + `GroupAndCondenseElements`, `TPZMatRedSolver` `EDarcyH1Hybrid` vs direct MKL/MUMPS path; error/VTK blocks commented out; `results_*.txt` opened cwd-relative, two in truncate mode.
- Open question (Phase 5 candidate): `TPZCreateApproximationSpace` copy ctor / `operator=` appear to copy only `fp[8]` + some bools, not the `HDivFamily/H1Family/HCurlFamily` flavor flags or `fStyle` (Pre/pzcreateapproxspace.h:58-75, read cut at 75 — **verify full file before claiming**).
- Wiki bootstrapped: index, atlas, 14 code pages, 16 concept stubs (status: draft; concepts get filled in Phase 3).

### Corrections ledger
| # | Claim | Source | Correction | Status |
|---|-------|--------|------------|--------|
| C1 | `pzcreateapproxspace.h` in `Mesh/` | agent (a) | Actually `Pre/pzcreateapproxspace.h` | fixed in atlas + [[approx-space-creators]] |

### Phase 2 — shallow vertical slices (this session)
- Traced 5 slices; pages created: [[flow-iter-elast]], [[flow-dupl-connects]], [[flow-mhm-hdivconstant]], [[flow-dfreebubbles-1el]], [[flow-unit-test-hdiv-creator]]; `EXECUTION_FLOWS.md` v1 written.
- **OQ2 RESOLVED** [repo]: installed NeoPZ stamps `PZ_BRANCH="SemiHybridElasticity"`, `PZ_REVISION="852a5116c"` (= origin/develop tip at last fetch; develop + 2 delta commits). Caveat: config stamp regenerates at cmake-configure, dylib rebuilt Jun 30 — stamp is a lower bound; delta to worktree HEAD is at most the `MultiplyByScalar`-virtual commit. Runtime label: `[run @ 852a5116c(+)]`.
- Verified [repo]: `HDivFamily {EHDivStandard, EHDivConstant, EHDivKernel, EHDivOptimized}`, `H1Family {EH1Standard, EH1WidePrism}`, `HCurlFamily {EHCurlStandard, EHCurlNoGrads}` (Shape/TPZEnumApproxFamily.h:5-11). hdiv-space page corrected (EHDivOptimized added).
- Verified [repo]: `divfree/TPZMatRedSolver.h:15` `ProblemOrigin {EDarcyHDiv, EElasticityHDiv, EDarcyH1Hybrid, EElasticityH1Hybrid}`.
- New observations (app-side unless noted):
  - iter_elast uses `EDarcyH1Hybrid` mode for an elasticity problem (`EElasticityH1Hybrid` exists) → Phase 4 must read `TPZMatRedSolver::Solve` before classifying (naming debt vs mis-selection).
  - Benchmark drivers have error/VTK legs commented out and internally inconsistent analytic lambdas (3D sinh solution in 2D sweeps; forcing ≠ −Δu of exact; dead assignments overwriting exact solutions in main_1element) — benchmarks validate performance only; correctness legs live in dFreeBubbles1el + unit tests.
  - Both benchmarks run `TPZLinearAnalysis(..., RenumType::ENone)` — renumbering off; why? → Phase 6.
  - Same-name class `TPZMHMHDivApproxCreator` exists in NeoPZ `Pre/` AND `divfreebubbles/divfree/` — include-path-dependent selection; migration-in-progress pattern → Phase 5 risk.
  - `TestDeRham` mechanics verified first-hand [repo:TestDeRham.cpp:49-120]: rank/kernel + inclusion checks via SVD, dims 2&3, k=1..3, pairs H1→HCurl, HCurl→L2/HDiv, HDiv(Const)→L2.
  - `TestHDivApproxSpaceCreator` grid verified first-hand [repo:152-216]: 3 HDiv families × Darcy/Elastic × 4 mesh types × p{1,2} × extra-p{0,1} × hybridization {ENone,EStandard,ESemi} × RB × condensed × refined (+lib-side MHM creator at :519-524).

### Phase 3 — domain & reference bootstrap (this session)
- Network research performed (WebSearch/WebFetch OK). 9 source pages created (see index §Sources). Key identifications, code-driven:
  - `EStandardSquared` ↔ Avancini–Shauer–Oliveira–Devloo CMAME 2025 primal *double-hybrid* elasticity (H(div)–L² displacements/pressure, weak tangential continuity via shear-traction multiplier). Mapping marked hypothesis-level.
  - `ESemi` ↔ Carvalho–Devloo IJNME 2024 *semi-hybrid-mixed* (strong normal, weak tangential continuity; duplicated connects).
  - Shape/ H(div)/H(curl) construction ↔ De Siqueira–Devloo–Gomes JCAM 240 (2013): geometry-based vectors × hierarchical H1 scalars — matches the code-structure hypothesis in [[shape-functions]].
  - MHM ↔ Araya–Harder–Paredes–Valentin SINUM 51(6) 2013 (+ Devloo-group polygonal-elasticity variant).
  - H(div) flavor/divergence-order variance is published (IJNME 2018 two-space paper; arXiv:1808.03625) → flavor surprises default to "intentional variant" pending trace.
- `DOMAIN_PRIMER.md` written (10 sections incl. reviewer caution list). Concept pages hybridization/mhm/piola/hdiv-space updated with reference-evidence blocks.
- Taraschi–Correa arXiv:2601.21635 (2026) fetched: primal hybrid (u,m,p) elasticity analysis — locking-free coercivity + inf-sup; related-community analysis anchor for the elasticity-hybrid line.
- Note: `Publications/` in-tree companions (hdiv2d/3d 2015, hdivCurved JCAM) still unmatched to exact papers — acceptable; not on critical path.

### Phase 4 — deep review (this session)
- **iter_elast executed** [run @ 852a5116c(+), scratch cwd, p=1 iterative sweep, exit 0]: condensation 168,200→19,600 eqs (50²); reduced split "High Order Flux" 4,900 + "Linear Flux" 14,700; K00 factorized; **CG = 19 iterations at 50² and 400²** (mesh-independent). Output-column semantics confirmed; memory column = KiB-on-macOS bug → [[finding-rusage-memory-units]].
- **Piola question RESOLVED** (agent trace, key lines re-verified first-hand): contravariant Piola in split factorization — master directions from Topology (+`fSideOrient` signs +permutation gather), `(1/|detJ|)·J` applied in `TPZCompElHDiv::ComputeShape` (`pzelchdiv.cpp:1032-1033`), FAD branch for curved derivatives. |detJ|-vs-signed-detJ composition left as expert-validation item. → [[piola-transformations]] status: reviewed/high.
- **H1 hybrid creator traced at develop** (agent, key lines re-verified): EStandardSquared = literal double hybridization (2nd interface/Lagrange layer; 1st-level flux absorbed into volume condensation groups; only 2nd-level skeleton global). Lagrange-level enum {EL2,EFlux,EDistFlux,EDelayDec,EAvSol,EHybFlux}; multipliers Darcy {1,1,1,−1} / Elastic {−1,−1,−1,1} (verified `TPZApproxCreator.cpp:780-795`).
- **HDiv creator + ESemi traced** (agent): duplicated-connect mechanics precise (even=constant-flux connect rebound to wrap on sideOrient=−1 side; odd=higher-order stays continuous); ESemi requires EHDivConstant/EHDivOptimized; elastic path adds rotation space (weak symmetry); "singular K00" guard documents why condensed HDivConstant elasticity needs ESemi or RB spaces.
- **Delta-file diffs characterized**: develop lacks TPZHybridElasticity2D body-force RHS (fixed 2 commits later upstream) → [[finding-hybridelasticity2d-missing-rhs-at-pin]] (major, confirmed, fixed-upstream); pzmatrix delta = virtual MultiplyByScalar; TPZSYSMPMatrix delta = self-assignment guard (same missing guard still present in `TPZMatRed::CopyFrom` [repo pzmatred.h:66-79]).
- Findings created: hybridelasticity2d-missing-rhs (major), hdivconstant-fad-index (minor, confirmed inconsistency REAL vs FAD branch, verified :191 vs :304), approxspace-copy-drops-families (minor), approx-creator-hygiene cluster (minor: tautology, mis-scoped if — both verified via `git show develop:`; stubs; dead Backup/UnitaryLagrange machinery), rusage-memory-units (app), voronoi-null-ganalytic (app).
- OQ4 RESOLVED (confirmed app bug). Remaining in-flight: MatRedSolver trace (agent), OQ6/OQ7.

### Phase 4→7 runtime evidence batch (this session)
- MatRedSolver trace complete [agent, key lines verified]: split by Lagrange level {1}→K00 (mode-independent); matrix-free Schur CG + block-diag(ELU) preconditioner; **EDarcyH1Hybrid vs EElasticityH1Hybrid differ only in preconditioner block size** ⇒ iter_elast mislabel is consequential for benchmarks → [[finding-matred-solver-mode-mislabel]]. OQ6 resolved (enum removed at app HEAD `fbe9696`; 3-arg ctor gone) → [[finding-mhm-target-uncompilable]]. OQ7 resolved: `RenumType::ENone` because `TPZSparseMatRed::ReorderEquations` imposes its own Lagrange-level-contiguous ordering.
- **Direct-vs-iterative sweep** [run, p=1]: 50² parity (81/69 vs 78/75 ms); 400²: direct solve 7,851 ms vs Schur-CG 4,718 ms; assembly ≈ equal (~5.3 s); memory similar. Direct-solve growth superlinear; iterative ≈ linear (mesh-independent CG). → Phase 6.
- **ctest on existing build** (Release, 40 tests) [run]: 35 pass / **5 crash** (TestCondensedHangingNodes SIGTRAP; TestReduced, TestErrorAnalysis, TestHangingNode, TestSBFem bus errors). Failing binaries same-session as libpz ⇒ not stale-ABI. Upstream CI **green at pin `6ffd38b` and at `852a511`** [web api]. Attribution: local working-tree state (build = 852a5116c + then-uncommitted edits) or machine config → [[finding-local-test-crashes-workingtree]] (insufficient-evidence; not counted against the pin).

### Phase 5 — C++ & architecture review (this session)
- Two sweeps (ownership/threads; API/build) + first-hand verification of all H-severity claims. `CPP_TECHNICAL_REVIEW.md` written.
- Verified H-items: DebugStop throws bad_exception unconditionally (pzerror.cpp:15-28, guards commented out); `gTolerance` header-static per-TU (TPZTopologyUtils.h:23) ⇒ SetTolerance silent no-op; asymmetric gmesh/cmesh destructor protection (pzgmesh.cpp:97-103 vs pzcmesh.cpp:178-187); `Contribute` non-const (TPZMatSingleSpace.h:112-114) while shared across assembly threads; build-config triad (CMakeLists.txt:82-90 — RelWithDebInfo gets neither macro; 14 dead `#ifdef DEBUG`; no -Wall, Xcode warnings silenced).
- Notable positives recorded: OR/OT parallel assembly design (coloring + condvar ordering, per-thread scratch), material mixin discipline, 3757 `override`, atomic TPZAutoPointer, good header-doc coverage, extensibility walkthroughs quantified (material ~335 lines; solver 3 overrides; element family = the weak point, ~500 lines with copy-paste dispatch across 14 files).
- Layering: single `pz` target makes CMake order decorative; 6 backward include edges; Mesh⇄Pre SCC. Counts: 963 `new TPZ*`, TPZAutoPointer 758, unique_ptr 0.
- 5 new finding pages (see index). OR/OT question resolved in [[structural-matrices]].

### Phases 6–8 (this session)
- Phase 6: `FINDINGS_AND_ROADMAP.md` written — measured Schur-vs-direct crossover table, K00-dominance analysis, shared-memory-only ceiling (no MPI; BDDC dormant), perf-infra gaps, tradeoff-annotated suggestions, consolidated findings register, roadmap v1.
- Phase 7: `TESTING_AND_VALIDATION_REVIEW.md` written — inventory + proves/doesn't-prove analysis (invariant-strong, rate-weak), live ctest results integrated, 8-point validation strategy topped by a manufactured-solution **rate** matrix with nonzero body forces (directly motivated by the at-pin RHS bug).
- Phase 8 lint: 60 wiki pages, 61 link targets, **0 broken links, 0 orphans** (scripted check). Duplicate scan: hp-adaptivity vs refinement-hanging-nodes intentionally split & cross-linked. Stale-claim scan: DOMAIN_PRIMER §4 softened — EStandardSquared is a *structurally matching H1-primal variant* of the Avancini 2025 method, not a verbatim implementation (matches Phase 4 trace + source-page caveat). Flows all link both code+concepts; all source pages linked from concepts/findings; index updated through Phase 5 findings.

### Phase 9 — final report (this session)
- `NEOPZ_TECHNICAL_ASSESSMENT.md` written (§1–§11 per spec): confidence-labeled ([HC]/[MC]/[LC]), every §5 tension classified 5-way, §10 = 10 expert questions, §11 = evidence-boundary map. Engagement complete; wiki is the durable knowledge base; log ends here for session 1.

### Corrections ledger (cont.)
| # | Claim | Source | Correction | Status |
|---|-------|--------|------------|--------|
| C2 | `TPZMatRedSolver` has modes `EDefault`/`EMHMSparse` | agent (b) report + older driver code | Current header has neither; only 4 ProblemOrigin values | fixed in [[matrix-and-solvers]]; spawned OQ6 |

## 2026-07-06 — Session 2: library-breadth rebalance

### Mandate & method
- User feedback: Session 1 over-weighted divfreebubbles/HDiv; rebalance the docs toward the library itself (Topology→Geom→Shape stack, all space types incl. discontinuous, TPZConnect restraints, TPZMultiphysicsCompMesh, TPZCondensedCompEl/TPZElementGroup/TPZSubCompMesh, materials, matrices, StrMatrix), and survey the 5 most recent NeoPZ application projects in `~/GitHub` (each embeds its own near-identical neopz).
- Method: 9 parallel read-only explorers (5 downstream apps + 4 library subsystems); all load-bearing claims spot-verified first-hand before entering the wiki. Same pin (`develop @ 6ffd38b12`); downstream claims carry app-repo evidence class at each repo's HEAD — never conflated with pin evidence.
- Downstream recency determined by the app repos' own git logs [repo]: Iterative-Saddle_Point (2026-03-27), GFEM (2025-12-18), ErrorEstimation (2025-11-05), wann (2025-09-23), MixedElasticity (2025-05-19). NeoPZ_masterResearch excluded (build/examples area, not an app repo).

### New wiki content
- `apps/` (new section): [[apps-overview]] + [[app-iterative-saddle-point]], [[app-gfem]], [[app-error-estimation]], [[app-wann]], [[app-mixed-elasticity]]. Highlights: wann ≠ electromagnetics (it is wellbore-Darcy + ANN); GFEM subclasses `TPZCompElH1::ComputeShape` for enrichment; ErrorEstimation implements 4 estimator families + closed hp loops; MixedElasticity runs 3–7-field tensor mixed methods; Iterative-Saddle_Point drives Uzawa loops over cloned Pardiso matrices. Cross-cutting observations (9) recorded in the overview.
- `code/` (new deep-dive pages, agent-traced + line-verified): [[element-families]] (H1/HCurl/disc/interfaces — incl. resolutions: H1Family is creation-time & prism-only; HCurlFamily is a live runtime switch; covariant Piola confirmed in `TPZCompElHCurl::TransformShape`; L² pressure = broken-H1 (p>0) vs TPZCompElDisc (p=0); no `EDisconnected` enum exists), [[TPZConnect]] (packed flags; dependency = L2 projection of coarse trace; complex-correct `ApplyConstraints`; condensed∧dependent illegal; `SaddlePermute` mechanics — resolves the old fBlock/fSolutionBlock OQ), [[multiphysics-composition]] (AddElements ancestor-walk; AddConnects dependency re-offset; datavec order = mesh-vector order — resolves material-system OQ), [[condensation-groups-submeshes]] (Resequence rules; K11Reduced/UGlobal; SubCompMesh dual inheritance + rigid-body modes; ordering constraints), [[geometry-refinement-maps]] (TPZGeoEl policy stack; .rpt format; genealogy→RestrainSide bridge; TPZGeoElMapped exact-map inheritance; pyramid→6 pyr + 4 tet; gRefDBase deserialization coupling).
- `concepts/`: new [[discontinuous-l2-dg]]; upgraded to reviewed: [[h1-space]], [[hcurl-space]], [[sbfem]] (Hamiltonian + direct `dgeev_`, bypasses TPZEigenSolver stack), [[hp-adaptivity]] (drivers confirmed downstream); targeted updates: [[mixed-methods]] (weak symmetry resolved), [[refinement-hanging-nodes]], [[geometric-mappings]], [[static-condensation]], [[piola-transformations]] (HCurl covariant note), [[error-estimation-convergence]].
- `code/` updates: [[material-system]] (taxonomy verified vs CMakeLists; CSTATE electromagnetics + eigen mixins; BC framework; TPZMatWithMem; post-processing seams), [[matrix-and-solvers]] (eigen-solver stack verified), [[TPZCompMesh]], [[TPZCompElHDiv]], [[divfree-support-lib]].

### Deliverables revised
- `NEOPZ_TECHNICAL_ASSESSMENT.md`: Session-2 method note; exec summary gains "How it is actually used" (six-repo evidence); §3 domain-map additions; §6 extensibility re-weighed (material layer absorbs ~all downstream variation; GFEM shows family *modification* is cheap); §11 evidence classes updated.
- `DOMAIN_PRIMER.md` §10–12 (families as code structures; eigen/complex/SBFem; downstream usage); `ALGORITHM_NOTES.md` §9–11; `CODEBASE_ATLAS.md` §7 downstream landscape; `EXECUTION_FLOWS.md` scope note; `CPP_TECHNICAL_REVIEW.md` M8 corrected + §6 corroboration.

### Corrections ledger (cont.)
| # | Claim | Source | Correction | Status |
|---|-------|--------|------------|--------|
| C3 | `needrefactor/` "still compiles into `pz`" | Session-1 CPP sweep (M8) + final report §6 | No `add_subdirectory(needrefactor)`; `libpz.dylib` has zero needrefactor symbols [verified by nm]. Risk is header/name shadowing only; severity M→L | fixed in CPP review, atlas (§9 ledger), final report §6+§9 |
| C4 | Agent cited wann sources under `sources/` | app explorer | Actual dir is `src/` | corrected in [[app-wann]] |

### Notable cross-repo observations (for future sessions)
- [repo] `TPZHybridElasticity2D` exists app-side in ErrorEstimation (`ErrorEstimation/Material/TPZHybridElasticity2D.h:25`, edited 2025-11-05) *and* in the library (a 5-delta file of this working tree) — the SemiHybridElasticity line spans both repos; same-name migration risk live right now.
- [repo] Embedded neopz branches: develop (Iterative, MixedElasticity), `Australia25` (GFEM, ErrorEstimation), develop@f3b4000 (wann) — recent commits in each are app-motivated library fixes (HDiv side-orient dependency check for wann; eigenvector-loading/subdivision fixes for GFEM).
- MatRed-pattern reimplementations: divfreebubbles `TPZSparseMatRed`, GFEM `TPZSSpMatRedStructMatrix`+`TPZSparseMatRed`, Iterative's manual Schur loop → candidate library feature.
- None of the five apps uses H(curl)/CSTATE — that coverage rests on Electromagnetics materials + unit tests (older WGMAResearch line downstream).

### Open questions (running list)
- OQ1: Copy semantics of `TPZCreateApproximationSpace` (see above) — Phase 5.
- ~~OQ2~~ RESOLVED (see Phase 2 notes): installed build = `852a5116c(+)`.
- OQ3: `[!shouldfail]`-tagged unit tests (SVD, some skyline ops) — known gaps; what do they imply for LAPACK-path reliability? Phase 7.
- OQ4: divfreebubbles `voronoi_mixed_elas` possible null `gAnalytic` dereference [agent] — confirm by reading source in Phase 2/4.
- OQ5: docs claim "open-source" but no LICENSE file in tree [agent, verified absent by agent search] — re-verify and note in Phase 7 hygiene.
- OQ6: `main_MHM_HDivConstant.cpp:184` references `TPZMatRedSolver<STATE>::EMHMSparse` which no longer exists in the header [repo] — does the current app tree still build this target? Check git log of `TPZMatRedSolver.h` / try understanding when enum changed (read-only). App-side finding candidate.
- OQ7: Why do benchmark drivers disable renumbering (`RenumType::ENone`)? Interaction with MatRedSolver block structure? — Phase 4/6.
