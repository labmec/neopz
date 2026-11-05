# NeoPZ Testing & Validation Review

**Phase 7 deliverable.** Question: does the repository give confidence in (a) software behavior and (b) mathematical correctness? Evidence: full test-suite inventory (Phase 1 explorer, key files read first-hand in Phases 2/4), a live `ctest` run on the existing Release build [run], CI-status checks via the GitHub API [web], and the validation-relevant findings from Phases 4–5.

---

## 1. What exists (inventory)

- **Framework**: Catch2 v3.3.2 (FetchContent), 33 suites / ~56 cpp / ~180+ TEST_CASEs, custom event listener distinguishing new failures from `[!shouldfail]`-marked known failures. 40 ctest entries at the build used here.
- **Mathematically substantive suites** (the library's real crown jewels):
  - `TestDeRham` — discrete exact-sequence checks at basis level: rank(op(left)) = ker(right) + range-inclusion via SVD, dims 2/3, k=1..3, H1→HCurl→HDiv→L2 incl. HDivConst (read first-hand, TestDeRham.cpp:49-120).
  - `TestMesh/TestHDiv` — De Rham on real meshes **including under face/node permutations** (`drham_permute_check`), side-shape continuity, shape order, bilinear reproduction.
  - `TestTopology` — constant divergence/curl reproduction per topology + face-orientation data structures.
  - `TestHDivApproxSpaceCreator` — the creator pipeline grid: 3 HDiv families × Darcy/Elastic × 4 mesh types × p{1,2} × extra-p × hybridization {ENone,EStandard,ESemi} × RB × condensation × refinement (+MHM creator), asserting exact linear/constant representation, domain integrals, condensed equation counts (read first-hand, :152-216).
  - `TestHCurl` (traces, permutations, curls), `TestH1ApproxSpaceCreator`, `TestHDivCollapsed` (fracture elements), `TestSBFem` (actual **convergence-rate** tests), `TestSolverComparison` (MUMPS vs Pardiso cross-backend equality), `TestMultithreading`/`TestStruct` (parallel == serial assembly), `TestBlend`/`TestGeometry` (curved maps vs exact geometry), `TestHangingNode`/`TestCondensedSpace` (constraints), `TestIntegNum` (quadrature exactness), plus matrix/FAD/plumbing suites.
- **CI**: 5 GitHub Actions workflows. macOS job = always-on gate (build + ctest); Linux + MKL jobs conditional on a prebuilt `neopz-deps` image; consumer smoke-test workflow builds NeoPZExamples against an install; docs workflow. **CI is green at the analysis pin** (`6ffd38b` success; also `852a511`) [web].

## 2. What the tests prove — and what they don't

**Proven (strong, unusual for a research FEM code):**
1. Space-conformity invariants *by construction and by permutation* — the orientation protocol that most FEM bugs hide in is directly exercised.
2. Discrete De Rham exactness at rank level for all family pairs, 2D & 3D, k≤3.
3. Exact reproduction of constants/linears through the full creator pipeline across an enormous config grid (incl. hybridization + condensation + refinement).
4. Cross-backend solver equality (MUMPS↔Pardiso) and parallel↔serial assembly equality.
5. Element-level quadrature exactness (partially — several cases commented out).

**Not proven (the honest gaps):**
1. **Convergence *rates*** — only SBFem tests assert rates; the core H1/H(div)/hybrid paths assert exact representation of low-order solutions, which catches wiring errors but not order-of-accuracy regressions (e.g., a Piola-scaling subtlety on curved meshes would pass every current test). The at-pin missing-RHS bug ([[finding-hybridelasticity2d-missing-rhs-at-pin]]) is the concrete demonstration: **no test drives a material with nonzero body force against a manufactured solution.**
2. **Curved geometry × vector spaces** — TestBlend validates maps, De Rham suites run on affine meshes; the combination (H(div) on curved elements, the topic of an in-tree Publications companion!) is untested; the FAD-branch inconsistency ([[finding-hdivconstant-fad-index]]) lives exactly in that shadow.
3. **hp/hanging-node × vector families** — hanging-node suites exist (H1-centric); constrained H(div)/H(curl) under nonuniform refinement is not visibly covered.
4. **Persistence round-trips** — one matrix round-trip; no gmesh/cmesh save-restore test despite the elaborate versioned-translator machinery.
5. **Post-processing correctness** — VTK output has no golden-file or invariant test (only indirect smoke via a parallel-error test and PostProcessVTK calls in the creator suite).
6. **Known-failing debt is institutionalized**: `[!shouldfail]` on SVD and several skyline ops; commented-out quadrature cases — visible, at least, but unburned-down.

## 3. Operational health signals

- **Live run** [run @ working-tree build]: 35/40 pass in 16.6 s (fast suite!); **5 crash** (TestCondensedHangingNodes SIGTRAP; TestReduced/TestErrorAnalysis/TestHangingNode/TestSBFem bus errors). Upstream CI green at the same stamped revision ⇒ attribution points at the then-uncommitted local edits or machine config — [[finding-local-test-crashes-workingtree]]. Either way: **the crash cluster sits in condensation/constraints, the area under active refactor** — exactly where a rate/regression net is thinnest.
- **Granularity**: `catch_discover_tests` disabled (log4cxx interference) ⇒ ctest sees one test per suite; a single crashing TEST_CASE takes down the whole suite's reporting. No ctest LABELS.
- **No coverage tooling** anywhere; no test-result artifacts in CI; no Windows CI despite Windows build support; Linux/MKL legs skippable silently.
- **App repo (divfreebubbles)**: Catch2 tests exist but `BUILD_TESTS=OFF` default, several commented out; **no CI at all** — three independent drift bugs found (uncompilable target, null-deref driver, stale README) are the predictable consequence.
- Reproducibility: versioned persistence exists but untested for meshes; `pz_config.h` stamps git revision (good); refpattern data is runtime-loaded from an absolute path baked at configure time (fragile for relocated installs).
- Meta: no LICENSE/CONTRIBUTING in-tree (legal/hygiene, affects external validation contributions).

## 4. Recommended validation strategy (prioritized)

1. **Manufactured-solution rate matrix** (the single highest-value addition): one parametrized suite driving {H1, mixed HDiv (Standard/Constant), hybrid (EStandard/Squared/ESemi)} × {Darcy, elasticity **with nonzero body force**} × {affine, curved} meshes for 3 refinement levels, asserting L2/energy **rates** within tolerance. Directly nets: the RHS-class bug, Piola-on-curved risks, order regressions from family edits. Cost: builds entirely on `TPZAnalyticSolution` + creators already in-tree.
2. **Un-crash the working tree** and add the 5 crashing suites to a required pre-merge gate for the research branches (they run in seconds).
3. **Enable `catch_discover_tests`** (fix the log4cxx detection issue or gate logging in tests) so CI failures name the TEST_CASE; add LABELS for `ctest -L math|plumbing|solver`.
4. **App-repo CI**: build all registered targets + run the fast Catch2 tests on push — would have caught all three drift findings. Mirror NeoPZ's consumer-smoke workflow.
5. **Persistence round-trip test** for gmesh+cmesh (+ one refined mesh with dependencies) — protects the translator machinery that reproducibility claims rest on.
6. **Golden-file VTK test** (tiny mesh, fixed fields, byte-compare modulo float formatting) + a pyramid/prism cell-type check.
7. **Coverage + sanitizers in CI** (one Debug+ASan job; one TSan job running TestMultithreading and the OR/OT paths — pairs with [[finding-thread-shared-materials]]).
8. Burn down `[!shouldfail]`/commented tests or convert them to tracked issues; they currently encode known-unknowns invisibly.

## 5. Verdict

Software-behavior confidence: **moderate-to-good** (fast, broad, CI-gated on two platforms; weakened by granularity, coverage-blindness, and platform gaps). Mathematical-correctness confidence: **good at the invariant level, weak at the accuracy level** — the suite is unusually strong on structural/conformity invariants (permutation-proof De Rham checks are genuinely rare) and unusually thin on convergence rates and curved-geometry combinations; the two confirmed at-pin math-adjacent defects both live precisely in the untested shadows.
