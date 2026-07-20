# NeoPZ Findings & Roadmap

**Phase 6 deliverable (performance & scalability) + consolidated roadmap (finalized in Phase 9).** Evidence at `develop @ 6ffd38b12`; runtime data `[run @ 852a5116c(+)]` from the instrumented iter_elast sweeps (this engagement) — see `EXECUTION_FLOWS.md` and `wiki/findings/`.

---

## A. Performance & scalability review

### A1. Measured behavior (first-party evidence)

2D hybrid-squared elasticity, p=1, quad meshes n×n [run, same machine, Release library build]:

| n | eqs (condensed) | Schur t1/t2 (ms) | Direct assemble/solve (ms) | CG iters |
|---|---|---|---|---|
| 50 | 19,600 | 78 / 75 | 81 / 69 | 19 |
| 100 | 79,200 | 250 / 247 | 278 / 472 | — |
| 200 | 318,400 | 1,043 / 992 | 1,097 / 1,769 | — |
| 300 | 717,600 | 2,746 / 2,449 | 2,671 / 3,989 | — |
| 400 | 1,276,800 | 4,949 / 4,718 | 5,302 / 7,851 | 19 |

Readings: (i) **CG iteration count is mesh-independent (19 at both extremes)** — the K00-reduction acts as a robust preconditioner for this family; (ii) iterative t2 scales ≈ linearly in DOFs while direct solve grows superlinearly (69→7,851 ms ≈ n^1.14 in DOFs) — crossover before 100²; (iii) assembly ≈ linear and nearly identical in both paths (it is the same assembly); (iv) inside t2 at 400², K00 Cholesky dominates (3.4 s of 4.7 s) → further gains live in the K00 factorization, not CG.
Caveats: single machine; p=1 only; memory column of the app's tables is platform-distorted ([[finding-rusage-memory-units]]); the elasticity sweep ran with the Darcy-shaped preconditioner ([[finding-matred-solver-mode-mislabel]]) — true elastic-block numbers may differ.

### A2. Architecture-level performance properties

- **Assembly parallelism is real and layered**: OR (single-writer consumer) vs OT (coloring, lock-free scatter) vs OMP/TBB variants ([[structural-matrices]]). Parallel==serial is unit-tested. Tradeoff: OR serializes scatter (bottleneck at high thread counts); OT pays coloring precompute + condvar ordering. No benchmark in-tree currently quantifies the crossover (PerfTests stale).
- **Hot-path virtual dispatch is consciously split**: shape functions are static non-virtual per-topology classes (documented as a performance decision); materials *are* virtual per integration point — inherent to the extensible-material design; element `CalcStiff` virtual per element. This is the conventional FEM tradeoff; no evidence it dominates (assembly ≈ linear and matches the direct path).
- **Memory layout**: chunked vectors (`TPZAdmChunkVector`) give pointer stability at the cost of contiguity; `TPZFMatrix` column-major dense; `TPZManVector`/`TPZFNMatrix` small-buffer optimizations pervasive in element code — good. Sparse = CSR (sym/nonsym) with Pardiso/MUMPS backends.
- **Solver stack**: in-house skyline/Cholesky/LDLt for small-mid problems; MKL Pardiso / MUMPS for sparse direct; CG/GMRES + Jacobi/block/element preconditioners; eigen via Arnoldi/Krylov + LAPACK. The app-side `TPZSparseMatRed` shows the intended scalable pattern (factor small K00 blocks, iterate on the skeleton complement).
- **Scalability ceiling: shared-memory only.** No MPI anywhere; SubStruct/BDDC (Dohrmann) exists but is dormant (perf tests self-described as needing revision; not in CI). METIS optional for ordering only. For the group's problem sizes (10⁶–10⁷ DOFs observed) shared-memory + MUMPS/Pardiso is adequate; beyond that there is no distributed path.
- **Renumbering**: Sloan/Cuthill-McKee/Metis available via `RenumType`; benchmarks run `ENone` deliberately because `TPZSparseMatRed::ReorderEquations` imposes its own Lagrange-level ordering — coherent, but means the *direct* baseline in those same runs is unrenumbered (its fill could improve; the comparison slightly favors the iterative path) — worth one control run [suggested].

### A3. Build/compile performance

Single 35 MB dylib; 96 explicit-instantiation TUs; heaviest headers concentrated in Plasticity/needrefactor but core `pzcompel.h`/`pzcmesh.h` (800+ lines each) sit inside the Mesh⇄Pre include SCC → wide recompiles on core edits. No PCH, no unity builds, no ccache integration (CI uses ccache actions externally). Low-risk wins available; matters for iteration speed more than runtime.

### A4. Performance-infrastructure gaps

- `PerfTests/` stale (2021/2024 touches; own README: "in need of a revision"); not wired to CI; no regression tracking of assembly/solve times. The CDash config is legacy.
- No profiling hooks beyond `TPZSimpleTimer`/`TPZTimer`; LIKWID/PAPI options exist but ungated by any maintained target.
- The benchmark drivers that *do* exist (divfreebubbles) have the two measurement bugs found in Phase 4 (memory units, preconditioner mislabel) — i.e., current published-table pipelines carry avoidable noise.

### A5. Suggested improvements (each with tradeoff)

1. **Revive one thin perf CI job** (assemble+solve a fixed 2D/3D Darcy + elasticity case, OR vs OT vs serial, JSON output): catches regressions; cost = CI minutes; risk ≈ 0. Highest value/effort ratio.
2. **Fix the two benchmark-instrumentation bugs** (units, mode) before the next paper run; trivial.
3. **One control run with renumbering on** for the direct baseline; if fill improves materially, report both.
4. **K00 factorization**: at 400², 72% of t2 — try MUMPS BLR / Pardiso two-level or reuse symbolic factorization across the sweep (same sparsity each idiv? no — mesh changes; but across pOrders it repeats). Medium effort, needs measurement first.
5. **PCH/unity for Mesh+Material** and breaking the Mesh⇄Pre SCC: build-time win, low runtime risk.
6. **Do not** micro-optimize material virtual calls or replace chunk vectors on speculation — no evidence they dominate; measure first (premature-optimization risk flagged).

---

## B. Consolidated findings register (running; final ranking in Phase 9)

**Library (NeoPZ) — correctness/API:** [[finding-hybridelasticity2d-missing-rhs-at-pin]] (major, fixed upstream) · [[finding-hdivconstant-fad-index]] · [[finding-approxspace-copy-drops-families]] · [[finding-approx-creator-hygiene]] · [[finding-local-test-crashes-workingtree]] (working tree, not pin).
**Library — C++/architecture:** [[finding-debugstop-throws-release]] · [[finding-global-state-cluster]] · [[finding-mesh-lifetime-ownership]] · [[finding-thread-shared-materials]] · [[finding-build-config-gaps]].
**Application (divfreebubbles):** [[finding-matred-solver-mode-mislabel]] · [[finding-rusage-memory-units]] · [[finding-mhm-target-uncompilable]] · [[finding-voronoi-null-ganalytic]].

## C. Roadmap (v1 — completed in Phase 9 §9 of the final report)

- **Low-risk tidy**: build-config triad; DebugStop typed error; DEBUG→PZDEBUG; self-assignment guards; approx-creator hygiene cluster; benchmark instrumentation fixes; ~days total.
- **Medium structural**: const-Contribute; mesh-lifetime symmetrization + deleted copies; gTolerance inline-variable fix; layering lint + Mesh⇄Pre decycle; creators prose docs + in-tree tutorial; app-repo CI (build-all-targets); persistence round-trip test for meshes; ~weeks.
- **High-risk/deep**: element-family extension API (de-duplicate the 14-file dispatch); needrefactor retirement (Session-2 correction C3: already out of the `pz` build — remaining work is header deletion/relocation + fixing SubStruct/PerfTests/Publications/unit-test includes, so risk drops); visibility/ABI macros; unique_ptr migration for unique ownership; distributed-memory strategy decision (revive BDDC vs external solver coupling); MatRed-style block reduction as a first-class library citizen (Session 2: three independent downstream reimplementations); each needs expert review + regression nets.
