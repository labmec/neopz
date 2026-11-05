# NeoPZ C++ Technical Review

**Phase 5 deliverable.** Objective C++/architecture review at `develop @ 6ffd38b12`. Sources: two targeted sweeps (ownership/lifetime/thread-safety; API/organization/build) with all high-severity claims re-verified first-hand, plus the Phase 1–4 traces. Every item: severity, evidence, why it matters, improvement + risk, and an **essential-vs-accidental** verdict. Finding pages under `wiki/findings/` carry full detail.

Legend: severity H/M/L; [✓] = verified first-hand; [agent] = sweep-cited (line-checked spot sample).

---

## 1. What is genuinely good (credit where due)

- **The domain architecture is coherent and intentional**: gmesh/cmesh split, sides/topology layer, static shape classes, connect-based conformity, materials as weak-form objects — all traceable to the founding design ([[devloo-1997-pz-environment]]) and still load-bearing. Complexity here is **essential**.
- **The material mixin design** (`TPZMatBase<TVar, Interfaces...>`, verified) is a modern, disciplined solution to the "N physics × M capabilities" matrix; boilerplate for a new material is ~335 lines with clear pure-virtual surface [agent, walkthrough].
- **Parallel assembly is real engineering**: OR = producer/consumer with single-writer scatter; OT = graph-coloring with atomic work index + condition-variable ordering and lock-free concurrent scatter (safe by coloring) [agent, lines cited]. Per-thread `TPZMaterialData` scratch is correctly isolated [✓ via pzinterpolationspace.cpp:480 pattern].
- **`override` discipline (3757 uses), `[[nodiscard]]` (150), `std::function` callbacks** (the factory table and SetExact are already modern) [agent counts].
- **Header-level Doxygen coverage on core public APIs is good (~70–90%)** [agent sampling]; docs pipeline (Doxygen+Sphinx→gh-pages) exists.
- **Persistence design** (ClassId + chunk translators) is elaborate and versioned — unusual for research codes.

## 2. High-severity findings (all accidental complexity)

| # | Finding | Evidence | Why it matters |
|---|---------|----------|----------------|
| H1 | **`DebugStop()` throws messageless `std::bad_exception` unconditionally — also in Release** (guard commented out) | [✓ Common/pzerror.cpp:15-28]; ~3029 call sites vs ~9 `catch` [agent] | The de-facto assertion mechanism gives production users a terminate with no diagnostic beyond a cerr line; guards can't be distinguished from fatal invariants; 51 `throw` vs 3029 DebugStop shows exceptions aren't the real strategy. → [[finding-debugstop-throws-release]] |
| H2 | **`pztopology::gTolerance` is a header-scope `static`** → one copy per TU; `SetTolerance()` mutates only TPZTopologyUtils.cpp's copy | [✓ Topology/TPZTopologyUtils.h:23]; readers in other TUs incl. header default args [agent] | Geometry point-location tolerance silently inconsistent across translation units; the commented-out singleton above it shows the intended fix was known. → [[finding-global-state-cluster]] |
| H3 | **GeoMesh→CompMesh dangling on the common path**: `~TPZCompMesh` clears the peer's back-pointer, `~TPZGeoMesh` does not; `SetReference(TPZGeoMesh*)` deliberately drops the owning autopointer | [✓ pzgmesh.cpp:97-103; pzcmesh.cpp:181-186; pzcmesh.h:772-773] | Destroy order becomes a correctness rule the compiler can't check; the co-ownership overload exists but is opt-in. → [[finding-mesh-lifetime-ownership]] |
| H4 | **Materials are shared mutable objects invoked concurrently through non-const `Contribute`** during OR/OT parallel assembly | [✓ TPZMatSingleSpace.h:112-114; agent: pzstrmatrixor.cpp:624, pzstrmatrixot.cpp:732, pzinterpolationspace.cpp:522] | Thread-safety rests on an unenforced "materials must be stateless in Contribute" convention; a single cached member = silent data race. Const-qualifying `Contribute/Solution` would make the existing invariant compiler-checked. → [[finding-thread-shared-materials]] |
| H5 | **Build-config macro gaps**: default `RelWithDebInfo` matches neither `$<CONFIG:RELEASE>` (PZNODEBUG) nor `$<CONFIG:DEBUG>` (PZDEBUG); 14 dead `#ifdef DEBUG` blocks (wrong macro); no `-Wall/-Wextra` anywhere; Xcode warnings explicitly silenced | [✓ CMakeLists.txt:84-90; ✓ TPZCompElHDivDuplConnects.cpp:62; ✓ grep counts] | The recommended default build silently disables both the debug checks *and* the release fast-paths; 14 real consistency checks never compile in any config; warnings that would catch the above are off. → [[finding-build-config-gaps]] |

## 3. Medium severity

- **M1 Layering is aspirational, not enforced** — single `pz` target; `add_subdirectory` order is decorative; 6 verified backward include edges (`Mesh→Pre` [✓ pzcmesh.h:18], `Geom→Mesh`, `Analysis→StrMatrix`, `Material(needrefactor,Plasticity)→Mesh`, `Refine→Mesh`, Solvers-before-Matrix order inversion) [agent]. Mesh⇄Pre is a genuine module-level cycle putting the two heaviest core headers in one recompilation SCC. *Essential?* No — the layer *concept* exists; enforcement is absent. Improvement: acyclic include lint in CI, forward-declare `TPZCreateApproximationSpace` in `pzcmesh.h`; cost low, risk low.
- **M2 Raw-owning members without copy control**: `TPZAnalysis::fSolver` raw owner + no deleted copy ⇒ double-delete on copy of any concrete analysis [agent, ctor/dtor lines]; contrast `fStructMatrix` already `TPZAutoPointer`. Improvement: delete copies or hold via autopointer; trivial, low risk.
- **M3 Ownership monoculture**: 963 `new TPZ*` sites; `TPZAutoPointer` (758, atomic refcount [✓ earlier]) vs `std::unique_ptr` **0**, `weak_ptr` 0 [agent counts]. Unique ownership is always raw `new`+manual `delete`; cycles (gmesh↔cmesh) can't be expressed safely. *Partly essential* (predates C++11; consistent house style), but new code keeps inheriting the hazard.
- **M4 `NConnectShapeF` per-family dispatch duplicated across 14 element .cpp files** with drift; formula centralized in Shape but the switch wrapper is copy-paste [agent]. Adding a family = copy 7 overrides + creator plumbing (~500 lines) rather than extend a hook. *Essential-ish* variability, *accidental* duplication.
- **M5 Global mutable state on active paths**: `gRefDBase` (refinement pattern DB, mutated by Initialize*/Insert/clear) blocks concurrent adaptivity and makes tests order-dependent; `gSinglePointMemory`; legacy per-material statics in `needrefactor/` (non-reentrant `gCurrentEq`) [agent]. → [[finding-global-state-cluster]]
- **M6 Installed `pz_config.h` bakes absolute build-machine paths** (`PZSOURCEDIR`, `PZ_REFPATTERN_DIR`) → non-relocatable installs; scalar-type macros are ABI-affecting globals [agent; consistent with observed install]. Improvement: runtime resolution/install-relative; medium effort.
- **M7 No symbol-visibility control on the shared lib** (everything exported; no PZ_API macro) [agent]. Larger ABI, slower links, zero encapsulation; standard `GenerateExportHeader` fix, medium effort (touching public headers).
- **M8 `Material/needrefactor/` legacy island — CORRECTED (Session 2, ledger C3)**: 19+108 files [✓ count], duplicate physics (two `TPZMixedDarcyFlow`s), old naming, the reentrancy statics. It does **not** compile into `pz` (no `add_subdirectory`; zero symbols in `libpz.dylib` [✓ nm]). Remaining risk is include-path shadowing of duplicate class names, plus out-of-library targets (SubStruct, PerfTests, Publications, TestCondensedElement) that include its headers. Severity drops M→L for the library proper; the retirement item in the roadmap becomes "delete/relocate headers + fix downstream includes" rather than a build-target carve-out.

## 4. Low severity / hygiene

- Four template-impl conventions coexisting (`.h.h` ×4 [✓], `_impl.h` ×5, inline-in-header, explicit-instantiation ×96 files) — standardize on `_impl.h`; L.
- Include-guard chaos (3 styles; generic `ANALYSISH` collision-prone guard) vs 7 `#pragma once`; L.
- Two naming eras `pz*`/`TPZ*` ~50/50 in Mesh/Matrix/Shape with filename↔class mismatch systemic in the old era (`pzcmesh.h`→`TPZCompMesh`); onboarding tax; L (rename churn risk high — recommend new-code-policy only).
- Self-assignment guards missing in `CopyFrom` family (`TPZMatRed` [✓]; patched post-pin in `TPZSYsmpMatrix`) ; L.
- `TPZCreateApproximationSpace` copy ops drop family flags ([[finding-approxspace-copy-drops-families]], verified) + `const void` setter signatures; L→M if a clone path is live.
- Approx-creator hygiene cluster (tautological guard, mis-scoped if, stubs, dead Backup/UnitaryLagrange machinery) — verified, [[finding-approx-creator-hygiene]]; L.
- Persistence coverage decay: elaborate machinery, one round-trip test (Phase 7 angle); L here.
- Docs miss the modern entry point: zero prose for `TPZ*ApproxCreator`; getting-started funnel leaves the tree (README → NeoPZExamples) [agent]; L→M for adoption.

## 5. Essential vs accidental — the honest split

**Essential (do not "simplify" without domain expertise):** sides/topology combinatorics; per-topology templates + explicit instantiation; hierarchical shape composition and its orientation protocol; multiphysics atomic-mesh + Lagrange-level condensation machinery; the variadic material mixins; swappable scalar types (`REAL/STATE`, FAD); refinement-pattern database concept.
**Accidental (fixable without touching the math):** everything in §2; the layering non-enforcement; ownership conventions; naming/guard/template-impl inconsistency; dead code (needrefactor, Backup paths, `#ifdef DEBUG`); build hygiene. Notably, the five H items are all in this category — the mathematical core is in better shape than the surrounding engineering.

## 6. Extensibility verdict (evidence-based)

Adding a **material**: well-designed path, moderate boilerplate (~335 lines), persistence registration only [agent walkthrough]. Adding a **solver**: small surface (3 overrides) [agent]. Adding an **element family**: the weak point — 7+ overrides plus copy-paste dispatch plus creator/factory plumbing across `Pre/pzcreateapproxspace` and creators (~500 lines, duplication-driven) [agent + Phase 4 traces]. The app repo (divfreebubbles) demonstrates the real extension workflow: derive creators, add materials, wrap solvers — workable, but same-name classes migrating app→lib (`TPZMHMHDivApproxCreator` in both) show the extension boundary is porous ([[divfree-support-lib]]).
**Session-2 corroboration across five more apps** ([[apps-overview]]): ~20 downstream material subclasses vs exactly one computational-element subclass (GFEM's `TPZCompElH1` override — which shows *modifying* a family via the `ComputeShape` seam is cheap even though *adding* one is not); two custom StrMatrix/matrix pairs and three independent MatRed reimplementations (the recurring gap); more same-name migrations observed (`TPZMixedElasticityND`, `TPZHybridElasticity2D`, `TPZSparseMatRed`) — the porous-boundary risk is systemic, not a divfreebubbles quirk.

## 7. Priority recommendations (cost/risk-annotated)

1. **H5 build-config triad** (RelWithDebInfo macro case, `DEBUG`→`PZDEBUG` rename, `-Wall -Wextra` in CI): hours of work, near-zero risk, immediately surfaces latent defects. Do first.
2. **H1 DebugStop**: introduce `TPZFatal(file,line,msg)`-style typed exception or abort-in-debug; mechanical replacement; low risk, high diagnostic payoff.
3. **H4 const-Contribute**: API-breaking but mechanical (derived overrides updated by compiler errors); pairs naturally with a thread-safety doc note. Medium churn, high invariant value.
4. **H2/H3/M2/M5**: small, local, high-value lifetime/global-state fixes.
5. **M1 layering lint + Mesh⇄Pre forward-decl**: cheap guardrail preventing further erosion.
6. Defer: renames (L), visibility macros (M7) until an ABI-break release window.
