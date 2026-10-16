---
type: code
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - mesh
  - dof
  - restraints
---

# TPZConnect — DOF bundles, restraints, Lagrange levels

Session-2 deep dive (agent trace, load-bearing lines re-verified [✓]). The DOF unit of the library: everything the user described as "shape function restraints" lives here.

## Anatomy (`Mesh/pzconnect.h`)
- `fSequenceNumber` (:32) — block number in `TPZCompMesh::fBlock`; −1 = unused. The single indirection connect→equations.
- `fNElConnected` (:35) — reference count of elements using this connect; **the** input to condensation decisions (`NElConnected()==1` ⇒ internal).
- Packed `union {fFlags; struct{fOrder, fNState, fLagrangeMultiplier, fIsCondensed}}` [✓ :42-60] — order, states/shape, **Lagrange level** (a small integer, not a bool: "level n multipliers need to be numbered after the multipliers of level n−1"), condensed flag, all in 4 bytes.
- `fNShape` (:60); `NDof() = fNShape*fNState` (:174-177).

## Dependency (restraint) machinery
- `fDependList` (:131): singly linked `TPZDependBase` list; `TPZDepend<TVar>` (:95-128) adds `fDepMatrix` — real and complex restraints coexist by templating.
- Semantics: dependent connect's DOFs = dep-matrix combination of the master's. Two interpretations at assembly (`Mesh/TPZElementMatrixT.cpp:327-343`): *shape* restraints (block-diagonal expansion by nstate) vs *algebraic* restraints (matrix as-is).
- **Creation** — the hanging-node path: `TPZInterpolatedElement::RestrainSideT<TVar>` (`Mesh/pzintel.cpp:872-1131`) builds a side mass matrix `M` (small side) + cross matrix `MSL` against the large neighbor's trace (pulled back through the geometric `SideTransform3`), solves `M⁻¹·MSL` by LU [✓ :964-967], and registers per-connect-pair blocks via `AddDependency` [✓ :1047-1050] (blocks with norm <1e-8 zeroed). So a hanging connect's restraint is literally the L2 projection of the coarse trace basis; the only geometric input is the side transform ([[geometry-refinement-maps]] §3). HCurl and HDiv override `RestrainSide` (`TPZCompElHCurl.cpp:622`, `pzelchdiv.cpp:1193-1196`).
- **Application** — per element in `TPZElementMatrixT::ApplyConstraints` (`TPZElementMatrixT.cpp:135-393`): `BuildConnectList` closes the set over masters; `BuildDependencyOrder` (`pzconnect.cpp:623-662`) fixes a topological order; then the congruence transform `Dᴴ·K·D` is applied through `TPZMatrixWindow::MultAdd` with **conjugate-transpose** flag [✓ `transp_a=2`, :359-366] — complex-correct. The shape-level analog `ExpandShape` (`pzconnect.cpp:387-444`) is real-only and hard-aborts on complex (:402-408).
- Downstream evidence that this is public API: wann glues 3D/2D/1D spaces with hand-built `AddDependency` calls ([[app-wann]]); GFEM keys enrichment functions to connect indices ([[app-gfem]]).

## Invariants (enforced by DebugStop, not types)
- **Condensed and dependent are mutually exclusive**: `CleanUpUnconnectedNodes` DebugStops on the combination [✓ `pzcmesh.cpp:618-624`, PZDEBUG only]; `TPZCondensedCompElT::Resequence` refuses to condense a dependent connect (`pzcondensedcompel.cpp:310-314`).
- `Reset()` DebugStops if `fDependList` non-null (`pzconnect.h:151-155`) — callers must `RemoveDepend()` first (e.g. `TPZMultiphysicsCompMesh.cpp:44-47`).

## Block structure & ordering (`Mesh/pzcmesh.cpp`)
- `fBlock` (target layout) vs `fSolutionBlock` (current layout of `fSolution`); `ExpandSolutionInternal` (:484-525) resequences and copies block-by-block, then `fSolutionBlock = fBlock`.
- `CleanUpUnconnectedNodes` (:588-798) renumbers into strata: **independent → condensed → dependent → freed** (:615-704).
- `SaddlePermute` (:2465-2723; an older version at :2372 is retired) orders connects globally by ascending Lagrange level so saddle-point factorizations eliminate in the right order — the mechanism behind "Lagrange levels order condensability" ([[static-condensation]], [[hybridization]]). Submesh variants subtract external connects.

Related: [[TPZCompMesh]] · [[refinement-hanging-nodes]] · [[condensation-groups-submeshes]] · [[multiphysics-composition]] · [[assembly]]
