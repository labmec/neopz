---
type: code
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12 (working-tree notes marked)
tags:
  - neopz
  - mesh
  - condensation
  - substructuring
---

# Element composition: TPZElementGroup, TPZCondensedCompEl, TPZSubCompMesh

Session-2 deep dive (agent trace, load-bearing lines re-verified [✓]). The three composition layers the user-visible concepts ([[static-condensation]], [[mhm]]) are built from.

## TPZElementGroup (`Mesh/pzelementgroup.{h,cpp}`)
A `TPZCompEl` owning a stack of elements + the union of their connects (`.h:24-25`). `AddElement` merges connects and **hides the element from the mesh** (nulls its `fElementVec` slot, `.cpp:72`); `Unwrap` restores. `CalcStiffInternal` (`.cpp:218-323`) sums member `ek/ef` into one dense block by connect map — pure summation, no elimination. `ReorderConnects` (`.cpp:116-145`) puts internal connects (`NElConnected()==1`, no dependency) first; `ExpandConnects` (`.cpp:431-448`) closes the connect list over dependency masters — called right before condensation wrapping. Downstream subclassing exists (`TPZSBFemElementGroupPostProcess` in [[app-error-estimation]]; `TPZSBFemElementGroup` in-library, [[sbfem]]).

## TPZCondensedCompEl (`Mesh/pzcondensedcompel.{h,cpp}`)
Decorator holding `fReferenceCompEl` + split connect lists; `NConnects()/ConnectIndex()` expose **only active connects** (`.h:52-65`) — that is the whole trick. `TPZCondensedCompElT<TVar>` owns a `TPZMatRed<TVar,TPZFMatrix<TVar>> fCondensed` (`.h:192`).
- `Resequence` (`.cpp:286-377`): any connect with `NElConnected()==1 && !HasDependency()` is force-condensed (:302-305); dependent connects can never be condensed [✓ DebugStop :310-314]; partition = condensed | active | dependent-but-kept.
- `CalcStiff` (`.cpp:384-693`): reference `CalcStiff` → `ApplyConstraints` **first** (:420-422, restraints before condensation) → reorder to `[condensed|active]` → `fCondensed.K11Reduced(K11,F1)` Schur complement → only the active block reaches the global matrix. `SetKeepMatrix(false)` frees internal blocks after condensation (:687-692) — the memory mode used by the HDiv creator. An experimental `USING_DGER` in-place LDLᵀ path (:468-579) DebugStops for non-double.
- `LoadSolution` (`.cpp:733-841`): gathers the active solution, `fCondensed.UGlobal(u1, elsol)` back-substitutes internals (:810), scatters into the mesh solution, recurses into the wrapped element — the implicit recovery step observed in every flow.

## TPZSubCompMesh (`Mesh/pzsubcmesh.{h,cpp}`, `Analysis/pzsmanal.cpp`)
Multiple-inherits **`TPZCompMesh` and `TPZCompEl`** [✓ `.h:32-35`] — a mesh that is an element of its father (the Schur-complement substructuring unit the user described).
- Bookkeeping: `fConnectIndex` (father indices of external connects), `fExternalLocIndex` (−1 = internal), father↔local maps (`.h:46-54`). `TransferElement`/`MakeAllInternal` demote connects to internal when `NElConnected()==1`, co-transferring dependency masters (`.cpp:484-571`).
- `SetAnalysis{Sparse,NonSymSparse,Skyline,Frontal}` install an internal `TPZSubMeshAnalysis` and — order matters — run `SaddlePermute()` **then** `PermuteExternalConnects()` before matrix creation (e.g. `.cpp:1385,1394`), giving the layout `[internal | external | constrained]`; the struct matrix's equation filter is restricted to `0..numinternal` (:1400). A GMRES-preconditioned option exists (:1408-1421).
- Schur exposure: father calls `CalcStiffInternal` (`.cpp:999-1225`) → `TPZSubMeshAnalysis::CondensedSolution` → `matred->K11Reduced(ek,ef)` [✓ `pzsmanal.cpp:152-158`] — condensed stiffness over external connects only; `AssembleInternal` builds the `TPZMatRed(numeq,numinternal)` and passes the **rigid-body-mode count** (`pzsmanal.cpp:108`). `LoadSolution` reverses via `UGlobal`.
- Floating substructures: `SetNumberRigidBodyModes(nrigid, lagrange)` (`.cpp:2124-2169`) allocates a tagged singular connect so `TPZMatRed` avoids factorizing a singular K00 — the library-level mechanism behind MHM's rigid-body coarse spaces.
- Drivers: `TPZCompMeshTools::PutinSubmeshes` (`Mesh/TPZCompMeshTools.cpp:462-529`) with `KeepOneLagrangian` (fixes the rigid mode); MHM controllers create submeshes directly (`Pre/TPZMHMeshControl.cpp:1413`).

## Composition pipeline & ordering constraints (fragile, DebugStop-enforced)
Canonical creator sequence (`TPZHDivApproxCreator::GroupAndCondenseElements`, `Pre/TPZHDivApproxCreator.cpp:658-710`): associate → `TPZElementGroup`s → **`ComputeNodElCon()` after grouping** (:691) → wrap in `TPZCondensedCompElT` with `SetKeepMatrix(false)` → `CleanUpUnconnectedNodes()`. Invariants: restraints resolve before condensation; condensed ∧ dependent is illegal (PZDEBUG DebugStop, `pzcmesh.cpp:618-624`); `SaddlePermute` before external-connect permutation before matrix creation; multiphysics `AddConnects` must re-offset dependency masters. *Working-tree note*: commit `852a5116c` (post-pin) split the H1 creator's `GroupElements`/`CondenseElements` and made the latter virtual — `TPZCompMeshTools::CondenseElements` keeps connects with `LagrangeMultiplier() >= LagrangeLevelNotCondensed` out by pre-incrementing `NElConnected` (`TPZCompMeshTools.cpp:607-618`), with a deliberately commented-out guard at :631.

Downstream: manual `TPZElementGroup`+`TPZCondensedCompElT` composition in [[app-iterative-saddle-point]] and [[app-mixed-elasticity]]; submesh-aware estimators in [[app-error-estimation]].

Related: [[static-condensation]] · [[TPZConnect]] · [[matrix-and-solvers]] (TPZMatRed) · [[mhm]] · [[multiphysics-composition]]
