---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: major
evidence-commit: 6ffd38b12
tags:
  - neopz
  - cpp
  - ownership
---

# Mesh lifetime: asymmetric dangling protection gmesh↔cmesh; raw-owner copy hazards

## Repository evidence (verified first-hand)
- `TPZCompMesh::~TPZCompMesh` (`Mesh/pzcmesh.cpp:178-187`) clears the geometric mesh's back-pointer if it points to `this` (`ref->ResetReference()`).
- `TPZGeoMesh::~TPZGeoMesh` (`Mesh/pzgmesh.cpp:97-103`) only runs `CleanUp()` — it **never nulls a surviving TPZCompMesh's `fReference`**, leaving the cmesh (and its whole `TPZGeoEl*` reference graph) dangling if the gmesh dies first.
- The co-ownership escape hatch exists but is opt-in and self-cancelling: `TPZCompMesh::SetReference(TPZGeoMesh*)` explicitly drops the owning `fGMesh` autopointer (`pzcmesh.h:772-773`) [agent, consistent with the dual-member design verified at pzcmesh.h:49-54].
- Ownership facts (verified/agent): gmesh deletes its `TPZGeoEl*`s (`pzgmesh.cpp:106-123`); cmesh deletes elements in a multi-pass dynamic_cast order (submesh→condensed→groups→interfaces→rest) *and* deletes materials (`pzcmesh.cpp:189-257`) — materials are mesh-owned raw pointers.
- Related copy hazard: `TPZAnalysis` owns raw `fSolver` (deleted in CleanUp) but declares no copy/move control ⇒ implicit copy double-deletes [agent: TPZAnalysis.h:73,111; TPZAnalysis.cpp:262-264]; contrast `fStructMatrix` (already `TPZAutoPointer`).

## Why it matters
Destruction order and no-copy rules are enforced only by convention. Downstream code (incl. divfreebubbles drivers, which routinely `new TPZGeoMesh` and hand raw pointers around, sometimes never deleting the analysis) relies on stack discipline to avoid UB. The multi-pass CleanUp shows real aggregation complexity handled manually — workable, but fragile against new element wrapper types.

## Classification
**Confirmed hazardous pattern (accidental)** — no observed crash traced to it in this engagement, but the asymmetry is objective and the double-delete path is reachable by a one-line user mistake.

## Suggested improvement / risk
(1) `~TPZGeoMesh` nulls the peer like the cmesh side already does — 3 lines, low risk. (2) Delete `TPZAnalysis` copy ops — trivial. (3) Longer term: prefer the autopointer `SetReference` overload / migrate unique ownership to `unique_ptr` gradually (963 raw `new TPZ*` sites, 0 unique_ptr today [agent counts]).

## Related
[[TPZGeoMesh]] · [[TPZCompMesh]] · [[TPZAnalysis]] · [[TPZAutoPointer]] · CPP_TECHNICAL_REVIEW §2-H3/§3-M2/M3
