---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - util
  - memory
---

# TPZAutoPointer — reference-counted smart pointer

## Responsibility
NeoPZ's own shared-ownership smart pointer (predates std::shared_ptr; authored by P. Devloo). Wraps `T*` in a heap `TPZReference` holding the pointer + `std::atomic_int` counter (Util/tpzautopointer.h:36-60 [repo]).

## Facts [repo]
- Thread-safe counting via atomics (commented-out legacy mutex machinery still visible, lines 23-30).
- Non-intrusive (external control block), no weak-pointer concept, no custom deleters (as far as read; full API pending).
- Pervasive in APIs alongside raw pointers: e.g. `TPZStructMatrix` accepts both `TPZCompMesh*` (non-owning) and `TPZAutoPointer<TPZCompMesh>` (StrMatrix/TPZStructMatrix.h:55-58), `TPZCompMesh` holds both `fReference` raw and `fGMesh` auto (pzcmesh.h:49-54). This dual convention is a recurring ownership-ambiguity theme → Phase 5.

## Related
[[TPZCompMesh]] · [[structural-matrices]] · [[matrix-and-solvers]]

## Open questions
- Conversion semantics between templated types (`TPZAutoPointerDynamicCast`?), aliasing, and cycles (gmesh↔cmesh back-pointers are raw, so no cycle leak — intentional?).
- Interaction with persistence (auto pointers in `Read/Write`).
