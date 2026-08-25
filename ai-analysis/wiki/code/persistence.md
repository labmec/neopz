---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - neopz
  - persistence
---

# Save/ — persistence (object serialization)

## Responsibility
Binary save/restore of NeoPZ object graphs (meshes, matrices, analyses): every serializable class derives from `TPZSavable`, has a registered `ClassId`, and implements `Read/Write` against a `TPZStream`. `TPZPersistenceManager` orchestrates whole-graph writes with pointer fixup; `TPZChunkTranslator`/`TPZChunkInTranslation` provide *versioned* backward-compatible reads [agent; core files verified to exist: `Save/TPZSavable.h`, `TPZStream.h`, `TPZPersistenceManager.h`].

## Notables
- MD5-checksum stream (`Save/pzmd5stream.h`, `USING_OPENSSL`) — also used by cmake regression file comparison [agent].
- Coverage is thin: single unit test (`TestPersistence`, one `TPZFMatrix<CSTATE>` round-trip) despite the elaborate translator machinery; **no gmesh/cmesh round-trip test** [agent] → Phase 7 gap.
- `ClassId` (`Hash/TPZHash` [repo include in TPZMatBase.h:13]) hashes class names for stable ids.

## Related
[[TPZCompMesh]] · [[TPZGeoMesh]] · [[material-system]] · [[refinement-hanging-nodes]] (refpatterns are also persisted `.rpt` data)

## Open questions
- Is mesh save/restore actually used in current workflows (LabMEC restarts?) or mostly legacy?
- Do all *new* classes (approx creators, VTK generator) implement persistence, or is coverage decaying? → Phase 5/7.
