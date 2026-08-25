---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: possible mathematical risk
severity: major
evidence-commit: 6ffd38b12
tags:
  - neopz
  - cpp
  - thread-safety
  - assembly
---

# Parallel assembly calls non-const Contribute on shared material objects — statelessness is an unenforced convention

## Repository evidence
- One material instance per id, shared by all elements/threads (`Mesh/pzcmesh.h:65` `std::map<int,TPZMaterial*>`) [verified].
- `Contribute(const TPZMaterialDataT<TVar>&, REAL, ek, ef)` is **not const-qualified** — data is const, `this` is mutable (`Material/TPZMatSingleSpace.h:112-114`, verified first-hand).
- Both parallel strategies invoke it concurrently: OR producer/consumer (`pzstrmatrixor.cpp:624`; global scatter single-threaded), OT graph-coloring (`pzstrmatrixot.cpp:732`; concurrent non-atomic scatter, safe by coloring) via `pzinterpolationspace.cpp:522` [agent, structure verified in Phase 4/5 sweeps].
- Mitigations present: `TPZMaterialData` scratch is a per-call stack local (`pzinterpolationspace.cpp:480`) [agent]; parallel==serial equality is unit-tested for sample materials (`TestMultithreading` [agent Phase 1]).

## The risk (and why it's "possible mathematical risk", not confirmed bug)
Nothing observed misbehaves today; the invariant "materials keep no mutable state inside Contribute" evidently holds for the tested materials. But the invariant is **not compiler-enforced**: any material caching a member (e.g., a lazily-computed constitutive matrix, a stored last-point, `TPZMatWithMem` interactions) becomes a silent data race with non-deterministic assembly — the worst failure mode for a numerics library (wrong numbers, not crashes). Forcing functions (`std::function` members) invoked inside Contribute are likewise shared.

## What would resolve
Const-qualify `Contribute/ContributeBC/Solution` across the hierarchy (compiler then proves statelessness or flags violations; `mutable`+synchronized escape hatch for legitimate caches); or document + add a TSAN CI job exercising OR/OT on all shipped materials. Expert/maintainer input useful for `TPZMatWithMem` (per-point memory is *by design* mutable — how is it synchronized during parallel assembly? → open question).

## Related
[[assembly]] · [[structural-matrices]] · [[material-system]] · [[finding-global-state-cluster]] · CPP_TECHNICAL_REVIEW §2-H4
