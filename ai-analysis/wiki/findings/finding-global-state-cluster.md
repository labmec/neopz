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
  - thread-safety
---

# Global mutable state cluster: per-TU gTolerance (broken setter), gRefDBase, legacy statics

## 1. `pztopology::gTolerance` — header-scope `static` ⇒ per-TU copies (verified first-hand)
`Topology/TPZTopologyUtils.h:23`: `static REAL gTolerance = pow(10, …);` at namespace scope in a header — internal linkage, one copy per translation unit. `SetTolerance/GetTolerance` are defined in `TPZTopologyUtils.cpp:13-16` and touch **only that TU's copy** [agent]; other readers (`tpzprism.cpp:990`, header default arguments in `tpzline.h:126,129`) each see their own copy. **`SetTolerance()` is a de-facto no-op for most of the library.** A commented-out `Settings` singleton sits directly above (:14-21, visible in the verified read) — the intended fix was known. Fix: C++17 `inline` variable or the singleton; trivial, low risk.
Classification: **confirmed bug** (silent misbehavior of a public API).

## 2. `gRefDBase` — global mutable refinement-pattern database
`Refine/TPZRefPatternDataBase.h:101` `extern TPZRefPatternDataBase gRefDBase`, mutated at runtime by `InitializeUniformRefPattern/InitializeAllUniformRefPatterns/InsertRefPattern/ReadRefPatternDBase/clear` [agent]. Read on the refinement path ⇒ concurrent adaptivity + any mutation is unsynchronized; tests become order-dependent through shared state. Classification: **essential-ish concept (pattern cache), accidental realization** — should be injectable/owned or internally synchronized. Severity M.

## 3. Legacy statics (contained in `needrefactor/`)
Non-reentrant `TPZCoupledTransportDarcy::gCurrentEq` (set via `SetCurrentMaterial`), `TPZBurger::gStabilizationScheme`, biharmonic coefficient globals [agent]. Prevent two problem instances coexisting; already self-labeled as refactor debt by directory name. Severity L (quarantined but still compiled into `pz`).

Also: `gSinglePointMemory` (`Mesh/pzcompelwithmem.h:36`) flips element memory allocation globally [agent]; `gPrintLevel` (`Common/pzreal.h:170`).

## Why it matters
These are the classic blockers for (a) thread-safe adaptive refinement, (b) reproducible test isolation, (c) two meshes/problems with different tolerances/settings in one process. The tolerance item is the sharpest: a *public setter that silently does nothing* for most consumers.

## Related
[[topology-module]] · [[refinement-hanging-nodes]] · [[finding-thread-shared-materials]] · CPP_TECHNICAL_REVIEW §2-H2/§3-M5
