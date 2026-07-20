---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - neopz
  - approx-creators
  - cpp
---

# Approx-creator hygiene cluster: dead guards, mis-scoped if, stubs, dead machinery

Small verified defects in the space-creator layer. None mis-computes on today's live paths; each is a trap for the next maintainer. Grouped because a single review/PR could clear them.

## Verified items (first-hand [repo], develop pin)
1. **Tautological guard, dead DebugStop** — `develop:Pre/TPZH1ApproxCreator.cpp` (`CreateBoundaryHDivSpace`): `if (fHybridType != EStandard || fHybridType != EStandardSquared) {…} else DebugStop();` — the `||` makes the condition always true; intended `&&`. Read via `git show develop:` (≈:642).
2. **Mis-scoped if (missing braces)** — same file, `CreateAtomicMeshes` (≈:104-106): `if(HybridType() != ENone)` guards only `meshvec[0]=CreateBoundaryHDivSpace();` — the equally-indented `meshvec[1]=CreateL2Space();` runs unconditionally. Benign today (ENone DebugStops earlier), latent under-fill of `meshvec` for any future non-hybrid path.
3. **Stub method** — `develop:Pre/TPZH1ApproxCreator.h:99-102`: `CreateRotationSpace` is a `DebugStop()` stub; elastic rotational modes are folded into `CreateConstantSpace` state counts instead [agent, consistent with :224-254].
4. **Stored-but-unused config** — `fH1Fam` (`develop:Pre/TPZH1ApproxCreator.h:19`) is never branched on in the build path [agent grep].
5. **Dead interface machinery (HDiv creator)** — `AddInterfaceComputationalElementsBackup` (`Pre/TPZHDivApproxCreator.cpp:844-928`) and the purpose-built `TPZCompElUnitaryLagrange` element are uncalled; the live path uses the generic `TPZMultiphysicsInterfaceElement` (`:568,800`) [agent, caller-grep]. Header comment "Check with Jeferson if this variable is indeed necessary" (`Mesh/TPZCompElHDivDuplConnects.h:19`) marks known uncertainty around `fConnDuplicated` lifetime.
6. (Cross-ref) copy-op member omission in the layer-1 factory: [[finding-approxspace-copy-drops-families]].

## Why it matters
The creator layer is the library's *current* front door (unit tests + app repos build through it). Dead guards and unbraced ifs in exactly the methods being actively refactored (the develop→HEAD delta moves code between `GroupElements`/`CondenseElements` here) raise the odds that a future edit activates a latent path. Essential-vs-accidental: **accidental complexity** — none of these is FEM-driven.

## Suggested improvement (low risk)
`&&` fix + braces + delete-or-wire the Backup/UnitaryLagrange pair + remove or use `fH1Fam` + either implement or remove the rotation stub; all under existing unit-test cover (`TestH1ApproxSpaceCreator`, `TestHDivApproxSpaceCreator`).

## Open questions
- EStandardSquared: what keeps `EAvSol`-level connects globally coupled (the explicit `IncrementElConnected` runs only for EStandard, `develop:…:758-767`)? Works today per tests; mechanism unclear → expert Q (§10 of final report).
- Elastic multiplier asymmetry: only the right interface is reset to +1 (`develop:…:210-212`) vs Darcy resetting both — physical intent not determinable from code alone → expert Q.

## Related
[[approx-space-creators]] · [[hybridization]] · [[flow-iter-elast]] · [[flow-dupl-connects]]
