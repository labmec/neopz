---
type: finding
status: draft
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - neopz
  - approximation
  - cpp
---

# TPZCreateApproximationSpace copy operations silently drop the space-family flags

## Repository evidence
`Pre/pzcreateapproxspace.h` [repo]: members include `HDivFamily fhdivfam`, `H1Family fh1fam`, `HCurlFamily fhcurlfam` (:44-46) and `MApproximationStyle fStyle` (:53). The copy constructor (:63-69) copies only `fp[8]` + `fCreateHybridMesh` + `fCreateLagrangeMultiplier` + `fCreateWithMemory`; `operator=` (:71-80) the same. **Family flags and style are not copied** → they reset to `DefaultFamily::…` (EHDivStandard/EH1Standard/EHCurlStandard) on the copy while the creation function pointers still reflect the source configuration.

## Failure scenario (mechanism-level)
Any path that copies a configured `TPZCompMesh::ApproxSpace()` (mesh clone/copy, or user code assigning one approx space from another) after `SetHDivFamily(EHDivConstant)`-style configuration yields a factory in a mixed state: fp[] creates H(div) elements, but elements consult the (now default) family flag → silently different space than intended (e.g. EHDivStandard instead of EHDivConstant), with no error.

## Impact qualifier (why severity is only *minor* pending Phase 5)
Impact depends on real call sites of the copy ops (TPZCompMesh copy ctor / Clone / persistence Read all candidates). If no live path copies a *family-configured* space, this stays a latent trap. → Phase 5: enumerate call sites; upgrade severity if a clone path is live in creators/tests.

## Adjacent hygiene (same header)
`const void SetHDivFamily(...)` etc. — `const void` return type; duplicated const/non-const getters both returning `const&` (:98-108). Cosmetic, but signals missing review on this header.

## Suggested improvement
Default the copy operations (`= default`) — all members are copyable — or copy the missing members explicitly; add a unit test cloning a `EHDivConstant`-configured mesh.

## Related
[[approx-space-creators]] · [[TPZCompMesh]] · [[hdiv-space]]
