---
type: source
status: reviewed
authority: canonical
scope: textbook
updated: 2026-07-02
confidence: high
tags:
  - fem
  - mixed-methods
---

# Mixed Finite Element Methods and Applications (Boffi–Brezzi–Fortin 2013)

## Bibliographic reference
D. Boffi, F. Brezzi, M. Fortin, *Mixed Finite Element Methods and Applications*, Springer Series in Computational Mathematics 44, 2013. DOI 10.1007/978-3-642-36519-5.

## Why this source matters
Canonical reference for every mixed-method expectation used in this assessment: inf-sup theory, H(div) conformity, RT/BDM families, Piola maps, hybridization basics.

## Claims extracted (the ones this assessment leans on)
- H(div) conformity ⇔ normal-trace continuity; H(curl) ⇔ tangential-trace continuity (Ch. 2).
- Contravariant Piola map preserves normal traces & the divergence pairing; covariant map preserves tangential traces — required on non-affine/curved maps for optimal rates (Ch. 2.1.3).
- Saddle-point well-posedness = ellipticity-on-kernel + inf-sup (Brezzi conditions) (Ch. 4-5).
- Discrete exactness/commuting diagrams underlie stability of compatible pairs (Ch. 2.5, FEEC-adjacent).
- Hybridization of mixed methods produces SPD condensed multiplier systems (Ch. 7 context; cf. [[cockburn-2009-unified-hybridization]]).

## Applicability to NeoPZ
Reference evidence for [[hdiv-space]], [[hcurl-space]], [[mixed-methods]], [[piola-transformations]], [[de-rham-complex]], [[hybridization]] — as *expected invariants*, not as prescriptions of NeoPZ's basis choices (NeoPZ uses its own hierarchical families, [[devloo-group-shape-construction]]).

## Limits of applicability
NeoPZ's families ≠ RT/BDM; dimension counts, DOF layouts and some stability proofs differ. Use for *properties* (traces, inf-sup, mapping requirements), not for family-specific facts.

## Related wiki pages
[[mixed-methods]] · [[hdiv-space]] · [[piola-transformations]] · [[de-rham-complex]] · [[hybridization]]

## Open questions
— none (background canon).
