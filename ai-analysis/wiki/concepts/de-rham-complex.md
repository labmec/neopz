---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
---

# De Rham complex (discrete exact sequence)

**Idea.** The sequence H1 →grad→ H(curl) →curl→ H(div) →div→ L² with image(left op) = kernel(right op) (on contractible domains). Discrete spaces that reproduce this exactness inherit stability/consistency for mixed problems (FEEC viewpoint).

**In NeoPZ (strong repo signal).** Dedicated test suite `UnitTest_PZ/TestDeRham/` compares "the dimension of the span of the differential operator of the left space against the kernel of the right space… rank(M_left) = ker(M_right)" via SVD, across pairs H1×HCurl, HCurl×HDiv, HDiv×L2 (needs LAPACK) [agent, header comment quoted]. Mesh-level checks: `TestMesh/TestHDiv.cpp` `CheckDRham(cel)` incl. under face permutations [agent]. Kernel-H(div) elements = explicit use of the complex (div-free fields as curls) → [[hdiv-space]], [[divfree-support-lib]].

**Why it matters for the review.** Exactness at the *basis* level is the library's own chosen correctness criterion for its space constructions — the assessment should trace exactly what property each test proves (rank equality ≠ full commuting-diagram property; clarify in Phase 3/4, mark what remains unproven, e.g. interpolation/commutativity, mesh-family uniformity).

**Reference anchors.** Arnold–Falk–Winther (FEEC); Demkowicz (exact sequences, projection-based interpolation); Devloo-group papers on compatible spaces.

Related: [[hdiv-space]] · [[hcurl-space]] · [[h1-space]] · [[mixed-methods]] · [[flow-unit-test-hdiv-creator]]
