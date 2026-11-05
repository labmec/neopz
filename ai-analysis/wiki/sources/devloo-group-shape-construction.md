---
type: source
status: reviewed
authority: canonical
scope: paper
updated: 2026-07-02
confidence: high
tags:
  - neopz
  - shape-functions
  - hdiv
  - hcurl
---

# Devloo-group shape-function construction papers (2009–2015)

## Bibliographic reference (tightly related group)
1. P.R.B. Devloo, C.M.A.A. Bravo, E.C. Rylo, "Systematic and generic construction of shape functions for p-adaptive meshes of multidimensional finite elements", CMAME 198 (2009) 1716–1725. *(H1 hierarchical construction; per-topology, per-side.)*
2. D. De Siqueira, P.R.B. Devloo, S.M. Gomes, "A new procedure for the construction of hierarchical high order Hdiv and Hcurl finite element spaces", J. Comput. Appl. Math. 240 (2013) 204–214. ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0377042712003998))
3. Related follow-ups: hierarchical H(div) bases on curved 2D manifolds (ResearchGate 282859729); "Two-Dimensional H(div)-Conforming Finite Element Spaces with hp-Adaptivity" (Springer, 10.1007/978-3-319-39929-4_9).

## Why this source matters
These papers *are* the published specification of what `Shape/` implements: NeoPZ's H(div)/H(curl) bases are not RT/BDM/Nédélec textbook families but the Devloo-group construction.

## Claims extracted
- **Construction principle** [2013 abstract]: choose vector fields based on each element's geometry, multiply them by hierarchical H1 scalar functions → vector basis with continuous normal (H(div)) or tangential (H(curl)) interface components.
- Bases are hierarchical, per-side organized → variable order per connect (hp) is by-construction.
- [2009] H1 shapes built systematically from topology side-closures + orientation rules.

## Applicability to NeoPZ
Direct: [[shape-functions]] (`TPZShapeHDiv*`, `TPZShapeHCurl*` = "vectors × H1 scalars" — matches the code-structure hypothesis recorded there), [[TPZCompElHDiv]], [[hdiv-space]], [[hcurl-space]]. Conformity-by-construction explains why unit tests focus on trace continuity + permutation invariance rather than comparing against RT/BDM.

## Limits of applicability
2D-centric in [2]; 3D families, `EHDivConstant`/`EHDivKernel`/`EHDivOptimized` flavors and later refinements are separate developments — do not over-apply. Textbook Piola expectations may not map 1:1 onto this construction ([[piola-transformations]]).

## Related wiki pages
[[shape-functions]] · [[hdiv-space]] · [[hcurl-space]] · [[topology-module]] · [[piola-transformations]]

## Open questions
- Which paper (if any) documents the 3D H(div) family as implemented at the pin? (Candidate: Devloo et al. IJNME 2018, see [[devloo-hdiv-variants-accuracy]].)
