---
type: source
status: reviewed
authority: canonical
scope: paper
updated: 2026-07-02
confidence: high
tags:
  - fem
  - mhm
---

# Multiscale Hybrid-Mixed method (Araya–Harder–Paredes–Valentin 2013)

## Bibliographic reference
R. Araya, C. Harder, D. Paredes, F. Valentin, "Multiscale Hybrid-Mixed Method", SIAM J. Numer. Anal. 51(6) (2013) 3505–3531. DOI [10.1137/120888223](https://epubs.siam.org/doi/abs/10.1137/120888223). (Origin; Harder–Valentin follow-ups 2015 extend.) Implementation-side companion in the same community: "On the Implementation of a Scalable Simulator for Multiscale Hybrid-Mixed Methods" ([arXiv:1703.10435](https://arxiv.org/pdf/1703.10435)).

## Why this source matters
Defines the MHM method NeoPZ's `TPZMHM*` controllers/creators implement (Devloo's group collaborates directly with the Valentin school).

## Claims extracted
- MHM relaxes continuity of the primal variable on a coarse skeleton via Lagrange multipliers while keeping **strong continuity of the normal flux component**; basis functions carry fine scales by solving **independent local problems** per coarse element (embarrassingly parallel); dual variable from postprocessing is **locally conservative** [abstracts].
- Well-posedness of local problems requires handling the constant/rigid-body kernel per subdomain (the multiplier system sees only skeleton unknowns + coarse constants).

## Applicability to NeoPZ
Grounds [[mhm]]: expected invariants for `TPZMHMHDivApproxCreator`/`PutinSubstructures`/`CondenseElements` — subdomain local solves = `TPZSubCompMesh` internal condensation; coarse constants = `IsRigidBodySpaces()=true` requirement observed in [[flow-mhm-hdivconstant]]; skeleton flux continuity = wrap/skeleton element structure.

## Limits of applicability
NeoPZ's realization (H(div) local solvers, EHDivConstant family, polygonal partitions) is a Devloo-group variant ([[devloo-mhm-elasticity-polygonal]]); the SINUM paper's exact spaces/estimates don't transfer verbatim.

## Related wiki pages
[[mhm]] · [[static-condensation]] · [[hybridization]] · [[flow-mhm-hdivconstant]]

## Open questions
- Which MHM generation (controllers vs creators) matches which paper generation — Phase 4/5 duplication check.
