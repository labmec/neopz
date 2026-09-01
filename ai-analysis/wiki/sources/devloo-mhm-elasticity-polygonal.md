---
type: source
status: draft
authority: strong
scope: paper
updated: 2026-07-02
confidence: medium
tags:
  - neopz
  - mhm
  - elasticity
---

# New H(div)-conforming MHM for elasticity on polygonal meshes (Devloo, Farias et al.)

## Bibliographic reference
P.R.B. Devloo, A.M. Farias, S.M. Gomes, et al., "New H(div)-conforming multiscale hybrid-mixed methods for the elasticity problem on polygonal meshes" ([Semantic Scholar](https://www.semanticscholar.org/paper/e7502db2ca2eea5cdea18dba1c5a0f11f38f7504); ESAIM:M2AN or similar venue — pin exact venue when needed).

## Why this source matters
The MHM+H(div)+polygonal-mesh+elasticity combination is precisely what the divfreebubbles MHM/voronoi drivers exercise ([[flow-mhm-hdivconstant]], `voronoi_mixed_elas`).

## Claims extracted
- Family of MHM methods for 2D linear elasticity on general polygonal meshes; approximate displacement and **stress divergence super-convergent in L²** [abstract via search].
- Local problems use H(div)-conforming stress spaces; rigid-body modes per subdomain are the coarse unknowns (consistent with `IsRigidBodySpaces`).

## Applicability to NeoPZ
Reference semantics for polygonal partitions (quadtree/Voronoi imports), scaling-center triangulation of polygons, and the role of `EHDivConstant` in making subdomain condensation exact. Feeds invariants for [[mhm]] and expected convergence for Phase 7 validation gaps.

## Limits of applicability
Paper-level; the at-pin code may implement a later/earlier variant. Venue/year still to pin down (marked draft).

## Related wiki pages
[[mhm]] · [[hdiv-space]] · [[flow-mhm-hdivconstant]] · [[mixed-methods]]

## Open questions
- Exact bibliographic detail + whether the divergence super-convergence is tested anywhere in-tree (Phase 7).
