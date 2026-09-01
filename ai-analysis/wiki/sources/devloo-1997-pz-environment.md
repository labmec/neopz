---
type: source
status: reviewed
authority: canonical
scope: paper
updated: 2026-07-02
confidence: high
tags:
  - neopz
  - architecture
---

# PZ: an object-oriented environment for scientific programming (Devloo 1997)

## Bibliographic reference
P.R.B. Devloo, "PZ: An object oriented environment for scientific programming", Computer Methods in Applied Mechanics and Engineering 150 (1997) 133–153. DOI 10.1016/s0045-7825(97)00097-2. (Cited as the reference publication in `README.md` [repo].)

## Why this source matters
The founding design paper of this exact library: states the original architectural intent (separation of geometry/topology/interpolation/algebra, extensibility via OO), against which today's structure can be honestly compared (what evolved vs. what ossified).

## Concepts covered
Geometric vs computational mesh split; element/side topology abstraction; hp interpolation; matrix abstraction; OO design for FEM.

## Claims extracted (relevant to the assessment)
- The gmesh/cmesh split and the "side" abstraction are *original, deliberate* design pillars — not accretions. → grounds [[TPZGeoMesh]], [[TPZCompMesh]], [[topology-module]].
- Extensibility via subclassing (materials, elements, matrices) is the intended extension mechanism → baseline for Phase 5 extensibility judgement.

## Applicability to NeoPZ
Architecture-intent evidence for `CODEBASE_ATLAS` §1-3 and the Phase 5 review ("essential vs accidental" complexity calls).

## Limits of applicability
1997 design predates C++11/17, multiphysics, HDiv/HCurl families, creators — do not treat as normative for those layers.

## Related wiki pages
[[TPZGeoMesh]] · [[TPZCompMesh]] · [[material-system]] · [[shape-functions]]

## Open questions
- Locate a copy for detail-level claims if Phase 5 needs direct quotes (paywalled; abstract-level use so far).
