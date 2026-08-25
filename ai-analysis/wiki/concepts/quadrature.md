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

# Numerical integration (quadrature)

**Idea.** Element integrals evaluated by quadrature on the master element; rule order must cover integrand order (2p + geometry effects), with special handling for singular integrands.

**In NeoPZ.** `Integral/` module: abstract `tpzintpoints.h`, Gauss rules `tpzgaussrule.h`, per-topology rules (`pzquad.h`: `TPZInt{Quad,Triang,Cube,Tetra3D,Pyram3D,Prism3D}`), tensor & collapsed constructions, `TPZIntQuadQuarterPoint` (singular quarter-point rule), adaptive `adapt.h`. Doxygen note: computations use/return **long double** internally [agent] — precision-vs-cost choice worth noting. Element integration order = f(connect orders) with user override (`SetIntegrationRule`); materials can request order bumps (`IntegrationRuleOrder`) [pattern; verify sites Phase 4]. `TestIntegNum` validates polynomial exactness per topology, but several cases commented out (2D quad, 3D cube/tetra/pyramid) [agent] → Phase 7 note.

**Invariants to check.** Order sufficiency for: nonconstant Jacobians (curved els), material coefficients, `fExtraInternalPOrder`-enriched spaces, and error integration (`SetExact` order param, e.g. iter_elast passes `solOrder=4` [repo:iter_elast.cpp:87,275]).

**Reference anchors.** Standard texts (Ern–Guermond); rule tables (Dunavant etc.) as needed only.

Related: [[assembly]] · [[geometric-mappings]] · [[shape-functions]]
