---
type: concept
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - adaptivity
---

# hp-adaptivity

**Idea.** Combine local mesh refinement (h) with local polynomial-order enrichment (p); with the right strategy gives exponential convergence for elliptic problems with singularities.

**In NeoPZ.** The *mechanisms* are all in-library: per-connect orders (hierarchical [[shape-functions]]), `TPZInterpolatedElement::PRefine`, refinement patterns + directional refinement ([[refinement-hanging-nodes]], [[geometry-refinement-maps]]), side-order compatibility rules ([[element-families]]), and hanging-node restraints resolved at assembly ([[TPZConnect]]).

**The adaptive *drivers* live downstream — confirmed (Session 2).** The ErrorEstimation app closes the loop: estimator → per-element refinement indicator → `Hrefinement`/`HPrefinement` (h vs hp selection per element, `ErrorNaca.cpp:487-489`) → re-solve, iterated (13-step NACA studies); helpers `Tools::hAdaptivity`, `RandomRefinement` ([[app-error-estimation]]). wann runs estimator-driven adaptive `Divide` loops with an H1-vs-mixed comparison estimator ([[app-wann]]). In-library remains: mechanisms + `pzmganalysis`/gradient reconstruction; no in-tree marking strategy.

**Invariants (unchanged).** Min/max order rules on shared sides; order propagation after PRefine; p-enrichment × H(div) flavors (`fExtraInternalPOrder`); H(div)/H(curl) restraints under nonuniform refinement remain the thin test area (Phase 7 gap).

**Reference anchors.** Devloo–Oden hp work; Szabó–Babuška; Demkowicz.

Related: [[refinement-hanging-nodes]] · [[error-estimation-convergence]] · [[shape-functions]] · [[app-error-estimation]]
