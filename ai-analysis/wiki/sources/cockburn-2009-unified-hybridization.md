---
type: source
status: reviewed
authority: canonical
scope: paper
updated: 2026-07-02
confidence: high
tags:
  - fem
  - hybridization
---

# Unified hybridization of DG, mixed and continuous Galerkin (Cockburn–Gopalakrishnan–Lazarov 2009)

## Bibliographic reference
B. Cockburn, J. Gopalakrishnan, R. Lazarov, "Unified hybridization of discontinuous Galerkin, mixed, and continuous Galerkin methods for second order elliptic problems", SIAM J. Numer. Anal. 47(2) (2009) 1319–1365. DOI 10.1137/070706616.

## Why this source matters
The standard modern frame for *what hybridization is* (break continuity, introduce skeleton multipliers, condense to a skeleton problem) — the vocabulary against which NeoPZ's `HybridizationType` taxonomy and wrap/interface/Lagrange geometry can be assessed as conventional-or-variant.

## Claims extracted
- Hybridized mixed methods: local solvers per element + a global equation only for the trace/multiplier unknown; the condensed system is SPD for symmetric elliptic problems.
- The multiplier is the trace of the primal variable (or numerical flux), and transmission conditions define the skeleton bilinear form.
- One frame covers mixed/DG/CG hybrid variants → variations (which continuity is broken, which multiplier space) are legitimate design choices, not errors.

## Applicability to NeoPZ
Baseline for [[hybridization]]: `EStandard` ≈ classic single-level hybridization; `EStandardSquared`/`ESemi` are Devloo-group extensions ([[avancini-2025-double-hybrid-elasticity]], [[carvalho-2024-semi-hybrid-stokes]]); SPD-ness expectation justifies LDLt/CG on condensed systems in the benchmark slices.

## Limits of applicability
Scalar elliptic focus; elasticity/Stokes variants need the specific papers. Static-condensation implementation details (grouping, connect levels) are NeoPZ-specific.

## Related wiki pages
[[hybridization]] · [[static-condensation]] · [[mixed-methods]] · [[flow-iter-elast]]

## Open questions
— none (background canon).
