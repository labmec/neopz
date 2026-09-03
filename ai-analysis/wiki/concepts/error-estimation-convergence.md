---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - validation
---

# Error computation & convergence validation

**Idea.** Manufactured/analytic solutions → compute norms of (u_h − u) per element and globally; convergence rate vs h (or p) validates implementation order. A posteriori estimators drive adaptivity.

**In NeoPZ.** `TPZAnalysis::SetExact(...)` + `PostProcessError(errors,...)` computes per-material norms via `TPZMatError*::Errors` interfaces ([[material-system]]); `Pre/TPZAnalyticSolution.h` supplies exact solutions+forcings (`TElasticity2DAnalytic` etc. [repo, used in iter_elast]). Error vector convention: typically [0]=energy?, [1]=L2?, per material — **norm indices are material-specific; document per slice** (iter_elast prints error[0..4] [repo, commented block]). Convergence *rate* tests in-library: `TestSBFem` explicitly; most suites test exact-representation (polynomial reproduction) instead of rates [agent] → Phase 7 theme: reproduction-tests vs rate-tests.
A posteriori: gradient reconstruction (Post/), dedicated ErrorEstimation work **confirmed downstream (Session 2)**: the ErrorEstimation app implements four estimator families (HDiv potential reconstruction, hybrid-H1 H1/HDiv reconstructions, MHM, partition-of-unity patch solves) and closed h/hp-adaptive loops, entirely from library primitives — see [[app-error-estimation]]. wann adds a two-discretization (H1-vs-mixed) comparison estimator ([[app-wann]]).

**Invariants to check.** Error integration order sufficiency; `ElementSolution` redim before PostProcessError (iter_elast commented code shows the required incantation [repo:391-401]); parallel post-process equality (tested [agent]).

**Reference anchors.** Ainsworth–Oden (a posteriori); standard a priori theory (Ern–Guermond) for expected rates per space/order.

Related: [[TPZAnalysis]] · [[material-system]] · [[mixed-methods]] · [[flow-dfreebubbles-1el]]
