---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - application
  - build
---

# App repo: MHM_HDivConstant target no longer compiles (enum + constructor removed at HEAD)

## Repository evidence [agent trace + first-hand grep]
- `targets/CMakeLists.txt:7-8` still registers `MHM_HDivConstant` [agent].
- `targets/main_MHM_HDivConstant.cpp:184` (live, compiled branch inside a runtime `if`) constructs `TPZMatRedSolver<STATE>(*Analysis, matBCAll, TPZMatRedSolver<STATE>::EMHMSparse)`.
- `divfree/TPZMatRedSolver.h` at HEAD `fbe9696` has neither `EMHMSparse` (enum now `ProblemOrigin{EDarcyHDiv,EElasticityHDiv,EDarcyH1Hybrid,EElasticityH1Hybrid}` — verified first-hand :15) nor any 3-argument constructor (only default + `(TPZLinearAnalysis&, ProblemOrigin)`, .h:17-19) [agent].
- Enum history (`git log -L`): `EDefault/ESparse` (ccd781c) → +`EMHMSparse` (6c74e3b) → replaced entirely at HEAD `fbe9696` ("Extending the functionality of the iterative solver to the 2d elastic equations") [agent].
- The existing binary predates the change (Mar 27 vs solver rework); other registered targets remain valid (iter_elast, dupl_connects2 use the 2-arg ctor; others only include the header) [agent].

## Assessment
**Confirmed bug — application repo build drift** (a registered target that cannot build at HEAD). Severity minor per target, but it is the third independent instance of the same systemic theme (with [[finding-voronoi-null-ganalytic]] and the stale README): **no CI/compile gate on the app repo**, so drivers rot silently as the shared `divfree/` library evolves. Feeds the Phase 7 recommendation: a build-all-targets CI job (mirroring NeoPZ's own `compile_externalprojects.yml` pattern) would catch all three classes of drift.

## What would resolve
`cmake --build` of the target (out of authorized scope) or user confirmation; fix = port the call to the 2-arg ctor with an appropriate mode or gate the ESemi branch out.

## Related
[[flow-mhm-hdivconstant]] · [[divfree-support-lib]] · [[finding-voronoi-null-ganalytic]]
