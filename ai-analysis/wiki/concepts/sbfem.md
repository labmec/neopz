---
type: concept
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
---

# SBFem — Scaled Boundary FEM

**Idea.** Semi-analytical method: discretize only the boundary of star-shaped subdomains; the radial direction is handled analytically via an eigenvalue problem (good for singularities/unbounded domains).

**In NeoPZ (Session-2 deep dive, agent-traced).** Classes: `Mesh/TPZSBFemElementGroup` (the eigenproblem owner — a [[condensation-groups-submeshes|TPZElementGroup]] subclass), `TPZSBFemVolume` (per-volume element carrying the eigen-basis), HDiv/L2/multiphysics variants (`TPZSBFemVolumeHdiv/L2/Multiphysics`, `TPZSBFemMultiphysicsElGroup`, `TPZCompElHDivSBFem`); builders `Pre/TPZBuildSBFem(Multiphysics)` (partition + per-partition center nodes + matid translation map).
**Eigen construction** (`TPZSBFemElementGroup.cpp`): assemble coefficient matrices E0,E1,E2 (+mass M0) from skeleton elements (`ComputeMatrices` :87-181); form the 2n×2n block **Hamiltonian** `[[E0⁻¹E1ᵀ, −E0⁻¹],[E1E0⁻¹E1ᵀ−E2, −E1E0⁻¹]]` with a ±(dim−2)/2 diagonal shift (:703-744); solve the **non-symmetric** eigenproblem — directly via LAPACK `dgeev_` (:1845-1864, verified) or an alternative blaze-lib path (`CalcStiffBlaze` :189-276); select modes with Re(λ)<0 as the radial basis; eigenpairs are `std::complex<double>` throughout; bubble modes get a second eigenproblem. Note it **bypasses** the library's `TPZEigenSolver` stack and calls LAPACK/blaze directly.
**Validation**: `TestSBFem`/`TestSBFemHdiv` — the only in-tree suites asserting actual **convergence rates** (Darcy + 3D elasticity).

**Downstream evidence (Session 2).** SBFem is alive as a *tool*, not just a method: GFEM extracts crack-tip singular enrichment modes from `TPZSBFemElementGroup::EigenValues()/LoadEigenVector()` ([[app-gfem]]); ErrorEstimation subclasses the builders/groups (`TPZBuildSBFemHybrid`, `TPZSBFemElementGroupPostProcess`) and offers SBFem as one of four mesh styles around singularities ([[app-error-estimation]]).

**Reference anchors.** Song & Wolf (SBFem origin); Devloo-group SBFem papers (sbfempaper branch exists [repo branch list]).

Related: [[hdiv-space]] · [[matrix-and-solvers]] · [[condensation-groups-submeshes]] · [[app-gfem]] · [[app-error-estimation]]
