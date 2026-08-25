---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - neopz
  - post-processing
  - vtk
---

# Post/ — post-processing & VTK output

## Responsibility
Turn FE solutions into visualization files. Two generations coexist:
1. **Legacy graph-mesh family** [agent]: `Post/pzgraphmesh.h` + per-format writers `pzvtkmesh.h` (VTK), `pzdxmesh.h` (OpenDX), `pzmvmesh.h` (MVGraphs), `pzv3dmesh.h`; graph elements subdivide each computational element for plotting; driven by `TPZAnalysis::DefineGraphMesh/PostProcess` with variable-name tables.
2. **Modern `TPZVTKGenerator`** (Post/TPZVTKGenerator.h [repo]): authored 2022, explicitly "Adapted from NGSolve's vtkoutput.hpp" (attribution in header, lines 1-6); writes legacy-format `.vtk` with cell types mapped in `TPZVTK::CellType` (point/line/tri/quad/tet/pyr/prism/hex, lines 30-56); resolution via uniform master-element subdivision (`vtkRes`).

Also: `Post/TPZVTKGeoMesh.h` — dump geometric meshes (+partition/materials) to VTK for debugging; `Post/pzpostprocanalysis.h` — L² projection of solutions onto a post-processing mesh (used for plasticity/state vars) [agent]; `Post/pzgradientreconstruction.h`.

## What the output means
Field names are resolved through the material's `VariableIndex/NSolutionVariables/Solution` interface ([[material-system]]) — i.e. output correctness depends on each material's `Solution()` implementation, per variable. → validation angle for Phase 7.

## Related
[[TPZAnalysis]] · [[material-system]] · [[vtk-output]] · [[flow-dfreebubbles-1el]]

## Open questions
- Legacy `.vtk` (ASCII legacy format) only, or also XML `.vtu`? (TPZVTKGenerator appears legacy-format; confirm + note ParaView implications in Phase 4.)
- How high-order fields are represented (subdivision only? no VTK Lagrange cells?) — matters for judging visualization fidelity of p>1 solutions.
- Pyramid handling in TPZVTKGenerator vs its `MAX_SUBEL{10}` comment (pyramid refinement) — check.
