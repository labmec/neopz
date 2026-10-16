---
type: source
status: draft
authority: strong
scope: paper
updated: 2026-07-02
confidence: medium
tags:
  - neopz
  - hdiv
---

# H(div) variants & divergence accuracy (Devloo group, 2018 + arXiv 1808.03625)

## Bibliographic reference (tightly related group)
1. P.R.B. Devloo et al., "Mixed finite element approximations based on 3-D hp-adaptive curved meshes with two types of H(div)-conforming spaces", Int. J. Numer. Methods Eng. (2018). DOI [10.1002/nme.5698](https://onlinelibrary.wiley.com/doi/10.1002/nme.5698).
2. Devloo group, "A remark concerning divergence accuracy order for H(div)-conforming finite element flux approximations", [arXiv:1808.03625](https://arxiv.org/pdf/1808.03625).

## Why this source matters
Documents (a) that NeoPZ deliberately ships **multiple H(div) space types** on 3D hp-adaptive *curved* meshes and (b) the group's own analysis of divergence accuracy orders — the published context for `HDivFamily` flavors (`EHDivStandard/EHDivConstant/...`) and the `fExtraInternalPOrder` (hdiv+/hdiv++) knob.

## Claims extracted
- Two H(div) space types with different internal enrichment yield different pressure/divergence accuracy on the same mesh [1, title/abstract level].
- Divergence order of H(div) flux approximations can differ from the flux order depending on family/geometry — a known, published subtlety, not an implementation accident [2].

## Applicability to NeoPZ
Direct context for [[hdiv-space]] flavors and the enriched-order options in [[approx-space-creators]]; supports classifying flavor-related surprises as *intentional variants*; feeds Phase 7 expectations (which orders should convergence tests show per family).

## Limits of applicability
Claims held at abstract level (paywalled/preprint skim pending); exact family definitions at the pin must come from code (Phase 4) — mark any mismatch as insufficient-evidence first.

## Related wiki pages
[[hdiv-space]] · [[approx-space-creators]] · [[piola-transformations]] · [[geometric-mappings]] · [[error-estimation-convergence]]

## Open questions
- Do the paper's "two types" correspond to `EHDivStandard` vs `EHDivConstant` at the pin, or to older families? What is `EHDivOptimized`?
