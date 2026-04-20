---
title: Printing Artifacts and Compensation
tags: [domain, lattice, am, compensation, woodward]
source: "Woodward & Fromen 2021, §4.3–4.4"
---

# Printing Artifacts and Compensation

The gap between **as-designed** and **as-printed** lattice dimensions, and the feedforward-compensation strategy Woodward & Fromen use to close that gap. This is an application concern the solver should support, not a first-wave requirement.

## The underlying phenomenon: overcure

In digital light processing (DLP) and continuous liquid interface production (CLIP), UV light polymerizes a photopolymer resin. Several effects cause **overcure** — more material solidifies than the digital design specifies:

- **Light scattering** inside the resin vat
- **Thermal gradients** during curing
- **Post-processing effects** — solvent exposure, resin bleed, thermal cure
- **Radial intensity variation** across the print plane

The result: printed windows (openings) are **smaller** than designed, and printed struts are **thicker** than designed. The effect is **non-uniform** — stronger at the part center, weaker at the periphery.

## Measured magnitude (Woodward §4.3)

For cubic lattices with different window sizes (printed in UMA 90 Cyan on Carbon M1):

| Designed window (mm) | Median A_w / A_design | Radial trend (slope × 1000) |
|---|---|---|
| 1.86 | 0.935 | 0.008 / mm  (~0.8%/mm) |
| 1.30 | 0.928 | 0.016 / mm  (~1.6%/mm) |
| 0.74 | 0.906 | 0.039 / mm  (~3.9%/mm) |

**Interpretation:** for the smallest windows, the central window area is ~85% of designed value; by 14 mm from center, it rises to ~52% (actually, the paper reports the theoretical central window area *falls* to ~52% at smallest L_c — the radial trend widens the windows as you move outward, but the baseline at center is dramatically shrunk). Smaller features suffer disproportionately.

Measured in UMA 90 Black, the radial trend is more consistent across configurations:
- A_w/A_design ≈ 0.012·R + 0.961 (1 part, centered)
- A_w/A_design ≈ 0.0125·R + 0.967 (3 parts, centered)
- A_w/A_design ≈ 0.0122·R + 0.942 (3 parts, off-center)

The slopes cluster around **12%/mm** radial growth; the intercepts (center value) cluster around ~95% of design.

## Resin-dependence

| Resin | Behavior |
|---|---|
| UMA 90 Cyan | As-printed window larger at periphery than at center; strong radial effect for small windows |
| UMA 90 Black | Darker pigment; starts closer to designed value at center; window area exceeds unity at periphery (slightly oversized) |
| EPU 40 | Elastomeric polyurethane; required a different print configuration (no flat base) to avoid delamination |
| SIL 30 | Silicone; most parts print successfully even below manufacturer-recommended feature sizes, but shows curing-related porosity deviations |

**Pigment matters.** Darker pigments absorb scattered UV and reduce overcure; lighter pigments let light penetrate further and cure larger regions.

## The compensation strategy (Woodward §4.4, Table 3)

A **two-stage feedforward** scheme, each stage a linear compensator inverting the observed trend.

### Stage 0: design input
- Target strut radius `r_0 = 0.183 mm` (in the paper's example)
- As-printed strut radius shows trend: `Δr_approx / r_0 = -0.015·R + 0.048`  (negative slope: shrinks away from center in this particular case)

### Stage 1: compensator C₁
- Apply compensator: `r_1(R) = 0.183 · 1/(1 - 0.015·R)`
- After printing: `Δr_approx / r_0 = 0.011·R + 0.037`  (slope now positive — overcorrected)

### Stage 2: compensator C₂
- Stack another compensator: `r_2(R) = 0.183 · 1/(1 - 0.015·R) · 1/(1 + 0.011·R)`
- After printing: `Δr_approx / r_0 = -0.002·R - 0.023`  (slope reduced ~85% compared to original, approximately flat)

The paper reports the residual slope after two iterations is **~85% smaller** than the uncompensated slope. This is a feedforward scheme — no closed-loop control, no in-situ feedback. Compensators are designed once per resin/geometry/hardware configuration and reused.

## When this strategy applies (and doesn't)

**Works well when:**
- Linear trend is a reasonable first-order approximation
- Print setup is reproducible (same printer, same resin, same timing, same post-processing)
- The number of compensation iterations is small (2 suffices here)

**Breaks down when:**
- The required compensation is larger than the printer's tolerance (fundamental resolution limit)
- Process variability exceeds the trend magnitude (noise > signal)
- Non-linear trends in higher-order terms matter
- A large number of iterations would be needed (polynomial compensation becomes impractical)
- Cross-configuration transfer (different materials, printers, orientations) — compensators do not transfer

## Project implications

1. **This is an optional pipeline stage.** Compensation is a post-design transformation of the strut-radius field, applied before mesh output. A function of the form `compensate(r_design: RadiusField, model: CompensationModel) -> RadiusField`.
2. **The compensation model is data-driven, not hardcoded.** Measured calibration data (per printer/resin/hardware) feeds into a fitted model (linear, for the paper's approach). Reading calibration data is an I/O concern — [[cli crate]] / [[lattice-gen crate]] boundary.
3. **Functional grading is the more general case.** Non-uniform strut radius as a function of position is how compensation is realized. The solver should support arbitrary `f(x,y,z) -> r` fields, not just per-part uniform `r`.
4. **Test oracles are tricky.** Compensation correctness requires measured-vs-designed data from a real printer — not a deterministic oracle. Bench-test compensator *construction* (given data, produce the inverse model); validate end-to-end against stored measurement datasets.
5. **Printability is separate from compensation.** The paper also produces Voronoi printability maps (Fig. 7) — regions of (r, porosity) space where parts fail outright. That is a **gating** function (is this design printable?), not a transformation.

## Data to carry forward

From the paper, for the Carbon M1 platform in UMA 90 Black:
- Manufacturer-recommended minimum feature size (XY): manufacturer-specified, printer-specific
- Minimum Z feature: tied to slice thickness, 100 µm in the paper
- Critical printable strut radius: ~0.065 mm (below this, lattices fail)
- Printable porosity bound: roughly > 0.25 (below a relative density of 0.75, lattices succeed; above it they close to solid)

These are starting values for a calibration table, not universal constants.

## Related

- [[Lattice Generation Pipeline]]
- [[Unit Cell Topologies]]
- [[Tests Should Make Programs Fail]]
- [[Reason About Correctness]]
