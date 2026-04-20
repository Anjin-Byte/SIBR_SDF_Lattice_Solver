---
title: AM Platforms
tags: [domain, additive-manufacturing, printers, calibration]
---

# AM Platforms

Per-printer/resin calibration numbers and constraints. The [[Lattice Generation Pipeline]] and [[Printing Artifacts and Compensation]] strategies are platform-agnostic in *framework*, platform-specific in *numbers*. This note keeps the numbers in one place.

**Numbers below are starting points**, not committed tolerances. Each entry should be verified against the printer's current firmware, the specific resin, and the user's own print results before being relied on.

## Formlabs (SLA / LFS / MSLA) — primary target

The project's primary target is Formlabs printers via PreForm as the slicer. The relevant printers and their defining parameters:

| Printer | Technology | XY resolution | Z slice thickness (typical) | Build volume |
|---|---|---|---|---|
| Form 3 / 3+ | LFS (laser) | 25 µm spot | 25–300 µm (typical 50–100) | 145 × 145 × 185 mm |
| Form 3B | LFS (laser) | 25 µm spot | 25–100 µm | 145 × 145 × 185 mm |
| Form 3L | LFS (dual laser) | 80 µm spot | 50–300 µm | 335 × 200 × 300 mm |
| Form 4 | MSLA | 50 µm pixel | 25–100 µm | 200 × 125 × 210 mm |
| Form 4B | MSLA | 50 µm pixel | 25–100 µm | 200 × 125 × 210 mm |

**Design rules (Formlabs-recommended, conservative):**

- Minimum wall thickness: 0.4 mm (supported), 0.6 mm (unsupported)
- Minimum unsupported strut: ~0.5 mm (down to ~0.3 mm with tuning per resin)
- Minimum hole diameter: 0.5 mm
- Minimum embossed/engraved feature: 0.4 mm wide × 0.4 mm deep
- Minimum clearance between moving parts: 0.5 mm

**PreForm-specific constraints:**

- Accepts STL, OBJ, 3MF
- Tolerates ~10 M triangles comfortably; slows noticeably past ~30 M
- Rejects non-manifold or non-watertight input (hard requirement)
- Self-supporting lattice interiors — **disable auto-supports** for the lattice body; keep supports only on external features and base
- Supports can be manually painted to exclude lattice regions

**Derived extraction-grid resolution `h`:**
- For minimum strut `r ≈ 0.2 mm`, `h = r/3 ≈ 0.07 mm` is the floor
- Typically `h = 0.1 mm` suffices (matches slice thickness)
- Decimation target edge length: 0.1 mm (see [[Mesh Decimation]])

## Carbon (CLIP / DLP) — Woodward paper's platform

Included for reference because the [[Lattice Generation Pipeline]] and [[Printing Artifacts and Compensation]] data come from Carbon M1.

| Parameter | Value (Carbon M1) |
|---|---|
| Technology | CLIP (continuous liquid interface production), DLP variant |
| XY resolution | ~75 µm |
| Z slice thickness | 100 µm |
| Critical printable strut radius | ~0.065 mm (UMA 90 Black) |
| Build volume | 189 × 118 × 326 mm |

Notable differences from Formlabs:
- **~3× finer minimum feature** due to the CLIP process's reduced O₂-inhibition dead zone.
- **Different overcure behavior** — Woodward's 2-stage linear compensation coefficients are Carbon-specific. The framework transfers; the numbers don't.
- **Resin library differs** — UMA 90, EPU 40, SIL 30 (Carbon) vs. Clear, Tough 1500, Durable, Flexible 80A, etc. (Formlabs).

## Nexa3D, Anycubic, Prusa SL, etc. — future targets

All MSLA/SLA printers have the same shape of constraints as Formlabs. Adding a new printer means filling in a row of the Formlabs table plus:
- Slicer software (Nexa xPRO, ChituBox, PrusaSlicer) — each has its own mesh-quality tolerances
- Design-rule recommendations (minimum features)
- Resin-specific overcure calibration (gather empirically)

Document per-printer data in the same structure as above. No architectural changes needed.

## FDM (not a target, for contrast)

Fused deposition modeling (e.g., Prusa MK4, Bambu X1C) has a different constraint shape — extruded plastic, layer adhesion, nozzle width. Minimum features are typically 0.4 mm nozzle × 0.2 mm layer height, and overhangs require supports or `<45°` geometry. **The SDF lattice approach is applicable**, but the compensation model and design rules are entirely different. Not a first-wave target.

## The platform-agnostic invariants

Regardless of printer:

1. **Output must be watertight + 2-manifold + self-intersection-free** — universal slicer requirement. MC33 delivers this; see [[Isosurface Extraction Methods]].
2. **Output resolution should match printer resolution** — finer is wasted. Decimation target edge length ≈ slice thickness; see [[Mesh Decimation]].
3. **Self-supporting designs reduce post-processing** — see [[Lattice Generation Pipeline]] "openness" invariant.
4. **Overcure/overgrowth compensation is resin-and-printer-specific** — framework in [[Printing Artifacts and Compensation]], numbers here.

## The platform-specific variables

| Variable | Changes between platforms |
|---|---|
| Minimum printable feature size | ✓ |
| Overcure compensation coefficients | ✓ |
| Slicer software + its mesh-quality tolerances | ✓ |
| Supported file formats | ✓ (most accept STL + 3MF) |
| Resin library | ✓ |
| Auto-support behavior | ✓ |
| Build volume | ✓ |

Keep a table per platform. Do **not** hardcode constants for one platform into code that should be platform-agnostic — route them through a `PlatformProfile` or similar configuration type.

## Related

- [[Lattice Generation Pipeline]]
- [[Printing Artifacts and Compensation]]
- [[Meshing Complexity Analysis]]
- [[Mesh Decimation]]
- [[cli crate]]
