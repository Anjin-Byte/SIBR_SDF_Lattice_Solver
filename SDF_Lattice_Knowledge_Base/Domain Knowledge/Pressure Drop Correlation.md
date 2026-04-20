---
title: Pressure Drop Correlation
tags: [domain, physics, pressure-drop, inayat]
source: "Inayat et al. 2016, Chem. Eng. J. 287, 704–719"
source-file: reference/pressureDropCorrelation.pdf
---

# Pressure Drop Correlation

A closed-form, theoretically-derived correlation for fluid pressure drop across open-cell foams, due to Inayat et al. (2016). The correlation has **no fitted constants** — it depends only on geometric parameters and fluid properties — and validates to <10% MRPE across a wide range of ceramic and metal foams. For the SDF Lattice Solver, it is the canonical **downstream consumer** of generated lattice geometry (Layer 3 of the [[Three-Layer Stack]]).

## The master equation

For flow through an open-cell medium of length L at superficial velocity V:

$$
\frac{\Delta P}{L} = \underbrace{32\,\tau^{2}\,\frac{\mu V}{\varepsilon_{o}\,d_{h}^{2}}}_{\text{viscous (Darcian)}} \;+\; \underbrace{\frac{\tau^{3}}{2}\,\frac{\rho V^{2}}{\varepsilon_{o}^{2}\,d_{h}}}_{\text{inertial (non-Darcian)}}
\quad(\text{Eq. 22})
$$

where:
- `μ` — dynamic viscosity of the fluid (Pa·s)
- `ρ` — fluid density (kg/m³)
- `V` — superficial fluid velocity (m/s)
- `ε_o` — **open porosity** (hydrodynamic porosity; fraction of the foam volume available to the fluid)
- `d_h` — **hydraulic diameter** of the foam (m)
- `τ` — **geometric tortuosity** of the flow path (dimensionless)

## The two auxiliary relations

### Hydraulic diameter

$$
d_{h} = \frac{4\,\varepsilon_{o}}{S_{v\text{-geo}}}
\quad(\text{Eq. 16})
$$

where `S_v-geo` is the geometric specific surface area (m² per m³ of geometric volume).

### Geometric tortuosity

$$
\tau = 1 + \phi\,\frac{1 - 0.971\,(1-\varepsilon_{o})^{0.5}}{4\,\varepsilon_{o}\,(1-\varepsilon_{o})^{0.5}}\,(1-\varepsilon_{o})
\quad(\text{Eq. 35})
$$

where `φ` is a shape factor depending on strut cross-section:

| Cross-section | φ |
|---|---|
| Cylindrical | 4.87 |
| Triangular | 5.62 |
| Concave triangular | 6.49 |

### Specific surface area (tetrakaidecahedron model)

$$
S_{v\text{-geo}} \cdot d_{w} = \phi\,\frac{1 - 0.971\,(1-\varepsilon_{o})^{0.5}}{(1-\varepsilon_{o})^{0.5}}\,(1-\varepsilon_{o})
\quad(\text{Eq. 34})
$$

where `d_w` is the **window diameter** — the average opening between adjacent cells. Rearranging gives `S_v-geo` from `ε_o` and `d_w` alone.

## The punch line

Substituting Eq. 34 → Eq. 16 → Eq. 22, the pressure drop is predictable from **just two geometric properties**:

- `ε_o` — open porosity
- `d_w` — window diameter

plus the strut-cross-section shape factor `φ`. This is the headline result of the paper.

## Geometry vocabulary (§3.2, Fig. 7)

The paper is strict about three often-conflated terms — the code should be equally strict:

| Term | Symbol | Meaning |
|---|---|---|
| **Cell diameter** | `d_c` | Diameter of the (approximately spherical) void enclosed by struts |
| **Window diameter** | `d_w` | Diameter of the opening between adjacent cells |
| **Strut diameter** | `d_s` | Diameter of the bar connecting nodes |
| **Pore size** | `D_p = (d_w + d_s)` | Combined length scale used in some correlations |

Empirical relation (foam literature): `d_w ≈ d_c / 2.3`, used when `d_w` isn't measured directly.

## Porosity taxonomy (§3.1)

The paper distinguishes multiple porosity definitions — code must not conflate them:

| Term | Symbol | Definition |
|---|---|---|
| **Open porosity** | `ε_o` | Void fraction available to fluid (relevant for fluid dynamics) |
| **Strut porosity** | `ε_s` | Internal void fraction inside struts (manufacturing artifact, e.g. hollow struts from polymer-replication ceramic foams) |
| **Total porosity** | `ε_t = ε_o + ε_s` | Sum |
| **Nominal porosity** | `ε_n` | Manufacturer-reported, typically ambiguous |

**For the SDF Lattice Solver**, strut porosity is **zero** — struts are modeled as solid SDF objects. Therefore `ε_o = ε_t` for our generated lattices, and the ambiguity the paper worries about doesn't apply. We still use `ε_o` explicitly to preserve traceability.

## Strut cross-section in practice (§3.3.2)

Foam struts change cross-section with porosity due to the manufacturing process:

| Foam type | Circular | Triangular | Concave-triangular |
|---|---|---|---|
| Ceramic | ε_o < 0.9 | ε_o > 0.9 | (insufficient data) |
| Metal | ε_o < 0.92 | 0.92–0.95 | ε_o > 0.95 |

**For lattices** (not foams), struts are **explicitly designed**. We choose the cross-section — most likely cylindrical (circular) in a first implementation, because it matches SDF capsule/cylinder primitives cleanly. This gives `φ = 4.87`.

**Important caveat:** the correlation was validated on foams, where cross-section is a consequence of manufacturing. Lattices with intentionally non-circular struts are outside the validated range. Use with care until we validate against lattice pressure-drop measurements directly.

## Dimensionless form

The paper also presents Eq. 22 in Reynolds/Hagen-number form:

$$
\mathrm{Re} = \frac{d_{h}\,V\,\rho}{\varepsilon_{o}\,\mu}
\qquad
\mathrm{Hg} = \frac{\Delta P}{L} \cdot \frac{d_{h}^{3}\,\rho}{\mu^{2}}
$$

$$
\mathrm{Hg} = (A)_{d_{h}}\,\mathrm{Re} \;+\; (B)_{d_{h}}\,\mathrm{Re}^{2}
\quad(\text{Eq. 27})
$$

with `(A)_{d_h} = 32τ²` and `(B)_{d_h} = τ³/2`. Useful when comparing against Forchheimer-style correlations or for range/validity checks.

## Validation range

The paper validates against measured pressure-drop data from 12+ sources covering:
- **Materials:** Aluminum, Al₂O₃ (alumina), Mullite, β-SiC, SSiC, Nickel
- **Fluids:** Air, water, glycerol, oil
- **PPI:** 5 to 45 (cells per inch)
- **Open porosity:** 0.75 to 0.99
- **Mean Relative Percentage Error (MRPE):** <10% on most sample sets; worst case ~26% on one anomalous sample

**For lattices (not foams)**, the validated range is empirically unknown. Lattices are *more regular* than foams, so first-principles reasoning suggests predictions should be at least as good. This is a hypothesis worth testing experimentally.

## Implementation shape in Rust

Keeping the correlation in [[lattice-gen crate]] (pure, deterministic, testable without a GPU — and logically belongs with the other lattice-geometric-property queries like porosity and window diameter):

```rust
pub struct FluidProperties {
    pub density: KilogramsPerCubicMeter,     // ρ
    pub viscosity: PascalSeconds,            // μ
}

pub struct GeometryProperties {
    pub open_porosity: f64,                  // ε_o ∈ (0, 1)
    pub window_diameter: Millimeters,        // d_w
    pub strut_cross_section: StrutCrossSection,
}

pub enum StrutCrossSection {
    Cylindrical,       // φ = 4.87
    Triangular,        // φ = 5.62
    ConcaveTriangular, // φ = 6.49
}

pub fn pressure_drop_per_length(
    fluid: &FluidProperties,
    geom: &GeometryProperties,
    velocity: MetersPerSecond,
) -> PascalsPerMeter { /* implement Eq. 22 via Eq. 16, 34, 35 */ }
```

Invariants (enforced by newtype constructors per [[Newtypes and Domain Types]]):
- `ε_o ∈ (0, 1)` strictly — neither zero nor one is physically valid
- All dimensions positive
- Velocity non-negative

## As an oracle for the solver

The correlation is a [[Sharp Oracles|sharp oracle]] at the Layer-2 / Layer-3 boundary of the [[Three-Layer Stack]]:

1. Generate a lattice with [[Lattice Generation Pipeline]].
2. Compute `ε_o` and `d_w` directly from the SDF representation (both are geometric queries on the generated mesh).
3. Predict ΔP/L from Eq. 22.
4. Compare against measured ΔP/L on the printed part.

Agreement within the paper's <10% MRPE band is the end-to-end verdict on the full stack. Disagreement is diagnostic — either the lattice is geometrically wrong, the pressure-drop model is inapplicable to the geometry, or the measurement is bad. Each is investigable.

## Open questions

1. Does `φ = 4.87` transfer from cylindrical-strut foams to cylindrical-strut lattices?
2. Does the tetrakaidecahedron-based Eq. 34 apply to non-Kelvin cell topologies, or only to Kelvin-like geometries?
3. For `ε_o < 0.75`, where the correlation is less validated, do our designed lattices still obey it?

These are testable hypotheses, not blocking concerns — but they are worth [[Regression Discipline|recording as adversarial test candidates]].

## Related

- [[Three-Layer Stack]]
- [[Lattice Generation Pipeline]]
- [[Unit Cell Topologies]]
- [[lattice-gen crate]]
- [[Newtypes and Domain Types]]
- [[Sharp Oracles]]
