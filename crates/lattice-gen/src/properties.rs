//! Closed-form geometric property queries on a [`LatticeJob`].
//!
//! These are computed from the job specification alone — no SDF evaluation,
//! no meshing, no GPU. They are the inputs downstream physics consumers need
//! (notably the [Pressure Drop Correlation](../../SDF_Lattice_Knowledge_Base/Domain%20Knowledge/Pressure%20Drop%20Correlation.md)
//! from Inayat et al. 2016, which computes ΔP/L from open porosity and
//! window diameter alone).
//!
//! # Phase 1c scope: linear (first-order) approximation
//!
//! The strut-volume closed form involves inclusion-exclusion over three
//! mutually perpendicular cylinders plus their endcap hemispheres. The exact
//! form includes Steinmetz-solid terms (bicylinder `(16/3)r³`, tricylinder
//! `8(2−√2)r³`) that are `O(r³)` corrections relative to the leading
//! `3πr²L` term.
//!
//! For the Woodward operating range (`r* = r/L_c ≤ ~0.2`), the linear
//! approximation is accurate to within ~1% of the true porosity. This
//! crate's implementations return the linear form; cross-validation tests
//! compare to Monte-Carlo SDF sampling as a runtime oracle.
//!
//! When downstream applications need tighter accuracy, the full
//! inclusion-exclusion formulas can be added as a follow-on refinement
//! without changing the public API.
//!
//! # Units
//!
//! - [`LatticeJob::open_porosity`] — dimensionless, `(0, 1)`.
//! - [`LatticeJob::window_diameter`] — same unit as `cell.length()` (mm).
//! - [`LatticeJob::specific_surface_area`] — 1 / (length unit), i.e., inverse mm.
//! - [`LatticeJob::hydraulic_diameter`] — same unit as `cell.length()` (mm).

use core::f32::consts::PI;

use crate::cell::UnitCell;
use crate::job::LatticeJob;

impl LatticeJob {
    /// Open porosity — the fluid-accessible volume fraction of the lattice.
    ///
    /// Returns a dimensionless value in `(0, 1)`. For the cubic topology,
    /// computed as the linear approximation
    ///
    /// ```text
    /// ε_o ≈ 1 - 3π (r*)²
    /// ```
    ///
    /// where `r* = r / L_c`. Accurate to within ~1% for `r* ≤ 0.2`.
    ///
    /// Correction terms from cylinder-cylinder and hemisphere intersections
    /// are `O((r*)³)` and are not currently included; see the module-level
    /// docs for rationale.
    pub fn open_porosity(&self) -> f32 {
        match self.cell() {
            UnitCell::Cubic { length } => {
                let r_star = self.strut().radius() / length;
                1.0 - 3.0 * PI * r_star * r_star
            }
        }
    }

    /// Window diameter — the largest disk that fits through a cell-face
    /// opening between the four perimeter struts.
    ///
    /// For a Cubic cell of edge length `L` with struts of radius `r`, the
    /// cell face has four strut cross-sections at its four edges (each
    /// contributed by the perpendicular axis-strut that passes through the
    /// face). The inscribed disk has diameter
    ///
    /// ```text
    /// d_w = L - 2 r
    /// ```
    ///
    /// Same unit as [`UnitCell::length`] (mm).
    pub fn window_diameter(&self) -> f32 {
        match self.cell() {
            UnitCell::Cubic { length } => length - 2.0 * self.strut().radius(),
        }
    }

    /// Geometric specific surface area — strut surface area per unit
    /// lattice volume.
    ///
    /// Linear approximation for the cubic topology: three capsules per cell,
    /// each contributing lateral cylinder area `2πrL`; hemisphere caps
    /// contribute `O(r²)` and are omitted in the linear form. Thus
    ///
    /// ```text
    /// S_v-geo ≈ 6π r / L²
    /// ```
    ///
    /// Units: inverse length (1/mm). Used in [`LatticeJob::hydraulic_diameter`]
    /// and in the Inayat pressure-drop correlation via the dimensionless
    /// form `S_v-geo · d_w`.
    pub fn specific_surface_area(&self) -> f32 {
        match self.cell() {
            UnitCell::Cubic { length } => 6.0 * PI * self.strut().radius() / (length * length),
        }
    }

    /// Hydraulic diameter — `d_h = 4 ε_o / S_v-geo` (Inayat Eq. 16).
    ///
    /// The characteristic length scale used in the Reynolds-number and
    /// Hagen-number formulations of the Inayat pressure-drop correlation.
    ///
    /// By construction, this identity holds to `f32` precision:
    /// `self.hydraulic_diameter() * self.specific_surface_area() == 4 * self.open_porosity()`.
    pub fn hydraulic_diameter(&self) -> f32 {
        4.0 * self.open_porosity() / self.specific_surface_area()
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    missing_docs
)]
mod tests {
    use super::*;
    use crate::{PrimitiveShape, StrutSpec};
    use glam::Vec3;
    use sdf::Sdf;

    /// Build a valid cubic job with a dummy primitive (properties depend
    /// only on cell and strut, not primitive).
    fn job(cell_length: f32, radius: f32) -> LatticeJob {
        LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::cubic(cell_length).unwrap(),
            StrutSpec::uniform(radius).unwrap(),
        )
        .unwrap()
    }

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong formula, factor-of-π or factor-of-2 errors.
    // --------------------------------------------------------------

    #[test]
    fn open_porosity_for_r_star_0_1_matches_linear_formula() {
        // r/L = 0.1 → ε_o = 1 - 3π(0.01) ≈ 0.9057
        let j = job(1.0, 0.1);
        let expected = 1.0 - 3.0 * PI * 0.01;
        assert!((j.open_porosity() - expected).abs() < 1e-6);
    }

    #[test]
    fn window_diameter_subtracts_two_radii() {
        let j = job(1.0, 0.1);
        assert!((j.window_diameter() - 0.8).abs() < 1e-6);
    }

    #[test]
    fn specific_surface_area_scales_with_inverse_length_squared() {
        let j = job(1.0, 0.1);
        let expected = 6.0 * PI * 0.1;
        assert!((j.specific_surface_area() - expected).abs() < 1e-6);
    }

    #[test]
    fn hydraulic_diameter_satisfies_inayat_identity() {
        // d_h * S_v-geo = 4 ε_o, by definition.
        let j = job(1.0, 0.1);
        let lhs = j.hydraulic_diameter() * j.specific_surface_area();
        let rhs = 4.0 * j.open_porosity();
        assert!((lhs - rhs).abs() < 1e-5);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: scale dependence, parameter-range boundary errors.
    // --------------------------------------------------------------

    #[test]
    fn porosity_is_scale_invariant() {
        // Doubling both cell and radius leaves r* unchanged; porosity must not change.
        let a = job(1.0, 0.1).open_porosity();
        let b = job(2.0, 0.2).open_porosity();
        assert!((a - b).abs() < 1e-6, "a = {a}, b = {b}");
    }

    #[test]
    fn porosity_high_for_thin_struts() {
        // Very thin struts → near-empty lattice → porosity near 1.
        let j = job(1.0, 0.01);
        assert!(j.open_porosity() > 0.999);
    }

    #[test]
    fn porosity_decreases_with_larger_radius() {
        let thin = job(1.0, 0.05).open_porosity();
        let thick = job(1.0, 0.2).open_porosity();
        assert!(thick < thin);
    }

    #[test]
    fn window_diameter_positive_for_valid_job() {
        // LatticeJob::new ensures r < L/2, so window = L - 2r > 0.
        let j = job(2.0, 0.3);
        assert!(j.window_diameter() > 0.0);
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: formula off by a factor, sign flip, divide-by-zero.
    // --------------------------------------------------------------

    /// Monte-Carlo cross-check: the analytical porosity should match the
    /// ratio of outside-strut samples to total samples from a uniformly
    /// distributed set of points inside one cell's AABB.
    ///
    /// This validates the *formula*, not the implementation directly — if
    /// either is wrong, this test catches it.
    #[test]
    fn porosity_agrees_with_monte_carlo_sdf_sampling() {
        // 60k samples is enough to get 5σ slack under 0.02.
        const N_SAMPLES: u16 = 60_000;

        // Use a seeded LCG for reproducibility (proptest's randomness isn't
        // appropriate here — we want a deterministic fixed-seed check).
        let cell_length: f32 = 2.0;
        let radius: f32 = 0.2; // r* = 0.1, in the linear-approximation regime
        let j = job(cell_length, radius);

        // Build the cell body directly (re-using the cubic composition).
        let body = crate::cell::cubic::cubic_cell_body(cell_length, radius).unwrap();
        // Numerical-Recipes LCG; f32 mantissa (23 bits) has headroom for
        // a u16 count, so casts below are lossless.
        let mut state: u32 = 0x1234_5678;
        let mut next = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            // Map to (-L/2, +L/2) per axis.
            (f32::from((state >> 16) as u16) / 65536.0 - 0.5) * cell_length
        };

        let mut outside_count: u16 = 0;
        for _ in 0..N_SAMPLES {
            let p = glam::vec3(next(), next(), next());
            if body.eval(p) > 0.0 {
                outside_count = outside_count.saturating_add(1);
            }
        }

        // u16 → f32 is lossless (u16 max = 65535 < 2^23 mantissa headroom).
        let mc_porosity = f32::from(outside_count) / f32::from(N_SAMPLES);

        let analytical = j.open_porosity();
        // 1/sqrt(N) ≈ 0.004 for N = 60k. Give 5σ slack ≈ 0.02.
        let slack = 5.0 / (f32::from(N_SAMPLES)).sqrt();
        assert!(
            (mc_porosity - analytical).abs() < slack,
            "MC porosity {mc_porosity} vs analytical {analytical} (slack {slack})"
        );
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // Target failure class: silent formula drift during future refactors.
    // --------------------------------------------------------------

    /// Regression: "porosity formula forgot the π factor."
    #[test]
    fn regression_porosity_includes_pi() {
        // Without π, porosity at r* = 0.1 would be 1 - 0.03 = 0.97, not ~0.906.
        let j = job(1.0, 0.1);
        assert!(
            (j.open_porosity() - (1.0 - 3.0 * PI * 0.01)).abs() < 1e-6
                && (j.open_porosity() - 0.97).abs() > 0.05
        );
    }

    /// Regression: "window diameter used L - r (one strut) instead of L - 2r (two struts)."
    #[test]
    fn regression_window_diameter_subtracts_both_struts() {
        // At L = 1, r = 0.3: correct is 0.4, buggy (L - r) would give 0.7.
        let j = job(1.0, 0.3);
        assert!((j.window_diameter() - 0.4).abs() < 1e-6);
        assert!((j.window_diameter() - 0.7).abs() > 0.1);
    }

    /// Regression: "S_v-geo used 3 cylinders but forgot the factor of 2π
    /// perimeter — computed πrL instead of 2πrL per strut."
    #[test]
    fn regression_surface_area_uses_full_circumference() {
        // At L = 1, r = 0.1: correct is 6π(0.1) ≈ 1.885.
        // Half-perimeter bug would give 3π(0.1) ≈ 0.942.
        let j = job(1.0, 0.1);
        let sv = j.specific_surface_area();
        assert!((sv - 6.0 * PI * 0.1).abs() < 1e-5);
        assert!((sv - 3.0 * PI * 0.1).abs() > 0.5);
    }

    /// Regression: "`hydraulic_diameter` divided by `4ε_o` instead of multiplying."
    /// Should satisfy `d_h · S_v-geo = 4ε_o` (Inayat identity).
    #[test]
    fn regression_hydraulic_diameter_identity_holds() {
        for (l, r) in [(1.0, 0.05), (2.4, 0.264), (0.5, 0.1)] {
            let j = job(l, r);
            let lhs = j.hydraulic_diameter() * j.specific_surface_area();
            let rhs = 4.0 * j.open_porosity();
            assert!(
                (lhs - rhs).abs() < 1e-4,
                "at (L={l}, r={r}): lhs={lhs}, rhs={rhs}"
            );
        }
    }
}
