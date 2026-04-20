//! `UniformScale` — isotropic scaling preserving exact distance.

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Uniformly scales the inner SDF by a factor `scale`. Evaluates
/// `inner.eval(p / scale) * scale`.
///
/// # Why `* scale` matters
///
/// Without the trailing multiplication, the returned value is in the inner's
/// scaled frame — a distance in "inner units," not world units. Multiplying
/// by `scale` converts back to world distance, preserving the SDF's
/// 1-Lipschitz property and keeping the function exact rather than a bound.
///
/// This is the only operator that changes metric, and the only one where
/// `* scale` is doing real work rather than just being notation.
///
/// # Invariant
///
/// `scale > 0` and finite. Non-uniform or negative scaling would break
/// the distance-preservation property. Established by [`UniformScale::new`].
///
/// # Example
///
/// ```
/// use sdf::{Sphere, UniformScale, Sdf};
/// use glam::vec3;
/// let s = Sphere::new(1.0).unwrap();
/// // Scale by 2 — effectively a sphere of radius 2.
/// let scaled = UniformScale::new(2.0, s).unwrap();
/// // The surface is at distance 2 from origin.
/// assert!((scaled.eval(vec3(2.0, 0.0, 0.0))).abs() < 1e-6);
/// // Inside at the center: distance -2.
/// assert_eq!(scaled.eval(vec3(0.0, 0.0, 0.0)), -2.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformScale<T> {
    scale: f32,
    inner: T,
}

impl<T> UniformScale<T> {
    /// Constructs a uniform scaling of `inner` by `scale`.
    ///
    /// # Errors
    ///
    /// - [`BuildError::NonFinite`] if `scale` is `NaN` or `±∞`.
    /// - [`BuildError::NonPositive`] if `scale <= 0`.
    pub fn new(scale: f32, inner: T) -> Result<Self, BuildError> {
        if !scale.is_finite() {
            return Err(BuildError::NonFinite {
                field: "scale",
                value: scale,
            });
        }
        if scale <= 0.0 {
            return Err(BuildError::NonPositive {
                field: "scale",
                value: scale,
            });
        }
        Ok(Self { scale, inner })
    }

    /// Returns the scale factor.
    pub fn scale(&self) -> f32 {
        self.scale
    }

    /// Returns a reference to the inner SDF.
    pub fn inner(&self) -> &T {
        &self.inner
    }
}

impl<T: Sdf> Sdf for UniformScale<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        self.inner.eval(p / self.scale) * self.scale
    }
}

impl<T: BoundSdf> BoundSdf for UniformScale<T> {}
impl<T: ExactSdf> ExactSdf for UniformScale<T> {}

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
    use crate::Sphere;
    use glam::vec3;
    use proptest::prelude::*;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: forgot * s, wrong scale direction.
    // --------------------------------------------------------------

    #[test]
    fn scale_by_two_doubles_effective_radius() {
        let s = Sphere::new(1.0).unwrap();
        let scaled = UniformScale::new(2.0, s).unwrap();
        // A sphere of radius 1 scaled by 2 has effective radius 2.
        // At the origin: distance = -2.
        assert_eq!(scaled.eval(Vec3::ZERO), -2.0);
        // On the surface at [2, 0, 0]: distance ≈ 0.
        assert!(scaled.eval(vec3(2.0, 0.0, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn scale_by_one_is_identity() {
        let s = Sphere::new(1.0).unwrap();
        let scaled = UniformScale::new(1.0, s).unwrap();
        for p in [Vec3::ZERO, vec3(0.5, 0.0, 0.0), vec3(2.0, 1.0, 0.0)] {
            assert_eq!(scaled.eval(p), s.eval(p));
        }
    }

    #[test]
    fn scale_by_half_halves_effective_radius() {
        let s = Sphere::new(1.0).unwrap();
        let scaled = UniformScale::new(0.5, s).unwrap();
        assert_eq!(scaled.eval(Vec3::ZERO), -0.5);
        assert!(scaled.eval(vec3(0.5, 0.0, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn accessors_return_construction_values() {
        let s = Sphere::new(1.0).unwrap();
        let scaled = UniformScale::new(3.0, s).unwrap();
        assert_eq!(scaled.scale(), 3.0);
        assert_eq!(*scaled.inner(), s);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_zero_scale() {
        let s = Sphere::new(1.0).unwrap();
        assert_eq!(
            UniformScale::new(0.0, s),
            Err(BuildError::NonPositive {
                field: "scale",
                value: 0.0,
            })
        );
    }

    #[test]
    fn new_rejects_negative_scale() {
        let s = Sphere::new(1.0).unwrap();
        assert_eq!(
            UniformScale::new(-1.0, s),
            Err(BuildError::NonPositive {
                field: "scale",
                value: -1.0,
            })
        );
    }

    #[test]
    fn new_rejects_nan_scale() {
        let s = Sphere::new(1.0).unwrap();
        assert!(matches!(
            UniformScale::new(f32::NAN, s),
            Err(BuildError::NonFinite { field: "scale", .. })
        ));
    }

    #[test]
    fn new_rejects_infinity_scale() {
        let s = Sphere::new(1.0).unwrap();
        assert!(UniformScale::new(f32::INFINITY, s).is_err());
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: Lipschitz violation at extreme scales.
    // --------------------------------------------------------------

    proptest! {
        /// Lipschitz-1 preserved under uniform scaling at any positive scale.
        #[test]
        fn prop_lipschitz_preserved_under_scale(
            r in 1e-2f32..1e2f32,
            scale in 1e-2f32..1e2f32,
            pc in (-50.0f32..50.0, -50.0f32..50.0, -50.0f32..50.0),
            qc in (-50.0f32..50.0, -50.0f32..50.0, -50.0f32..50.0),
        ) {
            let scaled = UniformScale::new(scale, Sphere::new(r).unwrap()).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            let q = vec3(qc.0, qc.1, qc.2);
            let df = (scaled.eval(p) - scaled.eval(q)).abs();
            let dpq = (p - q).length();
            // Relative slack scales with magnitude: 1e-5 * max eval.
            let slack = 1e-4 * (r * scale).max(1.0);
            prop_assert!(df <= dpq + slack,
                "df = {df}, dpq = {dpq}, slack = {slack}");
        }

        /// Scale composition: scale(s1, scale(s2, X)) should be equivalent to scale(s1*s2, X).
        #[test]
        fn prop_scale_composition(
            r in 1e-2f32..1e2f32,
            s1 in 1e-2f32..5.0,
            s2 in 1e-2f32..5.0,
            pc in (-20.0f32..20.0, -20.0f32..20.0, -20.0f32..20.0),
        ) {
            let sphere = Sphere::new(r).unwrap();
            let double = UniformScale::new(s1, UniformScale::new(s2, sphere).unwrap()).unwrap();
            let combined = UniformScale::new(s1 * s2, sphere).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            // Slack accounts for the chained division/multiplication rounding.
            prop_assert!((double.eval(p) - combined.eval(p)).abs() < 1e-3 * (r * s1 * s2).max(1.0));
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Forgot the `* scale` at the end — the function would
    /// return distances in the scaled frame rather than world units, making
    /// the SDF only a bound, not exact."
    /// Detection: distances from a scaled-by-2 sphere at the outer surface
    /// should be 0 (surface), not 0 / 2 = 0 (same value by coincidence),
    /// but at interior [1, 0, 0]: correct is -1 (world distance to surface
    /// at [2, 0, 0]), bugged would be -0.5 (scaled-frame distance).
    #[test]
    fn regression_multiply_by_scale_applied() {
        let s = Sphere::new(1.0).unwrap();
        let scaled = UniformScale::new(2.0, s).unwrap();
        // Point at [1, 0, 0] is inside the scaled-by-2 sphere (radius 2).
        // Correct world distance = 1 - 2 = -1.
        // Without the * 2, would be (1 - 2) / 2 = -0.5, wrong.
        assert_eq!(scaled.eval(vec3(1.0, 0.0, 0.0)), -1.0);
    }

    /// Regression: "Division direction wrong — scaled by `p * s` instead of
    /// `p / s`, producing a sphere that shrinks with increasing scale."
    #[test]
    fn regression_divide_not_multiply() {
        let s = Sphere::new(1.0).unwrap();
        let scaled = UniformScale::new(2.0, s).unwrap();
        // Scaled sphere's surface is at [2, 0, 0]. Query there → 0.
        // If inverted, surface would be at [0.5, 0, 0], and [2, 0, 0] would
        // return a large positive value.
        assert!(scaled.eval(vec3(2.0, 0.0, 0.0)).abs() < 1e-6);
    }

    /// Regression: "Zero scale silently accepted, causing divide-by-zero at eval."
    #[test]
    fn regression_zero_scale_rejected() {
        let s = Sphere::new(1.0).unwrap();
        assert!(UniformScale::new(0.0, s).is_err());
    }
}
