//! Sphere — the simplest exact SDF primitive.

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// A sphere of radius `r` centered at the origin.
///
/// # Closed form
///
/// ```text
/// f(p) = |p| - r
/// ```
///
/// This function is exact (returns the true Euclidean distance) and is
/// 1-Lipschitz — its gradient has unit magnitude everywhere except at the
/// origin, where it is undefined (interior pathological point).
///
/// # Invariant
///
/// `radius > 0` and finite. Established by [`Sphere::new`].
///
/// # Example
///
/// ```
/// use sdf::{Sdf, Sphere};
/// let s = Sphere::new(1.0).unwrap();
/// assert!(s.eval(glam::vec3(0.0, 0.0, 0.0)) < 0.0);   // inside
/// assert!((s.eval(glam::vec3(1.0, 0.0, 0.0))).abs() < 1e-6); // on surface
/// assert!(s.eval(glam::vec3(2.0, 0.0, 0.0)) > 0.0);   // outside
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    radius: f32,
}

impl Sphere {
    /// Constructs a sphere with the given radius.
    ///
    /// # Errors
    ///
    /// - [`BuildError::NonFinite`] if `radius` is `NaN` or `±∞`.
    /// - [`BuildError::NonPositive`] if `radius <= 0`.
    ///
    /// A zero-radius sphere is degenerate (it represents a single point, not
    /// a surface); we reject it to avoid propagating zero-volume geometry.
    pub fn new(radius: f32) -> Result<Self, BuildError> {
        if !radius.is_finite() {
            return Err(BuildError::NonFinite {
                field: "radius",
                value: radius,
            });
        }
        if radius <= 0.0 {
            return Err(BuildError::NonPositive {
                field: "radius",
                value: radius,
            });
        }
        Ok(Self { radius })
    }

    /// Returns the sphere's radius.
    pub fn radius(self) -> f32 {
        self.radius
    }
}

impl Sdf for Sphere {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        p.length() - self.radius
    }
}

impl BoundSdf for Sphere {}
impl ExactSdf for Sphere {}

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
    use proptest::prelude::*;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong formula, sign error, broken constructor.
    // --------------------------------------------------------------

    #[test]
    fn eval_at_origin_is_negative_radius() {
        let s = Sphere::new(1.0).unwrap();
        assert_eq!(s.eval(Vec3::ZERO), -1.0);
    }

    #[test]
    fn eval_outside_returns_distance_minus_radius() {
        let s = Sphere::new(1.0).unwrap();
        assert_eq!(s.eval(Vec3::new(2.0, 0.0, 0.0)), 1.0);
        assert_eq!(s.eval(Vec3::new(0.0, 3.0, 0.0)), 2.0);
        assert_eq!(s.eval(Vec3::new(0.0, 0.0, 5.0)), 4.0);
    }

    #[test]
    fn eval_on_surface_is_approximately_zero() {
        let s = Sphere::new(1.0).unwrap();
        assert_eq!(s.eval(Vec3::new(1.0, 0.0, 0.0)).abs(), 0.0);
        assert!(s.eval(Vec3::new(0.6, 0.8, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn radius_accessor_returns_construction_value() {
        let s = Sphere::new(2.5).unwrap();
        assert_eq!(s.radius(), 2.5);
    }

    #[test]
    fn sign_convention_negative_inside_positive_outside() {
        let s = Sphere::new(1.0).unwrap();
        // Inside points have negative SDF.
        assert!(s.eval(Vec3::new(0.5, 0.0, 0.0)) < 0.0);
        assert!(s.eval(Vec3::new(0.0, 0.9, 0.0)) < 0.0);
        // Outside points have positive SDF.
        assert!(s.eval(Vec3::new(1.1, 0.0, 0.0)) > 0.0);
        assert!(s.eval(Vec3::new(2.0, 2.0, 2.0)) > 0.0);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: degenerate-input panic, off-by-one at
    // threshold values, broken identity behavior.
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_zero_radius() {
        assert_eq!(
            Sphere::new(0.0),
            Err(BuildError::NonPositive {
                field: "radius",
                value: 0.0
            })
        );
    }

    #[test]
    fn new_rejects_negative_radius() {
        assert_eq!(
            Sphere::new(-1.0),
            Err(BuildError::NonPositive {
                field: "radius",
                value: -1.0
            })
        );
    }

    #[test]
    fn new_rejects_nan_radius() {
        assert!(matches!(
            Sphere::new(f32::NAN),
            Err(BuildError::NonFinite {
                field: "radius",
                ..
            })
        ));
    }

    #[test]
    fn new_rejects_positive_infinity_radius() {
        assert_eq!(
            Sphere::new(f32::INFINITY),
            Err(BuildError::NonFinite {
                field: "radius",
                value: f32::INFINITY
            })
        );
    }

    #[test]
    fn new_rejects_negative_infinity_radius() {
        assert_eq!(
            Sphere::new(f32::NEG_INFINITY),
            Err(BuildError::NonFinite {
                field: "radius",
                value: f32::NEG_INFINITY
            })
        );
    }

    #[test]
    fn new_accepts_smallest_normal_positive() {
        // Subnormals still pass the `> 0.0` check; this is a boundary case
        // we explicitly allow — callers with resolution requirements can
        // enforce minimums themselves.
        let s = Sphere::new(f32::MIN_POSITIVE).unwrap();
        assert_eq!(s.radius(), f32::MIN_POSITIVE);
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: numerical pathologies at scale extremes,
    // silent Lipschitz violations, sign inconsistencies.
    // --------------------------------------------------------------

    #[test]
    fn eval_tiny_radius_preserves_sign_structure() {
        // A 1-micrometer sphere: queries at 1e-7 should be inside, at 1e-5 outside.
        let s = Sphere::new(1e-6).unwrap();
        assert!(s.eval(Vec3::new(1e-7, 0.0, 0.0)) < 0.0);
        assert!(s.eval(Vec3::new(1e-5, 0.0, 0.0)) > 0.0);
    }

    #[test]
    fn eval_huge_radius_does_not_overflow() {
        // 1e4 is within f32 precision; sqrt(3 * 1e8) ~ 1.7e4 is fine.
        let s = Sphere::new(1e4).unwrap();
        let p = Vec3::splat(1e4);
        let d = s.eval(p);
        assert!(d.is_finite());
        // Outside: |p| = sqrt(3) * 1e4 ≈ 17320, minus 10000 = ~7320.
        assert!(d > 7000.0 && d < 7400.0);
    }

    proptest! {
        /// Lipschitz-1: |f(p) - f(q)| <= |p - q| + epsilon.
        ///
        /// The sphere SDF is mathematically exactly 1-Lipschitz. Floating-point
        /// rounding in `length()` and subtraction allows a small slack. If this
        /// property fails, sphere tracing that uses the returned distance as a
        /// safe step size will silently miss surfaces.
        #[test]
        fn prop_lipschitz_one(
            radius in 1e-3f32..1e3f32,
            pc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
            qc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
        ) {
            let s = Sphere::new(radius).unwrap();
            let p = Vec3::new(pc.0, pc.1, pc.2);
            let q = Vec3::new(qc.0, qc.1, qc.2);
            let df = (s.eval(p) - s.eval(q)).abs();
            let dpq = (p - q).length();
            // Slack of 1e-4 covers typical f32 rounding for magnitudes <= 100.
            prop_assert!(df <= dpq + 1e-4,
                "Lipschitz violated: |f(p) - f(q)| = {df}, |p - q| = {dpq}");
        }

        /// Sign consistency: interior points are negative, exterior points are
        /// positive. A point at fractional distance `t` from origin along any
        /// direction is inside iff `t * |dir| < r`.
        #[test]
        fn prop_sign_consistency(
            radius in 1e-2f32..1e2f32,
            dir in (-1.0f32..1.0, -1.0f32..1.0, -1.0f32..1.0).prop_filter(
                "non-zero direction",
                |(x, y, z)| (x.abs() + y.abs() + z.abs()) > 1e-3,
            ),
            t in 0.01f32..100.0,
        ) {
            let s = Sphere::new(radius).unwrap();
            let dir_v = Vec3::new(dir.0, dir.1, dir.2).normalize();
            let p = dir_v * t;
            let d = s.eval(p);
            if t < radius * 0.999 {
                prop_assert!(d < 0.0, "expected inside, got {d}");
            } else if t > radius * 1.001 {
                prop_assert!(d > 0.0, "expected outside, got {d}");
            }
        }

        /// Boundary tolerance: points sampled on the parametric surface
        /// (unit direction * radius) evaluate close to zero.
        #[test]
        fn prop_on_surface_is_near_zero(
            radius in 1e-2f32..1e2f32,
            dir in (-1.0f32..1.0, -1.0f32..1.0, -1.0f32..1.0).prop_filter(
                "non-zero direction",
                |(x, y, z)| (x.abs() + y.abs() + z.abs()) > 1e-3,
            ),
        ) {
            let s = Sphere::new(radius).unwrap();
            let p = Vec3::new(dir.0, dir.1, dir.2).normalize() * radius;
            let d = s.eval(p).abs();
            // Rounding is ~1 ULP of the radius magnitude.
            prop_assert!(d < radius * 1e-5 + 1e-6,
                "surface point has |eval| = {d}, radius = {radius}");
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests (seeded on plausible bugs)
    //
    // Per the project's Regression Discipline: these are permanent. Even
    // if no bug of this shape ever occurs, the cost of keeping the test
    // is trivial and the cost of losing one would be a silent correctness
    // regression during a future refactor.
    // --------------------------------------------------------------

    /// Regression: "Sphere computed as `|p|² - r²` instead of `|p| - r`".
    /// At `p = [2, 0, 0]` with `r = 1`, the correct value is 1.0, not 3.0.
    #[test]
    fn regression_not_squared_distance() {
        let s = Sphere::new(1.0).unwrap();
        assert_eq!(s.eval(Vec3::new(2.0, 0.0, 0.0)), 1.0);
    }

    /// Regression: "Sign convention flipped" — interior points should be
    /// negative, not positive.
    #[test]
    fn regression_interior_is_negative() {
        let s = Sphere::new(1.0).unwrap();
        assert!(s.eval(Vec3::ZERO) < 0.0);
    }

    /// Regression: "Constructor silently accepted zero radius, producing a
    /// degenerate point-only surface with undefined normals".
    #[test]
    fn regression_zero_radius_is_rejected() {
        assert!(Sphere::new(0.0).is_err());
    }

    /// Regression: "`eval` on NaN input silently returned a plausible-looking
    /// value instead of propagating." `p.length()` of a NaN vector is NaN;
    /// NaN - r = NaN. Behavior is documented as unspecified for NaN input,
    /// but we confirm NaN propagates rather than getting coerced.
    #[test]
    fn regression_nan_input_propagates() {
        let s = Sphere::new(1.0).unwrap();
        assert!(s.eval(Vec3::new(f32::NAN, 0.0, 0.0)).is_nan());
    }

    /// Regression: "Constructor rounded subnormal input to zero and then
    /// rejected it as non-positive". We specifically accept subnormals.
    #[test]
    fn regression_subnormal_radius_accepted() {
        let r = f32::MIN_POSITIVE;
        let s = Sphere::new(r).unwrap();
        assert_eq!(s.radius(), r);
    }
}
