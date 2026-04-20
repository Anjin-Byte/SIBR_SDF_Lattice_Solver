//! Capped cylinder with arbitrary-axis endpoints (flat endcaps).

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Minimum separation between cylinder endpoints, in world units.
const MIN_ENDPOINT_SEPARATION: f32 = 1e-6;

/// A capped cylinder defined by two endpoints `a` and `b` (the axis) and a
/// radius `r`. The cylinder has **flat** endcaps at the perpendicular planes
/// through `a` and `b` — distinguishing it from [`crate::Capsule`] which has
/// hemispherical endcaps.
///
/// # Closed form
///
/// From the IQ catalog (arbitrary-axis variant):
///
/// ```text
/// ba   = b - a
/// pa   = p - a
/// baba = ba · ba
/// paba = pa · ba
/// x    = |pa*baba - ba*paba| - r*baba
/// y    = |paba - baba/2| - baba/2
/// x2   = x²
/// y2   = y² * baba
/// d    = if max(x, y) < 0 { -min(x2, y2) } else { (if x > 0 { x2 } else { 0 }) + (if y > 0 { y2 } else { 0 }) }
/// f(p) = sign(d) * sqrt(|d|) / baba
/// ```
///
/// Geometrically: `x` measures perpendicular distance to the axis minus
/// radius; `y` measures axial position relative to the caps. The compound
/// formula handles all cases (inside, outside on side, outside past caps,
/// past caps and outside radially) with a single branch on signs.
///
/// # Invariants
///
/// - `radius > 0` and finite.
/// - `|b - a| >= MIN_ENDPOINT_SEPARATION` (endpoints are not coincident).
///
/// Both established by [`CappedCylinder::new`].
///
/// # Example
///
/// ```
/// use sdf::{CappedCylinder, Sdf};
/// use glam::{Vec3, vec3};
/// let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
/// // Midpoint of axis is inside.
/// assert!(c.eval(vec3(0.5, 0.0, 0.0)) < 0.0);
/// // Just past the +x cap, outside.
/// assert!(c.eval(vec3(1.1, 0.0, 0.0)) > 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CappedCylinder {
    a: Vec3,
    b: Vec3,
    radius: f32,
}

impl CappedCylinder {
    /// Constructs a capped cylinder with the given endpoints and radius.
    ///
    /// # Errors
    ///
    /// - [`BuildError::NonFinite`] if any input component is `NaN` or `±∞`.
    /// - [`BuildError::NonPositive`] if `radius <= 0`.
    /// - [`BuildError::Degenerate`] if `|b - a| < MIN_ENDPOINT_SEPARATION`.
    pub fn new(a: Vec3, b: Vec3, radius: f32) -> Result<Self, BuildError> {
        if !a.is_finite() {
            let (field, value) = non_finite_component("a", a);
            return Err(BuildError::NonFinite { field, value });
        }
        if !b.is_finite() {
            let (field, value) = non_finite_component("b", b);
            return Err(BuildError::NonFinite { field, value });
        }
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
        if (b - a).length() < MIN_ENDPOINT_SEPARATION {
            return Err(BuildError::Degenerate {
                reason: "cylinder endpoints are coincident or nearly so",
            });
        }
        Ok(Self { a, b, radius })
    }

    /// Returns the first endpoint (the flat cap at `a`).
    pub fn start(self) -> Vec3 {
        self.a
    }

    /// Returns the second endpoint (the flat cap at `b`).
    pub fn end(self) -> Vec3 {
        self.b
    }

    /// Returns the cylinder's radius.
    pub fn radius(self) -> f32 {
        self.radius
    }
}

/// See [`crate::primitives::capsule`] for the equivalent helper;
/// kept separate to avoid a cross-module dependency between primitives.
fn non_finite_component(prefix: &'static str, v: Vec3) -> (&'static str, f32) {
    if !v.x.is_finite() {
        match prefix {
            "a" => ("a.x", v.x),
            _ => ("b.x", v.x),
        }
    } else if !v.y.is_finite() {
        match prefix {
            "a" => ("a.y", v.y),
            _ => ("b.y", v.y),
        }
    } else {
        match prefix {
            "a" => ("a.z", v.z),
            _ => ("b.z", v.z),
        }
    }
}

impl Sdf for CappedCylinder {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let ba = self.b - self.a;
        let pa = p - self.a;
        let baba = ba.dot(ba);
        let paba = pa.dot(ba);
        let x = (pa * baba - ba * paba).length() - self.radius * baba;
        let y = (paba - baba * 0.5).abs() - baba * 0.5;
        let x2 = x * x;
        let y2 = y * y * baba;
        let d = if x.max(y) < 0.0 {
            -x2.min(y2)
        } else {
            let xs = if x > 0.0 { x2 } else { 0.0 };
            let ys = if y > 0.0 { y2 } else { 0.0 };
            xs + ys
        };
        d.signum() * d.abs().sqrt() / baba
    }
}

impl BoundSdf for CappedCylinder {}
impl ExactSdf for CappedCylinder {}

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
    use glam::vec3;
    use proptest::prelude::*;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong compound-formula signs, flat-cap vs
    // hemispherical-cap confusion.
    // --------------------------------------------------------------

    #[test]
    fn eval_on_axis_midpoint_is_negative_radius() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!((c.eval(vec3(0.5, 0.0, 0.0)) - (-0.1)).abs() < 1e-6);
    }

    #[test]
    fn eval_past_flat_cap_is_axial_distance() {
        // The flat cap at b = [1, 0, 0] means a point at [1.5, 0, 0] is
        // at perpendicular distance 0 from axis, 0.5 past the cap along axis.
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!((c.eval(vec3(1.5, 0.0, 0.0)) - 0.5).abs() < 1e-5);
    }

    #[test]
    fn eval_past_cap_and_off_axis_uses_edge_distance() {
        // Point past the cap AND outside the cylinder radius: distance is to
        // the circular edge of the cap.
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        // At [1.5, 0.3, 0]: past cap by 0.5, perp distance 0.3, radius 0.1.
        // Distance to edge (closest point on rim: [1, 0.1, 0]) = sqrt(0.5² + 0.2²)
        let expected = (0.5_f32.powi(2) + 0.2_f32.powi(2)).sqrt();
        let d = c.eval(vec3(1.5, 0.3, 0.0));
        assert!(
            (d - expected).abs() < 1e-4,
            "d = {d}, expected = {expected}"
        );
    }

    #[test]
    fn eval_on_side_at_radius_is_approximately_zero() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!(c.eval(vec3(0.5, 0.1, 0.0)).abs() < 1e-5);
    }

    #[test]
    fn eval_on_flat_cap_is_approximately_zero() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        // Exactly at the cap plane, inside the radial disk.
        assert!(c.eval(vec3(1.0, 0.05, 0.0)).abs() < 1e-5);
    }

    #[test]
    fn sign_convention_inside_is_negative() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!(c.eval(vec3(0.5, 0.05, 0.0)) < 0.0);
        assert!(c.eval(vec3(0.5, 0.2, 0.0)) > 0.0); // outside radially
        assert!(c.eval(vec3(1.5, 0.0, 0.0)) > 0.0); // outside past cap
        assert!(c.eval(vec3(-0.5, 0.0, 0.0)) > 0.0); // outside before cap
    }

    #[test]
    fn accessors_return_construction_values() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = CappedCylinder::new(a, b, 0.7).unwrap();
        assert_eq!(c.start(), a);
        assert_eq!(c.end(), b);
        assert_eq!(c.radius(), 0.7);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_coincident_endpoints() {
        let p = vec3(1.0, 2.0, 3.0);
        assert!(matches!(
            CappedCylinder::new(p, p, 0.1),
            Err(BuildError::Degenerate { .. })
        ));
    }

    #[test]
    fn new_rejects_zero_and_negative_radius() {
        assert!(CappedCylinder::new(Vec3::ZERO, Vec3::X, 0.0).is_err());
        assert!(CappedCylinder::new(Vec3::ZERO, Vec3::X, -0.5).is_err());
    }

    #[test]
    fn new_rejects_nan_in_any_field() {
        assert!(CappedCylinder::new(vec3(f32::NAN, 0.0, 0.0), Vec3::X, 0.1).is_err());
        assert!(CappedCylinder::new(Vec3::ZERO, vec3(0.0, f32::NAN, 0.0), 0.1).is_err());
        assert!(CappedCylinder::new(Vec3::ZERO, Vec3::X, f32::NAN).is_err());
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // --------------------------------------------------------------

    proptest! {
        /// Lipschitz-1: the capped cylinder SDF is exactly 1-Lipschitz in exact math.
        #[test]
        fn prop_lipschitz_one(
            ac in (-5.0f32..5.0, -5.0f32..5.0, -5.0f32..5.0),
            bc in (-5.0f32..5.0, -5.0f32..5.0, -5.0f32..5.0),
            r in 1e-3f32..1.0,
            pc in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
            qc in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
        ) {
            let a = vec3(ac.0, ac.1, ac.2);
            let b = vec3(bc.0, bc.1, bc.2);
            prop_assume!((b - a).length() >= 0.1);
            let c = CappedCylinder::new(a, b, r).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            let q = vec3(qc.0, qc.1, qc.2);
            let df = (c.eval(p) - c.eval(q)).abs();
            let dpq = (p - q).length();
            // CappedCylinder has sharper rounding near the rim; give a slightly
            // larger slack than other primitives.
            prop_assert!(df <= dpq + 1e-3,
                "Lipschitz violated: df = {df}, dpq = {dpq}");
        }

        /// Endpoint swap symmetry — the cylinder shape is the same
        /// regardless of which endpoint is labelled `a`.
        #[test]
        fn prop_endpoint_swap_symmetric(
            ac in (-5.0f32..5.0, -5.0f32..5.0, -5.0f32..5.0),
            bc in (-5.0f32..5.0, -5.0f32..5.0, -5.0f32..5.0),
            r in 1e-3f32..1.0,
            pc in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
        ) {
            let a = vec3(ac.0, ac.1, ac.2);
            let b = vec3(bc.0, bc.1, bc.2);
            prop_assume!((b - a).length() >= 0.1);
            let c1 = CappedCylinder::new(a, b, r).unwrap();
            let c2 = CappedCylinder::new(b, a, r).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            prop_assert!((c1.eval(p) - c2.eval(p)).abs() < 1e-4);
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Flat-cap vs hemispherical-cap confused — distance past cap
    /// was computed via capsule hemisphere instead of flat plane."
    /// For a cylinder, a point directly past the cap at [1.5, 0, 0] has
    /// distance 0.5 (to the flat cap), not something else.
    #[test]
    fn regression_flat_cap_not_hemispherical() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!((c.eval(vec3(1.5, 0.0, 0.0)) - 0.5).abs() < 1e-5);
    }

    /// Regression: "Sign computation on-axis midpoint returned positive."
    #[test]
    fn regression_on_axis_midpoint_is_inside() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!(c.eval(vec3(0.5, 0.0, 0.0)) < 0.0);
    }

    /// Regression: "`d.sqrt()` crashed or NaN-ed when `d` was slightly negative
    /// due to rounding." The formula uses `d.abs()` before sqrt precisely to
    /// avoid this.
    #[test]
    fn regression_no_nan_near_surface() {
        let c = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        // Sample many near-surface points; none should yield NaN.
        for t in 0..100_i32 {
            let theta = f32::from(i16::try_from(t).unwrap()) * 0.0628;
            let p = vec3(0.5, 0.1 * theta.cos(), 0.1 * theta.sin());
            let d = c.eval(p);
            assert!(d.is_finite(), "got non-finite at t={t}: d = {d}");
        }
    }

    /// Regression: "Coincident endpoints caused divide-by-zero → NaN."
    #[test]
    fn regression_coincident_endpoints_rejected() {
        assert!(CappedCylinder::new(Vec3::ZERO, Vec3::ZERO, 0.1).is_err());
    }
}
