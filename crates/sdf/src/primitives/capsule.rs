//! Capsule (line-segment swept sphere) — the canonical lattice strut primitive.

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Minimum separation between capsule endpoints, in world units.
/// Below this, the distance formula becomes numerically unstable.
const MIN_ENDPOINT_SEPARATION: f32 = 1e-6;

/// A capsule: a line segment of radius `r` with hemispherical endcaps.
///
/// The capsule is defined by two endpoints `a` and `b` (the axis of the
/// segment) and a radius `r`. It is the canonical strut primitive for
/// lattice generation.
///
/// # Closed form
///
/// ```text
/// pa = p - a
/// ba = b - a
/// h  = clamp(dot(pa, ba) / dot(ba, ba), 0, 1)
/// f(p) = length(pa - ba * h) - r
/// ```
///
/// Geometrically: `h` is the parameter along the axis nearest to `p`,
/// clamped to `[0, 1]`. The distance from `p` to that nearest-axis-point is
/// the core of the SDF; subtracting `r` inflates the line into a capsule.
/// The clamp produces hemispherical endcaps automatically.
///
/// # Invariants
///
/// - `radius > 0` and finite.
/// - `|b - a| >= MIN_ENDPOINT_SEPARATION` (endpoints are not coincident).
///
/// Both established by [`Capsule::new`].
///
/// # Example
///
/// ```
/// use sdf::{Capsule, Sdf};
/// use glam::vec3;
/// let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
/// # use glam::Vec3;
/// // Center of axis is inside at distance r from the cylinder surface.
/// assert!(c.eval(vec3(0.5, 0.0, 0.0)) < 0.0);
/// // A point at [0.5, 0.2, 0] is distance 0.2 from axis, minus r = 0.1.
/// assert!((c.eval(vec3(0.5, 0.2, 0.0)) - 0.1).abs() < 1e-6);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Capsule {
    a: Vec3,
    b: Vec3,
    radius: f32,
}

impl Capsule {
    /// Constructs a capsule with the given endpoints and radius.
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
                reason: "capsule endpoints are coincident or nearly so",
            });
        }
        Ok(Self { a, b, radius })
    }

    /// Returns the capsule's first endpoint.
    pub fn start(self) -> Vec3 {
        self.a
    }

    /// Returns the capsule's second endpoint.
    pub fn end(self) -> Vec3 {
        self.b
    }

    /// Returns the capsule's radius.
    pub fn radius(self) -> f32 {
        self.radius
    }
}

/// Extracts the first non-finite component name and value.
/// Callers invoke this only after they have established the input is
/// non-finite; at least one of `{x, y, z}` is guaranteed to be `NaN` or `±∞`.
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

impl Sdf for Capsule {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let pa = p - self.a;
        let ba = self.b - self.a;
        let h = (pa.dot(ba) / ba.dot(ba)).clamp(0.0, 1.0);
        (pa - ba * h).length() - self.radius
    }
}

impl BoundSdf for Capsule {}
impl ExactSdf for Capsule {}

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
    // Target failure class: wrong formula, wrong endcap behavior, axis mis-projection.
    // --------------------------------------------------------------

    #[test]
    fn eval_on_axis_midpoint_is_negative_radius() {
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert_eq!(c.eval(vec3(0.5, 0.0, 0.0)), -0.1);
    }

    #[test]
    fn eval_perpendicular_to_axis_middle() {
        // 0.2 away from the axis at the midpoint, radius 0.1 → distance 0.1.
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!((c.eval(vec3(0.5, 0.2, 0.0)) - 0.1).abs() < 1e-6);
    }

    #[test]
    fn eval_beyond_endpoint_a_uses_hemispherical_cap() {
        // Past the a endpoint, along -x: distance uses hemisphere around a.
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!((c.eval(vec3(-1.0, 0.0, 0.0)) - 0.9).abs() < 1e-6);
    }

    #[test]
    fn eval_beyond_endpoint_b_uses_hemispherical_cap() {
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!((c.eval(vec3(2.0, 0.0, 0.0)) - 0.9).abs() < 1e-6);
    }

    #[test]
    fn eval_at_endpoint_on_surface_is_approximately_zero() {
        // The hemisphere centered at `a` has its surface at `a + r * n` for
        // any unit `n` pointing away from `b`.
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!(c.eval(vec3(-0.1, 0.0, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn eval_on_cylindrical_surface_is_approximately_zero() {
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!(c.eval(vec3(0.5, 0.1, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn sign_convention_inside_is_negative() {
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        assert!(c.eval(vec3(0.5, 0.05, 0.0)) < 0.0);
        assert!(c.eval(vec3(0.5, 0.2, 0.0)) > 0.0);
    }

    #[test]
    fn eval_arbitrary_axis_capsule_works() {
        // A capsule along the diagonal of the unit cube.
        let a = Vec3::ZERO;
        let b = vec3(1.0, 1.0, 1.0);
        let c = Capsule::new(a, b, 0.1).unwrap();
        // The midpoint is inside, at distance r from a cylindrical surface.
        let mid = (a + b) * 0.5;
        assert_eq!(c.eval(mid), -0.1);
    }

    #[test]
    fn accessors_return_construction_values() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = Capsule::new(a, b, 0.7).unwrap();
        assert_eq!(c.start(), a);
        assert_eq!(c.end(), b);
        assert_eq!(c.radius(), 0.7);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: degenerate-endpoint divide-by-zero,
    // missing validation for any of {a, b, radius}.
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_coincident_endpoints() {
        let p = vec3(1.0, 2.0, 3.0);
        assert!(matches!(
            Capsule::new(p, p, 0.1),
            Err(BuildError::Degenerate { .. })
        ));
    }

    #[test]
    fn new_rejects_near_coincident_endpoints() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = a + Vec3::splat(1e-8); // below MIN_ENDPOINT_SEPARATION
        assert!(matches!(
            Capsule::new(a, b, 0.1),
            Err(BuildError::Degenerate { .. })
        ));
    }

    #[test]
    fn new_rejects_zero_radius() {
        assert_eq!(
            Capsule::new(Vec3::ZERO, Vec3::X, 0.0),
            Err(BuildError::NonPositive {
                field: "radius",
                value: 0.0
            })
        );
    }

    #[test]
    fn new_rejects_negative_radius() {
        assert_eq!(
            Capsule::new(Vec3::ZERO, Vec3::X, -0.5),
            Err(BuildError::NonPositive {
                field: "radius",
                value: -0.5
            })
        );
    }

    #[test]
    fn new_rejects_nan_in_any_field() {
        assert!(Capsule::new(vec3(f32::NAN, 0.0, 0.0), Vec3::X, 0.1).is_err());
        assert!(Capsule::new(Vec3::ZERO, vec3(f32::NAN, 0.0, 0.0), 0.1).is_err());
        assert!(Capsule::new(Vec3::ZERO, Vec3::X, f32::NAN).is_err());
    }

    #[test]
    fn new_accepts_minimum_separation() {
        // Exactly the threshold should succeed.
        let c = Capsule::new(Vec3::ZERO, vec3(MIN_ENDPOINT_SEPARATION, 0.0, 0.0), 0.1);
        assert!(c.is_ok());
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: Lipschitz violations, scale pathologies, axis-direction errors.
    // --------------------------------------------------------------

    proptest! {
        /// Lipschitz-1 for arbitrary-axis capsules.
        #[test]
        fn prop_lipschitz_one(
            ac in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
            bc in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
            r in 1e-3f32..1.0,
            pc in (-20.0f32..20.0, -20.0f32..20.0, -20.0f32..20.0),
            qc in (-20.0f32..20.0, -20.0f32..20.0, -20.0f32..20.0),
        ) {
            let a = vec3(ac.0, ac.1, ac.2);
            let b = vec3(bc.0, bc.1, bc.2);
            // Reject degenerate inputs from the test (the constructor would too).
            prop_assume!((b - a).length() >= 0.1);
            let c = Capsule::new(a, b, r).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            let q = vec3(qc.0, qc.1, qc.2);
            let df = (c.eval(p) - c.eval(q)).abs();
            let dpq = (p - q).length();
            prop_assert!(df <= dpq + 1e-4,
                "Lipschitz violated: df = {df}, dpq = {dpq}");
        }

        /// Symmetry: swapping endpoints must produce the same SDF value.
        #[test]
        fn prop_endpoint_swap_symmetric(
            ac in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
            bc in (-10.0f32..10.0, -10.0f32..10.0, -10.0f32..10.0),
            r in 1e-3f32..1.0,
            pc in (-20.0f32..20.0, -20.0f32..20.0, -20.0f32..20.0),
        ) {
            let a = vec3(ac.0, ac.1, ac.2);
            let b = vec3(bc.0, bc.1, bc.2);
            prop_assume!((b - a).length() >= 0.1);
            let c1 = Capsule::new(a, b, r).unwrap();
            let c2 = Capsule::new(b, a, r).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            prop_assert!((c1.eval(p) - c2.eval(p)).abs() < 1e-5);
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests (seeded on plausible bugs)
    // --------------------------------------------------------------

    /// Regression: "Axis clamp was missing, producing cylinder-extending-to-infinity
    /// behavior past the endpoints."
    #[test]
    fn regression_endpoints_produce_hemispherical_caps_not_infinite_cylinder() {
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        // Point far past the +x end: distance must grow along the axis
        // (hemisphere behavior), not stay at the perpendicular distance.
        let d_far = c.eval(vec3(10.0, 0.0, 0.0));
        // If clamp were missing, this would evaluate as distance to the infinite
        // axis (= 0) - 0.1 = -0.1 (wrong). Correct: 10 - 1 - 0.1 = 8.9.
        assert!((d_far - 8.9).abs() < 1e-5);
    }

    /// Regression: "Dot product used b-a instead of ba, producing wrong projection."
    #[test]
    fn regression_axis_projection_correct_direction() {
        // Near the 'a' endpoint, moving from a toward b, we should go from
        // hemisphere region (clamped at h=0) to cylindrical region (h>0).
        let c = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        // At x = -0.5 (past 'a' along -axis), h should clamp to 0, and distance
        // to point [0,0,0] is 0.5, minus 0.1 = 0.4.
        assert!((c.eval(vec3(-0.5, 0.0, 0.0)) - 0.4).abs() < 1e-6);
    }

    /// Regression: "Coincident endpoints produced NaN rather than a construction error."
    #[test]
    fn regression_coincident_endpoints_rejected() {
        let c = Capsule::new(Vec3::ZERO, Vec3::ZERO, 0.1);
        assert!(matches!(c, Err(BuildError::Degenerate { .. })));
    }

    /// Regression: "Sign of SDF depended on endpoint ordering."
    #[test]
    fn regression_endpoint_ordering_doesnt_affect_sign() {
        let c1 = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
        let c2 = Capsule::new(vec3(1.0, 0.0, 0.0), Vec3::ZERO, 0.1).unwrap();
        let q = vec3(0.5, 0.0, 0.0);
        assert_eq!(c1.eval(q), c2.eval(q));
    }
}
