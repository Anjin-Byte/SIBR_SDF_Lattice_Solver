//! Axis-aligned box centered at the origin.

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// An axis-aligned box centered at the origin, with half-extents `b`.
///
/// The box occupies the region `[-b.x, b.x] × [-b.y, b.y] × [-b.z, b.z]`.
///
/// # Closed form
///
/// ```text
/// q = abs(p) - b
/// f(p) = length(max(q, 0)) + min(max(q.x, q.y, q.z), 0)
/// ```
///
/// Geometrically: the first term is the distance from the query point to the
/// box when outside (distance from `q.max(0)` to origin); the second term is
/// the (negative) distance from the query point to the nearest face when
/// inside.
///
/// # Invariant
///
/// Every component of `half_extents` is strictly positive and finite.
/// Established by [`AxisBox::new`].
///
/// # Example
///
/// ```
/// use sdf::{AxisBox, Sdf};
/// use glam::vec3;
/// let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
/// assert!(b.eval(vec3(0.0, 0.0, 0.0)) < 0.0);       // inside
/// assert_eq!(b.eval(vec3(2.0, 0.0, 0.0)), 1.0);     // 1 unit past +x face
/// assert!(b.eval(vec3(1.0, 1.0, 1.0)).abs() < 1e-6); // on corner
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisBox {
    half_extents: Vec3,
}

impl AxisBox {
    /// Constructs a box with the given half-extents (along each axis).
    ///
    /// # Errors
    ///
    /// - [`BuildError::NonFinite`] if any component is `NaN` or `±∞`.
    /// - [`BuildError::NonPositive`] if any component is `≤ 0`.
    ///
    /// A zero-extent dimension produces a zero-volume box — a degenerate
    /// surface with undefined normals. We reject it rather than silently
    /// collapse to a lower-dimensional primitive.
    pub fn new(half_extents: Vec3) -> Result<Self, BuildError> {
        if !half_extents.is_finite() {
            // is_finite on Vec3 returns true only if all components are finite.
            // We pick the first non-finite component for the error field.
            let (field, value) = if !half_extents.x.is_finite() {
                ("half_extents.x", half_extents.x)
            } else if !half_extents.y.is_finite() {
                ("half_extents.y", half_extents.y)
            } else {
                ("half_extents.z", half_extents.z)
            };
            return Err(BuildError::NonFinite { field, value });
        }
        if half_extents.x <= 0.0 {
            return Err(BuildError::NonPositive {
                field: "half_extents.x",
                value: half_extents.x,
            });
        }
        if half_extents.y <= 0.0 {
            return Err(BuildError::NonPositive {
                field: "half_extents.y",
                value: half_extents.y,
            });
        }
        if half_extents.z <= 0.0 {
            return Err(BuildError::NonPositive {
                field: "half_extents.z",
                value: half_extents.z,
            });
        }
        Ok(Self { half_extents })
    }

    /// Returns the box's half-extents.
    pub fn half_extents(self) -> Vec3 {
        self.half_extents
    }
}

impl Sdf for AxisBox {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let q = p.abs() - self.half_extents;
        q.max(Vec3::ZERO).length() + q.max_element().min(0.0)
    }
}

impl BoundSdf for AxisBox {}
impl ExactSdf for AxisBox {}

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
    // Target failure class: wrong formula, sign error, wrong face/edge/corner distance.
    // --------------------------------------------------------------

    #[test]
    fn eval_at_origin_is_negative_min_half_extent() {
        // At the center, distance to nearest face = min half-extent, sign negative.
        let b = AxisBox::new(vec3(1.0, 2.0, 3.0)).unwrap();
        assert_eq!(b.eval(Vec3::ZERO), -1.0);
    }

    #[test]
    fn eval_past_face_returns_face_distance() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        assert_eq!(b.eval(vec3(2.0, 0.0, 0.0)), 1.0);
        assert_eq!(b.eval(vec3(0.0, 3.0, 0.0)), 2.0);
        assert_eq!(b.eval(vec3(0.0, 0.0, 5.0)), 4.0);
    }

    #[test]
    fn eval_at_corner_is_approximately_zero() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        assert!(b.eval(vec3(1.0, 1.0, 1.0)).abs() < 1e-6);
    }

    #[test]
    fn eval_at_edge_midpoint_is_approximately_zero() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        assert!(b.eval(vec3(1.0, 1.0, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn eval_at_face_midpoint_is_approximately_zero() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        assert!(b.eval(vec3(1.0, 0.0, 0.0)).abs() < 1e-6);
    }

    #[test]
    fn eval_diagonal_exterior_is_euclidean_distance_to_corner() {
        // At [2, 2, 2] with box of half-extent 1: closest point is corner
        // [1, 1, 1], distance = sqrt(3).
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        let d = b.eval(vec3(2.0, 2.0, 2.0));
        assert!((d - 3.0_f32.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn sign_convention_negative_inside_positive_outside() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        assert!(b.eval(vec3(0.5, 0.5, 0.5)) < 0.0);
        assert!(b.eval(vec3(1.5, 0.0, 0.0)) > 0.0);
    }

    #[test]
    fn half_extents_accessor_returns_construction_value() {
        let he = vec3(1.5, 2.5, 3.5);
        let b = AxisBox::new(he).unwrap();
        assert_eq!(b.half_extents(), he);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: degenerate-input panic, per-axis validation miss.
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_zero_component() {
        assert_eq!(
            AxisBox::new(vec3(0.0, 1.0, 1.0)),
            Err(BuildError::NonPositive {
                field: "half_extents.x",
                value: 0.0
            })
        );
        assert_eq!(
            AxisBox::new(vec3(1.0, 0.0, 1.0)),
            Err(BuildError::NonPositive {
                field: "half_extents.y",
                value: 0.0
            })
        );
        assert_eq!(
            AxisBox::new(vec3(1.0, 1.0, 0.0)),
            Err(BuildError::NonPositive {
                field: "half_extents.z",
                value: 0.0
            })
        );
    }

    #[test]
    fn new_rejects_negative_component() {
        assert_eq!(
            AxisBox::new(vec3(-1.0, 1.0, 1.0)),
            Err(BuildError::NonPositive {
                field: "half_extents.x",
                value: -1.0
            })
        );
    }

    #[test]
    fn new_rejects_nan_component() {
        assert!(matches!(
            AxisBox::new(vec3(f32::NAN, 1.0, 1.0)),
            Err(BuildError::NonFinite {
                field: "half_extents.x",
                ..
            })
        ));
        assert!(matches!(
            AxisBox::new(vec3(1.0, f32::NAN, 1.0)),
            Err(BuildError::NonFinite {
                field: "half_extents.y",
                ..
            })
        ));
        assert!(matches!(
            AxisBox::new(vec3(1.0, 1.0, f32::NAN)),
            Err(BuildError::NonFinite {
                field: "half_extents.z",
                ..
            })
        ));
    }

    #[test]
    fn new_rejects_infinite_component() {
        assert_eq!(
            AxisBox::new(vec3(f32::INFINITY, 1.0, 1.0)),
            Err(BuildError::NonFinite {
                field: "half_extents.x",
                value: f32::INFINITY
            })
        );
    }

    #[test]
    fn new_accepts_asymmetric_extents() {
        let b = AxisBox::new(vec3(1.0, 10.0, 100.0)).unwrap();
        // A point past the smallest half-extent in x returns x-distance;
        // past the largest in z returns z-distance.
        assert_eq!(b.eval(vec3(2.0, 0.0, 0.0)), 1.0);
        assert_eq!(b.eval(vec3(0.0, 0.0, 200.0)), 100.0);
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: numerical pathologies at extremes, Lipschitz violations.
    // --------------------------------------------------------------

    #[test]
    fn eval_tiny_box_preserves_sign_structure() {
        let b = AxisBox::new(Vec3::splat(1e-6)).unwrap();
        assert!(b.eval(Vec3::ZERO) < 0.0);
        assert!(b.eval(vec3(1e-5, 0.0, 0.0)) > 0.0);
    }

    #[test]
    fn eval_huge_box_does_not_overflow() {
        let b = AxisBox::new(Vec3::splat(1e4)).unwrap();
        let d = b.eval(Vec3::splat(1e5));
        assert!(d.is_finite());
        assert!(d > 0.0);
    }

    proptest! {
        /// Lipschitz-1 for the AxisBox SDF.
        #[test]
        fn prop_lipschitz_one(
            hx in 1e-3f32..1e3f32,
            hy in 1e-3f32..1e3f32,
            hz in 1e-3f32..1e3f32,
            pc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
            qc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
        ) {
            let b = AxisBox::new(vec3(hx, hy, hz)).unwrap();
            let p = vec3(pc.0, pc.1, pc.2);
            let q = vec3(qc.0, qc.1, qc.2);
            let df = (b.eval(p) - b.eval(q)).abs();
            let dpq = (p - q).length();
            prop_assert!(df <= dpq + 1e-4,
                "Lipschitz violated: |f(p) - f(q)| = {df}, |p - q| = {dpq}");
        }

        /// Sign consistency: points strictly inside the bounding region are
        /// negative; points outside the bounding region are positive.
        #[test]
        fn prop_sign_consistency(
            he in (1e-2f32..1e2, 1e-2f32..1e2, 1e-2f32..1e2),
            scale in 0.0f32..1.0,  // inside
            dir in (-1.0f32..1.0, -1.0f32..1.0, -1.0f32..1.0),
        ) {
            let he_v = vec3(he.0, he.1, he.2);
            let b = AxisBox::new(he_v).unwrap();
            // Interior: scale < 1 * he along each axis
            let inside = he_v * scale * 0.9;  // strictly inside
            prop_assert!(b.eval(inside) <= 0.0);
            // Exterior: |component| > half_extent by a multiplicative margin
            let outside = he_v * 2.0 + vec3(dir.0, dir.1, dir.2) * 0.01;
            prop_assert!(b.eval(outside) > 0.0);
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests (seeded on plausible bugs)
    // --------------------------------------------------------------

    /// Regression: "Per-axis validation only checked one component."
    #[test]
    fn regression_per_axis_validation() {
        // Each axis must be independently validated.
        assert!(AxisBox::new(vec3(1.0, -1.0, 1.0)).is_err());
        assert!(AxisBox::new(vec3(1.0, 1.0, -1.0)).is_err());
    }

    /// Regression: "Distance to corner used Chebyshev instead of Euclidean."
    #[test]
    fn regression_corner_distance_is_euclidean() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        // At [2, 2, 2], Euclidean distance to corner [1,1,1] is sqrt(3) ≈ 1.732.
        // Chebyshev would give 1.0 (wrong).
        let d = b.eval(vec3(2.0, 2.0, 2.0));
        assert!((d - 3.0_f32.sqrt()).abs() < 1e-5);
    }

    /// Regression: "Sign flipped for interior points."
    #[test]
    fn regression_interior_is_negative() {
        let b = AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap();
        assert!(b.eval(Vec3::ZERO) < 0.0);
        assert!(b.eval(vec3(0.5, 0.5, 0.5)) < 0.0);
    }

    /// Regression: "Zero half-extent silently accepted, producing degenerate
    /// surface with NaN normals."
    #[test]
    fn regression_zero_extent_rejected() {
        assert!(AxisBox::new(vec3(1.0, 0.0, 1.0)).is_err());
    }
}
