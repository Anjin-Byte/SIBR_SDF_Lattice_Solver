//! Axis-reflection symmetry operators — reflect the query point across a
//! coordinate plane before delegating.

use glam::Vec3;

use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Mirrors the inner SDF about the YZ plane (applies `p.x = |p.x|`).
///
/// The resulting SDF represents the union of the inner shape and its
/// reflection across the YZ plane. Reflection is an isometry, so distance
/// is preserved exactly.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, SymX, Translate, Sdf};
/// use glam::vec3;
/// // A sphere offset to +x, mirrored — creates a sphere on both sides of YZ.
/// let s = Sphere::new(1.0).unwrap();
/// let offset_sphere = Translate { offset: vec3(3.0, 0.0, 0.0), inner: s };
/// let mirrored = SymX { inner: offset_sphere };
/// // The inner sphere is at +3; the mirror is at -3.
/// assert_eq!(mirrored.eval(vec3(3.0, 0.0, 0.0)), -1.0);
/// assert_eq!(mirrored.eval(vec3(-3.0, 0.0, 0.0)), -1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymX<T> {
    /// The inner SDF being reflected.
    pub inner: T,
}

impl<T: Sdf> Sdf for SymX<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let mut q = p;
        q.x = q.x.abs();
        self.inner.eval(q)
    }
}

impl<T: BoundSdf> BoundSdf for SymX<T> {}
impl<T: ExactSdf> ExactSdf for SymX<T> {}

/// Mirrors the inner SDF about the XZ plane (applies `p.y = |p.y|`).
///
/// See [`SymX`] for the design rationale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymY<T> {
    /// The inner SDF being reflected.
    pub inner: T,
}

impl<T: Sdf> Sdf for SymY<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let mut q = p;
        q.y = q.y.abs();
        self.inner.eval(q)
    }
}

impl<T: BoundSdf> BoundSdf for SymY<T> {}
impl<T: ExactSdf> ExactSdf for SymY<T> {}

/// Mirrors the inner SDF about the XY plane (applies `p.z = |p.z|`).
///
/// See [`SymX`] for the design rationale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymZ<T> {
    /// The inner SDF being reflected.
    pub inner: T,
}

impl<T: Sdf> Sdf for SymZ<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let mut q = p;
        q.z = q.z.abs();
        self.inner.eval(q)
    }
}

impl<T: BoundSdf> BoundSdf for SymZ<T> {}
impl<T: ExactSdf> ExactSdf for SymZ<T> {}

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
    use crate::{Sphere, Translate};
    use glam::vec3;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: abs applied to wrong axis, reflection direction wrong.
    // --------------------------------------------------------------

    #[test]
    fn sym_x_mirrors_on_positive_and_negative_x() {
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(3.0, 0.0, 0.0),
            inner: s,
        };
        let m = SymX { inner: shifted };
        // Mirror image exists at -3 as well.
        assert_eq!(m.eval(vec3(3.0, 0.0, 0.0)), -1.0);
        assert_eq!(m.eval(vec3(-3.0, 0.0, 0.0)), -1.0);
    }

    #[test]
    fn sym_x_does_not_affect_other_axes() {
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(0.0, 3.0, 0.0),
            inner: s,
        };
        let m = SymX { inner: shifted };
        // Inner is on +y — NOT mirrored in y. Query at -y should NOT match.
        assert!(m.eval(vec3(0.0, 3.0, 0.0)) < 0.0); // matches inner
        assert!(m.eval(vec3(0.0, -3.0, 0.0)) > 0.0); // outside
    }

    #[test]
    fn sym_y_mirrors_on_y_axis() {
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(0.0, 3.0, 0.0),
            inner: s,
        };
        let m = SymY { inner: shifted };
        assert_eq!(m.eval(vec3(0.0, 3.0, 0.0)), -1.0);
        assert_eq!(m.eval(vec3(0.0, -3.0, 0.0)), -1.0);
    }

    #[test]
    fn sym_z_mirrors_on_z_axis() {
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(0.0, 0.0, 3.0),
            inner: s,
        };
        let m = SymZ { inner: shifted };
        assert_eq!(m.eval(vec3(0.0, 0.0, 3.0)), -1.0);
        assert_eq!(m.eval(vec3(0.0, 0.0, -3.0)), -1.0);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: degenerate behavior when inner is on the mirror plane.
    // --------------------------------------------------------------

    #[test]
    fn sym_of_origin_centered_inner_is_idempotent() {
        // If the inner is already centered, SymX has no effect.
        let s = Sphere::new(1.0).unwrap();
        let m = SymX { inner: s };
        for p in [
            Vec3::ZERO,
            vec3(0.5, 0.0, 0.0),
            vec3(2.0, 0.0, 0.0),
            vec3(-2.0, 0.0, 0.0),
        ] {
            assert_eq!(m.eval(p), s.eval(p));
        }
    }

    #[test]
    fn sym_is_idempotent() {
        // SymX(SymX(X)) should evaluate identically to SymX(X) at any point,
        // because p.x.abs().abs() == p.x.abs().
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(2.0, 1.0, 0.0),
            inner: s,
        };
        let single = SymX { inner: shifted };
        let double = SymX {
            inner: SymX { inner: shifted },
        };
        for p in [
            vec3(0.0, 0.0, 0.0),
            vec3(3.0, 1.0, 0.0),
            vec3(-3.0, 1.0, 0.0),
            vec3(5.0, 0.0, 0.5),
        ] {
            assert_eq!(single.eval(p), double.eval(p));
        }
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // --------------------------------------------------------------

    #[test]
    fn all_three_axes_together_produces_octant_symmetry() {
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(2.0, 3.0, 4.0),
            inner: s,
        };
        // Applying all three yields a sphere in every octant-equivalent position.
        let m = SymZ {
            inner: SymY {
                inner: SymX { inner: shifted },
            },
        };
        // The inner sphere exists at [±2, ±3, ±4] in all 8 octants.
        assert_eq!(m.eval(vec3(2.0, 3.0, 4.0)), -1.0);
        assert_eq!(m.eval(vec3(-2.0, -3.0, -4.0)), -1.0);
        assert_eq!(m.eval(vec3(-2.0, 3.0, 4.0)), -1.0);
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "`SymX` applied `-p.x` instead of `|p.x|`, which is
    /// rotation not reflection — doesn't produce a union of reflected copies."
    /// Detection: for an inner sphere at +3, with pure reflection the -3
    /// query returns inside; with negation (identity rotation by 180°),
    /// the result depends on which side the inner is, making only one query match.
    #[test]
    fn regression_abs_not_negate() {
        let s = Sphere::new(1.0).unwrap();
        let shifted = Translate {
            offset: vec3(3.0, 0.0, 0.0),
            inner: s,
        };
        let m = SymX { inner: shifted };
        // Both sides should be inside.
        assert!(m.eval(vec3(3.0, 0.0, 0.0)) < 0.0);
        assert!(m.eval(vec3(-3.0, 0.0, 0.0)) < 0.0);
    }

    /// Regression: "`SymX` used `p.abs()` on the whole vector, affecting y and z too."
    #[test]
    fn regression_only_affected_axis_mirrored() {
        let s = Sphere::new(1.0).unwrap();
        // Inner sphere deliberately offset only in +y.
        let shifted = Translate {
            offset: vec3(0.0, 3.0, 0.0),
            inner: s,
        };
        let m = SymX { inner: shifted };
        // Inner is at [0, 3, 0]. Query at [0, -3, 0] must NOT match (y isn't
        // mirrored by SymX). If bug applied abs to whole vector, it would.
        assert!(m.eval(vec3(0.0, -3.0, 0.0)) > 0.0);
    }
}
