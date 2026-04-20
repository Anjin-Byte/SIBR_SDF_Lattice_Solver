//! Primitive shapes — the outer boundary that trims the tiled lattice body.
//!
//! Phase 1a supports two primitive shapes: [`PrimitiveShape::Cube`] and
//! [`PrimitiveShape::Cylinder`]. Internally, each variant materializes as a
//! corresponding [`sdf`] primitive stored in a [`BoundaryShape`] enum that
//! itself implements the SDF trait hierarchy.
//!
//! The enum dispatch keeps `lattice_body`'s return type concrete (and thus
//! inlinable) while letting callers vary the primitive at runtime.

use glam::Vec3;
use sdf::{AxisBox, BoundSdf, CappedCylinder, ExactSdf, Sdf};

use crate::error::LatticeError;

/// A validated lattice-boundary shape.
///
/// Constructed via the `PrimitiveShape::cube` / `PrimitiveShape::cylinder`
/// helpers, each of which validates parameters and wraps construction errors
/// in [`LatticeError`]. The internal representation is the corresponding
/// `sdf` primitive — see [`BoundaryShape`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrimitiveShape {
    shape: BoundaryShape,
}

impl PrimitiveShape {
    /// Constructs an axis-aligned cube primitive with the given half-extents.
    ///
    /// # Errors
    ///
    /// Returns [`LatticeError::Sdf`] wrapping the underlying construction
    /// failure if any component of `half_extents` is non-finite or
    /// non-positive. See [`sdf::AxisBox::new`] for precise rules.
    pub fn cube(half_extents: Vec3) -> Result<Self, LatticeError> {
        Ok(Self {
            shape: BoundaryShape::Cube(AxisBox::new(half_extents)?),
        })
    }

    /// Constructs a capped-cylinder primitive with the given endpoints and radius.
    ///
    /// # Errors
    ///
    /// Returns [`LatticeError::Sdf`] wrapping the underlying construction
    /// failure if endpoints are coincident, the radius is non-positive, or
    /// any input is non-finite. See [`sdf::CappedCylinder::new`].
    pub fn cylinder(start: Vec3, end: Vec3, radius: f32) -> Result<Self, LatticeError> {
        Ok(Self {
            shape: BoundaryShape::Cylinder(CappedCylinder::new(start, end, radius)?),
        })
    }

    /// Returns the axis-aligned bounding box (min, max corners) of the primitive.
    ///
    /// Used internally to compute [`sdf::LimitedRepeat`] extents. The box is
    /// tight for `Cube` and a slight over-approximation for `Cylinder` (the
    /// axis-aligned bounding box of the cylinder's extent).
    pub(crate) fn aabb(self) -> (Vec3, Vec3) {
        match self.shape {
            BoundaryShape::Cube(b) => {
                let he = b.half_extents();
                (-he, he)
            }
            BoundaryShape::Cylinder(c) => {
                // Axis-aligned bounding box of a cylinder: expand the
                // start-end segment by `radius` in every axis. This is an
                // over-approximation, but only by at most one cell on each
                // side — acceptable for LimitedRepeat extent sizing.
                let lo = c.start().min(c.end()) - Vec3::splat(c.radius());
                let hi = c.start().max(c.end()) + Vec3::splat(c.radius());
                (lo, hi)
            }
        }
    }

    /// Returns the internal SDF representation.
    pub(crate) fn boundary(self) -> BoundaryShape {
        self.shape
    }
}

/// The concrete SDF used as the lattice's outer boundary.
///
/// This enum exists to keep `lattice_body`'s return type concrete and
/// inlinable. Each variant wraps the corresponding `sdf` primitive and
/// delegates evaluation; the `ExactSdf` impl propagates because each variant
/// is itself exact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum BoundaryShape {
    Cube(AxisBox),
    Cylinder(CappedCylinder),
}

impl Sdf for BoundaryShape {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        match self {
            Self::Cube(b) => b.eval(p),
            Self::Cylinder(c) => c.eval(p),
        }
    }
}

impl BoundSdf for BoundaryShape {}
impl ExactSdf for BoundaryShape {}

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

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong enum dispatch, missing variant.
    // --------------------------------------------------------------

    #[test]
    fn cube_eval_matches_underlying_axis_box() {
        let prim = PrimitiveShape::cube(Vec3::splat(1.0)).unwrap();
        let axis_box = AxisBox::new(Vec3::splat(1.0)).unwrap();
        for p in [Vec3::ZERO, vec3(0.5, 0.0, 0.0), vec3(2.0, 0.0, 0.0)] {
            assert_eq!(prim.boundary().eval(p), axis_box.eval(p));
        }
    }

    #[test]
    fn cylinder_eval_matches_underlying_capped_cylinder() {
        let prim = PrimitiveShape::cylinder(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.5).unwrap();
        let cyl = CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.5).unwrap();
        for p in [Vec3::ZERO, vec3(0.5, 0.0, 0.0), vec3(2.0, 0.0, 0.0)] {
            assert_eq!(prim.boundary().eval(p), cyl.eval(p));
        }
    }

    #[test]
    fn cube_aabb_matches_half_extents() {
        let he = vec3(1.0, 2.0, 3.0);
        let prim = PrimitiveShape::cube(he).unwrap();
        let (lo, hi) = prim.aabb();
        assert_eq!(lo, -he);
        assert_eq!(hi, he);
    }

    #[test]
    fn cylinder_aabb_contains_endpoints_plus_radius() {
        let prim = PrimitiveShape::cylinder(Vec3::ZERO, vec3(5.0, 0.0, 0.0), 0.5).unwrap();
        let (lo, hi) = prim.aabb();
        // Axis spans [0, 5] in x; radius inflates each axis by 0.5.
        assert_eq!(lo, vec3(-0.5, -0.5, -0.5));
        assert_eq!(hi, vec3(5.5, 0.5, 0.5));
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: construction error propagation.
    // --------------------------------------------------------------

    #[test]
    fn cube_rejects_zero_half_extent() {
        assert!(PrimitiveShape::cube(vec3(0.0, 1.0, 1.0)).is_err());
    }

    #[test]
    fn cube_rejects_nan_half_extent() {
        assert!(PrimitiveShape::cube(vec3(f32::NAN, 1.0, 1.0)).is_err());
    }

    #[test]
    fn cylinder_rejects_coincident_endpoints() {
        let p = vec3(1.0, 0.0, 0.0);
        assert!(PrimitiveShape::cylinder(p, p, 0.5).is_err());
    }

    #[test]
    fn cylinder_rejects_zero_radius() {
        assert!(PrimitiveShape::cylinder(Vec3::ZERO, Vec3::X, 0.0).is_err());
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "`BoundaryShape::Cube` and `BoundaryShape::Cylinder`
    /// variants accidentally delegated to each other, producing wrong SDF."
    #[test]
    fn regression_enum_dispatch_uses_correct_variant() {
        let cube = PrimitiveShape::cube(Vec3::splat(1.0)).unwrap();
        let cyl = PrimitiveShape::cylinder(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.5).unwrap();
        // A point far from the cylinder axis but inside the cube.
        // Cube returns negative (inside), cylinder returns positive (outside).
        let p = vec3(0.0, 0.0, 0.7);
        assert!(cube.boundary().eval(p) < 0.0);
        assert!(cyl.boundary().eval(p) > 0.0);
    }
}
