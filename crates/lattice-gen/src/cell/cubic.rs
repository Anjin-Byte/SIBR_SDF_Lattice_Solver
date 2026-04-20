//! Cubic cell topology — three axis-centered capsule struts per cell.
//!
//! # Cell geometry
//!
//! The fundamental domain of [`sdf::Repeat`] and [`sdf::LimitedRepeat`] is
//! the cube `[-L/2, L/2]^3` centered at the origin. For the struts to tile
//! correctly — joining end-to-end across cell boundaries — each strut is
//! centered on its axis, spanning `[-L/2, L/2]` along one axis with
//! `radius r < L/2` perpendicularly.
//!
//! The three struts meet at the origin (the cell's center node). When tiled,
//! cell-N's x-strut endpoint at `+L/2` joins cell-(N+1)'s x-strut endpoint
//! at `-L/2`, producing a continuous x-axis line. The full cubic lattice
//! skeleton emerges as three orthogonal families of such lines.
//!
//! # Why not edge-based struts?
//!
//! An alternative would place struts along the edges of the cube
//! `[0, L]^3`, but that frame doesn't match `Repeat`'s origin-centered
//! fundamental domain. Centering the cell at the origin aligns with the
//! SDF framework and keeps the composition exact.

use glam::vec3;
use sdf::{Capsule, Union};

use crate::error::LatticeError;

/// The concrete cell-body type for the Cubic topology.
///
/// Exposed internally as a type alias so `lattice_body`'s composition stays
/// visible to the compiler for inlining and precision-tracking purposes.
pub(crate) type CubicCellBody = Union<Capsule, Union<Capsule, Capsule>>;

/// Builds the cell-body SDF for the Cubic topology.
///
/// Constructs three axis-centered capsules — each spanning
/// `[-length/2, +length/2]` along its axis at radius `radius` — and unions
/// them. The result is a single unit cell's strut skeleton, ready for
/// tiling via [`sdf::LimitedRepeat`].
///
/// # Errors
///
/// Propagates any underlying [`sdf::BuildError`] from capsule construction.
/// Given `length > 0` and `radius > 0` (both validated upstream), the only
/// way this fails is if a caller bypasses `LatticeJob`'s invariants — so
/// this error path is defensive, not expected in practice.
pub(crate) fn cubic_cell_body(length: f32, radius: f32) -> Result<CubicCellBody, LatticeError> {
    let half = length * 0.5;
    let strut_x = Capsule::new(vec3(-half, 0.0, 0.0), vec3(half, 0.0, 0.0), radius)?;
    let strut_y = Capsule::new(vec3(0.0, -half, 0.0), vec3(0.0, half, 0.0), radius)?;
    let strut_z = Capsule::new(vec3(0.0, 0.0, -half), vec3(0.0, 0.0, half), radius)?;
    Ok(Union {
        a: strut_x,
        b: Union {
            a: strut_y,
            b: strut_z,
        },
    })
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
    use glam::Vec3;
    use sdf::Sdf;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong strut placement (e.g., [0,L] instead of
    // [-L/2, L/2]), wrong axis orientations, missing union arms.
    // --------------------------------------------------------------

    #[test]
    fn eval_at_origin_is_negative_radius() {
        // All three struts pass through the origin.
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        assert!((body.eval(Vec3::ZERO) - (-0.1)).abs() < 1e-6);
    }

    #[test]
    fn eval_on_x_strut_interior_is_negative_radius() {
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        // Anywhere on the x-axis within [-L/2, L/2] is inside the x-strut.
        assert!((body.eval(vec3(0.5, 0.0, 0.0)) - (-0.1)).abs() < 1e-6);
        assert!((body.eval(vec3(-0.5, 0.0, 0.0)) - (-0.1)).abs() < 1e-6);
    }

    #[test]
    fn eval_on_y_strut_interior_is_negative_radius() {
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        assert!((body.eval(vec3(0.0, 0.5, 0.0)) - (-0.1)).abs() < 1e-6);
    }

    #[test]
    fn eval_on_z_strut_interior_is_negative_radius() {
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        assert!((body.eval(vec3(0.0, 0.0, 0.5)) - (-0.1)).abs() < 1e-6);
    }

    #[test]
    fn eval_at_strut_endpoint_approximately_zero() {
        // x-strut endpoint at (+L/2, 0, 0) → on hemispherical cap at that end.
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        let d = body.eval(vec3(1.0, 0.0, 0.0));
        // Hemisphere at endpoint: distance 0 on axis = 0 - r = -r inside.
        // Actually the endpoint itself is the center of the hemisphere, so d = -r.
        assert!((d - (-0.1)).abs() < 1e-6);
    }

    #[test]
    fn eval_at_cell_interior_far_from_struts_is_positive() {
        // Near the cube corner [L/2, L/2, L/2]: all three struts pass far
        // from this point (struts are on the axes).
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        // At [1, 1, 1]: distance to each strut axis is sqrt(1 + 1) ≈ 1.414,
        // minus radius 0.1. All three give roughly the same value; min ≈ 1.314.
        let d = body.eval(vec3(1.0, 1.0, 1.0));
        assert!(d > 1.0, "expected far-from-strut positive, got {d}");
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: degenerate-input defensive path.
    // --------------------------------------------------------------

    #[test]
    fn construction_propagates_sdf_error_on_zero_radius() {
        assert!(cubic_cell_body(2.0, 0.0).is_err());
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Struts placed from `(0,0,0)` to `(L, 0, 0)` instead of
    /// centered on the origin." The `(0, L)` convention doesn't tile
    /// correctly under `Repeat` (whose fundamental domain is centered).
    /// Detection: with centered struts, a query at `(0, 0, 0)` is inside
    /// (center of all three struts). With `(0, L)` struts, `(0, 0, 0)` is
    /// still inside but only because it's an endpoint hemisphere, and more
    /// importantly `(L/2, 0, 0)` — the midpoint of the x-axis in that bad
    /// frame — would be interior to the x-strut. Both conventions coincide
    /// for this single query, but they diverge at other points:
    /// `(-L/4, 0, 0)` is inside the centered x-strut but would be OUTSIDE
    /// an edge-based `(0, L)` strut (which only covers positive x).
    #[test]
    fn regression_struts_are_axis_centered_not_edge_based() {
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        // With centered struts, a point on the -x half-axis is inside.
        let d = body.eval(vec3(-0.5, 0.0, 0.0));
        assert!(d < 0.0, "expected inside (centered strut), got {d}");
    }

    /// Regression: "Struts emitted along diagonal axes instead of pure
    /// coordinate axes — cubic lattice would look like a BCC or similar."
    #[test]
    fn regression_struts_follow_pure_coordinate_axes() {
        let body = cubic_cell_body(2.0, 0.1).unwrap();
        // A point equidistant from x and y axes but not on either:
        // if struts were on the xy-diagonal, this would report inside.
        let d = body.eval(vec3(0.7, 0.7, 0.0));
        assert!(
            d > 0.0,
            "expected point off coordinate axes to be outside, got {d}"
        );
    }
}
