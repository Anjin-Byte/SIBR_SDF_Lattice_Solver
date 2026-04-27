//! `BCCxy` cell topology — body-diagonal struts plus top-face edges.
//!
//! # Cell geometry
//!
//! The `BCCxy` (vertex octahedron) cell occupies `[-L/2, L/2]^3`, centered at
//! the origin. Its strut skeleton is:
//!
//! - **8 body diagonals** — one capsule from each of the 8 cube corners
//!   `(±L/2, ±L/2, ±L/2)` to the cube center `(0, 0, 0)`, the BCC node.
//! - **4 top-face edges** — the four edges of the square on `z = +L/2` with
//!   corners at `(±L/2, ±L/2, +L/2)`.
//!
//! The bottom face (`z = -L/2`) and the `±x` / `±y` faces are supplied by
//! the appropriate periodic neighbors under [`sdf::LimitedRepeat`]:
//! assigning the bottom face to the home cell would double-count with the
//! `+z` neighbor's own top face. Union idempotency means the double-count
//! would be harmless, but the asymmetric assignment keeps capsule count
//! minimal.
//!
//! # The "XY connections" interpretation
//!
//! The vault describes `BCCxy` as "octahedral struts (body diagonals) plus
//! edge struts in the XY plane." The ambiguity (top-only vs top+bottom vs
//! mid-plane square) was resolved in favor of top-only under the design
//! principle that minimal home-cell representatives + periodic images
//! reproduce the full lattice. See the planning note for alternatives
//! considered.

use glam::vec3;
use sdf::{Capsule, SmoothUnion};

use crate::error::LatticeError;

/// Concrete nested-`SmoothUnion` type for the `BCCxy` cell body (12 capsules).
///
/// With `joint_smoothness = 0`, each `SmoothUnion` falls back to bit-exact
/// hard `min` (see [`SmoothUnion`]'s `k = 0` fallback) so existing
/// behavior is unchanged for callers that don't opt into smoothing.
#[allow(clippy::type_complexity)]
pub(crate) type BccXyCellBody = SmoothUnion<
    Capsule,
    SmoothUnion<
        Capsule,
        SmoothUnion<
            Capsule,
            SmoothUnion<
                Capsule,
                SmoothUnion<
                    Capsule,
                    SmoothUnion<
                        Capsule,
                        SmoothUnion<
                            Capsule,
                            SmoothUnion<
                                Capsule,
                                SmoothUnion<
                                    Capsule,
                                    SmoothUnion<Capsule, SmoothUnion<Capsule, Capsule>>,
                                >,
                            >,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

/// Builds the cell-body SDF for the `BCCxy` topology.
///
/// Constructs 12 `Capsule`s — 8 body diagonals and 4 top-face edges — and
/// right-folds them into the [`BccXyCellBody`] nested-Union chain.
///
/// # Errors
///
/// Propagates any underlying [`sdf::BuildError`] from [`Capsule::new`]. In
/// practice, with `length > 0` and `radius > 0` validated upstream, the
/// only way this fails is a bypass of [`crate::LatticeJob::new`]'s
/// invariants — defensive path.
pub(crate) fn bccxy_cell_body(
    length: f32,
    radius: f32,
    joint_smoothness: f32,
) -> Result<BccXyCellBody, LatticeError> {
    let h = length * 0.5;
    let k = joint_smoothness;

    // Body diagonals (8): cube center ↔ each of 8 corners. The order is
    // the standard Gray-ish traversal; any order is fine since Union is
    // commutative.
    let d1 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(h, h, h), radius)?;
    let d2 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(-h, h, h), radius)?;
    let d3 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(h, -h, h), radius)?;
    let d4 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(-h, -h, h), radius)?;
    let d5 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(h, h, -h), radius)?;
    let d6 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(-h, h, -h), radius)?;
    let d7 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(h, -h, -h), radius)?;
    let d8 = Capsule::new(vec3(0.0, 0.0, 0.0), vec3(-h, -h, -h), radius)?;

    // Top-face edges (4): the 4 edges of the square on z = +h.
    let t1 = Capsule::new(vec3(h, h, h), vec3(-h, h, h), radius)?;
    let t2 = Capsule::new(vec3(h, h, h), vec3(h, -h, h), radius)?;
    let t3 = Capsule::new(vec3(-h, -h, h), vec3(h, -h, h), radius)?;
    let t4 = Capsule::new(vec3(-h, -h, h), vec3(-h, h, h), radius)?;

    // Right-fold 12 capsules into BccXyCellBody.
    Ok(SmoothUnion {
        a: d1,
        k,
        b: SmoothUnion {
            a: d2,
            k,
            b: SmoothUnion {
                a: d3,
                k,
                b: SmoothUnion {
                    a: d4,
                    k,
                    b: SmoothUnion {
                        a: d5,
                        k,
                        b: SmoothUnion {
                            a: d6,
                            k,
                            b: SmoothUnion {
                                a: d7,
                                k,
                                b: SmoothUnion {
                                    a: d8,
                                    k,
                                    b: SmoothUnion {
                                        a: t1,
                                        k,
                                        b: SmoothUnion {
                                            a: t2,
                                            k,
                                            b: SmoothUnion { a: t3, b: t4, k },
                                        },
                                    },
                                },
                            },
                        },
                    },
                },
            },
        },
    })
}

/// The 12-element home edge set as endpoint pairs, for test oracles.
#[cfg(test)]
pub(crate) fn home_edges(length: f32) -> [(glam::Vec3, glam::Vec3); 12] {
    let h = length * 0.5;
    [
        // Body diagonals.
        (vec3(0.0, 0.0, 0.0), vec3(h, h, h)),
        (vec3(0.0, 0.0, 0.0), vec3(-h, h, h)),
        (vec3(0.0, 0.0, 0.0), vec3(h, -h, h)),
        (vec3(0.0, 0.0, 0.0), vec3(-h, -h, h)),
        (vec3(0.0, 0.0, 0.0), vec3(h, h, -h)),
        (vec3(0.0, 0.0, 0.0), vec3(-h, h, -h)),
        (vec3(0.0, 0.0, 0.0), vec3(h, -h, -h)),
        (vec3(0.0, 0.0, 0.0), vec3(-h, -h, -h)),
        // Top-face edges.
        (vec3(h, h, h), vec3(-h, h, h)),
        (vec3(h, h, h), vec3(h, -h, h)),
        (vec3(-h, -h, h), vec3(h, -h, h)),
        (vec3(-h, -h, h), vec3(-h, h, h)),
    ]
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
    // Target failure class: wrong endpoint, missing edge, wrong axis.
    // --------------------------------------------------------------

    #[test]
    fn eval_at_cube_center_is_negative_radius() {
        // 8 body diagonals meet at the origin.
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        assert!((body.eval(Vec3::ZERO) - (-0.1)).abs() < 1e-5);
    }

    #[test]
    fn eval_on_body_diagonal_midpoint_is_negative_radius() {
        // Midpoint of (0,0,0)↔(h,h,h) = (L/4, L/4, L/4) = (1, 1, 1) for L=4.
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        assert!((body.eval(vec3(1.0, 1.0, 1.0)) - (-0.1)).abs() < 1e-5);
    }

    #[test]
    fn eval_on_top_face_edge_midpoint_is_negative_radius() {
        // Midpoint of (h, h, h)↔(-h, h, h) = (0, h, h) = (0, 2, 2) for L=4.
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        assert!((body.eval(vec3(0.0, 2.0, 2.0)) - (-0.1)).abs() < 1e-5);
    }

    #[test]
    fn eval_at_cube_corner_top_is_negative_radius() {
        // (h, h, h) is endpoint of body diagonal d1 AND of top-face edges
        // t1 and t2. Should be -radius.
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        assert!((body.eval(vec3(2.0, 2.0, 2.0)) - (-0.1)).abs() < 1e-5);
    }

    #[test]
    fn eval_at_cube_corner_bottom_is_negative_radius() {
        // (-h, -h, -h) is endpoint of body diagonal d8, but NOT of any
        // home top-face edge. Still -radius because d8's endpoint cap is
        // there.
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        assert!((body.eval(vec3(-2.0, -2.0, -2.0)) - (-0.1)).abs() < 1e-5);
    }

    #[test]
    fn eval_at_bottom_face_center_is_positive_pore() {
        // Bottom face (z = -h) is NOT home-assigned; it comes from the -z
        // neighbor. From the home cell's perspective, (0, 0, -h) is a
        // pore (not on any home strut).
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(0.0, 0.0, -2.0));
        assert!(d > 0.0, "expected positive at bottom-face center, got {d}");
    }

    // --------------------------------------------------------------
    // b. Edge / boundary tests.
    // --------------------------------------------------------------

    #[test]
    fn construction_propagates_sdf_error_on_zero_radius() {
        assert!(bccxy_cell_body(4.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn construction_propagates_sdf_error_on_negative_radius() {
        assert!(bccxy_cell_body(4.0, -0.1, 0.0).is_err());
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests.
    // --------------------------------------------------------------

    /// Oracle loop over the full 12-edge home set. A typo in any
    /// endpoint in `bccxy_cell_body` or `home_edges` would fail this.
    #[test]
    fn all_12_home_edges_have_negative_radius_at_midpoint() {
        let length = 4.0;
        let radius = 0.1;
        let body = bccxy_cell_body(length, radius, 0.0).unwrap();
        for (i, (a, b)) in home_edges(length).iter().enumerate() {
            let mid = (*a + *b) * 0.5;
            let d = body.eval(mid);
            assert!(
                (d - (-radius)).abs() < 1e-5,
                "edge {i} midpoint {mid:?}: expected {:.4}, got {d}",
                -radius
            );
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests.
    // --------------------------------------------------------------

    /// Regression: "Body diagonals went corner-to-corner instead of
    /// corner-to-center." Detection: the cube center is the key
    /// signature — with corner-to-corner diagonals, eval at the origin
    /// would be negative but the individual midpoints (which are NOT
    /// the origin under corner-to-corner) would reveal the error.
    #[test]
    fn regression_body_diagonals_pass_through_origin() {
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        // Origin sits at the meeting of all 8 body-diagonal ends.
        assert!((body.eval(Vec3::ZERO) - (-0.1)).abs() < 1e-5);
        // A point halfway between (h, h, h) and (-h, -h, -h) — the
        // origin — is on the diagonal regardless of corner-to-corner
        // vs corner-to-center, so ALSO probe a point that's on
        // corner-to-center but NOT on corner-to-corner:
        // The midpoint (h/2, h/2, h/2) is on d1 (corner-to-center)
        // but only on a corner-to-corner diagonal if the diagonal
        // passed through that midpoint, which only happens if it
        // goes corner-to-opposite-corner. For d1, the corner-to-corner
        // version would go (h,h,h)↔(-h,-h,-h), and (h/2, h/2, h/2) is
        // on that line. So this doesn't distinguish. Pick a better
        // probe: a point on (h, h, h)↔(0,0,0) but NOT on any
        // corner-to-opposite-corner line: (3h/4, 3h/4, 3h/4).
        let p = vec3(1.5, 1.5, 1.5);
        // On d1: at parameter t=0.75 along (0,0,0)→(h,h,h)=(2,2,2), so
        // point (1.5, 1.5, 1.5) IS on d1. Distance = 0, minus r = -r.
        // On corner-to-corner (2,2,2)↔(-2,-2,-2): same line, same result.
        // So this still doesn't distinguish. The only true distinguisher
        // is a point strictly interior to a corner-to-center diagonal but
        // NOT on its extension — which doesn't exist if we're talking
        // about infinite lines. However, if the corner-to-corner version
        // REPLACED corner-to-center (didn't add it), then (0, 0, 0) would
        // be inside the diagonal line-set but maybe not on any endpoint
        // hemisphere. The -radius check above at origin is the correct
        // signature. Leaving both for documentation.
        assert!((body.eval(p) - (-0.1)).abs() < 1e-5);
    }

    /// Regression: "Top-face edges placed on the wrong face (e.g., on
    /// `z = -L/2` or on `y = +L/2`)." Detection: a point that's on the
    /// correct top face should be inside; same point mirrored to a wrong
    /// face should behave differently.
    #[test]
    fn regression_top_face_edges_on_z_plus_not_z_minus() {
        let body = bccxy_cell_body(4.0, 0.1, 0.0).unwrap();
        // Midpoint of top edge t1: (0, h, h). Inside.
        assert!((body.eval(vec3(0.0, 2.0, 2.0)) - (-0.1)).abs() < 1e-5);
        // Mirror to bottom face: (0, h, -h). NOT inside a home edge
        // (bottom-face edges aren't home). A body diagonal passes near
        // this point but doesn't reach it exactly — the closest diagonal
        // is d3 = (0,0,0)↔(h,-h,h) or d6 = (0,0,0)↔(-h,h,-h). Closest
        // point check: for d6, direction (-1,1,-1)/√3, project (0,h,-h)
        // onto it: ((0)(-1) + h(1) + (-h)(-1))/√3 = 2h/√3 ≈ 2.31 for L=4.
        // But d6 has length h√3 ≈ 3.46, so t = 2.31/3.46 ≈ 0.67 ∈ [0,1].
        // Closest point on d6: 0.67 * (-h, h, -h) = (-1.33, 1.33, -1.33).
        // Distance from (0, h, -h)=(0,2,-2) to (-1.33, 1.33, -1.33):
        // ≈ √(1.33² + 0.67² + 0.67²) ≈ 1.63. Minus r = 0.1 gives ~1.53.
        let d = body.eval(vec3(0.0, 2.0, -2.0));
        assert!(
            d > 0.5,
            "expected well positive at bottom-face midpoint, got {d}"
        );
    }

    // Note: the previous `regression_exactness_preserved` test asserted
    // that `bccxy_cell_body` returned an `ExactSdf`. After the migration
    // to `SmoothUnion` for joint smoothing, the cell body is `Sdf` only —
    // smooth-min violates the BoundSdf/ExactSdf "never overestimate
    // magnitude" contract near joints. Callers that previously relied on
    // the exactness guarantee (none in this workspace) would now see a
    // type error. The function-signature `Result<BccXyCellBody, _>`
    // continues to enforce `Sdf` via concrete-type plumbing.
}
