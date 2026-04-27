//! Kelvin cell topology — truncated octahedron strut skeleton.
//!
//! # Cell geometry
//!
//! The truncated octahedron (TO) is inscribed in the origin-centered cube
//! `[-L/2, L/2]^3`. Its 24 vertices sit at every sign/axis permutation of
//! `(0, ±L/4, ±L/2)` — one coordinate zero, one at `±L/4`, one at `±L/2`.
//! Adjacency rule: two vertices are connected by an edge iff they differ
//! in exactly two of the three coordinates, each by `L/4`. Every edge has
//! the same length `L/(2·√2)`.
//!
//! # Home edge-set: all 36 edges
//!
//! Every cell contains the full 36-edge TO skeleton in its local frame:
//!
//! - **24 face-edges** — 4 edges forming the TO's square face embedded in
//!   each of the cube's 6 faces (`x, y, z = ±L/2`).
//! - **12 interior edges** — 4 edges per coordinate mid-plane (`x=0`,
//!   `y=0`, `z=0`), connecting vertices on perpendicular cube faces.
//!
//! Every face-edge is shared across the boundary with an adjacent cell,
//! but both cells represent the same edge in their respective local
//! frames (one as its `+face`, the other as its `-face`). Under
//! [`sdf::LimitedRepeat`], a query folds into exactly one cell frame and
//! evaluates once — the "duplicate" representation is never doubly
//! counted at evaluation time. Not doing this would produce missing
//! struts at cell boundaries because `round()`-based folding can send a
//! boundary point into a neighbor cell that didn't home-assign the edge.
//!
//! This mirrors the Cubic convention: its 3 axis struts span the full
//! period `[-L/2, L/2]` on their axis and appear identically in every
//! cell, for the same reason.

use glam::vec3;
use sdf::{Capsule, SmoothUnion};

use crate::error::LatticeError;

/// Concrete nested-`SmoothUnion` type for the Kelvin cell body (36 capsules).
///
/// Written out fully — 35 levels of right-folded `SmoothUnion` — so inlining
/// stays visible to the compiler. With `joint_smoothness = 0`, each
/// `SmoothUnion` falls back to bit-exact hard `min` (see [`SmoothUnion`]'s
/// `k = 0` fallback) so existing behavior is unchanged for callers that
/// don't opt into smoothing.
#[allow(clippy::type_complexity)]
pub(crate) type KelvinCellBody = SmoothUnion<
    Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule,
    SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule,
    SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule,
    SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule,
    SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule,
    SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, SmoothUnion<Capsule, Capsule>
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>;

/// Builds the cell-body SDF for the Kelvin topology.
///
/// Constructs 36 `Capsule`s per the edge-set listed in the module docs and
/// right-folds them into the [`KelvinCellBody`] nested-Union chain.
///
/// # Errors
///
/// Propagates any underlying [`sdf::BuildError`] from [`Capsule::new`]. In
/// practice, with `length > 0` and `radius > 0` validated upstream, the
/// only way this fails is a bypass of [`crate::LatticeJob::new`]'s
/// invariants — defensive path.
#[allow(clippy::too_many_lines)]
pub(crate) fn kelvin_cell_body(
    length: f32,
    radius: f32,
    joint_smoothness: f32,
) -> Result<KelvinCellBody, LatticeError> {
    let h = length * 0.5;
    let q = length * 0.25;
    let k = joint_smoothness;

    let edges = home_edges_inner(h, q);
    let mut caps: [Option<Capsule>; 36] = [None; 36];
    for (i, (a, b)) in edges.iter().enumerate() {
        caps[i] = Some(Capsule::new(*a, *b, radius)?);
    }
    // SAFETY-equivalent assertion (no unsafe used): every slot was just
    // populated by the loop above. The `expect` here is unreachable —
    // it documents the invariant without risking a silent panic path.
    #[allow(clippy::expect_used)]
    let c = caps.map(|x| x.expect("populated by loop above"));

    Ok(SmoothUnion {
        a: c[0], k,
        b: SmoothUnion {
            a: c[1], k,
            b: SmoothUnion {
                a: c[2], k,
                b: SmoothUnion {
                    a: c[3], k,
                    b: SmoothUnion {
                        a: c[4], k,
                        b: SmoothUnion {
                            a: c[5], k,
                            b: SmoothUnion {
                                a: c[6], k,
                                b: SmoothUnion {
                                    a: c[7], k,
                                    b: SmoothUnion {
                                        a: c[8], k,
                                        b: SmoothUnion {
                                            a: c[9], k,
                                            b: SmoothUnion {
                                                a: c[10], k,
                                                b: SmoothUnion {
                                                    a: c[11], k,
                                                    b: SmoothUnion {
                                                        a: c[12], k,
                                                        b: SmoothUnion {
                                                            a: c[13], k,
                                                            b: SmoothUnion {
                                                                a: c[14], k,
                                                                b: SmoothUnion {
                                                                    a: c[15], k,
                                                                    b: SmoothUnion {
                                                                        a: c[16], k,
                                                                        b: SmoothUnion {
                                                                            a: c[17], k,
                                                                            b: SmoothUnion {
                                                                                a: c[18], k,
                                                                                b: SmoothUnion {
                                                                                    a: c[19], k,
                                                                                    b: SmoothUnion {
                                                                                        a: c[20], k,
                                                                                        b: SmoothUnion {
                                                                                            a: c[21], k,
                                                                                            b: SmoothUnion {
                                                                                                a: c[22], k,
                                                                                                b: SmoothUnion {
                                                                                                    a: c[23], k,
                                                                                                    b: SmoothUnion {
                                                                                                        a: c[24], k,
                                                                                                        b: SmoothUnion {
                                                                                                            a: c[25], k,
                                                                                                            b: SmoothUnion {
                                                                                                                a: c[26], k,
                                                                                                                b: SmoothUnion {
                                                                                                                    a: c[27], k,
                                                                                                                    b: SmoothUnion {
                                                                                                                        a: c[28], k,
                                                                                                                        b: SmoothUnion {
                                                                                                                            a: c[29], k,
                                                                                                                            b: SmoothUnion {
                                                                                                                                a: c[30], k,
                                                                                                                                b: SmoothUnion {
                                                                                                                                    a: c[31], k,
                                                                                                                                    b: SmoothUnion {
                                                                                                                                        a: c[32], k,
                                                                                                                                        b: SmoothUnion {
                                                                                                                                            a: c[33], k,
                                                                                                                                            b: SmoothUnion {
                                                                                                                                                a: c[34], k,
                                                                                                                                                b: c[35],
                                                                                                                                            },
                                                                                                                                        },
                                                                                                                                    },
                                                                                                                                },
                                                                                                                            },
                                                                                                                        },
                                                                                                                    },
                                                                                                                },
                                                                                                            },
                                                                                                        },
                                                                                                    },
                                                                                                },
                                                                                            },
                                                                                        },
                                                                                    },
                                                                                },
                                                                            },
                                                                        },
                                                                    },
                                                                },
                                                            },
                                                        },
                                                    },
                                                },
                                            },
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

/// The 36 home-edge endpoints as `(a, b)` pairs, given already-split
/// half/quarter values. Factored out so both the builder and tests read
/// from the same table.
#[inline]
fn home_edges_inner(h: f32, q: f32) -> [(glam::Vec3, glam::Vec3); 36] {
    [
        // Face-edges on x = +h (4).
        (vec3(h, q, 0.0), vec3(h, 0.0, q)),
        (vec3(h, 0.0, q), vec3(h, -q, 0.0)),
        (vec3(h, -q, 0.0), vec3(h, 0.0, -q)),
        (vec3(h, 0.0, -q), vec3(h, q, 0.0)),
        // Face-edges on x = -h (4).
        (vec3(-h, q, 0.0), vec3(-h, 0.0, q)),
        (vec3(-h, 0.0, q), vec3(-h, -q, 0.0)),
        (vec3(-h, -q, 0.0), vec3(-h, 0.0, -q)),
        (vec3(-h, 0.0, -q), vec3(-h, q, 0.0)),
        // Face-edges on y = +h (4).
        (vec3(q, h, 0.0), vec3(0.0, h, q)),
        (vec3(0.0, h, q), vec3(-q, h, 0.0)),
        (vec3(-q, h, 0.0), vec3(0.0, h, -q)),
        (vec3(0.0, h, -q), vec3(q, h, 0.0)),
        // Face-edges on y = -h (4).
        (vec3(q, -h, 0.0), vec3(0.0, -h, q)),
        (vec3(0.0, -h, q), vec3(-q, -h, 0.0)),
        (vec3(-q, -h, 0.0), vec3(0.0, -h, -q)),
        (vec3(0.0, -h, -q), vec3(q, -h, 0.0)),
        // Face-edges on z = +h (4).
        (vec3(q, 0.0, h), vec3(0.0, q, h)),
        (vec3(0.0, q, h), vec3(-q, 0.0, h)),
        (vec3(-q, 0.0, h), vec3(0.0, -q, h)),
        (vec3(0.0, -q, h), vec3(q, 0.0, h)),
        // Face-edges on z = -h (4).
        (vec3(q, 0.0, -h), vec3(0.0, q, -h)),
        (vec3(0.0, q, -h), vec3(-q, 0.0, -h)),
        (vec3(-q, 0.0, -h), vec3(0.0, -q, -h)),
        (vec3(0.0, -q, -h), vec3(q, 0.0, -h)),
        // Interior edges at z = 0 (4).
        (vec3(h, q, 0.0), vec3(q, h, 0.0)),
        (vec3(h, -q, 0.0), vec3(q, -h, 0.0)),
        (vec3(-h, q, 0.0), vec3(-q, h, 0.0)),
        (vec3(-h, -q, 0.0), vec3(-q, -h, 0.0)),
        // Interior edges at y = 0 (4).
        (vec3(h, 0.0, q), vec3(q, 0.0, h)),
        (vec3(h, 0.0, -q), vec3(q, 0.0, -h)),
        (vec3(-h, 0.0, q), vec3(-q, 0.0, h)),
        (vec3(-h, 0.0, -q), vec3(-q, 0.0, -h)),
        // Interior edges at x = 0 (4).
        (vec3(0.0, h, q), vec3(0.0, q, h)),
        (vec3(0.0, h, -q), vec3(0.0, q, -h)),
        (vec3(0.0, -h, q), vec3(0.0, -q, h)),
        (vec3(0.0, -h, -q), vec3(0.0, -q, -h)),
    ]
}

/// The 36-element home edge set as endpoint pairs, for test oracles.
#[cfg(test)]
pub(crate) fn home_edges(length: f32) -> [(glam::Vec3, glam::Vec3); 36] {
    home_edges_inner(length * 0.5, length * 0.25)
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
    // Target failure class: wrong vertex coords, missing edges, wrong
    // edge endpoints.
    // --------------------------------------------------------------

    #[test]
    fn eval_at_vertex_plus_h_plus_q_zero_is_negative_radius() {
        // Vertex (L/2, L/4, 0) is the endpoint of three edges.
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(2.0, 1.0, 0.0));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn eval_at_vertex_plus_h_zero_plus_q_is_negative_radius() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(2.0, 0.0, 1.0));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn eval_at_vertex_plus_q_plus_h_zero_is_negative_radius() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(1.0, 2.0, 0.0));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn eval_on_face_edge_midpoint_is_negative_radius() {
        // Midpoint of (h, q, 0) ↔ (h, 0, q) = (h, q/2, q/2) = (2, 0.5, 0.5).
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(2.0, 0.5, 0.5));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn eval_on_negative_face_edge_midpoint_is_negative_radius() {
        // Midpoint of (-h, q, 0) ↔ (-h, 0, q) = (-h, q/2, q/2) = (-2, 0.5, 0.5).
        // This was previously not home-assigned; now it is.
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(-2.0, 0.5, 0.5));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn eval_on_interior_edge_midpoint_is_negative_radius() {
        // Midpoint of (h, q, 0) ↔ (q, h, 0) = (3L/8, 3L/8, 0) = (1.5, 1.5, 0).
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let d = body.eval(vec3(1.5, 1.5, 0.0));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn eval_at_cell_center_is_positive_pore() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        assert!(body.eval(Vec3::ZERO) > 0.0);
    }

    #[test]
    fn eval_at_hex_face_centroid_is_positive_pore() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        let p = vec3(1.5, 1.5, 1.5);
        assert!(
            body.eval(p) > 0.0,
            "expected positive, got {}",
            body.eval(p)
        );
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests.
    // --------------------------------------------------------------

    #[test]
    fn construction_propagates_sdf_error_on_zero_radius() {
        assert!(kelvin_cell_body(4.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn construction_propagates_sdf_error_on_negative_radius() {
        assert!(kelvin_cell_body(4.0, -0.1, 0.0).is_err());
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests.
    // --------------------------------------------------------------

    /// Every home-edge midpoint should evaluate to exactly -radius.
    /// One typo in any of the 36 edge entries would desynchronize and
    /// this test would fail for that index.
    #[test]
    fn all_36_home_edges_have_negative_radius_at_midpoint() {
        let length = 4.0;
        let radius = 0.1;
        let body = kelvin_cell_body(length, radius, 0.0).unwrap();
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

    /// Every home-edge endpoint should also evaluate to -radius (it's a
    /// node where multiple home edges meet).
    #[test]
    fn all_36_home_edge_endpoints_evaluate_negative_radius() {
        let length = 4.0;
        let radius = 0.1;
        let body = kelvin_cell_body(length, radius, 0.0).unwrap();
        for (i, (a, b)) in home_edges(length).iter().enumerate() {
            let da = body.eval(*a);
            let db = body.eval(*b);
            assert!(
                (da - (-radius)).abs() < 1e-5,
                "edge {i} start {a:?}: expected {:.4}, got {da}",
                -radius
            );
            assert!(
                (db - (-radius)).abs() < 1e-5,
                "edge {i} end {b:?}: expected {:.4}, got {db}",
                -radius
            );
        }
    }

    /// At the topology threshold r < L/(4√2), the bisector midpoint
    /// between two parallel square-face edges should still be positive
    /// (the binding non-collision case).
    #[test]
    fn parallel_square_face_edges_dont_collide_at_threshold() {
        let length = 4.0;
        // r_max = L/(4√2) ≈ 0.7071, take a bit under.
        let r = 0.7;
        let body = kelvin_cell_body(length, r, 0.0).unwrap();
        // Bisector between face-edge 1 and face-edge 3 on +x face = (+h, 0, 0).
        let bisector = vec3(length * 0.5, 0.0, 0.0);
        let d = body.eval(bisector);
        assert!(d > 0.0, "expected positive at bisector, got {d}");
    }

    // --------------------------------------------------------------
    // d. Regression tests.
    // --------------------------------------------------------------

    /// Regression: "Only positive-half faces were home-assigned; negative-
    /// half faces relied on periodic tiling, but `LimitedRepeat` folds via
    /// `round()` which can send cell-boundary queries to the neighbor that
    /// did NOT home-assign the edge — producing a spurious positive SDF
    /// right on a strut." Detection: probe a -x face-edge midpoint.
    #[test]
    fn regression_negative_half_faces_are_home_assigned() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        // (-h, q/2, q/2) is the midpoint of a -x face-edge.
        let d = body.eval(vec3(-2.0, 0.5, 0.5));
        assert!(
            (d - (-0.1)).abs() < 1e-5,
            "expected -0.1 at -x face-edge midpoint, got {d}"
        );
    }

    /// Regression: "Interior edge endpoints used corner-frame `[0, L]`
    /// coordinates instead of origin-centered `[-L/2, L/2]`."
    #[test]
    fn regression_interior_edges_use_signed_coordinates() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        // (-1.5, -1.5, 0) is the midpoint of an interior edge at z = 0.
        let d = body.eval(vec3(-1.5, -1.5, 0.0));
        assert!(
            (d - (-0.1)).abs() < 1e-5,
            "expected -0.1 at interior edge midpoint, got {d}"
        );
    }

    /// Regression: "Face-edge table used `±L/2` in place of `±L/4` for
    /// the in-face vertex coordinate — i.e., the small TO square was
    /// confused with the full cube face."
    #[test]
    fn regression_face_edges_use_quarter_not_half_for_in_face_coord() {
        let body = kelvin_cell_body(4.0, 0.1, 0.0).unwrap();
        // Cube corner — NOT a TO vertex.
        let d = body.eval(vec3(2.0, 2.0, 2.0));
        assert!(d > 0.0, "expected positive at cube corner, got {d}");
    }

    // Note: a previous regression test asserted `kelvin_cell_body`
    // returned an `ExactSdf`. That guarantee was deliberately dropped
    // when the cell body migrated from nested `Union` to nested
    // `SmoothUnion` for joint smoothing — `SmoothUnion` is `Sdf` only.
}
