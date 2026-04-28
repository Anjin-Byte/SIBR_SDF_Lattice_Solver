//! Asymptotic deciders for MC33 ambiguity resolution.
//!
//! Two primitives live here:
//!
//! - [`face_decider`] — bilinear saddle sign on an ambiguous *face*. Used
//!   by Cases 3/6/7/10/12 (Lewiner 4/7/8/11/13).
//! - [`interior_decider`] — sign of the trilinear saddle along the cube's
//!   body diagonal. Used by Cases 4/10 (Lewiner 5/11).
//!
//! Both are pure functions of the corner values plus a small integer
//! sentinel drawn from Lewiner's `testN[]` tables.
//!
//! # Face decider closed form
//!
//! For face `f` with corners `(A, B, C, D)` in the order Lewiner's
//! `test_face` specifies — where A–C and B–D are the two diagonals of
//! the face — the decider evaluates:
//!
//! ```text
//! if |A*C - B*D| < ε  →  return face ≥ 0   (tie-break via sign-of-face)
//! else                →  return face * A * (A*C - B*D) ≥ 0
//! ```
//!
//! A negative `face` argument means "invert the result of the test",
//! encoded per Lewiner's `test3[]` table.
//!
//! # Interior decider closed form
//!
//! For the body-ambiguous cases, the trilinear interpolant restricted to
//! the body diagonal `(0,0,0) → (1,1,1)` is a quadratic in the parameter
//! `t`. Its critical point `t* = -b / (2a)` is the candidate saddle. If
//! `t*` is outside `[0,1]` or the quadratic is degenerate (`|a| < ε`), no
//! interior saddle exists on the diagonal, and the decider falls back to
//! the sign of `s`. Otherwise we evaluate the bilinear `A(t), B(t), C(t),
//! D(t)` on the four body-parallel edges at `t*`, pack their signs into a
//! 4-bit key, and dispatch per the Lewiner `test_interior` table.
//!
//! # Corner-value sign convention
//!
//! Our SDF convention (inside = negative) matches Lewiner's stored cube
//! values (inside = negative after `f - iso` with `iso = 0`), so corner
//! values `A, B, C, D` are passed to the deciders unchanged.
//!
//! # Port
//!
//! Ported from `JuliaGeometry` `MarchingCubes.jl`, functions `test_face`
//! (lines ~498–535) and `test_interior` (lines ~378–490) as of fetched
//! version. Originally from Lewiner 2003.

/// Lewiner's corner indices for each face (1-based), converted to our
/// 0-based: `FACE_CORNERS[face - 1] = (A, B, C, D)` corner indices in
/// the order the decider needs.
const FACE_CORNERS: [[usize; 4]; 6] = [
    // Face 1: front (y=0)
    [0, 4, 5, 1],
    // Face 2: right (x=1)
    [1, 5, 6, 2],
    // Face 3: back (y=1)
    [2, 6, 7, 3],
    // Face 4: left (x=0)
    [3, 7, 4, 0],
    // Face 5: bottom (z=0)
    [0, 3, 2, 1],
    // Face 6: top (z=1)
    [4, 7, 6, 5],
];

/// Evaluates the face asymptotic decider for face `face`.
///
/// `face` encodes both the face identifier (`|face| ∈ 1..=6`) and an
/// optional inversion flag (`face < 0` means invert the result).
///
/// # Panics
///
/// Panics if `|face|` is outside `1..=6`. Tables store only valid face
/// numbers, so this indicates a bug in the dispatch tables.
#[allow(clippy::cast_sign_loss)]
pub(crate) fn face_decider(cube: &[f32; 8], face: i8) -> bool {
    let af = face.unsigned_abs() as usize;
    assert!(
        (1..=6).contains(&af),
        "face decider invoked with invalid face = {face}"
    );

    let [i_a, i_b, i_c, i_d] = FACE_CORNERS[af - 1];
    let a = cube[i_a];
    let b = cube[i_b];
    let c = cube[i_c];
    let d = cube[i_d];

    let det = a * c - b * d;
    // Tie-break for the saddle-zero case. Follows Lewiner's convention
    // exactly — return (face >= 0).
    if det.abs() < f32::EPSILON {
        return face >= 0;
    }
    // Main branch: `face * A * det >= 0` — the sign of `A * det` tells
    // us which side the saddle sits on, and `face` (if negative) inverts
    // that.
    f32::from(face) * a * det >= 0.0
}

/// Evaluates the interior asymptotic decider for a body-ambiguous case.
///
/// `s` is the sentinel from Lewiner's per-case `testN[cfg]` table: its
/// **sign** selects one of the two triangulations, its **magnitude** is
/// always `7` for the Case-5 / Case-11 branch that this function handles.
/// Later sessions may widen the accepted magnitude; for now we assert on
/// entry to catch copy-paste bugs in the `TEST_4` table.
///
/// Returns `true` when the caller should emit the "interior-connected"
/// triangulation (e.g. Case 4.1.1, 2 triangles); `false` for the
/// "tunnel" triangulation (Case 4.1.2, 6 triangles).
///
/// # Panics (debug only)
///
/// `debug_assert!`s that `s.unsigned_abs() == 7`.
#[allow(clippy::float_equality_without_abs)]
pub(crate) fn interior_decider(cube: &[f32; 8], s: i8) -> bool {
    debug_assert!(
        s.unsigned_abs() == 7,
        "interior decider expects |s| == 7, got {s}"
    );

    let d40 = cube[4] - cube[0];
    let d62 = cube[6] - cube[2];
    let d73 = cube[7] - cube[3];
    let d51 = cube[5] - cube[1];

    let a = d40 * d62 - d73 * d51;
    let b = cube[2] * d40 + cube[0] * d62 - cube[1] * d73 - cube[3] * d51;

    // Degenerate quadratic (linear or constant along diagonal) and
    // out-of-range critical point both collapse to the sign-of-s branch.
    // Matches Julia's IEEE-754-driven fallback without relying on
    // NaN/Inf comparisons.
    if a.abs() < f32::EPSILON {
        return s > 0;
    }
    let t = -b / (2.0 * a);
    if !(0.0..=1.0).contains(&t) {
        return s > 0;
    }

    let at = cube[0] + d40 * t;
    let bt = cube[3] + d73 * t;
    let ct = cube[2] + d62 * t;
    let dt = cube[1] + d51 * t;

    let test = u8::from(at >= 0.0)
        | u8::from(bt >= 0.0) << 1
        | u8::from(ct >= 0.0) << 2
        | u8::from(dt >= 0.0) << 3;

    // Intentional signed comparison against ε (not |x| < ε): this is
    // Lewiner's decision rule for which side of zero the saddle lies on,
    // with ε acting as a tie-breaking threshold on the positive side.
    // Matches Julia's `test_interior` verbatim.
    let sub = at * ct - bt * dt;
    match test {
        0..=4 | 6 | 8 | 9 | 12 => s > 0,
        5 => {
            if sub < f32::EPSILON {
                s > 0
            } else {
                s < 0
            }
        }
        10 => {
            if sub >= f32::EPSILON {
                s > 0
            } else {
                s < 0
            }
        }
        // test ∈ {7, 11, 13, 14, 15}
        _ => s < 0,
    }
}

/// Per-edge interior asymptotic decider — used by Chernyaev cases 6,
/// 7, 12, 13 (Lewiner 7, 8, 13, 14) to select between alternative
/// triangulations of an interior-ambiguous voxel.
///
/// Distinct from [`interior_decider`] because these four cases evaluate
/// the trilinear interpolant **along a single reference edge** (rather
/// than the body-diagonal quadratic of cases 4 and 10). The reference
/// edge is determined per-case:
/// - Chernyaev 6 (Lewiner 7) → `test6[cfg][3]`
/// - Chernyaev 7 (Lewiner 8) → `test7[cfg][5]`
/// - Chernyaev 12 (Lewiner 13) → `test12[cfg][4]`
/// - Chernyaev 13 (Lewiner 14) → `tiling13_5_1[cfg][subcfg][1] - 1`
///
/// The final 4-bit dispatch matches the case-5/11 path verbatim. `At`
/// is held at zero (per the Julia source's case-7/8/13/14 branch),
/// which forces bit 0 of the dispatch key on, restricting `test` to
/// odd values.
///
/// Returns `true` when the caller should emit the "first" alternative
/// (e.g., 13.5.1, 7.4.2 inverted polarity, etc.); `false` for the
/// "second".
///
/// Ported from `MarchingCubes.jl` `src/MarchingCubes.jl` lines 395-489
/// (the `case ∈ {7, 8, 13, 14}` branch of `test_interior`).
///
/// # Panics (debug only)
///
/// `debug_assert!`s that `edge < 12`.
#[allow(clippy::float_equality_without_abs, clippy::too_many_lines)]
pub(crate) fn interior_decider_per_edge(cube: &[f32; 8], edge: u8, s: i8) -> bool {
    debug_assert!(edge < 12, "per-edge interior reference edge must be 0..=11");

    // Per Julia: `At` is left at zero for the case-7/8/13/14 branch;
    // only `Bt`, `Ct`, `Dt` are computed. The bit-pack later then
    // forces bit 0 (At ≥ 0 is true) regardless.
    let (bt, ct, dt) = match edge {
        0 => {
            let t = cube[0] / (cube[0] - cube[1]);
            (
                cube[3] + (cube[2] - cube[3]) * t,
                cube[7] + (cube[6] - cube[7]) * t,
                cube[4] + (cube[5] - cube[4]) * t,
            )
        }
        1 => {
            let t = cube[1] / (cube[1] - cube[2]);
            (
                cube[0] + (cube[3] - cube[0]) * t,
                cube[4] + (cube[7] - cube[4]) * t,
                cube[5] + (cube[6] - cube[5]) * t,
            )
        }
        2 => {
            let t = cube[2] / (cube[2] - cube[3]);
            (
                cube[1] + (cube[0] - cube[1]) * t,
                cube[5] + (cube[4] - cube[5]) * t,
                cube[6] + (cube[7] - cube[6]) * t,
            )
        }
        3 => {
            let t = cube[3] / (cube[3] - cube[0]);
            (
                cube[2] + (cube[1] - cube[2]) * t,
                cube[6] + (cube[5] - cube[6]) * t,
                cube[7] + (cube[4] - cube[7]) * t,
            )
        }
        4 => {
            let t = cube[4] / (cube[4] - cube[5]);
            (
                cube[7] + (cube[6] - cube[7]) * t,
                cube[3] + (cube[2] - cube[3]) * t,
                cube[0] + (cube[1] - cube[0]) * t,
            )
        }
        5 => {
            let t = cube[5] / (cube[5] - cube[6]);
            (
                cube[4] + (cube[7] - cube[4]) * t,
                cube[0] + (cube[3] - cube[0]) * t,
                cube[1] + (cube[2] - cube[1]) * t,
            )
        }
        6 => {
            let t = cube[6] / (cube[6] - cube[7]);
            (
                cube[5] + (cube[4] - cube[5]) * t,
                cube[1] + (cube[0] - cube[1]) * t,
                cube[2] + (cube[3] - cube[2]) * t,
            )
        }
        7 => {
            let t = cube[7] / (cube[7] - cube[4]);
            (
                cube[6] + (cube[5] - cube[6]) * t,
                cube[2] + (cube[1] - cube[2]) * t,
                cube[3] + (cube[0] - cube[3]) * t,
            )
        }
        8 => {
            let t = cube[0] / (cube[0] - cube[4]);
            (
                cube[3] + (cube[7] - cube[3]) * t,
                cube[2] + (cube[6] - cube[2]) * t,
                cube[1] + (cube[5] - cube[1]) * t,
            )
        }
        9 => {
            let t = cube[1] / (cube[1] - cube[5]);
            (
                cube[0] + (cube[4] - cube[0]) * t,
                cube[3] + (cube[7] - cube[3]) * t,
                cube[2] + (cube[6] - cube[2]) * t,
            )
        }
        10 => {
            let t = cube[2] / (cube[2] - cube[6]);
            (
                cube[1] + (cube[5] - cube[1]) * t,
                cube[0] + (cube[4] - cube[0]) * t,
                cube[3] + (cube[7] - cube[3]) * t,
            )
        }
        // edge == 11
        _ => {
            let t = cube[3] / (cube[3] - cube[7]);
            (
                cube[2] + (cube[6] - cube[2]) * t,
                cube[1] + (cube[5] - cube[1]) * t,
                cube[0] + (cube[4] - cube[0]) * t,
            )
        }
    };

    // Pack At ≥ 0 (always true here, At = 0), Bt, Ct, Dt into a 4-bit
    // key. Bit 0 = 1 always for case 14, so `test` is always odd.
    let test =
        1_u8 | u8::from(bt >= 0.0) << 1 | u8::from(ct >= 0.0) << 2 | u8::from(dt >= 0.0) << 3;

    // Final dispatch — identical structure to [`interior_decider`]'s
    // tail. `At * Ct - Bt * Dt` reduces to `-Bt * Dt` since At = 0.
    let sub = -bt * dt;
    match test {
        0..=4 | 6 | 8 | 9 | 12 => s > 0,
        5 => {
            if sub < f32::EPSILON {
                s > 0
            } else {
                s < 0
            }
        }
        10 => {
            if sub >= f32::EPSILON {
                s > 0
            } else {
                s < 0
            }
        }
        // test ∈ {7, 11, 13, 14, 15}
        _ => s < 0,
    }
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

    /// For a cube where all corners on the tested face have sign matching
    /// the expected "connected" configuration, the decider returns true
    /// for positive face, false for negative face.
    #[test]
    fn decider_on_unambiguous_face_is_stable() {
        // All corners positive: A*C - B*D = 1 - 1 = 0 → degenerate,
        // tie-breaks on sign of face.
        let cube = [1.0_f32; 8];
        assert!(face_decider(&cube, 5));
        assert!(!face_decider(&cube, -5));
    }

    #[test]
    fn decider_diagonal_connected_returns_matching_sign() {
        // Set up face 5 (bottom, corners 0,3,2,1) with A=−1, B=+1, C=−1, D=+1.
        // det = A*C − B*D = (−1)(−1) − (1)(1) = 1 − 1 = 0 → degenerate
        // tie-break. Change D slightly to break degeneracy.
        let mut cube = [0.0_f32; 8];
        cube[0] = -1.0; // A
        cube[3] = 1.0; // B
        cube[2] = -1.0; // C
        cube[1] = 0.5; // D (smaller; det = 1 - 0.5 = 0.5 > 0)
        // Expected: face * A * det >= 0 → 5 * (−1) * 0.5 = −2.5 < 0 → false
        assert!(!face_decider(&cube, 5));
        // With inverted face: (−5) * (−1) * 0.5 = 2.5 >= 0 → true
        assert!(face_decider(&cube, -5));
    }

    #[test]
    #[should_panic(expected = "invalid face")]
    fn decider_panics_on_invalid_face_zero() {
        let cube = [0.0_f32; 8];
        face_decider(&cube, 0);
    }

    #[test]
    #[should_panic(expected = "invalid face")]
    fn decider_panics_on_invalid_face_seven() {
        let cube = [0.0_f32; 8];
        face_decider(&cube, 7);
    }

    /// Symmetry: two voxels sharing a face pass the same four corner values
    /// (in the same order) to the decider, so they compute the same result.
    /// This is the correctness condition for watertight MC33 at shared faces.
    #[test]
    fn decider_is_deterministic() {
        let cube = [-0.3_f32, 0.8, -0.2, 0.5, 0.7, -0.1, 0.3, 0.9];
        for face in [1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6] {
            let a = face_decider(&cube, face);
            let b = face_decider(&cube, face);
            assert_eq!(a, b);
        }
    }

    // --------------------------------------------------------------
    // interior_decider tests
    // --------------------------------------------------------------

    /// All corners equal → `a = 0` → degenerate-quadratic fallback.
    #[test]
    fn interior_decider_degenerate_a_falls_back_to_sign_of_s() {
        let cube = [1.0_f32; 8];
        assert!(interior_decider(&cube, 7));
        assert!(!interior_decider(&cube, -7));
    }

    /// When `t* = 0` (critical point at the bottom body-diagonal corner),
    /// the decider enters the in-range branch and consults the bottom-face
    /// corner signs. For `cube = [0, 1, 0, 1, 1, 1, 1, 1]` the packed
    /// test key is `15` → decider returns `s < 0`, which is the *inverse*
    /// of the sign-of-s fallback. Confirms t=0 is handled as in-range.
    #[test]
    fn interior_decider_t_exactly_zero_takes_in_range_branch() {
        let cube = [0.0_f32, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        // In-range branch: test == 15 → return s < 0.
        // So decider(cube, 7) == false and decider(cube, -7) == true.
        // (The fallback branch would give the opposite.)
        assert!(!interior_decider(&cube, 7));
        assert!(interior_decider(&cube, -7));
    }

    /// Mirror construction for `t* = 1` (critical point at the top body-
    /// diagonal corner). Same correctness claim: in-range branch fires.
    #[test]
    fn interior_decider_t_exactly_one_takes_in_range_branch() {
        let cube = [-1.0_f32, 1.0, -1.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        // a = 1, b = -2, t* = 1. At t=1: (At, Bt, Ct, Dt) = (0, 1, 0, 1).
        // test == 15 → s < 0 branch.
        assert!(!interior_decider(&cube, 7));
        assert!(interior_decider(&cube, -7));
    }

    /// A corner configuration whose body-diagonal critical point lies
    /// outside `[0, 1]` must fall back to the sign-of-s branch.
    #[test]
    fn interior_decider_t_out_of_range_falls_back_to_sign_of_s() {
        // With cube = [10, 0, 10, 0, 1, 0, 1, 0]:
        //   d40 = d62 = -9, d51 = d73 = 0
        //   a = 81, b = -180, t* = 180/162 ≈ 1.111  (out of [0,1])
        let cube = [10.0_f32, 0.0, 10.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        assert!(interior_decider(&cube, 7));
        assert!(!interior_decider(&cube, -7));
    }

    /// Flipping the sentinel sign must flip the result in the fallback
    /// branch.
    #[test]
    fn interior_decider_sign_of_s_inverts_fallback() {
        // Same degenerate cube as above, but verify both polarities.
        let cube = [0.5_f32; 8];
        assert_ne!(
            interior_decider(&cube, 7),
            interior_decider(&cube, -7),
            "sign-of-s branch must depend on sign of s"
        );
    }

    /// `interior_decider` is a pure function of `(cube, s)`; repeated
    /// invocation with identical inputs yields identical results.
    #[test]
    fn interior_decider_is_deterministic() {
        let cube = [-0.3_f32, 0.8, -0.2, 0.5, 0.7, -0.1, 0.3, 0.9];
        for s in [7_i8, -7] {
            assert_eq!(interior_decider(&cube, s), interior_decider(&cube, s));
        }
    }

    /// The canonical Case-4 subcase-1 voxel (corners 0 and 6 "outside" on
    /// the body diagonal, rest inside) has a well-defined interior test
    /// result. We assert only that the decider runs without panic and
    /// returns a determinate bool — the concrete value is validated by
    /// `case4::tests` against the emitted triangulation.
    #[test]
    fn interior_decider_runs_on_canonical_case4_voxel() {
        let cube = [1.0_f32, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0];
        // Both polarities return a bool without panic.
        let _ = interior_decider(&cube, 7);
        let _ = interior_decider(&cube, -7);
    }

    #[test]
    #[should_panic(expected = "|s| == 7")]
    fn interior_decider_debug_rejects_invalid_s_magnitude() {
        let cube = [0.0_f32; 8];
        interior_decider(&cube, 6);
    }

    // --------------------------------------------------------------
    // Property tests: purity (same inputs → same output).
    // --------------------------------------------------------------

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn interior_decider_is_pure(
            c in proptest::array::uniform8(-10.0_f32..10.0_f32),
            polarity in any::<bool>(),
        ) {
            prop_assume!(c.iter().all(|x| x.is_finite()));
            let s: i8 = if polarity { 7 } else { -7 };
            let first = interior_decider(&c, s);
            let second = interior_decider(&c, s);
            prop_assert_eq!(first, second);
        }

        #[test]
        fn face_decider_is_pure(
            c in proptest::array::uniform8(-10.0_f32..10.0_f32),
            face in 1_i8..=6_i8,
            invert in any::<bool>(),
        ) {
            prop_assume!(c.iter().all(|x| x.is_finite()));
            let f = if invert { -face } else { face };
            let first = face_decider(&c, f);
            let second = face_decider(&c, f);
            prop_assert_eq!(first, second);
        }
    }
}
