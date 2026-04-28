//! Case 10 (Chernyaev) / Lewiner base case 11 — "two opposite diagonal pairs".
//!
//! The configuration: four corners with the same sign, arranged as two
//! diagonal pairs on opposite faces. There are 6 voxel configurations
//! (`base_case == 11` rows in [`super::cases::CASES`]).
//!
//! # Dispatch
//!
//! 2 face tests + 1 body-diagonal interior test, branching to 5
//! sub-tilings. The interior test uses [`interior_decider`] (the
//! body-diagonal quadratic, same primitive as Chernyaev 4 / Lewiner 5).
//!
//! | (face1, face2) | Subcase | Triangles | Needs c-vertex? |
//! |----------------|---------|-----------|-----------------|
//! | (true, true)   | 10.1.1_ | 4         | no              |
//! | (true, false)  | 10.2    | 8         | yes             |
//! | (false, true)  | 10.2_   | 8         | yes             |
//! | (false, false) | interior decides 10.1.1 vs 10.1.2 | 4 or 8 | 10.1.2 only |
//!
//! # Sign / index conventions
//!
//! - **Tile edges** are stored 0-based 0..=11; edge index 12 (1-based 13)
//!   refers to the c-vertex computed by [`super::compute_c_vertex`].
//!   Used only by `TILING_10_2` and `TILING_10_2_`.
//! - **`TEST_10` columns**:
//!   - 0/1: face indices (1..=6) for [`face_decider`].
//!   - 2: body-diagonal sentinel `s` (always `7`).
//!
//! # Port
//!
//! Tables transcribed from `JuliaGeometry/MarchingCubes.jl` `src/lut.jl`:
//! - `test10` at lines 1198-1206
//! - `tiling10_1_1` at lines 1217-1234
//! - `tiling10_1_1_` at lines 1236-1253
//! - `tiling10_1_2` at lines 1255-1272
//! - `tiling10_2` at lines 1274-1291
//! - `tiling10_2_` at lines 1293-1310
//!
//! Dispatch mirrors `case == 11` in `src/MarchingCubes.jl`, lines 259-276.

use super::super::{CellCoord, Mesh};
use super::decider::{face_decider, interior_decider};
use super::emit_triangles;
use crate::mesh::grid::GridSpec;

/// Per-cfg 3-tuple: `[face_a, face_b, interior_s]`.
///
/// Source: Lewiner `test10[6][3]`. Cols 0/1 are face indices; col 2 is
/// the body-diagonal sentinel for [`interior_decider`] (always `7`).
#[rustfmt::skip]
const TEST_10: [[i8; 3]; 6] = [
    [  2,   4,   7],
    [  5,   6,   7],
    [  1,   3,   7],
    [  1,   3,   7],
    [  5,   6,   7],
    [  2,   4,   7],
];

/// Triangulation for subcase 10.1.1 (4 triangles, no c-vertex).
///
/// Source: Lewiner `tiling10_1_1[6]`.
#[rustfmt::skip]
const TILING_10_1_1: [[u8; 12]; 6] = [
    [  5,  10,   7,  11,   7,  10,   8,   1,   9,   1,   8,   3],
    [  1,   2,   5,   6,   5,   2,   4,   3,   0,   3,   4,   7],
    [ 11,   0,   8,   0,  11,   2,   4,   9,   6,  10,   6,   9],
    [  9,   0,  10,   2,  10,   0,   6,   8,   4,   8,   6,  11],
    [  7,   2,   3,   2,   7,   6,   0,   1,   4,   5,   4,   1],
    [  7,   9,   5,   9,   7,   8,  10,   1,  11,   3,  11,   1],
];

/// Triangulation for subcase 10.1.1 inverted polarity (4 triangles, no c-vertex).
///
/// Source: Lewiner `tiling10_1_1_[6]`.
#[rustfmt::skip]
const TILING_10_1_1_: [[u8; 12]; 6] = [
    [  5,   9,   7,   8,   7,   9,  11,   1,  10,   1,  11,   3],
    [  3,   2,   7,   6,   7,   2,   4,   1,   0,   1,   4,   5],
    [ 10,   0,   9,   0,  10,   2,   4,   8,   6,  11,   6,   8],
    [  8,   0,  11,   2,  11,   0,   6,   9,   4,   9,   6,  10],
    [  5,   2,   1,   2,   5,   6,   0,   3,   4,   7,   4,   3],
    [  7,  10,   5,  10,   7,  11,   9,   1,   8,   3,   8,   1],
];

/// Triangulation for subcase 10.1.2 (8 triangles, no c-vertex).
///
/// Source: Lewiner `tiling10_1_2[6]`.
#[rustfmt::skip]
const TILING_10_1_2: [[u8; 24]; 6] = [
    [  3,  11,   7,   3,   7,   8,   9,   8,   7,   5,   9,   7,   9,   5,  10,   9,  10,   1,   3,   1,  10,  11,   3,  10],
    [  7,   6,   5,   7,   5,   4,   0,   4,   5,   1,   0,   5,   0,   1,   2,   0,   2,   3,   7,   3,   2,   6,   7,   2],
    [ 11,   2,  10,   6,  11,  10,  11,   6,   4,  11,   4,   8,   0,   8,   4,   9,   0,   4,   0,   9,  10,   0,  10,   2],
    [ 11,   2,  10,  11,  10,   6,   4,   6,  10,   9,   4,  10,   4,   9,   0,   4,   0,   8,  11,   8,   0,   2,  11,   0],
    [  7,   6,   5,   4,   7,   5,   7,   4,   0,   7,   0,   3,   2,   3,   0,   1,   2,   0,   2,   1,   5,   2,   5,   6],
    [  7,   8,   3,  11,   7,   3,   7,  11,  10,   7,  10,   5,   9,   5,  10,   1,   9,  10,   9,   1,   3,   9,   3,   8],
];

/// Triangulation for subcase 10.2 (8 triangles, uses c-vertex slot 12).
///
/// Source: Lewiner `tiling10_2[6]`.
#[rustfmt::skip]
const TILING_10_2: [[u8; 24]; 6] = [
    [ 12,   5,   9,  12,   9,   8,  12,   8,   3,  12,   3,   1,  12,   1,  10,  12,  10,  11,  12,  11,   7,  12,   7,   5],
    [ 12,   1,   0,  12,   0,   4,  12,   4,   7,  12,   7,   3,  12,   3,   2,  12,   2,   6,  12,   6,   5,  12,   5,   1],
    [  4,   8,  12,   6,   4,  12,  10,   6,  12,   9,  10,  12,   0,   9,  12,   2,   0,  12,  11,   2,  12,   8,  11,  12],
    [ 12,   9,   4,  12,   4,   6,  12,   6,  11,  12,  11,   8,  12,   8,   0,  12,   0,   2,  12,   2,  10,  12,  10,   9],
    [  0,   3,  12,   4,   0,  12,   5,   4,  12,   1,   5,  12,   2,   1,  12,   6,   2,  12,   7,   6,  12,   3,   7,  12],
    [ 10,   5,  12,  11,  10,  12,   3,  11,  12,   1,   3,  12,   9,   1,  12,   8,   9,  12,   7,   8,  12,   5,   7,  12],
];

/// Triangulation for subcase 10.2 inverted polarity (8 triangles, c-vertex).
///
/// Source: Lewiner `tiling10_2_[6]`.
#[rustfmt::skip]
const TILING_10_2_: [[u8; 24]; 6] = [
    [  8,   7,  12,   9,   8,  12,   1,   9,  12,   3,   1,  12,  11,   3,  12,  10,  11,  12,   5,  10,  12,   7,   5,  12],
    [  4,   5,  12,   0,   4,  12,   3,   0,  12,   7,   3,  12,   6,   7,  12,   2,   6,  12,   1,   2,  12,   5,   1,  12],
    [ 12,  11,   6,  12,   6,   4,  12,   4,   9,  12,   9,  10,  12,  10,   2,  12,   2,   0,  12,   0,   8,  12,   8,  11],
    [  6,  10,  12,   4,   6,  12,   8,   4,  12,  11,   8,  12,   2,  11,  12,   0,   2,  12,   9,   0,  12,  10,   9,  12],
    [ 12,   7,   4,  12,   4,   0,  12,   0,   1,  12,   1,   5,  12,   5,   6,  12,   6,   2,  12,   2,   3,  12,   3,   7],
    [ 12,   7,  11,  12,  11,  10,  12,  10,   1,  12,   1,   3,  12,   3,   8,  12,   8,   9,  12,   9,   5,  12,   5,   7],
];

/// Emits Case-10 triangles for the given voxel.
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    subcase: i8,
    mesh: &mut Mesh,
) {
    debug_assert!(
        (1..=6).contains(&subcase),
        "Case 10 subcase out of range: {subcase}"
    );
    let cfg = usize::try_from(subcase - 1).unwrap_or(0);
    let f1 = face_decider(corners, TEST_10[cfg][0]);
    let f2 = face_decider(corners, TEST_10[cfg][1]);
    let tiling: &[u8] = match (f1, f2) {
        (true, true) => &TILING_10_1_1_[cfg],
        (true, false) => &TILING_10_2[cfg],
        (false, true) => &TILING_10_2_[cfg],
        (false, false) => {
            let s = TEST_10[cfg][2];
            if interior_decider(corners, s) {
                &TILING_10_1_1[cfg]
            } else {
                &TILING_10_1_2[cfg]
            }
        }
    };
    emit_triangles(grid, cell, corners, tiling, mesh);
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    missing_docs
)]
mod tests {
    use super::*;
    use crate::mesh::Mesh;
    use crate::mesh::grid::GridSpec;
    use glam::{UVec3, Vec3};

    fn unit_grid_with_corners(corners: &[f32; 8]) -> (GridSpec, Vec<f32>) {
        let grid = GridSpec::new(Vec3::ZERO, UVec3::splat(1), 1.0).unwrap();
        let mut field = vec![0.0_f32; grid.sample_count()];
        field[grid.sample_index(0, 0, 0)] = corners[0];
        field[grid.sample_index(1, 0, 0)] = corners[1];
        field[grid.sample_index(1, 1, 0)] = corners[2];
        field[grid.sample_index(0, 1, 0)] = corners[3];
        field[grid.sample_index(0, 0, 1)] = corners[4];
        field[grid.sample_index(1, 0, 1)] = corners[5];
        field[grid.sample_index(1, 1, 1)] = corners[6];
        field[grid.sample_index(0, 1, 1)] = corners[7];
        (grid, field)
    }

    fn case11_idx_subcase_pairs() -> Vec<(u8, i8)> {
        let mut hits = Vec::new();
        for our_idx in 0..=255_u8 {
            let lewiner_idx = !our_idx;
            let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
            if base_case == 11 {
                hits.push((our_idx, subcase));
            }
        }
        hits
    }

    #[test]
    fn case11_dispatch_has_six_configs() {
        assert_eq!(case11_idx_subcase_pairs().len(), 6);
    }

    /// Pin: `TEST_10[cfg][2] == 7` for every cfg (body-diagonal sentinel).
    #[test]
    fn test_10_column_2_is_seven() {
        for (cfg, row) in TEST_10.iter().enumerate() {
            assert_eq!(row[2], 7, "TEST_10[{cfg}][2] = {} (expected 7)", row[2]);
        }
    }

    /// Transcription oracle: `test10[1]` (Julia line 1199 → `(2, 4, 7)`).
    #[test]
    fn test_10_first_row_matches_lut() {
        assert_eq!(TEST_10[0], [2, 4, 7]);
    }

    /// Transcription oracle: `tiling10_1_1[1]` (Julia line 1218 →
    /// `(6, 11, 8, 12, 8, 11, 9, 2, 10, 2, 9, 4)` → 0-based).
    #[test]
    fn tiling_10_1_1_first_row_matches_lut() {
        assert_eq!(TILING_10_1_1[0], [5, 10, 7, 11, 7, 10, 8, 1, 9, 1, 8, 3]);
    }

    /// Sharp oracle: every entry in every tiling is `≤ 12`.
    #[test]
    fn tilings_use_only_legal_indices() {
        for cfg in 0..6 {
            for &e in &TILING_10_1_1[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_10_1_1_[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_10_1_2[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_10_2[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_10_2_[cfg] {
                assert!(e <= 12);
            }
        }
    }

    #[test]
    fn tiling_dimensions_match_lewiner() {
        assert_eq!(TEST_10.len(), 6);
        assert_eq!(TILING_10_1_1.len(), 6);
        assert_eq!(TILING_10_1_1[0].len(), 12);
        assert_eq!(TILING_10_1_2[0].len(), 24);
        assert_eq!(TILING_10_2[0].len(), 24);
    }

    /// Behavioral: every Lewiner-11 voxel under ±1 inputs runs `emit`
    /// to completion with 4 or 8 triangles.
    #[test]
    fn case10_dispatch_runs_on_every_case11_voxel() {
        for (our_idx, subcase) in case11_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            emit(
                &grid,
                CellCoord { x: 0, y: 0, z: 0 },
                &corners,
                subcase,
                &mut mesh,
            );
            let n = mesh.triangle_count();
            assert!(matches!(n, 4 | 8), "expected 4/8, got {n}");
            assert_eq!(mesh.vertices.len(), 3 * n);
            for v in &mesh.vertices {
                assert!(v.is_finite());
                assert!((0.0..=1.0).contains(&v.x));
                assert!((0.0..=1.0).contains(&v.y));
                assert!((0.0..=1.0).contains(&v.z));
            }
        }
    }

    /// Adversarial: every triangle non-degenerate across 6 × 32 sweep.
    #[test]
    fn case10_triangles_are_nondegenerate() {
        for (our_idx, subcase) in case11_idx_subcase_pairs() {
            for perturb_seed in 0..32_u32 {
                let mut corners = [0.0_f32; 8];
                for (bit, c) in corners.iter_mut().enumerate() {
                    let inside = our_idx & (1 << bit) != 0;
                    let mag = 0.3
                        + 1.4
                            * f32::from(
                                ((perturb_seed.wrapping_mul(7) ^ u32::try_from(bit).unwrap_or(0))
                                    & 0xF) as u8,
                            )
                            / 15.0;
                    *c = if inside { -mag } else { mag };
                }
                let (grid, _) = unit_grid_with_corners(&corners);
                let mut mesh = Mesh::default();
                emit(
                    &grid,
                    CellCoord { x: 0, y: 0, z: 0 },
                    &corners,
                    subcase,
                    &mut mesh,
                );
                for tri in &mesh.indices {
                    let v0 = mesh.vertices[tri[0] as usize];
                    let v1 = mesh.vertices[tri[1] as usize];
                    let v2 = mesh.vertices[tri[2] as usize];
                    let area2 = (v1 - v0).cross(v2 - v0).length_squared();
                    assert!(area2 > 1e-12);
                }
            }
        }
    }
}
