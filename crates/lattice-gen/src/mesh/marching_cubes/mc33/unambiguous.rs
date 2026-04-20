//! Lewiner 2003 face-consistent tables for the 7 unambiguous MC33 base cases.
//!
//! # Why this module exists
//!
//! Classic Lorensen-Cline 1987's tables satisfy *self-manifoldness* within a
//! single voxel but do **not** guarantee that two neighboring voxels produce
//! matching triangle edges on their shared face. Even on unambiguous cases,
//! two neighbors can independently pick different face diagonals → overlapping
//! triangles → edges shared by 3+ faces. Empirically this is ~2–3% of edges
//! on typical SDF workloads.
//!
//! Lewiner 2003 published replacement tables constructed to be face-consistent
//! by construction across all 256 `case_index` values. This module ports those
//! tables for the 7 unambiguous base cases (Lewiner numbering):
//!
//! | Base case | Chernyaev | Description             | Source table | Tris |
//! |-----------|-----------|-------------------------|--------------|------|
//! | 2         | 1         | Single corner           | `tiling1`    | 1    |
//! | 3         | 2         | Two adjacent corners    | `tiling2`    | 2    |
//! | 6         | 5         | Three adjacent corners  | `tiling5`    | 3    |
//! | 9         | 8         | Four corners on a face  | `tiling8`    | 2    |
//! | 10        | 9         | Four corners, zigzag    | `tiling9`    | 4    |
//! | 12        | 11        | Cross-face              | `tiling11`   | 4    |
//! | 15        | 14        | Four "flower" corners   | `tiling14`   | 4    |
//!
//! The 6 ambiguous base cases (4, 5 ported in Sessions 1/2; 7, 8, 11, 13, 14
//! pending) and the empty case (1) do not use this module.
//!
//! # Port
//!
//! Tables from Lewiner 2003 `LookUpTable.h` via `JuliaGeometry`
//! `MarchingCubes.jl` (`src/lut.jl`, `tiling1` / `tiling2` / `tiling5` /
//! `tiling8` / `tiling9` / `tiling11` / `tiling14` as of fetched version).
//! Edge indices converted from Lewiner's 1-based to our 0-based at port time.
//!
//! # Winding
//!
//! Same as `case3` / `case4`: emit `[base, base+1, base+2]` with no reversal.
//! Inherits Session 1's empirical finding that Lewiner-sourced tables produce
//! outward normals under our "inside = f < 0" convention.

use glam::Vec3;

use super::super::tables::EDGE_CORNERS;
use super::super::{CellCoord, Mesh, interpolate_edge};
use crate::mesh::grid::GridSpec;

/// Source: Lewiner `tiling1[]`, 16 × 3 entries. Base case 2 (single corner).
#[rustfmt::skip]
const TILING_1: [[u8; 3]; 16] = [
    [ 0,  8,  3],
    [ 0,  1,  9],
    [ 1,  2, 10],
    [ 3, 11,  2],
    [ 4,  7,  8],
    [ 9,  5,  4],
    [10,  6,  5],
    [ 7,  6, 11],
    [ 7, 11,  6],
    [10,  5,  6],
    [ 9,  4,  5],
    [ 4,  8,  7],
    [ 3,  2, 11],
    [ 1, 10,  2],
    [ 0,  9,  1],
    [ 0,  3,  8],
];

/// Source: Lewiner `tiling2[]`, 24 × 6 entries. Base case 3 (two adjacent).
#[rustfmt::skip]
const TILING_2: [[u8; 6]; 24] = [
    [ 1,  8,  3,  9,  8,  1],
    [ 0, 11,  2,  8, 11,  0],
    [ 4,  3,  0,  7,  3,  4],
    [ 9,  2, 10,  0,  2,  9],
    [ 0,  5,  4,  1,  5,  0],
    [ 3, 10,  1, 11, 10,  3],
    [ 1,  6,  5,  2,  6,  1],
    [ 7,  2,  3,  6,  2,  7],
    [ 9,  7,  8,  5,  7,  9],
    [ 6,  8,  4, 11,  8,  6],
    [10,  4,  9,  6,  4, 10],
    [11,  5, 10,  7,  5, 11],
    [11, 10,  5,  7, 11,  5],
    [10,  9,  4,  6, 10,  4],
    [ 6,  4,  8, 11,  6,  8],
    [ 9,  8,  7,  5,  9,  7],
    [ 7,  3,  2,  6,  7,  2],
    [ 1,  5,  6,  2,  1,  6],
    [ 3,  1, 10, 11,  3, 10],
    [ 0,  4,  5,  1,  0,  5],
    [ 9, 10,  2,  0,  9,  2],
    [ 4,  0,  3,  7,  4,  3],
    [ 0,  2, 11,  8,  0, 11],
    [ 1,  3,  8,  9,  1,  8],
];

/// Source: Lewiner `tiling5[]`, 48 × 9 entries. Base case 6 (three adjacent).
#[rustfmt::skip]
const TILING_5: [[u8; 9]; 48] = [
    [ 2,  8,  3,  2, 10,  8, 10,  9,  8],
    [ 1, 11,  2,  1,  9, 11,  9,  8, 11],
    [ 4,  1,  9,  4,  7,  1,  7,  3,  1],
    [ 8,  5,  4,  8,  3,  5,  3,  1,  5],
    [ 0, 10,  1,  0,  8, 10,  8, 11, 10],
    [11,  4,  7, 11,  2,  4,  2,  0,  4],
    [ 7,  0,  8,  7,  6,  0,  6,  2,  0],
    [ 9,  3,  0,  9,  5,  3,  5,  7,  3],
    [ 3,  6, 11,  3,  0,  6,  0,  4,  6],
    [ 3,  9,  0,  3, 11,  9, 11, 10,  9],
    [ 5,  2, 10,  5,  4,  2,  4,  0,  2],
    [ 9,  6,  5,  9,  0,  6,  0,  2,  6],
    [ 0,  7,  8,  0,  1,  7,  1,  5,  7],
    [10,  0,  1, 10,  6,  0,  6,  4,  0],
    [ 6,  3, 11,  6,  5,  3,  5,  1,  3],
    [10,  7,  6, 10,  1,  7,  1,  3,  7],
    [ 1,  4,  9,  1,  2,  4,  2,  6,  4],
    [11,  1,  2, 11,  7,  1,  7,  5,  1],
    [ 8,  2,  3,  8,  4,  2,  4,  6,  2],
    [ 2,  5, 10,  2,  3,  5,  3,  7,  5],
    [ 7, 10,  6,  7,  8, 10,  8,  9, 10],
    [ 6,  9,  5,  6, 11,  9, 11,  8,  9],
    [ 5,  8,  4,  5, 10,  8, 10, 11,  8],
    [ 4, 11,  7,  4,  9, 11,  9, 10, 11],
    [ 4,  7, 11,  4, 11,  9,  9, 11, 10],
    [ 5,  4,  8,  5,  8, 10, 10,  8, 11],
    [ 6,  5,  9,  6,  9, 11, 11,  9,  8],
    [ 7,  6, 10,  7, 10,  8,  8, 10,  9],
    [ 2, 10,  5,  2,  5,  3,  3,  5,  7],
    [ 8,  3,  2,  8,  2,  4,  4,  2,  6],
    [11,  2,  1, 11,  1,  7,  7,  1,  5],
    [ 1,  9,  4,  1,  4,  2,  2,  4,  6],
    [10,  6,  7, 10,  7,  1,  1,  7,  3],
    [ 6, 11,  3,  6,  3,  5,  5,  3,  1],
    [10,  1,  0, 10,  0,  6,  6,  0,  4],
    [ 0,  8,  7,  0,  7,  1,  1,  7,  5],
    [ 9,  5,  6,  9,  6,  0,  0,  6,  2],
    [ 5, 10,  2,  5,  2,  4,  4,  2,  0],
    [ 3,  0,  9,  3,  9, 11, 11,  9, 10],
    [ 3, 11,  6,  3,  6,  0,  0,  6,  4],
    [ 9,  0,  3,  9,  3,  5,  5,  3,  7],
    [ 7,  8,  0,  7,  0,  6,  6,  0,  2],
    [11,  7,  4, 11,  4,  2,  2,  4,  0],
    [ 0,  1, 10,  0, 10,  8,  8, 10, 11],
    [ 8,  4,  5,  8,  5,  3,  3,  5,  1],
    [ 4,  9,  1,  4,  1,  7,  7,  1,  3],
    [ 1,  2, 11,  1, 11,  9,  9, 11,  8],
    [ 2,  3,  8,  2,  8, 10, 10,  8,  9],
];

/// Source: Lewiner `tiling8[]`, 6 × 6 entries. Base case 9 (four on face).
#[rustfmt::skip]
const TILING_8: [[u8; 6]; 6] = [
    [ 9,  8, 10, 10,  8, 11],
    [ 1,  5,  3,  3,  5,  7],
    [ 0,  4,  2,  4,  6,  2],
    [ 0,  2,  4,  4,  2,  6],
    [ 1,  3,  5,  3,  7,  5],
    [ 9, 10,  8, 10, 11,  8],
];

/// Source: Lewiner `tiling9[]`, 8 × 12 entries. Base case 10 (zigzag).
#[rustfmt::skip]
const TILING_9: [[u8; 12]; 8] = [
    [ 2, 10,  5,  3,  2,  5,  3,  5,  4,  3,  4,  8],
    [ 4,  7, 11,  9,  4, 11,  9, 11,  2,  9,  2,  1],
    [10,  7,  6,  1,  7, 10,  1,  8,  7,  1,  0,  8],
    [ 3,  6, 11,  0,  6,  3,  0,  5,  6,  0,  9,  5],
    [ 3, 11,  6,  0,  3,  6,  0,  6,  5,  0,  5,  9],
    [10,  6,  7,  1, 10,  7,  1,  7,  8,  1,  8,  0],
    [ 4, 11,  7,  9, 11,  4,  9,  2, 11,  9,  1,  2],
    [ 2,  5, 10,  3,  5,  2,  3,  4,  5,  3,  8,  4],
];

/// Source: Lewiner `tiling11[]`, 12 × 12 entries. Base case 12 (cross-face).
#[rustfmt::skip]
const TILING_11: [[u8; 12]; 12] = [
    [ 2, 10,  9,  2,  9,  7,  2,  7,  3,  7,  9,  4],
    [ 1,  6,  2,  1,  8,  6,  1,  9,  8,  8,  7,  6],
    [ 8,  3,  1,  8,  1,  6,  8,  6,  4,  6,  1, 10],
    [ 0,  8, 11,  0, 11,  5,  0,  5,  1,  5, 11,  6],
    [ 9,  5,  7,  9,  7,  2,  9,  2,  0,  2,  7, 11],
    [ 5,  0,  4,  5, 11,  0,  5, 10, 11, 11,  3,  0],
    [ 5,  4,  0,  5,  0, 11,  5, 11, 10, 11,  0,  3],
    [ 9,  7,  5,  9,  2,  7,  9,  0,  2,  2, 11,  7],
    [ 0, 11,  8,  0,  5, 11,  0,  1,  5,  5,  6, 11],
    [ 8,  1,  3,  8,  6,  1,  8,  4,  6,  6, 10,  1],
    [ 1,  2,  6,  1,  6,  8,  1,  8,  9,  8,  6,  7],
    [ 2,  9, 10,  2,  7,  9,  2,  3,  7,  7,  4,  9],
];

/// Source: Lewiner `tiling14[]`, 12 × 12 entries. Base case 15 (flower).
#[rustfmt::skip]
const TILING_14: [[u8; 12]; 12] = [
    [ 5,  9,  8,  5,  8,  2,  5,  2,  6,  3,  2,  8],
    [ 2,  1,  5,  2,  5,  8,  2,  8, 11,  4,  8,  5],
    [ 9,  4,  6,  9,  6,  3,  9,  3,  1, 11,  3,  6],
    [ 1, 11, 10,  1,  4, 11,  1,  0,  4,  7, 11,  4],
    [ 8,  2,  0,  8,  5,  2,  8,  7,  5, 10,  2,  5],
    [ 0,  7,  3,  0, 10,  7,  0,  9, 10,  6,  7, 10],
    [ 0,  3,  7,  0,  7, 10,  0, 10,  9,  6, 10,  7],
    [ 8,  0,  2,  8,  2,  5,  8,  5,  7, 10,  5,  2],
    [ 1, 10, 11,  1, 11,  4,  1,  4,  0,  7,  4, 11],
    [ 9,  6,  4,  9,  3,  6,  9,  1,  3, 11,  6,  3],
    [ 2,  5,  1,  2,  8,  5,  2, 11,  8,  4,  5,  8],
    [ 5,  8,  9,  5,  2,  8,  5,  6,  2,  3,  8,  2],
];

/// Emits triangles for a voxel in one of the 7 unambiguous Lewiner base
/// cases via the corresponding face-consistent tiling table.
///
/// `base_case` must be in `{2, 3, 6, 9, 10, 12, 15}`; `subcase` is Lewiner's
/// 1-based configuration index from [`super::cases::CASES`].
///
/// # Panics
///
/// `debug_assert!`s that `base_case` is one of the expected unambiguous
/// values. An unexpected value indicates a dispatch bug in
/// [`super::process_cell`].
#[allow(clippy::cast_sign_loss, clippy::expect_used)]
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    base_case: i8,
    subcase: i8,
    mesh: &mut Mesh,
) {
    debug_assert!(subcase >= 1, "subcase must be 1-based: {subcase}");
    let cfg = (subcase - 1) as usize;

    let tiling: &[u8] = match base_case {
        2 => &TILING_1[cfg],
        3 => &TILING_2[cfg],
        6 => &TILING_5[cfg],
        9 => &TILING_8[cfg],
        10 => &TILING_9[cfg],
        12 => &TILING_11[cfg],
        15 => &TILING_14[cfg],
        _ => unreachable!("unambiguous::emit called with base_case {base_case}"),
    };

    let mut edge_verts: [Option<Vec3>; 12] = [None; 12];
    let mut edge_vertex = |e: u8, corners: &[f32; 8]| -> Vec3 {
        let idx = e as usize;
        if let Some(v) = edge_verts[idx] {
            v
        } else {
            let (a, b) = EDGE_CORNERS[idx];
            let v = interpolate_edge(grid, cell, a, b, corners[a as usize], corners[b as usize]);
            edge_verts[idx] = Some(v);
            v
        }
    };

    let mut chunks = tiling.chunks_exact(3);
    for tri in chunks.by_ref() {
        let v0 = edge_vertex(tri[0], corners);
        let v1 = edge_vertex(tri[1], corners);
        let v2 = edge_vertex(tri[2], corners);
        let base = u32::try_from(mesh.vertices.len()).expect("vertex count exceeds u32");
        mesh.vertices.push(v0);
        mesh.vertices.push(v1);
        mesh.vertices.push(v2);
        mesh.indices.push([base, base + 1, base + 2]);
    }
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
    use glam::UVec3;

    /// Every edge index in every tiling must be in `0..=11`.
    /// Catches transcription errors.
    #[test]
    fn all_tables_have_valid_edge_indices() {
        let flat1: &[u8] = TILING_1.as_flattened();
        let flat2: &[u8] = TILING_2.as_flattened();
        let flat5: &[u8] = TILING_5.as_flattened();
        let flat8: &[u8] = TILING_8.as_flattened();
        let flat9: &[u8] = TILING_9.as_flattened();
        let flat11: &[u8] = TILING_11.as_flattened();
        let flat14: &[u8] = TILING_14.as_flattened();

        for (name, flat) in [
            ("TILING_1", flat1),
            ("TILING_2", flat2),
            ("TILING_5", flat5),
            ("TILING_8", flat8),
            ("TILING_9", flat9),
            ("TILING_11", flat11),
            ("TILING_14", flat14),
        ] {
            for (i, &e) in flat.iter().enumerate() {
                assert!(
                    e <= 11,
                    "{name}[{i}] = {e} out of 0..=11 (transcription error)"
                );
            }
        }
    }

    /// Build a single-voxel grid whose corner values match `corners`.
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

    /// Returns corner values `[±1, ...]` for the given `our_idx`, where bit
    /// `i` set means corner `i` is inside (`-1.0`).
    fn corners_for(our_idx: u8) -> [f32; 8] {
        let mut out = [1.0_f32; 8];
        for (i, c) in out.iter_mut().enumerate() {
            if our_idx & (1 << i) != 0 {
                *c = -1.0;
            }
        }
        out
    }

    /// For each of the 7 unambiguous base cases we handle, find some
    /// `our_idx` that lands in that base case and verify `emit` produces
    /// non-zero triangles with in-cell, finite vertices.
    #[test]
    fn emit_one_voxel_per_base_case() {
        let targets = [2_i8, 3, 6, 9, 10, 12, 15];
        for &target in &targets {
            let mut hit: Option<(u8, i8)> = None;
            for our_idx in 0..=255_u8 {
                let lewiner_idx = !our_idx;
                let (b, s) = super::super::cases::CASES[lewiner_idx as usize];
                if b == target {
                    hit = Some((our_idx, s));
                    break;
                }
            }
            let (our_idx, subcase) = hit.expect("no case_index maps to this base case");
            let corners = corners_for(our_idx);
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            emit(&grid, cell, &corners, target, subcase, &mut mesh);
            assert!(
                mesh.triangle_count() > 0,
                "base_case {target} (our_idx={our_idx}, subcase={subcase}) \
                 produced no triangles"
            );
            for v in &mesh.vertices {
                assert!(
                    v.is_finite()
                        && (0.0..=1.0).contains(&v.x)
                        && (0.0..=1.0).contains(&v.y)
                        && (0.0..=1.0).contains(&v.z),
                    "base_case {target}: vertex out of unit cube: {v:?}"
                );
            }
        }
    }

    /// Complementary bit patterns (global sign flip) must produce equal
    /// triangle counts — topology of an SDF is invariant under f ↦ -f.
    #[test]
    fn emit_sign_flip_preserves_triangle_count() {
        let targets = [2_i8, 3, 6, 9, 10, 12, 15];
        for &target in &targets {
            let mut hit: Option<(u8, i8)> = None;
            for our_idx in 0..=255_u8 {
                let lewiner_idx = !our_idx;
                let (b, s) = super::super::cases::CASES[lewiner_idx as usize];
                if b == target {
                    hit = Some((our_idx, s));
                    break;
                }
            }
            let (our_idx, subcase_a) = hit.unwrap();
            let flipped = !our_idx;
            let subcase_b = super::super::cases::CASES[(!flipped) as usize].1;

            let corners_a = corners_for(our_idx);
            let corners_b = corners_for(flipped);

            let run = |corners: &[f32; 8], s: i8| -> usize {
                let (grid, _) = unit_grid_with_corners(corners);
                let mut mesh = Mesh::default();
                emit(
                    &grid,
                    CellCoord { x: 0, y: 0, z: 0 },
                    corners,
                    target,
                    s,
                    &mut mesh,
                );
                mesh.triangle_count()
            };

            assert_eq!(
                run(&corners_a, subcase_a),
                run(&corners_b, subcase_b),
                "base_case {target}: sign-flipped configurations produce \
                 different triangle counts"
            );
        }
    }

    #[test]
    #[should_panic(expected = "unambiguous::emit called with base_case 4")]
    fn emit_panics_on_ambiguous_base_case_4() {
        let corners = [0.0_f32; 8];
        let (grid, _) = unit_grid_with_corners(&corners);
        let mut mesh = Mesh::default();
        emit(
            &grid,
            CellCoord { x: 0, y: 0, z: 0 },
            &corners,
            4,
            1,
            &mut mesh,
        );
    }

    #[test]
    #[should_panic(expected = "unambiguous::emit called with base_case 7")]
    fn emit_panics_on_ambiguous_base_case_7() {
        let corners = [0.0_f32; 8];
        let (grid, _) = unit_grid_with_corners(&corners);
        let mut mesh = Mesh::default();
        emit(
            &grid,
            CellCoord { x: 0, y: 0, z: 0 },
            &corners,
            7,
            1,
            &mut mesh,
        );
    }
}
