//! Case 4 (Chernyaev) / Case 5 (Lewiner) — body-ambiguous disambiguation.
//!
//! The configuration: two corners with the same sign sit at opposite ends
//! of one of the cube's body diagonals; the other six have the opposite
//! sign. The trilinear interpolant has a saddle *inside* the cube, so the
//! face decider cannot resolve the ambiguity — we need the interior
//! decider (see [`super::decider::interior_decider`]).
//!
//! # Subcases
//!
//! - **4.1.1**: interior-connected (2 triangles). Emitted when the
//!   decider returns `true` — the surface threads through the cube as two
//!   simply-connected caps.
//! - **4.1.2**: tunnel (6 triangles). Emitted when the decider returns
//!   `false` — the two same-sign corners are joined by a thin bridge
//!   through the cube's interior.
//!
//! # Port
//!
//! Tables from Lewiner 2003 `LookUpTable.h` via `JuliaGeometry`'s
//! `MarchingCubes.jl` (`src/lut.jl`, `test4` / `tiling4_1` / `tiling4_2`,
//! lines ~490–541 as of fetched version). Edge indices converted from
//! Lewiner's 1-based to our 0-based convention at port time.
//!
//! # Winding
//!
//! Same convention as [`super::case3`]: emit `[base, base+1, base+2]`
//! without reversal. The empirical winding check in Session 1
//! established that Lewiner-sourced tables produce outward normals under
//! our "inside = f < 0" convention; this session inherits that finding.

use glam::Vec3;

use super::super::tables::EDGE_CORNERS;
use super::super::{CellCoord, Mesh, interpolate_edge};
use super::decider::interior_decider;
use crate::grid::GridSpec;

/// For each Case-4 subcase (indexed by `subcase - 1`), the sentinel
/// passed to `interior_decider`. All entries have `|s| == 7` — the
/// interior-test sentinel. Positive vs negative sign selects which
/// polarity of the test to apply, per the Lewiner test-table convention.
///
/// Source: Lewiner `test4[]`, 8 entries.
const TEST_4: [i8; 8] = [7, 7, 7, 7, -7, -7, -7, -7];

/// Triangulation for Case 4.1.1 (interior-connected): 2 triangles =
/// 6 edge indices per subcase. Edge indices are 0-based (converted from
/// Lewiner's 1-based at port time).
///
/// Source: Lewiner `tiling4_1[]`, 8 × 6 entries.
#[rustfmt::skip]
const TILING_4_1: [[u8; 6]; 8] = [
    [ 0,  8,  3,  5, 10,  6],
    [ 0,  1,  9, 11,  7,  6],
    [ 1,  2, 10,  8,  4,  7],
    [ 9,  5,  4,  2,  3, 11],
    [ 4,  5,  9, 11,  3,  2],
    [10,  2,  1,  7,  4,  8],
    [ 9,  1,  0,  6,  7, 11],
    [ 3,  8,  0,  6, 10,  5],
];

/// Triangulation for Case 4.1.2 (interior tunnel): 6 triangles =
/// 18 edge indices per subcase.
///
/// Source: Lewiner `tiling4_2[]`, 8 × 18 entries.
#[rustfmt::skip]
const TILING_4_2: [[u8; 18]; 8] = [
    [ 8,  5,  0,  5,  8,  6,  3,  6,  8,  6,  3, 10,  0, 10,  3, 10,  0,  5],
    [ 9,  6,  1,  6,  9,  7,  0,  7,  9,  7,  0, 11,  1, 11,  0, 11,  1,  6],
    [10,  7,  2,  7, 10,  4,  1,  4, 10,  4,  1,  8,  2,  8,  1,  8,  2,  7],
    [11,  4,  3,  4, 11,  5,  2,  5, 11,  5,  2,  9,  3,  9,  2,  9,  3,  4],
    [ 3,  4, 11,  5, 11,  4, 11,  5,  2,  9,  2,  5,  2,  9,  3,  4,  3,  9],
    [ 2,  7, 10,  4, 10,  7, 10,  4,  1,  8,  1,  4,  1,  8,  2,  7,  2,  8],
    [ 1,  6,  9,  7,  9,  6,  9,  7,  0, 11,  0,  7,  0, 11,  1,  6,  1, 11],
    [ 0,  5,  8,  6,  8,  5,  8,  6,  3, 10,  3,  6,  3, 10,  0,  5,  0, 10],
];

/// Emits Case-4 triangles for the given voxel, dispatching to subcase
/// 4.1.1 or 4.1.2 via the interior decider.
///
/// `subcase` is Lewiner's 1-based subcase index in `1..=8` from
/// [`super::cases::CASES`].
#[allow(clippy::cast_sign_loss, clippy::expect_used)]
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    subcase: i8,
    mesh: &mut Mesh,
) {
    debug_assert!(
        (1..=8).contains(&subcase),
        "Case 4 subcase out of range: {subcase}"
    );
    let cfg = (subcase - 1) as usize;

    let tiling: &[u8] = if interior_decider(corners, TEST_4[cfg]) {
        &TILING_4_1[cfg]
    } else {
        &TILING_4_2[cfg]
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
    use crate::Mesh;
    use crate::grid::GridSpec;
    use glam::UVec3;

    /// Build a single-cell grid whose corner values match `corners` in
    /// our standard sample ordering (x fastest, then y, then z).
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

    /// Emits for the canonical Case-4 subcase-1 voxel: corners 0 and 6
    /// (opposite ends of the main body diagonal) outside, rest inside.
    ///
    /// Our `case_index` has bits 1, 2, 3, 4, 5, 7 set (corners inside) =
    /// `0b10111110 = 190`. Lewiner's index = `!190 & 0xFF = 65`, and
    /// `CASES[65] = (5, 1)` — Lewiner base case 5 (Chernyaev 4),
    /// subcase 1.
    #[test]
    fn case4_emits_for_subcase_1_canonical() {
        let corners: [f32; 8] = [1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0];

        // Confirm the dispatch index is what we claim.
        let our_idx: u8 = 0b1011_1110;
        let lewiner_idx = !our_idx;
        let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
        assert_eq!(
            base_case, 5,
            "our_idx {our_idx:#010b} must map to Lewiner base case 5"
        );
        assert_eq!(subcase, 1, "subcase must be 1 for this canonical voxel");

        let (grid, _field) = unit_grid_with_corners(&corners);
        let mut mesh = Mesh::default();
        let cell = CellCoord { x: 0, y: 0, z: 0 };
        emit(&grid, cell, &corners, subcase, &mut mesh);

        let n = mesh.triangle_count();
        assert!(
            n == 2 || n == 6,
            "Case 4 must produce 2 or 6 triangles, got {n}"
        );

        for v in &mesh.vertices {
            assert!(v.is_finite());
            assert!((0.0..=1.0).contains(&v.x));
            assert!((0.0..=1.0).contains(&v.y));
            assert!((0.0..=1.0).contains(&v.z));
        }
    }

    /// Every one of Case 4's 8 subcases must emit a legal triangulation
    /// for a voxel that lands in it. We drive each subcase via its own
    /// `case_index` value — the 8 `case_index` values that map to
    /// Lewiner base case 5 in the CASES table.
    #[test]
    fn case4_covers_all_eight_subcases() {
        // Scan the CASES table for all `our_idx` values that map to
        // Lewiner base case 5. There are exactly 8.
        let mut hits: Vec<(u8, i8)> = Vec::new();
        for our_idx in 0..=255_u8 {
            let lewiner_idx = !our_idx;
            let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
            if base_case == 5 {
                hits.push((our_idx, subcase));
            }
        }
        assert_eq!(
            hits.len(),
            8,
            "expected exactly 8 Case-4 configurations, got {}",
            hits.len()
        );

        for (our_idx, subcase) in hits {
            // Materialize corners from the bit pattern (bit set = inside
            // = f < 0).
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }

            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            emit(&grid, cell, &corners, subcase, &mut mesh);

            let n = mesh.triangle_count();
            assert!(
                n == 2 || n == 6,
                "our_idx={our_idx:#010b} subcase={subcase}: expected 2 or 6 triangles, got {n}"
            );
            for v in &mesh.vertices {
                assert!(v.is_finite(), "non-finite vertex for our_idx {our_idx}");
                assert!((0.0..=1.0).contains(&v.x));
                assert!((0.0..=1.0).contains(&v.y));
                assert!((0.0..=1.0).contains(&v.z));
            }
        }
    }

    /// Winding check on the canonical subcase-1 voxel: corners 0 and 6
    /// are the two *outside* (f > 0) corners; the rest are inside.
    /// Outward normals therefore point *toward* whichever outside
    /// corner is nearer to the triangle (and away from the surrounding
    /// inside bulk).
    #[test]
    fn case4_triangle_normals_point_outward() {
        let corners: [f32; 8] = [1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0];
        let (grid, _field) = unit_grid_with_corners(&corners);
        let mut mesh = Mesh::default();
        let cell = CellCoord { x: 0, y: 0, z: 0 };
        let (_, subcase) = super::super::cases::CASES[!0b1011_1110_u8 as usize];
        emit(&grid, cell, &corners, subcase, &mut mesh);
        assert!(!mesh.indices.is_empty(), "expected Case-4 triangles");

        let outside_a = Vec3::ZERO; // corner 0
        let outside_b = Vec3::ONE; // corner 6
        let mut checked = 0_u32;
        for (t_idx, tri) in mesh.indices.iter().enumerate() {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            let centroid = (v0 + v1 + v2) / 3.0;
            let normal = (v1 - v0).cross(v2 - v0);
            // Outward = toward the closer outside corner.
            let outside_corner = if (centroid - outside_a).length_squared()
                < (centroid - outside_b).length_squared()
            {
                outside_a
            } else {
                outside_b
            };
            let toward_outside = outside_corner - centroid;
            let dot = normal.dot(toward_outside);
            assert!(
                dot >= 0.0,
                "triangle {t_idx} has inward normal: centroid={centroid:?}, \
                 normal={normal:?}, toward_outside={toward_outside:?}, dot={dot} \
                 (Case 4 winding regression)"
            );
            checked += 1;
        }
        assert!(checked > 0, "no Case-4 triangles checked");
    }

    /// Sign-flip invariant (topology, not boolean): a voxel whose corner
    /// signs are globally inverted must produce the same triangle count.
    /// Matches the plan's invariant #2.
    #[test]
    fn case4_sign_flip_preserves_triangle_count() {
        let corners_a: [f32; 8] = [1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0];
        let corners_b: [f32; 8] = [-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0];

        let run = |corners: &[f32; 8]| -> usize {
            let (grid, _field) = unit_grid_with_corners(corners);
            let our_idx: u8 = {
                let mut idx = 0_u8;
                for (i, c) in corners.iter().enumerate() {
                    if *c < 0.0 {
                        idx |= 1 << i;
                    }
                }
                idx
            };
            let (_, subcase) = super::super::cases::CASES[!our_idx as usize];
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            emit(&grid, cell, corners, subcase, &mut mesh);
            mesh.triangle_count()
        };

        assert_eq!(
            run(&corners_a),
            run(&corners_b),
            "sign-flipped Case-4 configurations must yield the same \
             triangle count (topology is sign-symmetric even when the \
             interior decider's boolean is not)"
        );
    }
}
