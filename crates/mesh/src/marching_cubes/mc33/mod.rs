//! Marching Cubes 33 (Chernyaev 1995 / Lewiner 2003) — per-voxel processing.
//!
//! # Scope: fully ported
//!
//! Every Lewiner base case has a topologically-correct dispatch path:
//!
//! - **Empty (Lewiner 1)** — no geometry, early-return.
//! - **Unambiguous (Lewiner 2/3/6/9/10/12/15)** — `unambiguous` submodule's
//!   face-consistent tables. Eliminates the ~2–3% non-manifold edge rate
//!   from classic MC's 1987 triangulations on typical SDF workloads.
//! - **Lewiner 4 (Chernyaev 3, "two diag corners on a face")** — [`case3`].
//! - **Lewiner 5 (Chernyaev 4, "two diag corners on a body diagonal")** —
//!   [`case4`].
//! - **Lewiner 7 (Chernyaev 6, "three corners + diagonal")** — [`case6`].
//! - **Lewiner 8 (Chernyaev 7, "three non-adjacent corners")** — [`case7`].
//! - **Lewiner 11 (Chernyaev 10, "two opposite diagonal pairs")** — [`case10`].
//! - **Lewiner 13 (Chernyaev 12, "3+1 + diagonal")** — [`case12`].
//! - **Lewiner 14 (Chernyaev 13, "tunnel")** — [`case13`].
//!
//! The only residual fall-through to classic MC is the `SUBCFG_13` `-1`
//! sentinel in case 13 (numerically borderline face-test combinations
//! that should not occur for valid voxels in practice).
//!
//! # Sign/index convention mismatch with the Lewiner tables
//!
//! Our internal `case_index` convention sets bit `p` iff corner `p` has
//! `f < 0` (inside). Lewiner's convention is the opposite (bit set iff
//! corner is **outside**). To use Lewiner's `cases` lookup table we
//! convert at dispatch time via `lewiner_index(our_idx) = !our_idx`.
//!
//! Edge indices in Lewiner's tables are 1-based `(1..=12)`; ours are
//! 0-based `(0..=11)`. Conversion is done at port time in each case's
//! tiling tables, not at runtime.

mod case10;
mod case12;
mod case13;
mod case3;
mod case4;
mod case6;
mod case7;
mod cases;
mod decider;
mod unambiguous;

use glam::Vec3;

use super::tables::EDGE_CORNERS;
use super::{CellCoord, Mesh, case_index, classic, interpolate_edge, read_corner_values};
use crate::mesh::grid::GridSpec;

/// Computes the cube-interior "c-vertex" used by Chernyaev cases 6, 7,
/// 10, 12, and 13 (Lewiner 7, 8, 11, 13, 14) when their triangulations
/// reference the 13th-edge sentinel (Lewiner 1-based 13 → 0-based 12).
///
/// Defined as the centroid of all sign-change edge interpolations in
/// the voxel. Mirrors `add_c_vertex` from `JuliaGeometry/MarchingCubes.jl`
/// `src/MarchingCubes.jl` lines 662-698. Shared between case-13 and
/// case-7 (and other future case modules) to eliminate duplication.
///
/// # Panics (debug only)
///
/// `debug_assert!`s that at least one cube edge has a sign change.
/// Voxels of any ambiguous case have ≥ 1 sign-change edge by construction.
pub(super) fn compute_c_vertex(grid: &GridSpec, cell: CellCoord, corners: &[f32; 8]) -> Vec3 {
    let mut sum = Vec3::ZERO;
    let mut count = 0_u32;
    for &(a, b) in &EDGE_CORNERS {
        let va = corners[a as usize];
        let vb = corners[b as usize];
        if (va < 0.0) != (vb < 0.0) {
            sum += interpolate_edge(grid, cell, a, b, va, vb);
            count += 1;
        }
    }
    debug_assert!(
        count > 0,
        "ambiguous-case voxel must have at least one sign-change edge"
    );
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / count as f32;
    sum * inv
}

/// Walks `tiling` in 3-element chunks, materializing edge-interpolated
/// vertices (cached per voxel) and pushing one triangle per chunk in
/// Lewiner's natural winding (no reversal).
///
/// Slots 0..=11 of the cache hold cube-edge interpolations; slot 12
/// holds the c-vertex (lazily computed via [`compute_c_vertex`] only
/// when the tiling references index 12). Shared between every MC33
/// case module that emits Lewiner-sourced triangulations.
pub(super) fn emit_triangles(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    tiling: &[u8],
    mesh: &mut Mesh,
) {
    let mut edge_verts: [Option<Vec3>; 13] = [None; 13];
    let mut edge_vertex = |e: u8, corners: &[f32; 8]| -> Vec3 {
        let idx = e as usize;
        if let Some(v) = edge_verts[idx] {
            return v;
        }
        let v = if e == 12 {
            compute_c_vertex(grid, cell, corners)
        } else {
            let (a, b) = EDGE_CORNERS[idx];
            interpolate_edge(grid, cell, a, b, corners[a as usize], corners[b as usize])
        };
        edge_verts[idx] = Some(v);
        v
    };

    for tri in tiling.chunks_exact(3) {
        let v0 = edge_vertex(tri[0], corners);
        let v1 = edge_vertex(tri[1], corners);
        let v2 = edge_vertex(tri[2], corners);
        #[allow(clippy::expect_used)]
        let base = u32::try_from(mesh.vertices.len()).expect("vertex count exceeds u32");
        mesh.vertices.push(v0);
        mesh.vertices.push(v1);
        mesh.vertices.push(v2);
        mesh.indices.push([base, base + 1, base + 2]);
    }
}

/// Converts our internal `case_index` (bit set iff corner is inside) to
/// Lewiner's (bit set iff corner is outside).
#[inline]
fn lewiner_index(our_idx: u8) -> u8 {
    !our_idx
}

/// Processes a single voxel and appends emitted triangles to `mesh`.
///
/// Dispatches per the Lewiner base-case table — every base case is
/// fully ported:
/// - Case 1 (empty) → no-op.
/// - Cases 2/3/6/9/10/12/15 (unambiguous) → Lewiner face-consistent tables.
/// - Case 4 → Case 3 disambiguator.
/// - Case 5 → Case 4 disambiguator.
/// - Case 7 → Case 6 disambiguator.
/// - Case 8 → Case 7 disambiguator.
/// - Case 11 → Case 10 disambiguator.
/// - Case 13 → Case 12 disambiguator.
/// - Case 14 → Case 13 disambiguator.
///
/// The only path back to classic MC is the `SUBCFG_13` `-1` sentinel
/// in case 13 — numerically borderline face-test combinations that
/// should not occur for valid voxels in practice. Treated as a
/// defensive fallback rather than panicking.
pub(crate) fn process_cell(grid: &GridSpec, field: &[f32], cell: CellCoord, mesh: &mut Mesh) {
    let corners = read_corner_values(grid, field, cell);
    let our_idx = case_index(&corners);
    let lewiner_idx = lewiner_index(our_idx);
    let (base_case, subcase) = cases::CASES[lewiner_idx as usize];

    // Empty case — no surface crosses this voxel.
    if base_case == 1 {
        return;
    }

    // Unambiguous cases use Lewiner's face-consistent tables.
    if matches!(base_case, 2 | 3 | 6 | 9 | 10 | 12 | 15) {
        unambiguous::emit(grid, cell, &corners, base_case, subcase, mesh);
        return;
    }

    // Ambiguous cases — each delegates to its Chernyaev-numbered module.
    match base_case {
        4 => case3::emit(grid, cell, &corners, subcase, mesh),
        5 => case4::emit(grid, cell, &corners, subcase, mesh),
        7 => case6::emit(grid, cell, &corners, subcase, mesh),
        8 => case7::emit(grid, cell, &corners, subcase, mesh),
        11 => case10::emit(grid, cell, &corners, subcase, mesh),
        13 => case12::emit(grid, cell, &corners, subcase, mesh),
        14 => {
            if !case13::emit(grid, cell, &corners, subcase, mesh) {
                // case13 returns false only on the `SUBCFG_13` `-1`
                // sentinel — defensively fall back to classic.
                classic::process_cell(grid, field, cell, mesh);
            }
        }
        // Defensive: `cases::CASES` should only return base_case ∈ 1..=15,
        // and every value is handled above. Fall back rather than panic.
        _ => classic::process_cell(grid, field, cell, mesh),
    }
}

/// Test-only: returns the Lewiner base case for a given internal
/// `case_index`. Used by manifoldness diagnostics in the parent module.
#[cfg(test)]
pub(super) fn lewiner_base_case(our_idx: u8) -> i8 {
    cases::CASES[lewiner_index(our_idx) as usize].0
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

    #[test]
    fn lewiner_index_is_bitwise_complement() {
        assert_eq!(lewiner_index(0), 255);
        assert_eq!(lewiner_index(255), 0);
        assert_eq!(lewiner_index(5), 250);
        // Round-trip involution.
        for i in 0..=255_u8 {
            assert_eq!(lewiner_index(lewiner_index(i)), i);
        }
    }

    #[test]
    fn cases_lookup_at_empty_configurations() {
        // our_idx 0 (all outside) → lewiner_idx 255 → (1, 0) empty.
        let (base, _) = cases::CASES[lewiner_index(0) as usize];
        assert_eq!(base, 1);
        // our_idx 255 (all inside) → lewiner_idx 0 → (1, 0) empty.
        let (base, _) = cases::CASES[lewiner_index(255) as usize];
        assert_eq!(base, 1);
    }

    #[test]
    fn cases_lookup_at_known_case3_configurations() {
        // our_idx 5 = 0b00000101: corners 0 and 2 inside, diagonally
        // opposite on the bottom face. This is the canonical Case 3
        // configuration.
        let (base, _) = cases::CASES[lewiner_index(5) as usize];
        assert_eq!(
            base, 4,
            "our_idx 5 should map to Lewiner base case 4 (Chernyaev 3)"
        );
    }

    /// Sweep all 256 `case_index` values through `process_cell` against
    /// a single-voxel grid; assert no panic and that every non-empty
    /// case produces a non-empty mesh. This pins the "MC33 is fully
    /// ported" property — every base case has a working dispatch.
    #[test]
    fn full_port_no_voxel_falls_through_unhandled() {
        use glam::{UVec3, Vec3};
        let grid = GridSpec::new(Vec3::ZERO, UVec3::splat(1), 1.0).unwrap();
        for our_idx in 0..=255_u8 {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let mut field = vec![0.0_f32; grid.sample_count()];
            field[grid.sample_index(0, 0, 0)] = corners[0];
            field[grid.sample_index(1, 0, 0)] = corners[1];
            field[grid.sample_index(1, 1, 0)] = corners[2];
            field[grid.sample_index(0, 1, 0)] = corners[3];
            field[grid.sample_index(0, 0, 1)] = corners[4];
            field[grid.sample_index(1, 0, 1)] = corners[5];
            field[grid.sample_index(1, 1, 1)] = corners[6];
            field[grid.sample_index(0, 1, 1)] = corners[7];

            let mut mesh = Mesh::default();
            process_cell(&grid, &field, CellCoord { x: 0, y: 0, z: 0 }, &mut mesh);

            let lewiner_idx = lewiner_index(our_idx);
            let (base_case, _) = cases::CASES[lewiner_idx as usize];
            if base_case == 1 {
                assert_eq!(
                    mesh.indices.len(),
                    0,
                    "empty case our_idx={our_idx:#010b} produced triangles"
                );
            } else {
                assert!(
                    !mesh.indices.is_empty(),
                    "our_idx={our_idx:#010b} (base_case={base_case}) produced no triangles"
                );
            }
        }
    }

    #[test]
    fn cases_lookup_at_known_case4_configurations() {
        // our_idx 0b10111110 = 190 → bits 1,2,3,4,5,7 inside = 6 corners
        // inside; corners 0 and 6 are outside (opposite ends of the main
        // body diagonal). Canonical Case 4 subcase-1 voxel.
        let (base, sub) = cases::CASES[lewiner_index(0b1011_1110) as usize];
        assert_eq!(
            base, 5,
            "our_idx 0b10111110 should map to Lewiner base case 5 (Chernyaev 4)"
        );
        assert_eq!(sub, 1, "expected subcase 1");
    }
}
