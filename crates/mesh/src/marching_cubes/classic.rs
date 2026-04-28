//! Classic Marching Cubes (Lorensen & Cline, 1987) per-voxel processing.
//!
//! The 15 unique cases (256 after sign/rotation) are encoded in the
//! private `tables::classic` module. Known limitation: at the ~14
//! ambiguous configurations, this implementation picks the
//! Lorensen-default triangulation, which can produce holes and
//! non-manifold edges at ~5% of edges on smooth SDFs. See
//! [`super::mc33`] for the topologically correct alternative.

use glam::Vec3;

use super::tables::EDGE_CORNERS;
use super::tables::classic::{EDGE_TABLE, TRI_TABLE};
use super::{CellCoord, Mesh, case_index, interpolate_edge, read_corner_values};
use crate::grid::GridSpec;

/// Processes a single voxel and appends emitted triangles to `mesh`.
///
/// # Winding convention
///
/// The Bourke-style tables produce triangles whose normals point inward
/// under our "inside = f < 0" sign convention. To match the project-wide
/// "outward normals" convention (verified by the signed-volume integration
/// test in `tests/meshing.rs`), we reverse the emitted triangle winding:
/// `[base, base+2, base+1]` instead of `[base, base+1, base+2]`.
///
/// # Casts and `expect`
///
/// - `tri[i] as usize` casts `i8` to `usize`; we verify `tri[i] >= 0`
///   immediately before each cast, so this is safe.
/// - `u32::try_from(mesh.vertices.len()).expect(...)` asserts the mesh has
///   not overflowed `u32` vertex count. For the grid sizes we support
///   this is unreachable in practice; the message documents the invariant.
#[allow(clippy::cast_sign_loss, clippy::expect_used)]
pub(crate) fn process_cell(grid: &GridSpec, field: &[f32], cell: CellCoord, mesh: &mut Mesh) {
    let corner_values = read_corner_values(grid, field, cell);
    let case_idx = case_index(&corner_values);

    let edge_mask = EDGE_TABLE[case_idx as usize];
    if edge_mask == 0 {
        return;
    }

    // Compute surface vertex positions on each active edge.
    let mut edge_vertices = [Vec3::ZERO; 12];
    for e in 0..12 {
        if edge_mask & (1 << e) != 0 {
            let (a, b) = EDGE_CORNERS[e];
            edge_vertices[e] = interpolate_edge(
                grid,
                cell,
                a,
                b,
                corner_values[a as usize],
                corner_values[b as usize],
            );
        }
    }

    // Emit triangles per the lookup table.
    let tri = TRI_TABLE[case_idx as usize];
    let mut i = 0;
    while i < tri.len() && tri[i] >= 0 {
        let e0 = tri[i] as usize;
        let e1 = tri[i + 1] as usize;
        let e2 = tri[i + 2] as usize;
        let base = u32::try_from(mesh.vertices.len()).expect("vertex count exceeds u32");
        mesh.vertices.push(edge_vertices[e0]);
        mesh.vertices.push(edge_vertices[e1]);
        mesh.vertices.push(edge_vertices[e2]);
        // Reversed winding — see doc comment above.
        mesh.indices.push([base, base + 2, base + 1]);
        i += 3;
    }
}
