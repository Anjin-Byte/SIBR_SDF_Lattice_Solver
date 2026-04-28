//! CPU-reference triangle meshing via Marching Cubes.
//!
//! This module implements the [CPU Reference Path](../../../SDF_Lattice_Knowledge_Base/Architecture/Patterns/CPU%20Reference%20Path.md)
//! for lattice meshing: a slow-but-correct isosurface extraction that serves
//! as the oracle for the later GPU implementation.
//!
//! # Extraction methods
//!
//! Two methods are available, selectable at runtime via [`ExtractionMethod`]:
//!
//! - **Classic Marching Cubes** (Lorensen & Cline, 1987). The simplest,
//!   well-understood algorithm. Produces holes and non-manifold edges at
//!   ~14 ambiguous configurations (~5% of edges for smooth SDFs).
//! - **Marching Cubes 33** (Chernyaev, 1995). Topologically correct.
//!   Currently scaffolded — the disambiguator tables are a bounded
//!   follow-on session; until then, MC33 produces output identical to
//!   classic MC.
//!
//! See [`marching_cubes`] for the implementation and the Domain Knowledge
//! note "Isosurface Extraction Methods" for background.
//!
//! # Pipeline stages
//!
//! - **Extraction** ([`marching_cubes`]) emits an **unwelded** mesh:
//!   every voxel's triangles reference fresh vertex indices, so
//!   `mesh.vertices.len() ≈ 3 * mesh.triangle_count()`. The unwelded
//!   form is the natural output of Marching Cubes and is the one the
//!   library's public extraction API returns.
//! - **Welding** ([`weld`]) deduplicates vertices that share a quantum
//!   bucket and drops triangles that become degenerate after dedup.
//!   CLI callers apply welding before STL export so downstream
//!   consumers see a properly indexed mesh. The library's `mesh_with` /
//!   `run` stay unwelded so callers retain explicit control.
//!
//! # Other limitations (unchanged from earlier phases)
//!
//! - **Dense grid only, single-threaded.** This is the reference
//!   implementation — the `gpu` crate will provide the sparse, parallel
//!   production version. For Woodward Fig. 4B-scale parts (~1M voxels),
//!   this runs in seconds; for Fig. 4A-scale (~100M voxels), it is not
//!   performant.

pub mod export;
pub mod grid;
pub mod marching_cubes;
pub mod smooth;
pub mod subdivide;
pub mod weld;

pub use export::Format;
pub use grid::GridSpec;
pub use marching_cubes::ExtractionMethod;
pub use smooth::{SmoothError, TaubinParams, taubin, taubin_with_progress};
pub use subdivide::{ButterflyParams, butterfly, butterfly_with_progress};
pub use weld::weld_by_position;

use glam::Vec3;

/// A triangle mesh produced by [`mesh`].
///
/// Triangles are stored as index triples into `vertices`. No vertex welding
/// is performed — see the module-level docs.
#[derive(Debug, Clone, Default)]
pub struct Mesh {
    /// Vertex positions.
    pub vertices: Vec<Vec3>,
    /// Triangle indices (each triple references three `vertices`).
    pub indices: Vec<[u32; 3]>,
}

impl Mesh {
    /// Returns the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len()
    }

    /// Returns the number of vertices (pre-welding; may include duplicates).
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Returns true if the mesh has no triangles.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}

/// Meshes the lattice body for a validated job at the given grid resolution,
/// using the caller-selected [`ExtractionMethod`].
///
/// Runs marching cubes on a dense voxel grid, evaluating the composed
/// lattice-body SDF at every grid point. The output is an unwelded triangle
/// list — see the module-level docs for limitations.
///
/// # Performance
///
/// This is the CPU reference implementation, single-threaded. Cost is
/// `O(k · N³)` for an `N³`-resolution grid with `k` struts per cell. For
/// Woodward Fig. 4B-scale parts (~200³ voxels), runs in seconds; for
/// Fig. 4A-scale (~500³+ voxels), consider the production GPU path.
pub fn mesh_with(job: &crate::LatticeJob, grid: &GridSpec, method: ExtractionMethod) -> Mesh {
    mesh_with_progress(job, grid, method, &mut ())
}

/// Like [`mesh_with`], but reports progress to `progress`. See
/// [`marching_cubes::run_with_progress`] for the tick schedule.
pub fn mesh_with_progress(
    job: &crate::LatticeJob,
    grid: &GridSpec,
    method: ExtractionMethod,
    progress: &mut impl crate::Progress,
) -> Mesh {
    let body = crate::lattice_body(job);
    marching_cubes::run_with_progress(&body, grid, method, progress)
}

/// Convenience wrapper around [`mesh_with`] that uses
/// [`ExtractionMethod::default()`] (classic MC).
///
/// Preserves the signature used by earlier phases before MC33 selection
/// landed.
pub fn mesh(job: &crate::LatticeJob, grid: &GridSpec) -> Mesh {
    mesh_with(job, grid, ExtractionMethod::default())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, missing_docs)]
mod tests {
    use super::*;

    #[test]
    fn mesh_default_is_empty() {
        let m = Mesh::default();
        assert!(m.is_empty());
        assert_eq!(m.triangle_count(), 0);
        assert_eq!(m.vertex_count(), 0);
    }
}
