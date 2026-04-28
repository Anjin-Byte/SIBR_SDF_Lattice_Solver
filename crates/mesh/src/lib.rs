//! SDF → triangle mesh: CPU-reference Marching Cubes, vertex welding, Taubin
//! smoothing, butterfly subdivision, and feature-gated format encoders.
//!
//! This crate is an SDF-generic [Consumer Crate](../../../SDF_Lattice_Knowledge_Base/Architecture/Patterns/Consumer%20Crates.md):
//! it depends only on [`sdf`] and operates on any `impl Sdf`. The lattice
//! domain is one consumer; CSG trees, procedural fields, and any other SDF
//! source can drive these algorithms equally well.
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
//!   Callers apply welding before export so downstream consumers see a
//!   properly indexed mesh.
//! - **Post-processing** ([`postprocess`]) groups the optional
//!   mesh-quality stages: Taubin smoothing reduces noise without
//!   shrinkage; butterfly subdivision densifies the mesh while
//!   preserving original vertex positions. Each is independent — a
//!   caller skips a stage simply by not invoking it.
//! - **Export** ([`export`]) writes the welded mesh to a printer-friendly
//!   format. STL and OBJ are gated behind the `stl` and `obj` Cargo
//!   features (both default-on); 3MF is reserved for a future feature.
//!
//! # Other limitations (unchanged from earlier phases)
//!
//! - **Dense grid only, single-threaded.** This is the reference
//!   implementation — the `gpu` crate will provide the sparse, parallel
//!   production version. For Woodward Fig. 4B-scale parts (~1M voxels),
//!   this runs in seconds; for Fig. 4A-scale (~100M voxels), it is not
//!   performant.
//!
//! # Precision contract
//!
//! Algorithms in this crate accept any `impl Sdf` — they require only
//! `Sdf::eval`, never `BoundSdf` or `ExactSdf`. That is deliberate: the
//! lattice-body composition propagates `Sdf`-only upward through smooth
//! combinators, and consumers downstream of those combinators (Marching
//! Cubes, slicing, rendering) only need point evaluation. See the vault's
//! "SDF Primitives Catalog" bound-contamination section.

pub mod export;
pub mod grid;
pub mod marching_cubes;
pub mod postprocess;
pub mod progress;
pub mod weld;

pub use export::Format;
pub use grid::GridSpec;
pub use marching_cubes::ExtractionMethod;
pub use postprocess::smooth::{SmoothError, TaubinParams, taubin, taubin_with_progress};
pub use postprocess::subdivide::{ButterflyParams, butterfly, butterfly_with_progress};
pub use progress::Progress;
pub use weld::weld_by_position;

use glam::Vec3;
use sdf::Sdf;

/// A triangle mesh produced by [`mesh`].
///
/// Triangles are stored as index triples into `vertices`. No vertex welding
/// is performed by extraction — see the module-level docs.
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

/// Meshes any SDF at the given grid resolution using the caller-selected
/// [`ExtractionMethod`].
///
/// Runs marching cubes on a dense voxel grid, evaluating `sdf` at every
/// grid point. The output is an unwelded triangle list — see the
/// module-level docs for limitations.
///
/// # Performance
///
/// This is the CPU reference implementation, single-threaded. Cost is
/// `O(N³ · sdf_eval)` for an `N³`-resolution grid. For Woodward Fig. 4B-scale
/// parts (~200³ voxels), runs in seconds; for Fig. 4A-scale (~500³+ voxels),
/// consider the production GPU path.
pub fn mesh_with<S: Sdf>(sdf: &S, grid: &GridSpec, method: ExtractionMethod) -> Mesh {
    mesh_with_progress(sdf, grid, method, &mut ())
}

/// Like [`mesh_with`], but reports progress to `progress`. See
/// [`marching_cubes::run_with_progress`] for the tick schedule.
pub fn mesh_with_progress<S: Sdf, P: Progress>(
    sdf: &S,
    grid: &GridSpec,
    method: ExtractionMethod,
    progress: &mut P,
) -> Mesh {
    marching_cubes::run_with_progress(sdf, grid, method, progress)
}

/// Convenience wrapper around [`mesh_with`] that uses
/// [`ExtractionMethod::default()`] (classic MC).
///
/// Preserves the signature used by earlier phases before MC33 selection
/// landed.
pub fn mesh<S: Sdf>(sdf: &S, grid: &GridSpec) -> Mesh {
    mesh_with(sdf, grid, ExtractionMethod::default())
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
