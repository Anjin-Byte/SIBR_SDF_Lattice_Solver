//! Lattice-aware grid factory.
//!
//! [`grid_spec_for_job`] sizes a [`mesh::GridSpec`] to a `LatticeJob`'s
//! primitive AABB. This lives in `lattice-gen` (not `mesh`) because it
//! depends on the lattice-specific `LatticeJob` type. See
//! [Consumer Crates](../../../SDF_Lattice_Knowledge_Base/Architecture/Patterns/Consumer%20Crates.md)
//! for the architectural rationale: the `mesh` crate is SDF-generic and
//! must not know about lattice types; the lattice-aware glue lives here.

use glam::{UVec3, Vec3};

use crate::error::LatticeError;
use crate::job::LatticeJob;

/// Builds a [`mesh::GridSpec`] tight to a [`LatticeJob`]'s primitive at the
/// given cell size.
///
/// The grid is sized to just contain the primitive's AABB, rounded up to
/// the nearest whole cell on each axis and padded by one cell on each side
/// (so the marching-cubes algorithm sees both sign changes at the boundary).
///
/// # Errors
///
/// Propagates [`LatticeError`] from underlying [`mesh::GridSpec::new`]
/// validation. Given a valid `job` and `cell_size > 0`, this does not fail.
// Cast notes:
// - `cells.{x,y,z}.max(1.0) as u32`: `cells` is `ceil(positive_finite)`,
//   clamped to ≥ 1, so it is a finite non-negative f32. Values up to
//   roughly 10^7 (typical production upper bound) fit in u32 without
//   truncation. Pathological jobs with >10^7 cells per axis would
//   exceed the grid-memory budget long before hitting the u32 limit.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn grid_spec_for_job(job: &LatticeJob, cell_size: f32) -> Result<mesh::GridSpec, LatticeError> {
    let (lo, hi) = job.primitive().aabb();
    // Pad by one cell on each side so the boundary is cleanly captured.
    let padding = Vec3::splat(cell_size);
    let padded_lo = lo - padding;
    let padded_hi = hi + padding;
    let extent = padded_hi - padded_lo;
    // Resolution = ceil(extent / cell_size) on each axis.
    let cells = (extent / cell_size).ceil();
    // Clamp to at least 1 cell per axis; `mesh::GridSpec::new` would reject 0.
    let res = UVec3::new(
        cells.x.max(1.0) as u32,
        cells.y.max(1.0) as u32,
        cells.z.max(1.0) as u32,
    );
    Ok(mesh::GridSpec::new(padded_lo, res, cell_size)?)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp, missing_docs)]
mod tests {
    use super::*;
    use crate::{PrimitiveShape, StrutSpec, UnitCell};

    #[test]
    fn for_job_produces_valid_grid() {
        let job = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
            UnitCell::cubic(1.0).unwrap(),
            StrutSpec::uniform(0.1).unwrap(),
        )
        .unwrap();
        let g = grid_spec_for_job(&job, 0.1).unwrap();
        // Primitive half-extent 2 → extent 4, padded to 4.2, cells = 42.
        assert!(g.resolution().x >= 40);
        assert!(g.cell_size() == 0.1);
    }
}
