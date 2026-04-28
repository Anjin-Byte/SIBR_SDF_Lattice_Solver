//! `GridSpec` — the sampling grid for marching cubes.
//!
//! A grid is an axis-aligned box subdivided into `resolution` cubic cells
//! per axis. Field evaluation occurs at `(resolution + 1)^3` grid *points*
//! (the corners of those cells). Each cell is inspected by the MC algorithm
//! to produce 0-5 triangles.

use glam::{UVec3, Vec3};

use crate::error::LatticeError;
use crate::job::LatticeJob;

/// A sampling grid for marching cubes.
///
/// Constructed via [`GridSpec::new`] or [`GridSpec::for_job`]. All
/// constructors validate that `cell_size > 0`, `resolution` is non-zero on
/// every axis, and `origin` is finite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridSpec {
    origin: Vec3,
    resolution: UVec3,
    cell_size: f32,
}

impl GridSpec {
    /// Constructs a `GridSpec` with the given origin, resolution, and cell size.
    ///
    /// The grid's sampled region extends from `origin` to
    /// `origin + resolution.as_vec3() * cell_size` on each axis.
    ///
    /// # Errors
    ///
    /// - [`LatticeError::Sdf`] wrapping [`sdf::BuildError::NonFinite`] if
    ///   any component of `origin` is `NaN` or `±∞`.
    /// - [`LatticeError::Sdf`] wrapping [`sdf::BuildError::NonPositive`] if
    ///   `cell_size <= 0`.
    /// - [`LatticeError::Sdf`] wrapping [`sdf::BuildError::Degenerate`] if
    ///   any component of `resolution` is zero.
    pub fn new(origin: Vec3, resolution: UVec3, cell_size: f32) -> Result<Self, LatticeError> {
        if !origin.is_finite() {
            return Err(sdf::BuildError::NonFinite {
                field: "grid.origin",
                value: 0.0,
            }
            .into());
        }
        if !cell_size.is_finite() {
            return Err(sdf::BuildError::NonFinite {
                field: "grid.cell_size",
                value: cell_size,
            }
            .into());
        }
        if cell_size <= 0.0 {
            return Err(sdf::BuildError::NonPositive {
                field: "grid.cell_size",
                value: cell_size,
            }
            .into());
        }
        if resolution.x == 0 || resolution.y == 0 || resolution.z == 0 {
            return Err(sdf::BuildError::Degenerate {
                reason: "grid resolution must be non-zero on every axis",
            }
            .into());
        }
        Ok(Self {
            origin,
            resolution,
            cell_size,
        })
    }

    /// Builds a grid tight to a `LatticeJob`'s primitive at the given cell size.
    ///
    /// The grid is sized to just contain the primitive's AABB, rounded up to
    /// the nearest whole cell on each axis and padded by one cell on each
    /// side (so the MC algorithm sees both sign changes at the boundary).
    ///
    /// # Errors
    ///
    /// Propagates [`LatticeError`] from underlying validation. Given a
    /// valid `job` and `cell_size > 0`, this does not fail.
    // Cast notes:
    // - `cells.{x,y,z}.max(1.0) as u32`: `cells` is `ceil(positive_finite)`,
    //   clamped to ≥ 1, so it is a finite non-negative f32. Values up to
    //   roughly 10^7 (typical production upper bound) fit in u32 without
    //   truncation. Pathological jobs with >10^7 cells per axis would
    //   exceed the grid-memory budget long before hitting the u32 limit.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub fn for_job(job: &LatticeJob, cell_size: f32) -> Result<Self, LatticeError> {
        let (lo, hi) = job.primitive().aabb();
        // Pad by one cell on each side so the boundary is cleanly captured.
        let padding = Vec3::splat(cell_size);
        let padded_lo = lo - padding;
        let padded_hi = hi + padding;
        let extent = padded_hi - padded_lo;
        // Resolution = ceil(extent / cell_size) on each axis.
        let cells = (extent / cell_size).ceil();
        // Clamp to at least 1 cell per axis; `new` would reject 0.
        let res = UVec3::new(
            cells.x.max(1.0) as u32,
            cells.y.max(1.0) as u32,
            cells.z.max(1.0) as u32,
        );
        Self::new(padded_lo, res, cell_size)
    }

    /// Returns the grid origin (the lower corner).
    pub fn origin(&self) -> Vec3 {
        self.origin
    }

    /// Returns the per-axis cell resolution.
    pub fn resolution(&self) -> UVec3 {
        self.resolution
    }

    /// Returns the cubic cell size.
    pub fn cell_size(&self) -> f32 {
        self.cell_size
    }

    /// Returns the world-space position of the grid point at integer
    /// indices `(sx, sy, sz)`. Valid for `sx <= resolution.x`, etc.
    ///
    /// Casts `u32 → f32` lose precision for values > 2²³ ≈ 8.4M, which
    /// exceeds any reasonable dense-grid resolution (the grid memory would
    /// be many terabytes long before this). Accepted as a documented
    /// limitation of dense grids.
    #[inline]
    #[allow(clippy::cast_precision_loss)]
    pub fn sample_point(&self, sx: u32, sy: u32, sz: u32) -> Vec3 {
        self.origin + Vec3::new(sx as f32, sy as f32, sz as f32) * self.cell_size
    }

    /// Returns the total number of grid sample points:
    /// `(resolution.x + 1) * (resolution.y + 1) * (resolution.z + 1)`.
    pub fn sample_count(&self) -> usize {
        ((self.resolution.x + 1) as usize)
            * ((self.resolution.y + 1) as usize)
            * ((self.resolution.z + 1) as usize)
    }

    /// Returns the linear index into a flat field buffer for grid point
    /// `(sx, sy, sz)`. Row-major: x varies fastest, then y, then z.
    #[inline]
    pub fn sample_index(&self, sx: u32, sy: u32, sz: u32) -> usize {
        let stride_y = (self.resolution.x + 1) as usize;
        let stride_z = stride_y * ((self.resolution.y + 1) as usize);
        (sx as usize) + (sy as usize) * stride_y + (sz as usize) * stride_z
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp, missing_docs)]
mod tests {
    use super::*;
    use crate::{PrimitiveShape, StrutSpec, UnitCell};

    #[test]
    fn new_accepts_valid_input() {
        let g = GridSpec::new(Vec3::ZERO, UVec3::splat(4), 0.5).unwrap();
        assert_eq!(g.origin(), Vec3::ZERO);
        assert_eq!(g.resolution(), UVec3::splat(4));
        assert_eq!(g.cell_size(), 0.5);
    }

    #[test]
    fn new_rejects_zero_resolution_component() {
        assert!(GridSpec::new(Vec3::ZERO, UVec3::new(0, 4, 4), 0.5).is_err());
        assert!(GridSpec::new(Vec3::ZERO, UVec3::new(4, 0, 4), 0.5).is_err());
        assert!(GridSpec::new(Vec3::ZERO, UVec3::new(4, 4, 0), 0.5).is_err());
    }

    #[test]
    fn new_rejects_non_positive_cell_size() {
        assert!(GridSpec::new(Vec3::ZERO, UVec3::splat(4), 0.0).is_err());
        assert!(GridSpec::new(Vec3::ZERO, UVec3::splat(4), -1.0).is_err());
    }

    #[test]
    fn new_rejects_nan_origin() {
        assert!(GridSpec::new(Vec3::new(f32::NAN, 0.0, 0.0), UVec3::splat(4), 0.5).is_err());
    }

    #[test]
    fn new_rejects_nan_cell_size() {
        assert!(GridSpec::new(Vec3::ZERO, UVec3::splat(4), f32::NAN).is_err());
    }

    #[test]
    fn sample_point_is_origin_at_zero_index() {
        let origin = Vec3::new(1.0, 2.0, 3.0);
        let g = GridSpec::new(origin, UVec3::splat(4), 0.5).unwrap();
        assert_eq!(g.sample_point(0, 0, 0), origin);
    }

    #[test]
    fn sample_point_at_last_index_is_origin_plus_extent() {
        let origin = Vec3::ZERO;
        let g = GridSpec::new(origin, UVec3::new(4, 4, 4), 0.5).unwrap();
        // Grid spans 4 cells × 0.5 per cell = 2.0 per axis.
        assert_eq!(g.sample_point(4, 4, 4), Vec3::splat(2.0));
    }

    #[test]
    fn sample_count_is_nplus1_cubed() {
        let g = GridSpec::new(Vec3::ZERO, UVec3::splat(4), 0.5).unwrap();
        // 5 × 5 × 5 = 125.
        assert_eq!(g.sample_count(), 125);
    }

    #[test]
    fn sample_index_is_row_major() {
        let g = GridSpec::new(Vec3::ZERO, UVec3::new(4, 4, 4), 0.5).unwrap();
        assert_eq!(g.sample_index(0, 0, 0), 0);
        assert_eq!(g.sample_index(1, 0, 0), 1);
        assert_eq!(g.sample_index(0, 1, 0), 5); // stride_y = 5
        assert_eq!(g.sample_index(0, 0, 1), 25); // stride_z = 25
        assert_eq!(g.sample_index(4, 4, 4), 124);
    }

    #[test]
    fn for_job_produces_valid_grid() {
        let job = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
            UnitCell::cubic(1.0).unwrap(),
            StrutSpec::uniform(0.1).unwrap(),
        )
        .unwrap();
        let g = GridSpec::for_job(&job, 0.1).unwrap();
        // Primitive half-extent 2 → extent 4, padded to 4.2, cells = 42.
        assert!(g.resolution().x >= 40);
        assert!(g.cell_size() == 0.1);
    }
}
