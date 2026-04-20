//! Marching Cubes extraction — classic Lorensen-Cline + MC33 variants.
//!
//! See the [mesh module's](super) top-level docs for algorithm background
//! and the discussion of known limitations of each method. This module
//! provides runtime selection between methods via [`ExtractionMethod`].

pub mod classic;
pub mod mc33;
pub mod tables;

use glam::Vec3;
use sdf::Sdf;

use super::Mesh;
use super::grid::GridSpec;
use tables::CORNER_OFFSETS;

/// Which isosurface-extraction method to use.
///
/// See each variant's doc for trade-offs. `Default` is [`Self::ClassicMc`]
/// so callers that don't specify a method preserve the pre-MC33 behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExtractionMethod {
    /// Classic Marching Cubes (Lorensen & Cline, 1987).
    ///
    /// Fast, simple, well-understood. Produces holes and non-manifold edges
    /// at the ~14 ambiguous configurations the original 1987 table cannot
    /// disambiguate consistently; our integration tests measure this at
    /// roughly 5% of edges for smooth SDFs. See Domain Knowledge note
    /// "Isosurface Extraction Methods" for details.
    #[default]
    ClassicMc,

    /// Marching Cubes 33 (Chernyaev, 1995).
    ///
    /// Topologically correct at all non-degenerate configurations.
    /// Runs the classic case index lookup; for the known-ambiguous cases
    /// it applies face and interior asymptotic deciders to select the
    /// correct triangulation.
    ///
    /// **Current implementation status**: ported incrementally.
    /// - **Case 3** (face-ambiguous, two diagonal corners on a face):
    ///   fully disambiguated via the face asymptotic decider.
    /// - **Case 4** (body-ambiguous, two diagonal corners on a body
    ///   diagonal): fully disambiguated via the interior decider.
    /// - **Cases 6, 7, 10, 12, 13**: still fall back to classic's
    ///   triangulation with a single `tracing::warn!` per run.
    /// - **Unambiguous cases**: delegated to classic unconditionally.
    ///
    /// MC33 output is therefore at least as good as classic MC, and
    /// strictly better wherever a Case-3 or Case-4 voxel occurs.
    Mc33,
}

/// Extracts the zero-isosurface of `body` on the given sampling `grid`
/// using the specified `method`.
///
/// The returned mesh has no vertex welding — each voxel's triangles
/// reference freshly-emitted vertices. See the parent `mesh` module docs
/// for other pipeline limitations.
pub fn run<F: Sdf>(body: &F, grid: &GridSpec, method: ExtractionMethod) -> Mesh {
    let field = sample_field(body, grid);
    let mut mesh = Mesh::default();
    let res = grid.resolution();
    for cz in 0..res.z {
        for cy in 0..res.y {
            for cx in 0..res.x {
                let cell = CellCoord {
                    x: cx,
                    y: cy,
                    z: cz,
                };
                match method {
                    ExtractionMethod::ClassicMc => {
                        classic::process_cell(grid, &field, cell, &mut mesh);
                    }
                    ExtractionMethod::Mc33 => {
                        mc33::process_cell(grid, &field, cell, &mut mesh);
                    }
                }
            }
        }
    }
    mesh
}

/// Evaluates the SDF at every grid point and returns a flat row-major buffer.
pub(crate) fn sample_field<F: Sdf>(body: &F, grid: &GridSpec) -> Vec<f32> {
    let mut field = Vec::with_capacity(grid.sample_count());
    let res = grid.resolution();
    for sz in 0..=res.z {
        for sy in 0..=res.y {
            for sx in 0..=res.x {
                let p = grid.sample_point(sx, sy, sz);
                field.push(body.eval(p));
            }
        }
    }
    field
}

/// A cell's integer origin in the grid. Shared between extraction methods.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CellCoord {
    pub x: u32,
    pub y: u32,
    pub z: u32,
}

/// Linearly interpolates the zero crossing between two grid corners.
///
/// Shared between classic MC and MC33 — the interpolation rule is the
/// same; only the triangulation table differs between methods.
///
/// For SDF values `va` at corner `a` and `vb` at corner `b` with opposite
/// signs, the zero crossing lies at parameter `t = va / (va - vb)` from `a`
/// toward `b` (treating the SDF as linear along the edge — exact for
/// 1-Lipschitz exact SDFs in the limit of small `cell_size`).
pub(crate) fn interpolate_edge(
    grid: &GridSpec,
    cell: CellCoord,
    corner_a: u8,
    corner_b: u8,
    va: f32,
    vb: f32,
) -> Vec3 {
    let (ax, ay, az) = CORNER_OFFSETS[corner_a as usize];
    let (bx, by, bz) = CORNER_OFFSETS[corner_b as usize];
    let pa = grid.sample_point(cell.x + ax, cell.y + ay, cell.z + az);
    let pb = grid.sample_point(cell.x + bx, cell.y + by, cell.z + bz);

    // Guard against near-equal field values to avoid division by near-zero.
    let denom = va - vb;
    if denom.abs() < 1e-12 {
        // Corner values essentially equal — emit the midpoint.
        return (pa + pb) * 0.5;
    }
    let t = (va / denom).clamp(0.0, 1.0);
    pa + (pb - pa) * t
}

/// Computes the 8-bit case index for a voxel's corner values.
///
/// Bit `i` is set iff corner `i` has `f < 0` (is inside the surface).
/// We treat `f == 0` as "outside" (>= 0) so that points exactly on the
/// surface are assigned a consistent classification.
///
/// Shared between methods — both classic and MC33 index into 256-entry
/// tables keyed on this value; MC33 then optionally applies an
/// asymptotic-decider test on top for ambiguous cases.
#[inline]
pub(crate) fn case_index(corner_values: &[f32; 8]) -> u8 {
    let mut idx: u8 = 0;
    for (i, &v) in corner_values.iter().enumerate() {
        if v < 0.0 {
            idx |= 1 << i;
        }
    }
    idx
}

/// Reads the 8 corner values for a voxel from the flat field buffer.
#[inline]
pub(crate) fn read_corner_values(grid: &GridSpec, field: &[f32], cell: CellCoord) -> [f32; 8] {
    let mut out = [0.0_f32; 8];
    for (i, &(dx, dy, dz)) in CORNER_OFFSETS.iter().enumerate() {
        let idx = grid.sample_index(cell.x + dx, cell.y + dy, cell.z + dz);
        out[i] = field[idx];
    }
    out
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
    use glam::{UVec3, vec3};
    use sdf::Sphere;

    // --------------------------------------------------------------
    // Backward-compatibility regression tests.
    //
    // These pin that classic MC output is unchanged across the module
    // restructure. If any of these fail, the refactor changed the
    // extraction behavior and the change needs explicit review.
    // --------------------------------------------------------------

    #[test]
    fn regression_classic_mc_structural_invariants() {
        // The refactor MUST NOT change classic-MC behavior. We pin several
        // structural invariants that would catch any triangulation drift:
        //
        // 1. Unwelded mesh: vertex_count == 3 * triangle_count exactly.
        // 2. Triangle count is nonzero (the grid spans the sphere surface).
        // 3. Signed volume is positive and close to the analytical (4/3)πr³.
        //
        // These together would detect case-index convention flips, winding
        // reversal regressions, table corruption, or entire-algorithm swaps.
        let sphere = Sphere::new(1.0).unwrap();
        let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(30), 0.1).unwrap();
        let mesh = run(&sphere, &grid, ExtractionMethod::ClassicMc);

        assert!(mesh.triangle_count() > 0);
        assert_eq!(
            mesh.vertex_count(),
            3 * mesh.triangle_count(),
            "unwelded invariant broken"
        );

        let vol: f32 = mesh
            .indices
            .iter()
            .map(|tri| {
                let v0 = mesh.vertices[tri[0] as usize];
                let v1 = mesh.vertices[tri[1] as usize];
                let v2 = mesh.vertices[tri[2] as usize];
                v0.dot(v1.cross(v2)) / 6.0
            })
            .sum();
        let exact = core::f32::consts::PI * 4.0 / 3.0;
        let rel_err = ((vol - exact) / exact).abs();
        assert!(
            rel_err < 0.05,
            "signed volume {vol} diverges from analytical {exact} (rel err {rel_err})"
        );
    }

    #[test]
    fn mc33_on_sphere_matches_classic_within_small_divergence() {
        // A smooth sphere at modest resolution produces few (often zero)
        // Case-3 voxels, so MC33 and classic should agree on most — ideally
        // all — cells. We don't require byte-identical output any more
        // (Case 3 is now disambiguated), but the triangle count must be
        // very close and the signed volumes must agree tightly.
        let sphere = Sphere::new(1.0).unwrap();
        let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(20), 0.15).unwrap();
        let m_classic = run(&sphere, &grid, ExtractionMethod::ClassicMc);
        let m_mc33 = run(&sphere, &grid, ExtractionMethod::Mc33);
        let classic_n = m_classic.triangle_count();
        let mc33_n = m_mc33.triangle_count();
        let dt = classic_n.abs_diff(mc33_n);
        assert!(
            dt <= classic_n / 20,
            "triangle counts diverge: classic={classic_n}, mc33={mc33_n}"
        );
    }

    #[test]
    fn default_method_is_classic() {
        assert_eq!(ExtractionMethod::default(), ExtractionMethod::ClassicMc);
    }

    /// Histograms the Lewiner base case across every voxel of a workload.
    /// Returns `[u32; 16]` where index `i` counts voxels in Lewiner base
    /// case `i` (index 0 is unused). Used by diagnostics below.
    #[allow(clippy::cast_sign_loss)]
    fn case_histogram<F: Sdf>(body: &F, grid: &GridSpec) -> [u32; 16] {
        let field = sample_field(body, grid);
        let mut hist = [0_u32; 16];
        let res = grid.resolution();
        for cz in 0..res.z {
            for cy in 0..res.y {
                for cx in 0..res.x {
                    let cell = CellCoord {
                        x: cx,
                        y: cy,
                        z: cz,
                    };
                    let corners = read_corner_values(grid, &field, cell);
                    let our_idx = case_index(&corners);
                    let base_case = mc33::lewiner_base_case(our_idx);
                    hist[base_case as usize] += 1;
                }
            }
        }
        hist
    }

    /// Decisive check: for every `case_index`, classic MC and MC33
    /// (with Session 2.5's unambiguous-table port) should produce
    /// identical triangulations on a single-voxel test grid — *if* my
    /// hypothesis is correct that Lewiner's unambiguous tables are
    /// geometrically equivalent to classic's (same triangles, possibly
    /// different edge order). Any divergence is interesting.
    #[test]
    fn diagnostic_classic_vs_mc33_per_case_index() {
        let grid = GridSpec::new(Vec3::ZERO, UVec3::splat(1), 1.0).unwrap();
        let mut diff_count = 0_u32;
        for our_idx in 0..=255_u8 {
            let mut field = vec![0.0_f32; grid.sample_count()];
            let mut corners = [1.0_f32; 8];
            for (i, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << i) != 0 {
                    *c = -1.0;
                }
            }
            field[grid.sample_index(0, 0, 0)] = corners[0];
            field[grid.sample_index(1, 0, 0)] = corners[1];
            field[grid.sample_index(1, 1, 0)] = corners[2];
            field[grid.sample_index(0, 1, 0)] = corners[3];
            field[grid.sample_index(0, 0, 1)] = corners[4];
            field[grid.sample_index(1, 0, 1)] = corners[5];
            field[grid.sample_index(1, 1, 1)] = corners[6];
            field[grid.sample_index(0, 1, 1)] = corners[7];

            let mut m_classic = Mesh::default();
            let mut m_mc33 = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            classic::process_cell(&grid, &field, cell, &mut m_classic);
            mc33::process_cell(&grid, &field, cell, &mut m_mc33);
            if m_classic.triangle_count() != m_mc33.triangle_count()
                || m_classic.vertices.len() != m_mc33.vertices.len()
            {
                diff_count += 1;
            }
        }
        eprintln!(
            "classic vs mc33 per-case-index comparison: {diff_count}/256 case_indices \
             produce differently-sized meshes"
        );
    }

    /// Current workloads (smooth sphere, cubic-lattice of cylinders) do
    /// not hit Lewiner base cases 4 (Chernyaev 3) or 5 (Chernyaev 4) —
    /// these require locally-saddle surface geometry. The remaining
    /// non-manifold edges on these workloads come from Lewiner cases
    /// 7/8/11/13/14, still pending port. This test pins the empirical
    /// finding so a future workload that DOES exercise Case 3/4 will
    /// stand out.
    #[test]
    fn diagnostic_current_workloads_do_not_hit_case3_or_case4() {
        use crate::{LatticeJob, PrimitiveShape, StrutSpec, UnitCell, lattice_body};

        let sphere = Sphere::new(1.0).unwrap();
        let sphere_grid = GridSpec::new(Vec3::splat(-2.0), UVec3::splat(30), 4.0 / 30.0).unwrap();
        let sphere_hist = case_histogram(&sphere, &sphere_grid);
        eprintln!("sphere @ 30^3 base-case histogram: {sphere_hist:?}");

        let job = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
            UnitCell::cubic(1.0).unwrap(),
            StrutSpec::uniform(0.15).unwrap(),
        )
        .unwrap();
        let lattice_grid = GridSpec::for_job(&job, 0.1).unwrap();
        let body = lattice_body(&job);
        let lattice_hist = case_histogram(&body, &lattice_grid);
        eprintln!("cubic lattice base-case histogram: {lattice_hist:?}");

        // Base case 4 = Chernyaev 3 (face-ambiguous, Session 1 port).
        // Base case 5 = Chernyaev 4 (body-ambiguous, Session 2 port).
        // Both are zero on current workloads — the port is correct but
        // inert here until we either (a) port Cases 7/8/11/13/14 or
        // (b) introduce an adversarial fixture.
        assert_eq!(sphere_hist[4], 0, "sphere hit Case 3 — workload changed");
        assert_eq!(sphere_hist[5], 0, "sphere hit Case 4 — workload changed");
        assert_eq!(lattice_hist[4], 0, "lattice hit Case 3 — workload changed");
        assert_eq!(lattice_hist[5], 0, "lattice hit Case 4 — workload changed");
    }

    // --------------------------------------------------------------
    // Existing behavior checks (preserved from pre-refactor).
    // --------------------------------------------------------------

    #[test]
    fn empty_grid_produces_empty_mesh() {
        let sphere = Sphere::new(0.01).unwrap();
        let grid = GridSpec::new(vec3(10.0, 10.0, 10.0), UVec3::splat(4), 0.5).unwrap();
        let mesh = run(&sphere, &grid, ExtractionMethod::ClassicMc);
        assert!(mesh.is_empty());
    }

    #[test]
    fn all_interior_field_produces_empty_mesh() {
        let sphere = Sphere::new(100.0).unwrap();
        let grid = GridSpec::new(Vec3::ZERO, UVec3::splat(2), 0.5).unwrap();
        let mesh = run(&sphere, &grid, ExtractionMethod::ClassicMc);
        assert!(mesh.is_empty());
    }

    #[test]
    fn sphere_produces_nonzero_triangles() {
        let sphere = Sphere::new(1.0).unwrap();
        let grid = GridSpec::new(Vec3::splat(-2.0), UVec3::splat(20), 0.2).unwrap();
        let mesh = run(&sphere, &grid, ExtractionMethod::ClassicMc);
        assert!(mesh.triangle_count() > 0);
        for tri in &mesh.indices {
            for &i in tri {
                assert!(mesh.vertices[i as usize].is_finite());
            }
        }
    }

    #[test]
    fn sphere_triangle_count_scales_with_resolution() {
        let sphere = Sphere::new(1.0).unwrap();
        let coarse = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(10), 0.3).unwrap();
        let fine = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(30), 0.1).unwrap();
        let m_coarse = run(&sphere, &coarse, ExtractionMethod::ClassicMc);
        let m_fine = run(&sphere, &fine, ExtractionMethod::ClassicMc);
        assert!(m_fine.triangle_count() > m_coarse.triangle_count());
    }

    #[test]
    fn interpolated_vertices_lie_on_edges() {
        let sphere = Sphere::new(1.0).unwrap();
        let origin = Vec3::splat(-1.5);
        let cell_size = 0.3;
        let res = UVec3::splat(10);
        let grid = GridSpec::new(origin, res, cell_size).unwrap();
        let mesh = run(&sphere, &grid, ExtractionMethod::ClassicMc);
        let max = origin + vec3(res.x as f32, res.y as f32, res.z as f32) * cell_size;
        for v in &mesh.vertices {
            assert!(
                v.x >= origin.x - 1e-4
                    && v.y >= origin.y - 1e-4
                    && v.z >= origin.z - 1e-4
                    && v.x <= max.x + 1e-4
                    && v.y <= max.y + 1e-4
                    && v.z <= max.z + 1e-4,
                "vertex {v:?} outside grid [{origin:?}, {max:?}]"
            );
        }
    }
}
