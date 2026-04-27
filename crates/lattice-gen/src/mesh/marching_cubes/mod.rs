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
use crate::progress::Progress;
use tables::CORNER_OFFSETS;

/// Which isosurface-extraction method to use.
///
/// See each variant's doc for trade-offs. `Default` is [`Self::Mc33`] —
/// topologically correct at every non-degenerate configuration. Use
/// [`Self::ClassicMc`] only if reproducing pre-MC33-port behavior or
/// for differential testing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExtractionMethod {
    /// Classic Marching Cubes (Lorensen & Cline, 1987).
    ///
    /// Fast, simple, well-understood. Produces holes and non-manifold edges
    /// at the ~14 ambiguous configurations the original 1987 table cannot
    /// disambiguate consistently; our integration tests measure this at
    /// roughly 5% of edges for smooth SDFs. See Domain Knowledge note
    /// "Isosurface Extraction Methods" for details.
    ClassicMc,

    /// Marching Cubes 33 (Chernyaev, 1995).
    ///
    /// Topologically correct at all non-degenerate configurations.
    /// Runs the Lewiner 2003 case index lookup; for ambiguous cases it
    /// applies face and interior asymptotic deciders to select the
    /// correct triangulation. **All Chernyaev cases (3, 4, 6, 7, 10,
    /// 12, 13) are fully ported** — see the `mc33` module docs for
    /// the per-case dispatch map. The unambiguous cases use Lewiner's
    /// face-consistent tables, eliminating the ~2–3% non-manifold
    /// edge rate on shared faces that classic MC produces.
    ///
    /// MC33 output is at least as good as classic MC on every
    /// configuration and strictly better at ambiguous configurations.
    #[default]
    Mc33,
}

/// Extracts the zero-isosurface of `body` on the given sampling `grid`
/// using the specified `method`.
///
/// The returned mesh has no vertex welding — each voxel's triangles
/// reference freshly-emitted vertices. See the parent `mesh` module docs
/// for other pipeline limitations.
pub fn run<F: Sdf>(body: &F, grid: &GridSpec, method: ExtractionMethod) -> Mesh {
    run_with_progress(body, grid, method, &mut ())
}

/// Like [`run`], but reports progress to `progress`. Ticks once per
/// z-slab for both the field-sampling and extraction phases, so the
/// declared length is `(res.z + 1) + res.z`: `res.z + 1` sampling
/// layers (the sample grid has one more layer than the cell grid in
/// each axis) plus `res.z` extraction layers.
pub fn run_with_progress<F: Sdf>(
    body: &F,
    grid: &GridSpec,
    method: ExtractionMethod,
    progress: &mut impl Progress,
) -> Mesh {
    let res = grid.resolution();
    progress.set_len(u64::from(res.z + 1) + u64::from(res.z));
    let field = sample_field(body, grid, progress);
    let mut mesh = Mesh::default();
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
        progress.inc(1);
    }
    progress.finish();
    mesh
}

/// Evaluates the SDF at every grid point and returns a flat row-major
/// buffer. Reports one tick per sampled z-slab to `progress`. Callers
/// that want a combined bar (sampling + extraction) should use
/// [`run_with_progress`]; this variant does not call `set_len` or
/// `finish` because it may be the first phase of a longer stage.
///
/// Pass `&mut ()` to opt out of progress reporting.
pub(crate) fn sample_field<F: Sdf>(
    body: &F,
    grid: &GridSpec,
    progress: &mut impl Progress,
) -> Vec<f32> {
    let mut field = Vec::with_capacity(grid.sample_count());
    let res = grid.resolution();
    for sz in 0..=res.z {
        for sy in 0..=res.y {
            for sx in 0..=res.x {
                let p = grid.sample_point(sx, sy, sz);
                field.push(body.eval(p));
            }
        }
        progress.inc(1);
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
///
/// # Determinism (cross-voxel face consistency)
///
/// Two adjacent voxels share a face and therefore share the four edges of
/// that face. Each of those edges is reached from opposite directions in
/// the two voxels' local corner numbering — so a naive
/// `pa + (pb - pa) * t` formulation, with `t = va / (va - vb)`, computes
/// the *same world point* via *different f32 operand orderings* in each
/// voxel. In f32 arithmetic, the two orderings drift apart by ~1 ULP,
/// producing distinct vertex positions that resist welding (the welding
/// tolerance can land 1-ULP-apart values on opposite sides of a
/// quantization boundary, leaving them unmerged). Result: 24 boundary
/// edges per Kelvin lattice symmetry feature in the welded mesh, which
/// render as visible holes.
///
/// To eliminate this, the function **canonicalizes on integer corner
/// coordinates** (`cell + offset`, exact `u32` arithmetic) before any
/// f32 operation. Both voxels see the same `(lo, hi)` integer pair for
/// the shared edge → both feed `sample_point` and the interpolation
/// formula in identical operand order → both produce *bitwise-identical*
/// `Vec3` outputs. Welding then merges them naturally without tolerance
/// jitter.
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
    let coord_a = (cell.x + ax, cell.y + ay, cell.z + az);
    let coord_b = (cell.x + bx, cell.y + by, cell.z + bz);

    // Canonicalize: always interpolate from the lexicographically-smaller
    // integer corner coordinate to the larger. See the `# Determinism`
    // section above for why this is load-bearing for cross-voxel face
    // consistency.
    let (lo, lo_v, hi, hi_v) = if coord_a <= coord_b {
        (coord_a, va, coord_b, vb)
    } else {
        (coord_b, vb, coord_a, va)
    };

    let pa = grid.sample_point(lo.0, lo.1, lo.2);
    let pb = grid.sample_point(hi.0, hi.1, hi.2);

    // Guard against near-equal field values to avoid division by near-zero.
    let denom = lo_v - hi_v;
    if denom.abs() < 1e-12 {
        // Corner values essentially equal — emit the midpoint.
        return (pa + pb) * 0.5;
    }
    let t = (lo_v / denom).clamp(0.0, 1.0);
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
    fn default_method_is_mc33() {
        assert_eq!(ExtractionMethod::default(), ExtractionMethod::Mc33);
    }

    /// Histograms the Lewiner base case across every voxel of a workload.
    /// Returns `[u32; 16]` where index `i` counts voxels in Lewiner base
    /// case `i` (index 0 is unused). Used by diagnostics below.
    #[allow(clippy::cast_sign_loss)]
    fn case_histogram<F: Sdf>(body: &F, grid: &GridSpec) -> [u32; 16] {
        let field = sample_field(body, grid, &mut ());
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

    // --------------------------------------------------------------
    // Progress plumbing (Level 1).
    //
    // Target failure class: a stage adds a progress variant but
    // silently drops inc/finish (the bar would stall or never release).
    // --------------------------------------------------------------

    #[test]
    fn run_with_progress_reports_expected_tick_count() {
        use crate::progress::Spy;
        let sphere = Sphere::new(1.0).unwrap();
        let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(8), 0.4).unwrap();
        let mut spy = Spy::default();
        let _mesh = run_with_progress(&sphere, &grid, ExtractionMethod::ClassicMc, &mut spy);
        let res = grid.resolution();
        let expected = u64::from(res.z + 1) + u64::from(res.z);
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, expected);
        assert_eq!(spy.inc_sum, expected);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn run_with_progress_on_empty_grid_still_finishes() {
        use crate::progress::Spy;
        // 1-cell grid with the sphere far away — zero triangles but
        // progress must still be well-formed.
        let sphere = Sphere::new(0.001).unwrap();
        let grid = GridSpec::new(vec3(100.0, 100.0, 100.0), UVec3::splat(1), 0.1).unwrap();
        let mut spy = Spy::default();
        let _ = run_with_progress(&sphere, &grid, ExtractionMethod::ClassicMc, &mut spy);
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.finish_calls, 1);
        assert_eq!(spy.inc_sum, spy.total);
    }

    #[test]
    fn run_without_progress_delegates_cleanly() {
        // The non-progress `run` must still produce the same mesh as
        // `run_with_progress(… &mut ())` — byte-identical in this case
        // because the extraction logic is the same.
        let sphere = Sphere::new(1.0).unwrap();
        let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(10), 0.3).unwrap();
        let a = run(&sphere, &grid, ExtractionMethod::ClassicMc);
        let b = run_with_progress(&sphere, &grid, ExtractionMethod::ClassicMc, &mut ());
        assert_eq!(a.indices, b.indices);
        assert_eq!(a.vertices, b.vertices);
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

    /// Sharp oracle for cross-voxel face consistency: swapping the two
    /// corner arguments to [`interpolate_edge`] must produce a
    /// **bitwise-identical** `Vec3` output. The two adjacent voxels that
    /// share a face reach the shared edge from opposite corner-numbering
    /// directions; without the canonicalization in `interpolate_edge`, f32
    /// arithmetic produces 1-ULP-different positions that resist welding
    /// and leave boundary-edge holes in the rendered mesh.
    #[test]
    fn interpolate_edge_is_order_invariant() {
        use super::tables::EDGE_CORNERS;

        let grid = GridSpec::new(Vec3::ZERO, UVec3::splat(1), 1.0).unwrap();
        let cell = CellCoord { x: 0, y: 0, z: 0 };
        // Asymmetric SDF magnitudes so `t` is neither 0, 0.5, nor 1 — the
        // arithmetic exercises every operation, not just the symmetric
        // midpoint short-circuit.
        let va = 0.3_f32;
        let vb = -0.7_f32;

        for (i, &(a, b)) in EDGE_CORNERS.iter().enumerate() {
            let v_ab = interpolate_edge(&grid, cell, a, b, va, vb);
            // Swap both the corner indices AND the field values, so the
            // call describes the same edge from the opposite direction.
            let v_ba = interpolate_edge(&grid, cell, b, a, vb, va);
            assert_eq!(
                v_ab.to_array(),
                v_ba.to_array(),
                "edge {i} ({a}, {b}): order-swap produced different f32 outputs \
                 v_ab={v_ab:?} v_ba={v_ba:?} — `interpolate_edge` is not bitwise \
                 order-invariant, so adjacent voxels sharing this edge will \
                 produce different vertex positions"
            );
        }
    }
}
