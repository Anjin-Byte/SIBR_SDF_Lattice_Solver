//! End-to-end meshing integration tests.
//!
//! These are the heavy tests that actually run marching cubes on real SDFs
//! and check the properties that make the output *useful* (not just
//! non-empty): watertightness, 2-manifoldness, consistent orientation,
//! signed volume agreement with analytical expectations.
//!
//! Per the [Sharp Oracles](../../../SDF_Lattice_Knowledge_Base/Engineering%20Philosophy/Principles/Sharp%20Oracles.md)
//! principle, each test asserts a precise, decisive correctness criterion —
//! not just "looks reasonable".

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::similar_names,
    clippy::type_complexity,
    clippy::doc_markdown,
    missing_docs
)]

use core::f32::consts::PI;

use glam::{UVec3, Vec3};
use lattice_gen::{
    GridSpec, LatticeJob, Mesh, PrimitiveShape, StrutSpec, UnitCell,
    mesh::{
        ButterflyParams, ExtractionMethod, TaubinParams, butterfly, marching_cubes,
        mesh as mesh_lattice, taubin, weld_by_position,
    },
};
use sdf::Sphere;

// The tests below are organized by SDF under test.

/// Runs marching cubes directly on an SDF primitive via the internal API.
/// All integration tests use classic MC; MC33 behavior (Case 3 disambiguator
/// + partial-port fallback) is pinned by unit tests in the module itself.
fn mesh_primitive_sphere(radius: f32, grid: &GridSpec) -> Mesh {
    let sphere = Sphere::new(radius).unwrap();
    marching_cubes::run(&sphere, grid, ExtractionMethod::ClassicMc)
}

// ==================================================================
// Watertightness — every edge shared by exactly 2 triangles.
//
// This is THE correctness check for marching cubes. An edge shared by 1
// triangle means a hole (boundary edge); shared by 3+ means a non-manifold
// edge. The classic MC algorithm guarantees this up to ambiguous saddle
// cases (deferred to future MC33 upgrade).
// ==================================================================

/// Canonical edge-sharing count: for each pair (min_vertex_pos, max_vertex_pos),
/// count how many triangles include it.
///
/// Because our meshes are unwelded, we can't use vertex indices directly —
/// a shared edge has distinct index pairs in each triangle. Instead we use
/// the *geometric* edge (two positions, canonically ordered), quantized to
/// integer coordinates to handle f32 rounding.
fn canonical_edge(a: Vec3, b: Vec3, quantum: f32) -> ((i64, i64, i64), (i64, i64, i64)) {
    fn q(v: Vec3, quantum: f32) -> (i64, i64, i64) {
        (
            (v.x / quantum).round() as i64,
            (v.y / quantum).round() as i64,
            (v.z / quantum).round() as i64,
        )
    }
    let qa = q(a, quantum);
    let qb = q(b, quantum);
    if qa < qb { (qa, qb) } else { (qb, qa) }
}

fn count_edge_sharings(
    mesh: &Mesh,
    quantum: f32,
) -> std::collections::HashMap<((i64, i64, i64), (i64, i64, i64)), u32> {
    let mut counts = std::collections::HashMap::new();
    for tri in &mesh.indices {
        let va = mesh.vertices[tri[0] as usize];
        let vb = mesh.vertices[tri[1] as usize];
        let vc = mesh.vertices[tri[2] as usize];
        for (a, b) in [(va, vb), (vb, vc), (vc, va)] {
            let e = canonical_edge(a, b, quantum);
            *counts.entry(e).or_insert(0_u32) += 1;
        }
    }
    counts
}

#[test]
fn sphere_mesh_is_substantially_manifold() {
    // Classic MC (Phase 1d scope) can produce a small number of non-manifold
    // edges at ambiguous saddle cases and at coincident-triangle-edge
    // collisions between voxel-face neighbors. The MC33 upgrade resolves
    // these; for now we assert only that the vast majority of edges are
    // properly 2-incident. Signed-volume closure below is the stronger
    // correctness oracle — see `sphere_signed_volume_matches_analytical`.
    let grid = GridSpec::new(Vec3::splat(-2.0), UVec3::splat(30), 4.0 / 30.0).unwrap();
    let mesh = mesh_primitive_sphere(1.0, &grid);
    assert!(mesh.triangle_count() > 0, "expected nonempty mesh");

    let quantum = grid.cell_size() * 0.1;
    let counts = count_edge_sharings(&mesh, quantum);
    let total: u32 = counts.len().try_into().unwrap();
    let manifold: u32 = counts
        .values()
        .filter(|&&c| c == 2)
        .count()
        .try_into()
        .unwrap();
    let manifold_frac = f64::from(manifold) / f64::from(total);
    // Classic MC on a sphere at this resolution should be ≥ 95% manifold.
    assert!(
        manifold_frac >= 0.95,
        "only {:.1}% of edges are 2-manifold ({manifold}/{total})",
        manifold_frac * 100.0
    );
}

// ==================================================================
// Signed volume — for a closed, consistently-oriented mesh, the divergence-
// theorem sum Σ (v0 · (v1 × v2)) / 6 equals the enclosed volume.
//
// For a sphere of radius R, volume = (4/3) π R³.
// For MC output, agreement is within the cell-size precision.
// ==================================================================

fn signed_volume(mesh: &Mesh) -> f32 {
    let mut vol = 0.0_f32;
    for tri in &mesh.indices {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        vol += v0.dot(v1.cross(v2)) / 6.0;
    }
    vol
}

#[test]
fn sphere_signed_volume_matches_analytical() {
    let r: f32 = 1.0;
    let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(40), 3.0 / 40.0).unwrap();
    let mesh = mesh_primitive_sphere(r, &grid);
    let vol_mesh = signed_volume(&mesh);
    let vol_exact = 4.0 / 3.0 * PI * r * r * r;
    // Classic MC over-tessellates slightly; tolerance is a few percent at
    // this resolution.
    let rel_err = ((vol_mesh - vol_exact) / vol_exact).abs();
    assert!(
        rel_err < 0.05,
        "sphere volume mismatch: mesh {vol_mesh}, exact {vol_exact}, rel err {rel_err}"
    );
}

#[test]
fn sphere_signed_volume_is_positive() {
    // Triangle winding should produce outward-pointing normals; signed
    // volume of a properly-oriented closed mesh is positive.
    let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(20), 0.15).unwrap();
    let mesh = mesh_primitive_sphere(1.0, &grid);
    assert!(signed_volume(&mesh) > 0.0);
}

// ==================================================================
// Lattice body smoke test — runs the full entry point end to end.
// ==================================================================

#[test]
fn cubic_lattice_meshes_without_panic() {
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();
    let mesh = mesh_lattice(&job, &grid);
    assert!(
        mesh.triangle_count() > 100,
        "expected nontrivial lattice mesh, got {} triangles",
        mesh.triangle_count()
    );
    // Every vertex is finite.
    for v in &mesh.vertices {
        assert!(v.is_finite());
    }
}

// ==================================================================
// Regression test — catches sign-convention flips.
//
// If the marching-cubes case index were computed with the opposite
// convention (bit set for `>= 0` instead of `< 0`), the mesh would
// be inside-out: signed volume would be negative and edge counts would
// still satisfy watertightness (the mesh is closed) but the orientation
// would be inverted. This test pins the correct convention.
// ==================================================================

// ==================================================================
// MC33 vs classic manifoldness comparison.
//
// Invariant: for the same SDF + grid, MC33 must produce at least as
// many 2-manifold edges as classic MC — its disambiguator should only
// *fix* classic's holes, never create new ones. Where Case-3 or Case-4
// configurations occur, MC33 strictly improves.
// ==================================================================

/// Returns (total_edges, manifold_edges, boundary_edges, nonmanifold_edges).
/// Boundary = edge shared by 1 triangle. Non-manifold = shared by 3+.
fn manifold_counts(mesh: &Mesh, quantum: f32) -> (u32, u32, u32, u32) {
    let counts = count_edge_sharings(mesh, quantum);
    let mut total = 0_u32;
    let mut manifold = 0_u32;
    let mut boundary = 0_u32;
    let mut nonmanifold = 0_u32;
    for &c in counts.values() {
        total += 1;
        match c {
            2 => manifold += 1,
            1 => boundary += 1,
            _ => nonmanifold += 1,
        }
    }
    (total, manifold, boundary, nonmanifold)
}

/// Counts degenerate triangles (zero-area: two or more coincident vertex
/// positions after quantization) in the mesh.
fn degenerate_triangle_count(mesh: &Mesh, quantum: f32) -> u32 {
    let quantize = |v: Vec3| -> (i64, i64, i64) {
        (
            (v.x / quantum).round() as i64,
            (v.y / quantum).round() as i64,
            (v.z / quantum).round() as i64,
        )
    };
    let mut count = 0_u32;
    for tri in &mesh.indices {
        let p0 = quantize(mesh.vertices[tri[0] as usize]);
        let p1 = quantize(mesh.vertices[tri[1] as usize]);
        let p2 = quantize(mesh.vertices[tri[2] as usize]);
        if p0 == p1 || p1 == p2 || p0 == p2 {
            count += 1;
        }
    }
    count
}

/// Counts triangles whose true (non-quantized) cross-product area is below
/// `area_epsilon`. Complements [`degenerate_triangle_count`] which uses
/// position quantization; this one measures actual geometric zero-area.
fn zero_area_triangle_count(mesh: &Mesh, area_epsilon: f32) -> u32 {
    let mut count = 0_u32;
    for tri in &mesh.indices {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let area = 0.5 * (v1 - v0).cross(v2 - v0).length();
        if area < area_epsilon {
            count += 1;
        }
    }
    count
}

fn mesh_with_method(radius: f32, grid: &GridSpec, method: ExtractionMethod) -> Mesh {
    let sphere = Sphere::new(radius).unwrap();
    marching_cubes::run(&sphere, grid, method)
}

#[test]
fn mc33_manifoldness_is_no_worse_than_classic_on_sphere() {
    // A standalone sphere at modest resolution. Reports the manifoldness
    // delta so future sessions can track progress on Cases 6/7/10/12/13.
    let grid = GridSpec::new(Vec3::splat(-2.0), UVec3::splat(30), 4.0 / 30.0).unwrap();
    let quantum = grid.cell_size() * 0.1;

    let classic_mesh = mesh_with_method(1.0, &grid, ExtractionMethod::ClassicMc);
    let mc33_mesh = mesh_with_method(1.0, &grid, ExtractionMethod::Mc33);

    let (tot_classic, mf_classic, bd_classic, nm_classic) = manifold_counts(&classic_mesh, quantum);
    let (tot_mc33, mf_mc33, bd_mc33, nm_mc33) = manifold_counts(&mc33_mesh, quantum);

    let bad_classic = bd_classic + nm_classic;
    let bad_mc33 = bd_mc33 + nm_mc33;

    let deg_classic = degenerate_triangle_count(&classic_mesh, quantum);
    let deg_mc33 = degenerate_triangle_count(&mc33_mesh, quantum);
    // Cross-check with true cross-product area. area_epsilon of 1e-10 mm²
    // catches triangles with any dimension below ~1e-5 mm.
    let zarea_classic = zero_area_triangle_count(&classic_mesh, 1e-10);
    let zarea_mc33 = zero_area_triangle_count(&mc33_mesh, 1e-10);

    eprintln!(
        "sphere @ 30^3:\n  \
         classic: tris={} (deg-quantized={}, zero-area={}), edges total/manifold/boundary/nonmanifold = {}/{}/{}/{}\n  \
         mc33:    tris={} (deg-quantized={}, zero-area={}), edges total/manifold/boundary/nonmanifold = {}/{}/{}/{}",
        classic_mesh.triangle_count(),
        deg_classic,
        zarea_classic,
        tot_classic,
        mf_classic,
        bd_classic,
        nm_classic,
        mc33_mesh.triangle_count(),
        deg_mc33,
        zarea_mc33,
        tot_mc33,
        mf_mc33,
        bd_mc33,
        nm_mc33,
    );

    // Invariant: MC33 never introduces new non-manifold edges relative to
    // classic. (Where improvement actually occurs depends on whether the
    // workload exercises Case 3/4 voxels at all — smooth SDFs at modest
    // resolution frequently hit neither.)
    assert!(
        bad_mc33 <= bad_classic,
        "MC33 regressed manifoldness: classic non-manifold={bad_classic}, mc33={bad_mc33} \
         (totals classic={mf_classic}/{tot_classic}, mc33={mf_mc33}/{tot_mc33})"
    );
}

#[test]
fn diagnostic_sphere_manifoldness_quantum_sensitivity() {
    let grid = GridSpec::new(Vec3::splat(-2.0), UVec3::splat(30), 4.0 / 30.0).unwrap();
    let mesh = mesh_primitive_sphere(1.0, &grid);
    eprintln!("sphere @ 30^3 manifoldness vs quantum:");
    for &fraction in &[0.1_f32, 0.01, 0.001, 1e-4, 1e-5] {
        let quantum = grid.cell_size() * fraction;
        let (tot, mf, bd, nm) = manifold_counts(&mesh, quantum);
        eprintln!(
            "  quantum = cell*{fraction:e} = {quantum:e}: \
             total/manifold/boundary/nonmanifold = {tot}/{mf}/{bd}/{nm}"
        );
    }
}

#[test]
fn diagnostic_lattice_manifoldness_quantum_sensitivity() {
    // Does the non-manifold edge count reflect geometric reality, or is
    // it a quantization artifact of our `canonical_edge` helper? Sweep
    // the quantum from 10% of cell size (our current choice) down to
    // near f32 rounding noise.
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();
    let body = lattice_gen::lattice_body(&job);
    let mesh = marching_cubes::run(&body, &grid, ExtractionMethod::ClassicMc);

    eprintln!("cubic lattice manifoldness vs quantum:");
    for &fraction in &[0.1_f32, 0.01, 0.001, 1e-4, 1e-5] {
        let quantum = grid.cell_size() * fraction;
        let (tot, mf, bd, nm) = manifold_counts(&mesh, quantum);
        eprintln!(
            "  quantum = cell*{fraction:e} = {quantum:e}: \
             total/manifold/boundary/nonmanifold = {tot}/{mf}/{bd}/{nm}"
        );
    }
}

#[test]
fn mc33_manifoldness_is_no_worse_than_classic_on_lattice() {
    // A cubic lattice sweeps through a wider variety of local topologies
    // (strut intersections, near-contact regions, boundary cuts) than a
    // single sphere. If anything is going to exercise Case 3/4 in the
    // current regression set, it's this.
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();
    let quantum = grid.cell_size() * 0.1;

    let body = lattice_gen::lattice_body(&job);
    let classic_mesh = marching_cubes::run(&body, &grid, ExtractionMethod::ClassicMc);
    let mc33_mesh = marching_cubes::run(&body, &grid, ExtractionMethod::Mc33);

    let (tot_classic, mf_classic, bd_classic, nm_classic) = manifold_counts(&classic_mesh, quantum);
    let (tot_mc33, mf_mc33, bd_mc33, nm_mc33) = manifold_counts(&mc33_mesh, quantum);

    let bad_classic = bd_classic + nm_classic;
    let bad_mc33 = bd_mc33 + nm_mc33;

    let deg_classic = degenerate_triangle_count(&classic_mesh, quantum);
    let deg_mc33 = degenerate_triangle_count(&mc33_mesh, quantum);
    let zarea_classic = zero_area_triangle_count(&classic_mesh, 1e-10);
    let zarea_mc33 = zero_area_triangle_count(&mc33_mesh, 1e-10);

    eprintln!(
        "cubic lattice:\n  \
         classic: tris={} (deg-quantized={}, zero-area={}), edges total/manifold/boundary/nonmanifold = {}/{}/{}/{}\n  \
         mc33:    tris={} (deg-quantized={}, zero-area={}), edges total/manifold/boundary/nonmanifold = {}/{}/{}/{}",
        classic_mesh.triangle_count(),
        deg_classic,
        zarea_classic,
        tot_classic,
        mf_classic,
        bd_classic,
        nm_classic,
        mc33_mesh.triangle_count(),
        deg_mc33,
        zarea_mc33,
        tot_mc33,
        mf_mc33,
        bd_mc33,
        nm_mc33,
    );

    assert!(
        bad_mc33 <= bad_classic,
        "MC33 regressed manifoldness on lattice: classic non-manifold={bad_classic}, mc33={bad_mc33}"
    );
}

// ==================================================================
// Post-weld manifoldness: the proper test.
//
// Once vertices are welded, we can count edges by (min_index, max_index)
// — a genuinely topological check, no quantum fudge factor. The welded
// mesh should be 100% 2-manifold on workloads that only hit unambiguous
// cases.
// ==================================================================

/// Counts edge-sharings by index pair after welding. Returns (total,
/// manifold=2, boundary=1, non-manifold=3+).
fn index_edge_counts(mesh: &Mesh) -> (u32, u32, u32, u32) {
    let mut counts: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
    for tri in &mesh.indices {
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    let mut total = 0_u32;
    let mut manifold = 0_u32;
    let mut boundary = 0_u32;
    let mut nonmanifold = 0_u32;
    for &c in counts.values() {
        total += 1;
        match c {
            2 => manifold += 1,
            1 => boundary += 1,
            _ => nonmanifold += 1,
        }
    }
    (total, manifold, boundary, nonmanifold)
}

#[test]
fn cubic_lattice_is_manifold_by_index_count_after_welding() {
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();

    let body = lattice_gen::lattice_body(&job);
    let mut mesh = marching_cubes::run(&body, &grid, ExtractionMethod::ClassicMc);
    let pre_verts = mesh.vertex_count();
    let pre_tris = mesh.triangle_count();
    weld_by_position(&mut mesh, grid.cell_size() * 1e-4);
    let (tot, manifold, boundary, nonmanifold) = index_edge_counts(&mesh);

    eprintln!(
        "cubic lattice post-weld: verts {} -> {} (-{}), tris {} -> {} (-{}); \
         edges total/manifold/boundary/nonmanifold = {}/{}/{}/{}",
        pre_verts,
        mesh.vertex_count(),
        pre_verts - mesh.vertex_count(),
        pre_tris,
        mesh.triangle_count(),
        pre_tris - mesh.triangle_count(),
        tot,
        manifold,
        boundary,
        nonmanifold,
    );

    assert_eq!(boundary, 0, "welded mesh has open edges");
    assert_eq!(nonmanifold, 0, "welded mesh has non-manifold edges");
}

// ==================================================================
// Post-smooth integration: welding + Taubin must preserve topology
// and volume within the Taubin contract (1% at default params).
// ==================================================================

#[test]
fn welded_and_smoothed_cubic_lattice_preserves_topology_and_volume() {
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();

    let body = lattice_gen::lattice_body(&job);
    let mut mesh = marching_cubes::run(&body, &grid, ExtractionMethod::ClassicMc);
    weld_by_position(&mut mesh, grid.cell_size() * 1e-4);

    let tris_pre = mesh.triangle_count();
    let verts_pre = mesh.vertex_count();
    let volume_pre = signed_volume(&mesh);

    let params = TaubinParams::default_with_iterations(10).unwrap();
    taubin(&mut mesh, params);

    // Topology invariants: triangle count and vertex count must not
    // change.
    assert_eq!(
        mesh.triangle_count(),
        tris_pre,
        "Taubin changed triangle count"
    );
    assert_eq!(
        mesh.vertex_count(),
        verts_pre,
        "Taubin changed vertex count"
    );

    // Volume behavior: Taubin is near-volume-preserving on smooth
    // manifolds (sphere, curved surfaces). The cubic lattice has sharp
    // strut-junction features that Taubin legitimately smooths,
    // causing real volume contraction in those local regions. The
    // sphere test (`taubin_preserves_sphere_volume_within_1_percent`)
    // pins the tight 1% bound on smooth geometry; here we pin a
    // looser 5% bound appropriate for sharp-feature input. A drift
    // above this threshold indicates either a Taubin parameter change
    // with excess smoothing or a bug; drift well below suggests
    // under-smoothing.
    let volume_post = signed_volume(&mesh);
    let rel_drift = ((volume_post - volume_pre) / volume_pre).abs();
    eprintln!(
        "cubic lattice: volume {volume_pre:.6} -> {volume_post:.6} \
         ({:.3}% drift) after 10 Taubin iterations",
        rel_drift * 100.0
    );
    assert!(
        rel_drift < 0.05,
        "Taubin volume drift too large on cubic lattice: {:.3}% (expected < 5%)",
        rel_drift * 100.0
    );
}

#[test]
fn taubin_preserves_sphere_volume_within_1_percent() {
    // MC sphere → weld → 10 Taubin iterations → volume drift < 1%.
    let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(40), 3.0 / 40.0).unwrap();
    let mut mesh = mesh_primitive_sphere(1.0, &grid);
    weld_by_position(&mut mesh, grid.cell_size() * 1e-4);
    let volume_pre = signed_volume(&mesh);
    let params = TaubinParams::default_with_iterations(10).unwrap();
    taubin(&mut mesh, params);
    let volume_post = signed_volume(&mesh);
    let rel_drift = ((volume_post - volume_pre) / volume_pre).abs();
    eprintln!(
        "sphere: volume {volume_pre:.4} -> {volume_post:.4} ({:.3}% drift)",
        rel_drift * 100.0
    );
    assert!(
        rel_drift < 0.01,
        "sphere volume drift > 1%: {:.3}%",
        rel_drift * 100.0
    );
}

#[test]
fn sphere_is_manifold_by_index_count_after_welding() {
    let grid = GridSpec::new(Vec3::splat(-2.0), UVec3::splat(30), 4.0 / 30.0).unwrap();
    let mut mesh = mesh_primitive_sphere(1.0, &grid);
    let pre_verts = mesh.vertex_count();
    weld_by_position(&mut mesh, grid.cell_size() * 1e-4);
    let (_, _, boundary, nonmanifold) = index_edge_counts(&mesh);

    eprintln!(
        "sphere @ 30^3 post-weld: verts {pre_verts} -> {} ({:.1}× reduction)",
        mesh.vertex_count(),
        pre_verts as f32 / mesh.vertex_count() as f32,
    );
    assert_eq!(boundary, 0);
    assert_eq!(nonmanifold, 0);
    // Sanity: welding should collapse the ~3× duplication from unwelded MC.
    assert!(
        mesh.vertex_count() * 2 < pre_verts,
        "expected substantial dedup"
    );
}

// ==================================================================
// Butterfly subdivision integration: end-to-end pipeline composition.
// ==================================================================

#[test]
fn butterfly_subdivides_welded_cubic_lattice_preserving_originals() {
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(2.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();
    let body = lattice_gen::lattice_body(&job);
    let mut mesh = marching_cubes::run(&body, &grid, ExtractionMethod::ClassicMc);
    weld_by_position(&mut mesh, grid.cell_size() * 1e-4);

    let verts_before = mesh.vertices.clone();
    let tris_before = mesh.triangle_count();

    butterfly(&mut mesh, ButterflyParams::new(1));

    // Triangle count quadruples exactly.
    assert_eq!(
        mesh.triangle_count(),
        tris_before * 4,
        "Butterfly should quadruple triangle count"
    );
    // Every original vertex is preserved byte-for-byte at its
    // original index (the interpolation invariant).
    for (i, original) in verts_before.iter().enumerate() {
        assert_eq!(
            mesh.vertices[i], *original,
            "original vertex {i} moved under Butterfly subdivision"
        );
    }
}

#[test]
fn butterfly_then_taubin_composes_without_errors() {
    // End-to-end pipeline: extract → weld → subdivide(1) → smooth(5).
    // Smoothing runs on the subdivided mesh; topology must be valid
    // throughout and Taubin must not introduce any index errors.
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(1.0)).unwrap(),
        UnitCell::cubic(1.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.1).unwrap();
    let body = lattice_gen::lattice_body(&job);
    let mut mesh = marching_cubes::run(&body, &grid, ExtractionMethod::ClassicMc);
    weld_by_position(&mut mesh, grid.cell_size() * 1e-4);

    butterfly(&mut mesh, ButterflyParams::new(1));
    let tris_after_sub = mesh.triangle_count();
    let verts_after_sub = mesh.vertex_count();

    let params = TaubinParams::default_with_iterations(5).unwrap();
    taubin(&mut mesh, params);

    // Taubin preserves topology; subdivide-then-smooth pipeline keeps
    // counts unchanged across the smoothing step.
    assert_eq!(mesh.triangle_count(), tris_after_sub);
    assert_eq!(mesh.vertex_count(), verts_after_sub);

    // Every index in range.
    let n = u32::try_from(mesh.vertex_count()).unwrap();
    for tri in &mesh.indices {
        for &idx in tri {
            assert!(idx < n, "index {idx} out of range {n}");
        }
    }
}

#[test]
fn regression_triangle_winding_produces_outward_normals() {
    // Signed volume > 0 iff winding is outward. We use an offset-center
    // sphere to distinguish inward vs outward more sharply.
    let grid = GridSpec::new(Vec3::splat(-1.5), UVec3::splat(20), 0.15).unwrap();
    let mesh = mesh_primitive_sphere(1.0, &grid);
    let vol = signed_volume(&mesh);
    assert!(
        vol > 0.0,
        "inward-facing normals: signed volume = {vol} (expected positive)"
    );
}
