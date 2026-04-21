//! End-to-end integration test for the progress-reporting seam.
//!
//! Runs a small but realistic pipeline — mesh → weld → Taubin → STL
//! export — with a counting `Progress` spy on every stage, and asserts
//! that each stage declared a length, ticked up to it, and finished.
//!
//! This complements the Level-1 per-stage tests in each module: those
//! verify the plumbing of one stage in isolation with synthetic input;
//! this test verifies the whole pipeline using a real lattice job, so a
//! regression that breaks the `Progress` contract at the boundary
//! between library and binary still fails a test that lives inside
//! `lattice-gen` (where the trait is defined).

#![allow(clippy::unwrap_used, missing_docs)]

use glam::Vec3;
use lattice_gen::{
    GridSpec, LatticeJob, PrimitiveShape, Progress, StrutSpec, UnitCell,
    mesh::{
        ExtractionMethod, TaubinParams, export, mesh_with_progress, taubin_with_progress,
        weld_by_position,
    },
};

/// Counting spy mirroring the one in `src/progress.rs` tests.
/// Duplicated here because the crate's test module is private and the
/// integration-test target cannot see it.
#[derive(Debug, Default)]
struct Spy {
    set_len_calls: u32,
    total: u64,
    inc_sum: u64,
    finish_calls: u32,
}

impl Progress for Spy {
    fn set_len(&mut self, total: u64) {
        self.set_len_calls += 1;
        self.total = total;
    }
    fn inc(&mut self, delta: u64) {
        self.inc_sum += delta;
    }
    fn finish(&mut self) {
        self.finish_calls += 1;
    }
}

fn assert_well_formed(spy: &Spy, stage: &str) {
    assert_eq!(spy.set_len_calls, 1, "{stage}: set_len not called exactly once");
    assert_eq!(
        spy.finish_calls, 1,
        "{stage}: finish not called exactly once"
    );
    assert_eq!(
        spy.inc_sum, spy.total,
        "{stage}: tick sum ({}) != declared total ({})",
        spy.inc_sum, spy.total
    );
}

#[test]
fn full_pipeline_progress_is_well_formed_on_small_cubic_lattice() {
    // Small workload — fast enough for the test suite, real enough to
    // exercise each stage meaningfully.
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(3.0)).unwrap(),
        UnitCell::cubic(2.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::for_job(&job, 0.15).unwrap();

    let mut mesh_spy = Spy::default();
    let mut mesh = mesh_with_progress(&job, &grid, ExtractionMethod::ClassicMc, &mut mesh_spy);
    assert_well_formed(&mesh_spy, "mesh");
    assert!(mesh_spy.total > 0, "mesh stage declared zero work");

    // Welding has no progress signal — the CLI shows a spinner for it.
    weld_by_position(&mut mesh, 0.15 * 1e-4);

    let mut smooth_spy = Spy::default();
    let params = TaubinParams::default_with_iterations(3).unwrap();
    taubin_with_progress(&mut mesh, params, &mut smooth_spy);
    assert_well_formed(&smooth_spy, "taubin");
    assert_eq!(smooth_spy.total, 3);

    let mut export_spy = Spy::default();
    let mut out = Vec::<u8>::new();
    export::stl::write_with_progress(&mesh, &mut out, &mut export_spy).unwrap();
    assert_well_formed(&export_spy, "stl_export");
    assert!(!out.is_empty());
}

#[test]
fn empty_mesh_through_export_still_produces_well_formed_progress() {
    // Far-away primitive + tiny grid → no triangles. Every stage must
    // still finish cleanly.
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(0.5)).unwrap(),
        UnitCell::cubic(2.0).unwrap(),
        StrutSpec::uniform(0.15).unwrap(),
    )
    .unwrap();
    let grid = GridSpec::new(Vec3::splat(100.0), glam::UVec3::splat(2), 0.1).unwrap();

    let mut mesh_spy = Spy::default();
    let mut mesh = mesh_with_progress(&job, &grid, ExtractionMethod::ClassicMc, &mut mesh_spy);
    assert_well_formed(&mesh_spy, "mesh(empty)");

    let mut export_spy = Spy::default();
    let mut out = Vec::<u8>::new();
    export::stl::write_with_progress(&mesh, &mut out, &mut export_spy).unwrap();
    assert_well_formed(&export_spy, "stl_export(empty)");

    // Sanity: welding an empty mesh doesn't panic.
    weld_by_position(&mut mesh, 1e-4);
}
