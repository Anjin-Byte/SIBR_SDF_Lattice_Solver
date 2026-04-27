//! Orchestration — the actual action the CLI performs.
//!
//! This is the only file that wires the `sdf` / `lattice-gen` libraries
//! to file I/O. It validates args, constructs a `LatticeJob`, runs the
//! mesher, opens the output file, and invokes the chosen exporter.

use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

use anyhow::{Context, Result, anyhow};
use lattice_gen::mesh::{
    ButterflyParams, ExtractionMethod, Format, TaubinParams, butterfly_with_progress, export,
    mesh_with_progress, taubin_with_progress, weld_by_position,
};
use lattice_gen::{GridSpec, LatticeError, LatticeJob, PrimitiveShape, StrutSpec, UnitCell};

use crate::args::{Args, CellTopologyArg, ExtractionMethodArg, OutputFormat, PrimitiveKind};
use crate::progress::Pipeline;

/// Executes the CLI action described by `args`, rendering progress via
/// the supplied [`Pipeline`].
pub fn run(args: &Args, pipeline: &Pipeline) -> Result<()> {
    let primitive = build_primitive(args).context("while constructing the primitive")?;
    let cell = build_cell(args.cell_topology, args.cell_length)
        .with_context(|| format!("cell_length = {}", args.cell_length))?;
    let strut = StrutSpec::uniform(args.strut_radius)
        .with_context(|| format!("strut_radius = {}", args.strut_radius))?
        .with_joint_smoothness(args.strut_smoothness)
        .with_context(|| format!("strut_smoothness = {}", args.strut_smoothness))?;

    let job = LatticeJob::new(primitive, cell, strut)
        .context("validating lattice job")?
        .with_boundary_smoothness(args.boundary_smoothness)
        .with_context(|| format!("boundary_smoothness = {}", args.boundary_smoothness))?;

    log_job_properties(args, &job);

    let grid_cell = resolve_grid_cell(args)?;
    let grid =
        GridSpec::for_job(&job, grid_cell).with_context(|| format!("grid_cell = {grid_cell}"))?;
    tracing::info!(
        "grid: resolution = {:?}, cell = {} mm ({} sample points)",
        grid.resolution(),
        grid.cell_size(),
        grid.sample_count(),
    );
    log_derived_quantities(args, grid.cell_size());

    let method = match args.extraction_method {
        ExtractionMethodArg::Classic => ExtractionMethod::ClassicMc,
        ExtractionMethodArg::Mc33 => ExtractionMethod::Mc33,
    };
    tracing::info!("extraction method: {method:?}");

    let t0 = Instant::now();
    let mut mesh_bar = pipeline.stage("meshing");
    let mut m = mesh_with_progress(&job, &grid, method, &mut mesh_bar);
    let elapsed = t0.elapsed();
    let unwelded_tris = m.triangle_count();
    let unwelded_verts = m.vertex_count();
    tracing::info!(
        "meshed in {elapsed:.2?}: {unwelded_tris} triangles, {unwelded_verts} vertices (unwelded)"
    );

    if m.is_empty() {
        return Err(anyhow!(
            "mesh has zero triangles — the grid may not intersect the lattice body"
        ));
    }

    // Weld vertices before export so downstream consumers (PreForm, GPU
    // paths, etc.) see a properly indexed mesh. Tolerance is chosen
    // empirically — see `weld_by_position`'s doc comment for the
    // rationale. No inner progress signal — welding is fast; the CLI
    // shows a spinner so the user sees the stage is active.
    let mut weld_bar = pipeline.spinner("welding");
    let tolerance = grid.cell_size() * 1e-4;
    weld_by_position(&mut m, tolerance);
    lattice_gen::Progress::finish(&mut weld_bar);
    tracing::info!(
        "welded: {} triangles, {} vertices ({} duplicates merged, {} degenerates dropped)",
        m.triangle_count(),
        m.vertex_count(),
        unwelded_verts - m.vertex_count(),
        unwelded_tris - m.triangle_count(),
    );

    // Optional Butterfly subdivision. Interpolating scheme — original
    // vertex positions are preserved. Each iteration quadruples
    // triangle count and adds edge-count new vertices. Applied before
    // smoothing so Taubin can polish any high-frequency artifacts the
    // stencil introduces.
    if args.subdivide_iterations > 0 {
        let predicted_tris = m.triangle_count() * 4_usize.pow(args.subdivide_iterations);
        tracing::info!(
            "subdividing: {} Butterfly iteration(s) — predicted output: {} triangles",
            args.subdivide_iterations,
            predicted_tris,
        );
        let t_sub = Instant::now();
        let mut sub_bar = pipeline.stage("subdivide");
        butterfly_with_progress(
            &mut m,
            ButterflyParams::new(args.subdivide_iterations),
            &mut sub_bar,
        );
        tracing::info!(
            "subdivided: {} triangles, {} vertices in {:.2?}",
            m.triangle_count(),
            m.vertex_count(),
            t_sub.elapsed(),
        );
    }

    // Optional Taubin smoothing. Uses Taubin 1995's recommended
    // parameters; topology is preserved and vertex motion is bounded
    // well below printer tolerance at reasonable iteration counts.
    if args.smooth_iterations > 0 {
        let params = TaubinParams::default_with_iterations(args.smooth_iterations)
            .context("building Taubin parameters")?;
        let t_smooth = Instant::now();
        let mut smooth_bar = pipeline.stage("smoothing");
        taubin_with_progress(&mut m, params, &mut smooth_bar);
        tracing::info!(
            "smoothed: {} Taubin iterations in {:.2?}",
            args.smooth_iterations,
            t_smooth.elapsed(),
        );
    }

    let format = resolve_format(args)?;
    let output_path = &args.output;
    let file = File::create(output_path)
        .with_context(|| format!("creating output file {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);
    let mut export_bar = pipeline.stage("writing");
    match format {
        Format::Stl => export::stl::write_with_progress(&m, &mut writer, &mut export_bar)
            .context("writing STL")?,
        Format::Obj => export::obj::write_with_progress(&m, &mut writer, &mut export_bar)
            .context("writing OBJ")?,
    }
    tracing::info!("wrote {} ({:?})", output_path.display(), format);
    Ok(())
}

/// Resolves the mesh `grid_cell` (mm) from whichever of the three
/// mutually-exclusive CLI flags the caller provided. The clap
/// `ArgGroup` guarantees exactly one is `Some`; this function enforces
/// per-flag range validity and emits user-friendly errors.
///
/// Formulas (see the Plan file's Correctness model for derivation):
/// - `--grid-cell h` → `h`
/// - `--grid-ratio N` → `strut_radius / N`
/// - `--chord-accuracy ε` → `strut_radius * sqrt(8 * ε)`
fn resolve_grid_cell(args: &Args) -> Result<f32> {
    if let Some(h) = args.grid_cell {
        if !h.is_finite() || h <= 0.0 {
            return Err(anyhow!("--grid-cell must be > 0 and finite, got {h}"));
        }
        return Ok(h);
    }
    if let Some(n) = args.grid_ratio {
        if !n.is_finite() || n <= 0.0 {
            return Err(anyhow!("--grid-ratio must be > 0 and finite, got {n}"));
        }
        return Ok(args.strut_radius / n);
    }
    // ArgGroup guarantees at least one flag is set; if we reach here,
    // it must be chord_accuracy.
    let eps = args
        .chord_accuracy
        .ok_or_else(|| anyhow!("clap ArgGroup invariant violated — no resolution flag set"))?;
    if !eps.is_finite() || eps <= 0.0 || eps >= 1.0 {
        return Err(anyhow!(
            "--chord-accuracy must be in (0, 1) and finite, got {eps}"
        ));
    }
    Ok(args.strut_radius * (8.0 * eps).sqrt())
}

/// Logs the three derived mesh-resolution quantities plus Woodward's
/// lattice ratio `r*`, so every run shows the user the full picture
/// regardless of which input dialect they chose.
fn log_derived_quantities(args: &Args, resolved_grid_cell: f32) {
    let n = args.strut_radius / resolved_grid_cell;
    let chord_err = 1.0 / (8.0 * n * n);
    let r_star = args.strut_radius / args.cell_length;
    tracing::info!(
        "mesh: grid_cell = {resolved_grid_cell:.4} mm, r/N = {n:.2}, \
         chord error = {:.3}% of strut radius, r* = {r_star:.3}",
        chord_err * 100.0,
    );
}

/// Dispatches to the right [`UnitCell`] constructor based on the CLI
/// topology choice.
///
/// Returns [`LatticeError`] (not `anyhow::Error`) so callers can attach
/// their own context with `.with_context(...)`.
fn build_cell(topology: CellTopologyArg, length: f32) -> Result<UnitCell, LatticeError> {
    match topology {
        CellTopologyArg::Cubic => UnitCell::cubic(length),
        CellTopologyArg::Kelvin => UnitCell::kelvin(length),
        CellTopologyArg::Bccxy => UnitCell::bccxy(length),
    }
}

/// Logs the closed-form job-level properties (open porosity, window
/// diameter, hydraulic diameter) when the topology supports them;
/// otherwise emits a one-line "not yet implemented" note.
///
/// The non-cubic property formulas are deferred to a follow-up. Calling
/// `job.open_porosity()` et al. on Kelvin / `BccXy` jobs panics via
/// `todo!()` — this gate is the thin shim keeping the CLI usable for
/// those topologies in the meantime.
fn log_job_properties(args: &Args, job: &LatticeJob) {
    match args.cell_topology {
        CellTopologyArg::Cubic => {
            tracing::info!(
                "job: open_porosity = {:.4}, window_diameter = {:.4} mm, hydraulic_diameter = {:.4} mm",
                job.open_porosity(),
                job.window_diameter(),
                job.hydraulic_diameter(),
            );
        }
        CellTopologyArg::Kelvin | CellTopologyArg::Bccxy => {
            tracing::info!(
                "job: closed-form properties (open_porosity, window_diameter, \
                 hydraulic_diameter) are not yet implemented for {:?}; \
                 see crates/lattice-gen/src/properties.rs",
                args.cell_topology,
            );
        }
    }
}

/// Constructs a `PrimitiveShape` from the parsed args, validating required
/// fields per primitive kind.
fn build_primitive(args: &Args) -> Result<PrimitiveShape> {
    match args.primitive {
        PrimitiveKind::Cube => {
            let he = args
                .half_extents
                .ok_or_else(|| anyhow!("--half-extents is required for --primitive cube"))?;
            Ok(PrimitiveShape::cube(he)?)
        }
        PrimitiveKind::Cylinder => {
            let start = args
                .cylinder_start
                .ok_or_else(|| anyhow!("--cylinder-start is required for --primitive cylinder"))?;
            let end = args
                .cylinder_end
                .ok_or_else(|| anyhow!("--cylinder-end is required for --primitive cylinder"))?;
            let radius = args
                .cylinder_radius
                .ok_or_else(|| anyhow!("--cylinder-radius is required for --primitive cylinder"))?;
            Ok(PrimitiveShape::cylinder(start, end, radius)?)
        }
    }
}

/// Determines the output format from `--format` if set, otherwise from the
/// output path's file extension.
fn resolve_format(args: &Args) -> Result<Format> {
    if let Some(explicit) = args.format {
        return Ok(match explicit {
            OutputFormat::Stl => Format::Stl,
            OutputFormat::Obj => Format::Obj,
        });
    }
    let ext = args
        .output
        .extension()
        .and_then(|e| e.to_str())
        .ok_or_else(|| {
            anyhow!(
                "cannot infer output format from {}: no file extension. Use --format.",
                args.output.display()
            )
        })?;
    Format::from_extension(ext).ok_or_else(|| {
        anyhow!(
            "unknown output extension '.{ext}' (supported: stl, obj). Use --format to override."
        )
    })
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
    use std::path::PathBuf;

    fn args_with_output(output: &str, format: Option<OutputFormat>) -> Args {
        Args {
            primitive: PrimitiveKind::Cube,
            half_extents: Some(glam::Vec3::splat(1.0)),
            cylinder_start: None,
            cylinder_end: None,
            cylinder_radius: None,
            cell_topology: CellTopologyArg::Cubic,
            cell_length: 1.0,
            strut_radius: 0.1,
            strut_smoothness: 0.0,
            boundary_smoothness: 0.0,
            grid_cell: Some(0.1),
            grid_ratio: None,
            chord_accuracy: None,
            output: PathBuf::from(output),
            format,
            extraction_method: ExtractionMethodArg::Classic,
            smooth_iterations: 0,
            subdivide_iterations: 0,
        }
    }

    #[test]
    fn resolve_format_infers_stl_from_extension() {
        let args = args_with_output("out.stl", None);
        assert_eq!(resolve_format(&args).unwrap(), Format::Stl);
    }

    #[test]
    fn resolve_format_infers_obj_from_extension() {
        let args = args_with_output("out.obj", None);
        assert_eq!(resolve_format(&args).unwrap(), Format::Obj);
    }

    #[test]
    fn resolve_format_explicit_overrides_extension() {
        let args = args_with_output("out.stl", Some(OutputFormat::Obj));
        assert_eq!(resolve_format(&args).unwrap(), Format::Obj);
    }

    #[test]
    fn resolve_format_rejects_no_extension() {
        let args = args_with_output("out", None);
        assert!(resolve_format(&args).is_err());
    }

    #[test]
    fn resolve_format_rejects_unknown_extension() {
        let args = args_with_output("out.xyz", None);
        assert!(resolve_format(&args).is_err());
    }

    // --------------------------------------------------------------
    // resolve_grid_cell: the three-dialect resolver.
    // --------------------------------------------------------------

    fn args_with_grid_cell(h: Option<f32>, n: Option<f32>, eps: Option<f32>) -> Args {
        let mut args = args_with_output("out.stl", None);
        args.strut_radius = 0.3;
        args.grid_cell = h;
        args.grid_ratio = n;
        args.chord_accuracy = eps;
        args
    }

    // Ordinary.

    #[test]
    fn resolve_grid_cell_with_explicit_h_passes_through() {
        let args = args_with_grid_cell(Some(0.05), None, None);
        assert_eq!(resolve_grid_cell(&args).unwrap(), 0.05);
    }

    #[test]
    fn resolve_grid_cell_from_ratio_divides() {
        let args = args_with_grid_cell(None, Some(3.0), None);
        // strut_radius 0.3 / 3 = 0.1
        let resolved = resolve_grid_cell(&args).unwrap();
        assert!((resolved - 0.1).abs() < 1e-6, "got {resolved}");
    }

    #[test]
    fn resolve_grid_cell_from_chord_accuracy_uses_correct_formula() {
        // For N = 3: chord_accuracy = 1/(8*9) = 1/72.
        // Reverse formula: h = r * sqrt(8 * eps) = 0.3 * sqrt(8/72) = 0.3 * sqrt(1/9) = 0.1.
        let eps = 1.0 / 72.0;
        let args = args_with_grid_cell(None, None, Some(eps));
        let resolved = resolve_grid_cell(&args).unwrap();
        assert!(
            (resolved - 0.1).abs() < 1e-6,
            "expected 0.1, got {resolved}"
        );
    }

    #[test]
    fn resolve_grid_cell_ratio_and_explicit_match_numerically() {
        // Same geometry via two dialects should produce the same h.
        let r = 0.3_f32;
        let n = 5.0_f32;
        let args_ratio = args_with_grid_cell(None, Some(n), None);
        let args_explicit = args_with_grid_cell(Some(r / n), None, None);
        assert_eq!(
            resolve_grid_cell(&args_ratio).unwrap(),
            resolve_grid_cell(&args_explicit).unwrap(),
        );
    }

    // Edge.

    #[test]
    fn resolve_grid_cell_rejects_zero_h() {
        let args = args_with_grid_cell(Some(0.0), None, None);
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_negative_h() {
        let args = args_with_grid_cell(Some(-0.1), None, None);
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_nan_h() {
        let args = args_with_grid_cell(Some(f32::NAN), None, None);
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_infinite_h() {
        let args = args_with_grid_cell(Some(f32::INFINITY), None, None);
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_zero_ratio() {
        let args = args_with_grid_cell(None, Some(0.0), None);
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_negative_ratio() {
        let args = args_with_grid_cell(None, Some(-5.0), None);
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_chord_accuracy_at_zero() {
        let args = args_with_grid_cell(None, None, Some(0.0));
        assert!(resolve_grid_cell(&args).is_err());
    }

    #[test]
    fn resolve_grid_cell_rejects_chord_accuracy_at_or_above_one() {
        for eps in [1.0_f32, 1.5, 100.0] {
            let args = args_with_grid_cell(None, None, Some(eps));
            assert!(
                resolve_grid_cell(&args).is_err(),
                "chord_accuracy = {eps} should be rejected"
            );
        }
    }

    #[test]
    fn resolve_grid_cell_rejects_negative_chord_accuracy() {
        let args = args_with_grid_cell(None, None, Some(-0.001));
        assert!(resolve_grid_cell(&args).is_err());
    }

    // Adversarial.

    #[test]
    fn resolve_grid_cell_handles_extreme_chord_accuracy_without_overflow() {
        // Extremely tight.
        let tight = args_with_grid_cell(None, None, Some(1e-12));
        let h = resolve_grid_cell(&tight).unwrap();
        assert!(h.is_finite() && h > 0.0);
        // Extremely loose (just under 1).
        let loose = args_with_grid_cell(None, None, Some(0.999));
        let h = resolve_grid_cell(&loose).unwrap();
        assert!(h.is_finite() && h > 0.0);
    }

    #[test]
    fn resolve_grid_cell_handles_extreme_ratio_without_overflow() {
        let args = args_with_grid_cell(None, Some(1e6), None);
        let h = resolve_grid_cell(&args).unwrap();
        assert!(h.is_finite() && h > 0.0);
    }

    // Round-trip (property-style, multiple-value).

    #[test]
    fn resolve_grid_cell_ratio_round_trips() {
        // For several r and N, confirm h = r/N then N' = r/h recovers N
        // to within f32 precision.
        for r in [0.01_f32, 0.1, 0.3, 1.5, 10.0] {
            for n in [1.0_f32, 2.5, 3.0, 7.5, 10.0, 100.0] {
                let mut args = args_with_grid_cell(None, Some(n), None);
                args.strut_radius = r;
                let h = resolve_grid_cell(&args).unwrap();
                let n_recovered = r / h;
                assert!(
                    (n_recovered - n).abs() / n < 1e-5,
                    "round-trip failed: r={r}, n={n}, h={h}, recovered={n_recovered}"
                );
            }
        }
    }

    #[test]
    fn resolve_grid_cell_chord_accuracy_round_trips() {
        // For several r and ε, confirm h = r·√(8ε) then ε' = h²/(8r²)
        // recovers ε.
        for r in [0.01_f32, 0.1, 0.3, 1.5, 10.0] {
            for eps in [0.001_f32, 0.01, 0.05, 0.1, 0.5, 0.9] {
                let mut args = args_with_grid_cell(None, None, Some(eps));
                args.strut_radius = r;
                let h = resolve_grid_cell(&args).unwrap();
                let eps_recovered = (h * h) / (8.0 * r * r);
                assert!(
                    (eps_recovered - eps).abs() / eps < 1e-5,
                    "round-trip failed: r={r}, eps={eps}, h={h}, recovered={eps_recovered}"
                );
            }
        }
    }

    // --------------------------------------------------------------
    // build_cell — topology dispatch.
    // --------------------------------------------------------------

    // Ordinary.

    #[test]
    fn build_cell_dispatches_to_cubic() {
        let c = build_cell(CellTopologyArg::Cubic, 2.0).unwrap();
        assert!(matches!(c, UnitCell::Cubic { .. }));
    }

    #[test]
    fn build_cell_dispatches_to_kelvin() {
        let c = build_cell(CellTopologyArg::Kelvin, 2.0).unwrap();
        assert!(matches!(c, UnitCell::Kelvin { .. }));
    }

    #[test]
    fn build_cell_dispatches_to_bccxy() {
        let c = build_cell(CellTopologyArg::Bccxy, 2.0).unwrap();
        assert!(matches!(c, UnitCell::BccXy { .. }));
    }

    // Edge.

    #[test]
    fn build_cell_propagates_lattice_error_on_negative_length() {
        assert!(build_cell(CellTopologyArg::Kelvin, -1.0).is_err());
        assert!(build_cell(CellTopologyArg::Bccxy, 0.0).is_err());
        assert!(build_cell(CellTopologyArg::Cubic, f32::NAN).is_err());
    }

    // --------------------------------------------------------------
    // Default topology back-compat and parse tests.
    // --------------------------------------------------------------

    use clap::Parser;

    /// Regression: "--cell-topology default changed; existing invocations
    /// that omit the flag started producing a non-cubic lattice silently."
    /// Detection: parse a minimal arg vector without --cell-topology and
    /// assert the default is Cubic.
    #[test]
    fn regression_default_topology_is_cubic() {
        let args = Args::try_parse_from([
            "sibr-lattice",
            "--primitive",
            "cube",
            "--half-extents",
            "1,1,1",
            "--cell-length",
            "1",
            "--strut-radius",
            "0.1",
            "--grid-cell",
            "0.1",
            "-o",
            "x.stl",
        ])
        .unwrap();
        assert_eq!(args.cell_topology, CellTopologyArg::Cubic);
    }

    #[test]
    fn parse_cell_topology_kelvin_spelling() {
        let args = Args::try_parse_from([
            "sibr-lattice",
            "--primitive",
            "cube",
            "--half-extents",
            "1,1,1",
            "--cell-topology",
            "kelvin",
            "--cell-length",
            "1",
            "--strut-radius",
            "0.1",
            "--grid-cell",
            "0.1",
            "-o",
            "x.stl",
        ])
        .unwrap();
        assert_eq!(args.cell_topology, CellTopologyArg::Kelvin);
    }

    #[test]
    fn parse_cell_topology_bccxy_spelling() {
        let args = Args::try_parse_from([
            "sibr-lattice",
            "--primitive",
            "cube",
            "--half-extents",
            "1,1,1",
            "--cell-topology",
            "bccxy",
            "--cell-length",
            "1",
            "--strut-radius",
            "0.1",
            "--grid-cell",
            "0.1",
            "-o",
            "x.stl",
        ])
        .unwrap();
        assert_eq!(args.cell_topology, CellTopologyArg::Bccxy);
    }

    // --------------------------------------------------------------
    // run() — end-to-end smoke tests for every topology.
    // Target failure class: a future refactor accidentally calls a
    // todo!() property method on a non-cubic job, panicking the CLI.
    // --------------------------------------------------------------

    /// Build a smoke-test `Args` with a small cube primitive, cell length
    /// 2.0, and a comfortably small strut radius that satisfies `r_max` for
    /// every topology (`r_max = L/(4√2) ≈ 0.354` for Kelvin, the tightest).
    fn smoke_args(topology: CellTopologyArg, out: &std::path::Path) -> Args {
        let mut args = args_with_output(out.to_str().unwrap(), None);
        args.half_extents = Some(glam::Vec3::splat(3.0));
        args.cell_topology = topology;
        args.cell_length = 2.0;
        args.strut_radius = 0.15;
        args.grid_cell = Some(0.15);
        args
    }

    fn smoke_output_path(tag: &str) -> std::path::PathBuf {
        // Unique-enough per process; parallel tests in the same crate
        // would collide without a nonce, so include the tag.
        let mut p = std::env::temp_dir();
        p.push(format!("sibr-cli-smoke-{}-{}.stl", tag, std::process::id()));
        p
    }

    #[test]
    fn run_cubic_smoke_completes() {
        let out = smoke_output_path("cubic");
        let args = smoke_args(CellTopologyArg::Cubic, &out);
        run(&args, &Pipeline::new()).expect("cubic run should succeed");
        assert!(out.exists(), "output file not written");
        let _ = std::fs::remove_file(&out);
    }

    /// Regression: Kelvin via the CLI must not panic on the properties
    /// log. Before gating, `job.open_porosity()` was called
    /// unconditionally and hit `todo!()` for non-cubic topologies.
    #[test]
    fn regression_cli_kelvin_does_not_panic_on_properties_log() {
        let out = smoke_output_path("kelvin");
        let args = smoke_args(CellTopologyArg::Kelvin, &out);
        run(&args, &Pipeline::new()).expect("kelvin run should succeed without panic");
        assert!(out.exists(), "output file not written");
        let _ = std::fs::remove_file(&out);
    }

    /// Regression: `BccXy` via the CLI must not panic on the properties log.
    #[test]
    fn regression_cli_bccxy_does_not_panic_on_properties_log() {
        let out = smoke_output_path("bccxy");
        let args = smoke_args(CellTopologyArg::Bccxy, &out);
        run(&args, &Pipeline::new()).expect("bccxy run should succeed without panic");
        assert!(out.exists(), "output file not written");
        let _ = std::fs::remove_file(&out);
    }
}
