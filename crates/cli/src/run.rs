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
    ButterflyParams, ExtractionMethod, Format, TaubinParams, butterfly, export, mesh_with, taubin,
    weld_by_position,
};
use lattice_gen::{GridSpec, LatticeJob, PrimitiveShape, StrutSpec, UnitCell};

use crate::args::{Args, ExtractionMethodArg, OutputFormat, PrimitiveKind};

/// Executes the CLI action described by `args`.
pub fn run(args: &Args) -> Result<()> {
    let primitive = build_primitive(args).context("while constructing the primitive")?;
    let cell = UnitCell::cubic(args.cell_length)
        .with_context(|| format!("cell_length = {}", args.cell_length))?;
    let strut = StrutSpec::uniform(args.strut_radius)
        .with_context(|| format!("strut_radius = {}", args.strut_radius))?;

    let job = LatticeJob::new(primitive, cell, strut).context("validating lattice job")?;

    tracing::info!(
        "job: open_porosity = {:.4}, window_diameter = {:.4} mm, hydraulic_diameter = {:.4} mm",
        job.open_porosity(),
        job.window_diameter(),
        job.hydraulic_diameter(),
    );

    let grid = GridSpec::for_job(&job, args.grid_cell)
        .with_context(|| format!("grid_cell = {}", args.grid_cell))?;
    tracing::info!(
        "grid: resolution = {:?}, cell = {} mm ({} sample points)",
        grid.resolution(),
        grid.cell_size(),
        grid.sample_count(),
    );

    let method = match args.extraction_method {
        ExtractionMethodArg::Classic => ExtractionMethod::ClassicMc,
        ExtractionMethodArg::Mc33 => ExtractionMethod::Mc33,
    };
    tracing::info!("extraction method: {method:?}");

    let t0 = Instant::now();
    let mut m = mesh_with(&job, &grid, method);
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
    // rationale.
    let tolerance = grid.cell_size() * 1e-4;
    weld_by_position(&mut m, tolerance);
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
        butterfly(&mut m, ButterflyParams::new(args.subdivide_iterations));
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
        taubin(&mut m, params);
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
    match format {
        Format::Stl => export::stl::write(&m, &mut writer).context("writing STL")?,
        Format::Obj => export::obj::write(&m, &mut writer).context("writing OBJ")?,
    }
    tracing::info!("wrote {} ({:?})", output_path.display(), format);
    Ok(())
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
#[allow(clippy::unwrap_used, clippy::panic, missing_docs)]
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
            cell_length: 1.0,
            strut_radius: 0.1,
            grid_cell: 0.1,
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
}
