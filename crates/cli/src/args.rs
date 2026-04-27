//! Argument parsing via [clap].
//!
//! The CLI accepts flat flags — no subcommands. See the README / `--help`
//! output for usage.

use std::path::PathBuf;

use clap::{ArgGroup, Parser, ValueEnum};
use glam::Vec3;

/// Generate an SDF-based lattice mesh and write it to a file.
#[derive(Debug, Parser)]
#[command(name = "sibr-lattice", version, about)]
#[command(group(
    ArgGroup::new("mesh_resolution")
        .required(true)
        .args(["grid_cell", "grid_ratio", "chord_accuracy"]),
))]
pub struct Args {
    /// Primitive boundary shape.
    #[arg(long, value_enum)]
    pub primitive: PrimitiveKind,

    /// Cube half-extents as `x,y,z` in mm (required if `--primitive cube`).
    #[arg(long, value_parser = parse_vec3, required_if_eq("primitive", "cube"))]
    pub half_extents: Option<Vec3>,

    /// Cylinder start endpoint `x,y,z` in mm (required if `--primitive cylinder`).
    #[arg(long, value_parser = parse_vec3, required_if_eq("primitive", "cylinder"))]
    pub cylinder_start: Option<Vec3>,

    /// Cylinder end endpoint `x,y,z` in mm (required if `--primitive cylinder`).
    #[arg(long, value_parser = parse_vec3, required_if_eq("primitive", "cylinder"))]
    pub cylinder_end: Option<Vec3>,

    /// Cylinder radius in mm (required if `--primitive cylinder`).
    #[arg(long, required_if_eq("primitive", "cylinder"))]
    pub cylinder_radius: Option<f32>,

    /// Unit cell topology. Default: cubic. Kelvin and `BccXy` are newly
    /// supported; closed-form properties (`open_porosity`,
    /// `window_diameter`, `hydraulic_diameter`) are cubic-only pending a
    /// follow-up.
    #[arg(long, value_enum, default_value_t = CellTopologyArg::Cubic)]
    pub cell_topology: CellTopologyArg,

    /// Unit cell edge length in mm.
    #[arg(long)]
    pub cell_length: f32,

    /// Strut radius in mm. Must be strictly less than `cell_length / 2`.
    #[arg(long)]
    pub strut_radius: f32,

    /// Smoothing radius (in mm) for strut-to-strut joints (cosmetic).
    /// Default `0.0` gives sharp (hard-`min`) joints — bit-identical to
    /// pre-flag behavior. Positive values fillet the joints with that
    /// approximate radius, producing C1-smooth joins. Recommended
    /// starting value is `strut_radius * 0.3` (e.g. `0.09` for
    /// `strut_radius = 0.3`). Larger values visibly round the joints and
    /// can mitigate stress concentrations in printed parts.
    ///
    /// Note: this does **not** eliminate Marching-Cubes "noise islands"
    /// at the lattice-primitive boundary — for that, use
    /// `--boundary-smoothness`.
    #[arg(long, default_value_t = 0.0)]
    pub strut_smoothness: f32,

    /// Smoothing radius (in mm) for the lattice-primitive boundary
    /// intersection. Default `0.0` gives a sharp (hard-`max`) boundary —
    /// bit-identical to pre-flag behavior. Positive values fillet the
    /// lattice-primitive interface with a C1-smooth blend, eliminating
    /// Marching-Cubes "noise islands" that arise where lattice struts
    /// cross curved primitive boundaries (e.g., cylinder walls) at
    /// grazing angles. A few µm (`0.001` mm) is enough to suppress the
    /// islands; the geometric impact is invisible at that scale. Larger
    /// values (`0.05–0.1` mm) visibly round the lattice-primitive corner.
    #[arg(long, default_value_t = 0.0)]
    pub boundary_smoothness: f32,

    /// Marching-cubes grid cell size in mm (absolute). Smaller is finer;
    /// recommended starting point is `strut_radius / 3` per the Meshing
    /// Complexity Analysis note (equivalent to `--grid-ratio 3`). Exactly
    /// one of `--grid-cell`, `--grid-ratio`, or `--chord-accuracy` must
    /// be specified.
    #[arg(long, group = "mesh_resolution")]
    pub grid_cell: Option<f32>,

    /// Mesh sampling density: number of grid cells per strut radius.
    /// Convenience alternative to `--grid-cell` — the CLI computes
    /// `grid_cell = strut_radius / grid_ratio`. Higher is finer. Typical
    /// useful range 3–10; below 3 surfaces look faceted, above 10 the
    /// triangle count outgrows printer tolerance. Exactly one of
    /// `--grid-cell`, `--grid-ratio`, or `--chord-accuracy` must be
    /// specified.
    #[arg(long, group = "mesh_resolution")]
    pub grid_ratio: Option<f32>,

    /// Mesh chord-error target: maximum distance between the extracted
    /// surface and the true smooth iso-surface, expressed as a
    /// dimensionless fraction of strut radius. Convenience alternative
    /// to `--grid-cell` — the CLI computes
    /// `grid_cell = strut_radius * sqrt(8 * chord_accuracy)`. Must be
    /// in `(0, 1)`; lower is more accurate (e.g., `0.005` = 0.5% of
    /// strut radius). Exactly one of `--grid-cell`, `--grid-ratio`, or
    /// `--chord-accuracy` must be specified.
    #[arg(long, group = "mesh_resolution")]
    pub chord_accuracy: Option<f32>,

    /// Output file path. Format is inferred from the extension (`.stl`,
    /// `.obj`), or can be overridden with `--format`.
    #[arg(long, short)]
    pub output: PathBuf,

    /// Explicit output format (overrides the file-extension inference).
    #[arg(long, value_enum)]
    pub format: Option<OutputFormat>,

    /// Isosurface-extraction method. Default: MC33 (Chernyaev 1995) —
    /// topologically correct at every non-degenerate configuration.
    /// Pass `--extraction-method classic` for the legacy Lorensen-Cline
    /// 1987 algorithm (faster but produces holes / non-manifold edges
    /// at ambiguous configurations). See `ExtractionMethod::Mc33` docs.
    #[arg(long, value_enum, default_value_t = ExtractionMethodArg::Mc33)]
    pub extraction_method: ExtractionMethodArg,

    /// Number of Taubin smoothing iterations to apply after welding.
    /// Default: 0 (no smoothing). Typical useful range is 5–20; higher
    /// values produce smoother surfaces but accumulate more per-vertex
    /// motion (bounded by printer tolerance at realistic values). Uses
    /// Taubin 1995's recommended λ=0.5, μ=-0.5263 parameters.
    #[arg(long, default_value_t = 0)]
    pub smooth_iterations: u32,

    /// Number of Butterfly subdivision iterations to apply after
    /// welding, before smoothing. Default: 0 (no subdivision). Each
    /// iteration quadruples triangle count and adds edge-count many
    /// new vertices — typical useful range is 1–2. Butterfly is
    /// interpolating, so original vertex positions are preserved.
    #[arg(long, default_value_t = 0)]
    pub subdivide_iterations: u32,
}

/// Which primitive boundary to use.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum PrimitiveKind {
    /// Axis-aligned cube centered at the origin.
    Cube,
    /// Capped cylinder between two endpoints.
    Cylinder,
}

/// Which unit-cell topology to tile.
///
/// Deliberately a CLI-side enum separate from [`lattice_gen::UnitCell`]:
/// clap's `ValueEnum` derive needs a local type, and keeping them
/// distinct lets us rename CLI spellings without touching the library API.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum CellTopologyArg {
    /// Simple cubic — 3 axis-centered struts per cell. The only topology
    /// with closed-form properties currently implemented.
    Cubic,
    /// Kelvin (truncated octahedron) — 36-edge TO skeleton per cell.
    Kelvin,
    /// `BCCxy` (vertex octahedron) — 8 body diagonals + 4 top-face edges
    /// per cell.
    Bccxy,
}

/// Which file format to write.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OutputFormat {
    /// Binary STL (for `PreForm` / 3D printers).
    Stl,
    /// Wavefront OBJ ASCII (for 3D viewers).
    Obj,
}

/// Which isosurface extraction method to use.
///
/// Deliberately a CLI-side enum separate from
/// [`lattice_gen::mesh::ExtractionMethod`]: clap's `ValueEnum` derive needs a
/// local type, and keeping them distinct lets us rename CLI spellings
/// without touching the library API.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ExtractionMethodArg {
    /// Classic Marching Cubes (Lorensen & Cline, 1987).
    Classic,
    /// Marching Cubes 33 (Chernyaev, 1995). Cases 3 and 4 are fully
    /// disambiguated; remaining ambiguous cases fall back to classic.
    /// See `ExtractionMethod::Mc33` library docs.
    Mc33,
}

/// Parses a comma-separated `x,y,z` string into a [`Vec3`].
fn parse_vec3(s: &str) -> Result<Vec3, String> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 3 {
        return Err(format!(
            "expected 3 comma-separated values (x,y,z), got {}: {s:?}",
            parts.len()
        ));
    }
    let mut out = [0.0_f32; 3];
    for (i, p) in parts.iter().enumerate() {
        out[i] = p
            .trim()
            .parse::<f32>()
            .map_err(|e| format!("component {i} ({p:?}): {e}"))?;
    }
    Ok(Vec3::new(out[0], out[1], out[2]))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp, missing_docs)]
mod tests {
    use super::*;

    #[test]
    fn parse_vec3_accepts_three_values() {
        assert_eq!(parse_vec3("1,2,3").unwrap(), Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(
            parse_vec3("-1.5,0.0,2.5").unwrap(),
            Vec3::new(-1.5, 0.0, 2.5)
        );
    }

    #[test]
    fn parse_vec3_tolerates_whitespace() {
        assert_eq!(parse_vec3("1, 2, 3").unwrap(), Vec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn parse_vec3_rejects_wrong_count() {
        assert!(parse_vec3("1,2").is_err());
        assert!(parse_vec3("1,2,3,4").is_err());
        assert!(parse_vec3("").is_err());
    }

    #[test]
    fn parse_vec3_rejects_non_numeric() {
        assert!(parse_vec3("a,b,c").is_err());
        assert!(parse_vec3("1,b,3").is_err());
    }
}
