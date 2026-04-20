//! Argument parsing via [clap].
//!
//! The CLI accepts flat flags — no subcommands. See the README / `--help`
//! output for usage.

use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use glam::Vec3;

/// Generate an SDF-based lattice mesh and write it to a file.
#[derive(Debug, Parser)]
#[command(name = "sibr-solver", version, about)]
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

    /// Unit cell edge length in mm.
    #[arg(long)]
    pub cell_length: f32,

    /// Strut radius in mm. Must be strictly less than `cell_length / 2`.
    #[arg(long)]
    pub strut_radius: f32,

    /// Marching-cubes grid cell size in mm. Smaller is finer; recommended
    /// value is `strut_radius / 3` per the Meshing Complexity Analysis note.
    #[arg(long)]
    pub grid_cell: f32,

    /// Output file path. Format is inferred from the extension (`.stl`,
    /// `.obj`), or can be overridden with `--format`.
    #[arg(long, short)]
    pub output: PathBuf,

    /// Explicit output format (overrides the file-extension inference).
    #[arg(long, value_enum)]
    pub format: Option<OutputFormat>,

    /// Isosurface-extraction method. Default: classic (Lorensen-Cline 1987).
    ///
    /// MC33 (Chernyaev 1995) is being ported incrementally. Cases 3 and 4
    /// (face-ambiguous and body-ambiguous diagonals) are fully
    /// disambiguated; remaining ambiguous cases fall back to classic with
    /// a one-time warning. See `ExtractionMethod::Mc33` docs.
    #[arg(long, value_enum, default_value_t = ExtractionMethodArg::Classic)]
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
/// [`lattice_gen::ExtractionMethod`]: clap's `ValueEnum` derive needs a
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
