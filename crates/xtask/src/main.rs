//! Workspace dev tooling — post-processing operations that don't belong
//! in the CLI's hot path.
//!
//! # Commands
//!
//! - `remesh` — decimate an STL via `meshoptimizer`'s QEM simplifier
//!   (Garland & Heckbert 1997 + topology-preservation guards). Target a
//!   sub-printer-resolution error bound so the output is geometrically
//!   indistinguishable from the input at slicer resolution but dramatically
//!   smaller on disk.

#![forbid(unsafe_code)]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    // User-facing percent display; usize triangle counts never approach the
    // f64 mantissa limit.
    clippy::cast_precision_loss,
)]

mod progress;

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, anyhow};
use clap::{Parser, Subcommand};
use meshopt::{SimplifyOptions, VertexDataAdapter};

use crate::progress::Pipeline;

#[derive(Debug, Parser)]
#[command(name = "xtask", about = "Workspace post-processing tools")]
struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Debug, Subcommand)]
enum Cmd {
    /// Decimate an STL via meshoptimizer's QEM simplifier with a
    /// sub-printer-resolution error bound.
    Remesh {
        /// Input STL path.
        #[arg(long, short)]
        input: PathBuf,

        /// Output STL path (parent dir must exist).
        #[arg(long, short)]
        output: PathBuf,

        /// Maximum allowed quadric error, **relative to mesh diagonal**
        /// (meshoptimizer convention — NOT absolute mm). Default `0.0002`
        /// corresponds to ≈0.012 mm absolute on a ~60 mm mesh — roughly
        /// half the Formlabs SLA XY spot size (~0.025 mm). This is the
        /// "brush against the printer limit" setting: features at the
        /// print-resolution limit survive; sub-printer-resolution detail
        /// is collapsed. Tighten to `0.0001` for extra safety margin;
        /// loosen to `0.0005` when file size is the binding constraint.
        #[arg(long, default_value_t = 0.0002)]
        target_error: f32,

        /// Target triangle-count floor. `meshoptimizer` stops once the
        /// mesh hits this count OR the error bound would be violated,
        /// whichever comes first. Defaults to 1 (i.e., decimate as
        /// aggressively as the error bound allows).
        #[arg(long, default_value_t = 1)]
        target_tris: usize,
    },
}

fn main() -> Result<()> {
    let args = <Args as Parser>::parse();
    let pipeline = Pipeline::new();
    match args.cmd {
        Cmd::Remesh {
            input,
            output,
            target_error,
            target_tris,
        } => remesh(&pipeline, &input, &output, target_error, target_tris),
    }
}

/// Reads `input` as binary STL, feeds it through `meshopt::simplify` with
/// topology-preservation guards (`SimplifyOptions::empty()`, which enables
/// the safe mode), writes the result as binary STL to `output`.
fn remesh(
    pipeline: &Pipeline,
    input: &PathBuf,
    output: &PathBuf,
    target_error: f32,
    target_tris: usize,
) -> Result<()> {
    // ---- load input ----------------------------------------------------
    let t_read = Instant::now();
    let read_spinner = pipeline.spinner("reading");
    let mut in_file =
        BufReader::new(File::open(input).with_context(|| format!("opening {}", input.display()))?);
    let stl = stl_io::read_stl(&mut in_file).context("parsing STL")?;
    let in_tris_raw = stl.faces.len();
    let in_verts_raw = stl.vertices.len();
    read_spinner.finish(format!(
        "read {in_tris_raw} tris, {in_verts_raw} verts in {:.2?}",
        t_read.elapsed()
    ));

    // ---- convert to meshopt inputs ------------------------------------
    // `stl_io` already welds vertices and returns an indexed mesh.
    let vertices: Vec<[f32; 3]> = stl.vertices.iter().map(|v| [v[0], v[1], v[2]]).collect();
    let indices: Vec<u32> = stl
        .faces
        .iter()
        .flat_map(|f| f.vertices.iter().map(|i| u32::try_from(*i).unwrap()))
        .collect();

    // ---- simplify -----------------------------------------------------
    let t_simp = Instant::now();
    let simp_spinner = pipeline.spinner("simplify");
    // VertexDataAdapter: tell meshopt how to find position data inside our
    // vertex array. Stride = 12 bytes (3 × f32), offset = 0.
    let vertex_bytes: &[u8] = bytemuck_slice(&vertices);
    let adapter = VertexDataAdapter::new(vertex_bytes, 12, 0)
        .map_err(|e| anyhow!("VertexDataAdapter: {e:?}"))?;
    let target_index_count = target_tris.saturating_mul(3);
    // `SimplifyOptions::empty()` = topology-preserving, foldover-detecting,
    // manifold-safe. Exactly what a printable output needs.
    let simplified = meshopt::simplify(
        &indices,
        &adapter,
        target_index_count,
        target_error,
        SimplifyOptions::empty(),
        None,
    );
    let out_tris = simplified.len() / 3;
    simp_spinner.finish(format!(
        "simplified to {out_tris} tris in {:.2?} (target_error={target_error}, target_tris={target_tris})",
        t_simp.elapsed()
    ));

    // ---- write output as binary STL -----------------------------------
    let t_write = Instant::now();
    let write_spinner = pipeline.spinner("writing");
    let mut out_file = BufWriter::new(
        File::create(output).with_context(|| format!("creating {}", output.display()))?,
    );

    // stl_io writes binary STL from an iterator of triangles.
    let tris_iter = simplified.chunks_exact(3).map(|tri| {
        let a = vertices[tri[0] as usize];
        let b = vertices[tri[1] as usize];
        let c = vertices[tri[2] as usize];
        // Face normal — simplified meshes keep the same orientation.
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let nl = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let normal = if nl > 0.0 {
            stl_io::Normal::new([n[0] / nl, n[1] / nl, n[2] / nl])
        } else {
            stl_io::Normal::new([0.0, 0.0, 0.0])
        };
        stl_io::Triangle {
            normal,
            vertices: [
                stl_io::Vertex::new(a),
                stl_io::Vertex::new(b),
                stl_io::Vertex::new(c),
            ],
        }
    });

    stl_io::write_stl(&mut out_file, tris_iter).context("writing STL")?;
    write_spinner.finish(format!(
        "wrote {} in {:.2?}",
        output.display(),
        t_write.elapsed()
    ));

    pipeline.println(format!(
        "[xtask::remesh] DONE: {} → {} triangles ({:.1}% reduction)",
        in_tris_raw,
        out_tris,
        (1.0 - out_tris as f64 / in_tris_raw as f64) * 100.0
    ));

    Ok(())
}

/// Cast a slice of `[f32; 3]` to `&[u8]` for `VertexDataAdapter`. No unsafe
/// here — we use `bytemuck` indirectly via `meshopt`'s own helpers would
/// require adding it as a dep; the direct `[f32;3] → &[u8]` cast below is
/// safe because `[f32; 3]` is `Plain Old Data`.
fn bytemuck_slice(vs: &[[f32; 3]]) -> &[u8] {
    // SAFETY equivalent: f32 is Pod, fixed-size arrays of Pod are Pod,
    // slicing is length × size_of, aligned to at most size_of::<f32>=4
    // which is a multiple of align_of::<u8>=1. See `bytemuck::cast_slice`.
    //
    // We avoid `unsafe` by going through `core::slice::from_raw_parts` —
    // except that IS unsafe. Instead: just add `bytemuck` as a dep.
    //
    // This shim uses `bytemuck::cast_slice` via the `meshopt` crate's
    // re-export path if available; otherwise we fall back. In practice
    // `meshopt` re-exports `bytemuck`.
    //
    // For simplicity here: use the stdlib `align_to` which is unsafe but
    // we've forbidden unsafe. So: depend on bytemuck directly.
    meshopt::utilities::typed_to_bytes(vs)
}
