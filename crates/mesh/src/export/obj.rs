//! Wavefront OBJ ASCII writer.
//!
//! # Format
//!
//! OBJ is a line-based text format. Our writer emits:
//!
//! - A single `#` comment header.
//! - One `v x y z` line per vertex.
//! - One `f i j k` line per triangle (**1-indexed** per OBJ spec).
//!
//! We omit normals, texture coords, groups, and materials — the MC output
//! has no meaningful per-vertex normals or texture parameterization.

use std::io::{self, Write};

use crate::mesh::Mesh;
use crate::progress::Progress;

/// Lines per progress tick (combined vertex-line and face-line chunks).
const LINES_PER_TICK: usize = 1024;

/// Writes `mesh` to `writer` as ASCII Wavefront OBJ.
///
/// # Errors
///
/// Returns [`io::Error`] if any write fails.
///
/// # Index convention
///
/// OBJ faces are **1-indexed**. We add 1 to every index we emit. A regression
/// test pins this convention.
pub fn write<W: Write>(mesh: &Mesh, writer: &mut W) -> io::Result<()> {
    write_with_progress(mesh, writer, &mut ())
}

/// Like [`write`], but reports progress to `progress`. Ticks once per
/// [`LINES_PER_TICK`] lines; declared length combines the vertex pass
/// and face pass.
///
/// # Errors
///
/// Same as [`write`].
pub fn write_with_progress<W: Write>(
    mesh: &Mesh,
    writer: &mut W,
    progress: &mut impl Progress,
) -> io::Result<()> {
    let v_ticks = mesh.vertices.len().div_ceil(LINES_PER_TICK.max(1));
    let f_ticks = mesh.indices.len().div_ceil(LINES_PER_TICK.max(1));
    progress.set_len((v_ticks + f_ticks) as u64);

    writeln!(writer, "# SIBR SDF Lattice Generator")?;
    writeln!(
        writer,
        "# {} vertices, {} triangles",
        mesh.vertices.len(),
        mesh.indices.len()
    )?;

    for chunk in mesh.vertices.chunks(LINES_PER_TICK) {
        for v in chunk {
            writeln!(writer, "v {} {} {}", v.x, v.y, v.z)?;
        }
        progress.inc(1);
    }
    for chunk in mesh.indices.chunks(LINES_PER_TICK) {
        for tri in chunk {
            writeln!(writer, "f {} {} {}", tri[0] + 1, tri[1] + 1, tri[2] + 1)?;
        }
        progress.inc(1);
    }
    progress.finish();
    Ok(())
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
    use glam::vec3;

    fn single_triangle_mesh() -> Mesh {
        Mesh {
            vertices: vec![
                vec3(0.0, 0.0, 0.0),
                vec3(1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
            ],
            indices: vec![[0, 1, 2]],
        }
    }

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong field ordering, missing newlines.
    // --------------------------------------------------------------

    #[test]
    fn obj_contains_header_and_vertex_and_face_lines() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let v_lines: Vec<_> = text.lines().filter(|l| l.starts_with("v ")).collect();
        let f_lines: Vec<_> = text.lines().filter(|l| l.starts_with("f ")).collect();
        assert_eq!(v_lines.len(), 3);
        assert_eq!(f_lines.len(), 1);
        assert!(text.starts_with("# SIBR"));
    }

    #[test]
    fn obj_empty_mesh_writes_header_and_nothing_else() {
        let mesh = Mesh::default();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("# SIBR"));
        assert!(!text.contains("\nv "));
        assert!(!text.contains("\nf "));
    }

    #[test]
    fn obj_vertex_lines_have_three_coords() {
        let mesh = Mesh {
            vertices: vec![vec3(1.5, -2.5, 0.0)],
            indices: vec![],
        };
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("v 1.5 -2.5 0\n") || text.contains("v 1.5 -2.5 0 "));
    }

    // --------------------------------------------------------------
    // c. Progress plumbing (Level 1).
    // --------------------------------------------------------------

    #[test]
    fn obj_write_with_progress_ticks_per_chunk() {
        use crate::progress::Spy;
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        let mut spy = Spy::default();
        write_with_progress(&mesh, &mut buf, &mut spy).unwrap();
        // 3 vertices + 1 triangle → ceil(3/1024)=1 + ceil(1/1024)=1 = 2 ticks.
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, 2);
        assert_eq!(spy.inc_sum, 2);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn obj_write_with_progress_empty_mesh_still_finishes() {
        use crate::progress::Spy;
        let mesh = Mesh::default();
        let mut buf = Vec::new();
        let mut spy = Spy::default();
        write_with_progress(&mesh, &mut buf, &mut spy).unwrap();
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, 0);
        assert_eq!(spy.inc_sum, 0);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn obj_write_and_write_with_progress_produce_identical_bytes() {
        let mesh = single_triangle_mesh();
        let mut a = Vec::new();
        let mut b = Vec::new();
        write(&mesh, &mut a).unwrap();
        write_with_progress(&mesh, &mut b, &mut ()).unwrap();
        assert_eq!(a, b);
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "OBJ face indices written as 0-indexed." The spec requires
    /// 1-indexed. Any 3D viewer will either silently skip or visually scramble
    /// the mesh for 0-indexed input.
    #[test]
    fn regression_face_indices_are_one_indexed() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Exactly one face line, indices 1 2 3 (not 0 1 2).
        let f_line = text.lines().find(|l| l.starts_with("f ")).unwrap();
        assert_eq!(f_line, "f 1 2 3");
    }

    /// Regression: "Face indices referred to the wrong vertex — produced
    /// a scrambled mesh." A mesh with two triangles at distinct vertices
    /// must emit face lines that reference each triangle's actual indices,
    /// not a shared or reused set.
    #[test]
    fn regression_face_indices_match_input() {
        let mesh = Mesh {
            vertices: vec![
                vec3(0.0, 0.0, 0.0),
                vec3(1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
                vec3(0.0, 0.0, 1.0),
            ],
            indices: vec![[0, 1, 2], [0, 2, 3]],
        };
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let f_lines: Vec<_> = text.lines().filter(|l| l.starts_with("f ")).collect();
        assert_eq!(f_lines[0], "f 1 2 3");
        assert_eq!(f_lines[1], "f 1 3 4");
    }
}
