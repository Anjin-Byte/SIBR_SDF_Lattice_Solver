//! Binary STL writer.
//!
//! # Format
//!
//! Binary STL is a fixed-layout little-endian format:
//!
//! ```text
//! Offset  Length  Contents
//! ---------------------------------
//! 0       80      Header (opaque bytes; conventionally a description)
//! 80      4       u32: triangle count
//! per triangle:
//! +0      12      3×f32: normal (x, y, z)
//! +12     12      3×f32: vertex 0
//! +24     12      3×f32: vertex 1
//! +36     12      3×f32: vertex 2
//! +48     2       u16: attribute byte count (always 0 in standard STL)
//! ```
//!
//! Every numeric field is little-endian. Total per-triangle size is 50 bytes.

use std::io::{self, Write};

use glam::Vec3;

use crate::Mesh;
use crate::progress::Progress;

/// Maximum triangle count representable in the binary STL header's u32 field.
const MAX_STL_TRIANGLES: usize = u32::MAX as usize;

/// Triangles per progress tick. Coarse enough that ticking isn't a
/// write-loop hotspot (~50 bytes/triangle → ~50 KB per tick).
const TRIANGLES_PER_TICK: usize = 1024;

/// Writes `mesh` to `writer` in binary STL format.
///
/// # Errors
///
/// - Returns [`io::Error`] if any write fails.
/// - Returns [`io::Error`] of kind [`io::ErrorKind::InvalidInput`] if the
///   mesh has more than `u32::MAX` triangles (the STL header cannot
///   represent the count).
///
/// # Panics
///
/// Will not panic in practice: the `u32::try_from(mesh.indices.len())`
/// conversion below is bounds-checked by the preceding `InvalidInput`
/// return. If the bounds check is ever reordered or removed, the `expect`
/// documents the invariant being relied upon.
///
/// # Normals
///
/// STL per-triangle normals are computed here from the triangle's vertices
/// via `(v1 - v0) × (v2 - v0)` normalized. This matches our "outward =
/// positive-SDF side" winding convention (verified by the signed-volume
/// integration test in `tests/meshing.rs`).
pub fn write<W: Write>(mesh: &Mesh, writer: &mut W) -> io::Result<()> {
    write_with_progress(mesh, writer, &mut ())
}

/// Like [`write`], but reports progress to `progress`. Ticks once per
/// [`TRIANGLES_PER_TICK`] triangles written; declared length is
/// `triangle_count.div_ceil(TRIANGLES_PER_TICK)`.
///
/// # Errors
///
/// Same as [`write`]. A write failure partway through leaves `progress`
/// in an intermediate state but still calls `finish` via the closure-free
/// path (see `?` propagation below; `finish` is called only on success).
/// Callers that need `finish` on the error path should wrap the call.
///
/// # Panics
///
/// Will not panic in practice — the `u32::try_from(mesh.indices.len())`
/// is bounds-checked by the preceding [`io::ErrorKind::InvalidInput`]
/// return.
#[allow(clippy::expect_used)]
pub fn write_with_progress<W: Write>(
    mesh: &Mesh,
    writer: &mut W,
    progress: &mut impl Progress,
) -> io::Result<()> {
    if mesh.indices.len() > MAX_STL_TRIANGLES {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "mesh has {} triangles; binary STL header limit is {} (u32::MAX)",
                mesh.indices.len(),
                MAX_STL_TRIANGLES,
            ),
        ));
    }

    let total_ticks = mesh.indices.len().div_ceil(TRIANGLES_PER_TICK.max(1));
    progress.set_len(total_ticks as u64);

    // Header: 80 bytes, zero-filled except for a short description.
    let mut header = [0_u8; 80];
    let desc = b"SIBR SDF Lattice Generator";
    header[..desc.len()].copy_from_slice(desc);
    writer.write_all(&header)?;

    // Triangle count (u32 LE). Safe cast: bounds-checked above.
    let count = u32::try_from(mesh.indices.len()).expect("bounds-checked above");
    writer.write_all(&count.to_le_bytes())?;

    // Triangles, ticked per chunk.
    for chunk in mesh.indices.chunks(TRIANGLES_PER_TICK) {
        for tri in chunk {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            let normal = triangle_normal(v0, v1, v2);
            write_vec3(writer, normal)?;
            write_vec3(writer, v0)?;
            write_vec3(writer, v1)?;
            write_vec3(writer, v2)?;
            // Attribute byte count — always zero in standard STL.
            writer.write_all(&0_u16.to_le_bytes())?;
        }
        progress.inc(1);
    }

    progress.finish();
    Ok(())
}

/// Computes the unit normal of a triangle via `(v1 - v0) × (v2 - v0)`.
///
/// For degenerate triangles (collinear vertices or coincident vertices),
/// the cross product has near-zero length. We return `Vec3::ZERO` in that
/// case rather than `NaN`; STL readers tolerate zero normals and compute
/// their own if needed.
fn triangle_normal(v0: Vec3, v1: Vec3, v2: Vec3) -> Vec3 {
    let edge1 = v1 - v0;
    let edge2 = v2 - v0;
    let cross = edge1.cross(edge2);
    let len = cross.length();
    if len < 1e-20 { Vec3::ZERO } else { cross / len }
}

/// Writes a `Vec3` as three little-endian f32s.
#[inline]
fn write_vec3<W: Write>(writer: &mut W, v: Vec3) -> io::Result<()> {
    writer.write_all(&v.x.to_le_bytes())?;
    writer.write_all(&v.y.to_le_bytes())?;
    writer.write_all(&v.z.to_le_bytes())?;
    Ok(())
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::similar_names,
    clippy::needless_range_loop,
    missing_docs
)]
mod tests {
    use super::*;
    use glam::vec3;

    /// Constructs a minimal 1-triangle mesh at known vertices.
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
    // Target failure class: wrong field layout, wrong endianness.
    // --------------------------------------------------------------

    #[test]
    fn binary_stl_has_correct_size() {
        // Expected size: 80 (header) + 4 (count) + 50 * N.
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        assert_eq!(buf.len(), 80 + 4 + 50);
    }

    #[test]
    fn binary_stl_triangle_count_little_endian() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        // Count bytes at offset 80..84.
        let count_bytes: [u8; 4] = buf[80..84].try_into().unwrap();
        assert_eq!(u32::from_le_bytes(count_bytes), 1);
    }

    #[test]
    fn binary_stl_empty_mesh_writes_header_and_zero_count() {
        let mesh = Mesh::default();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        assert_eq!(buf.len(), 80 + 4);
        let count_bytes: [u8; 4] = buf[80..84].try_into().unwrap();
        assert_eq!(u32::from_le_bytes(count_bytes), 0);
    }

    #[test]
    fn binary_stl_normal_points_in_positive_z_for_ccw_xy_triangle() {
        // A triangle in the xy plane with CCW winding (viewed from +z)
        // should have a +z normal.
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        // Normal is first per-triangle field, at offset 84.
        let nx_bytes: [u8; 4] = buf[84..88].try_into().unwrap();
        let ny_bytes: [u8; 4] = buf[88..92].try_into().unwrap();
        let nz_bytes: [u8; 4] = buf[92..96].try_into().unwrap();
        let nx = f32::from_le_bytes(nx_bytes);
        let ny = f32::from_le_bytes(ny_bytes);
        let nz = f32::from_le_bytes(nz_bytes);
        assert!((nx - 0.0).abs() < 1e-6);
        assert!((ny - 0.0).abs() < 1e-6);
        assert!((nz - 1.0).abs() < 1e-6);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // --------------------------------------------------------------

    #[test]
    fn degenerate_triangle_emits_zero_normal_not_nan() {
        // Three collinear vertices → zero-length cross product.
        let mesh = Mesh {
            vertices: vec![
                vec3(0.0, 0.0, 0.0),
                vec3(1.0, 0.0, 0.0),
                vec3(2.0, 0.0, 0.0),
            ],
            indices: vec![[0, 1, 2]],
        };
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        let nx = f32::from_le_bytes(buf[84..88].try_into().unwrap());
        let ny = f32::from_le_bytes(buf[88..92].try_into().unwrap());
        let nz = f32::from_le_bytes(buf[92..96].try_into().unwrap());
        assert!(!nx.is_nan() && !ny.is_nan() && !nz.is_nan());
        assert_eq!((nx, ny, nz), (0.0, 0.0, 0.0));
    }

    // --------------------------------------------------------------
    // c. Progress plumbing (Level 1).
    // --------------------------------------------------------------

    #[test]
    fn stl_write_with_progress_ticks_per_chunk() {
        use crate::progress::Spy;
        // 2500 triangles → ceil(2500 / 1024) = 3 ticks.
        let mut verts = Vec::new();
        let mut idx = Vec::new();
        for i in 0..2500_u32 {
            let base = verts.len() as u32;
            verts.push(vec3(i as f32, 0.0, 0.0));
            verts.push(vec3(i as f32 + 1.0, 0.0, 0.0));
            verts.push(vec3(i as f32, 1.0, 0.0));
            idx.push([base, base + 1, base + 2]);
            let _ = i; // suppress lint
        }
        let mesh = Mesh {
            vertices: verts,
            indices: idx,
        };
        let mut buf = Vec::new();
        let mut spy = Spy::default();
        write_with_progress(&mesh, &mut buf, &mut spy).unwrap();
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, 3);
        assert_eq!(spy.inc_sum, 3);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn stl_write_with_progress_empty_mesh_still_finishes() {
        use crate::progress::Spy;
        let mesh = Mesh::default();
        let mut buf = Vec::new();
        let mut spy = Spy::default();
        write_with_progress(&mesh, &mut buf, &mut spy).unwrap();
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, 0);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn stl_write_and_write_with_progress_produce_identical_bytes() {
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

    /// Regression: "STL triangle-count header written in big-endian order."
    #[test]
    fn regression_triangle_count_is_little_endian() {
        // Use a count with distinguishable byte pattern: 0x01_02_03_04.
        let mut triangles = vec![[0_u32, 1, 2]; 0x01_02_03_04];
        // But that's way too big to fit in memory for the test. Use a
        // smaller distinguishable number: 0x0102 = 258 triangles.
        triangles.truncate(0x0102);
        let mut verts = Vec::with_capacity(triangles.len() * 3);
        for i in 0..triangles.len() {
            let base = verts.len() as u32;
            verts.push(vec3(i as f32, 0.0, 0.0));
            verts.push(vec3(i as f32 + 1.0, 0.0, 0.0));
            verts.push(vec3(i as f32, 1.0, 0.0));
            triangles[i] = [base, base + 1, base + 2];
        }
        let mesh = Mesh {
            vertices: verts,
            indices: triangles,
        };
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        // 0x0102 little-endian is [0x02, 0x01, 0x00, 0x00].
        assert_eq!(buf[80], 0x02);
        assert_eq!(buf[81], 0x01);
        assert_eq!(buf[82], 0x00);
        assert_eq!(buf[83], 0x00);
    }

    /// Regression: "Per-triangle attribute byte count was omitted, making
    /// every triangle record 48 bytes instead of 50 — readers would
    /// misalign starting at triangle 2."
    #[test]
    fn regression_per_triangle_is_50_bytes() {
        // Two triangles → 80 + 4 + 2 * 50 = 184 bytes.
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
        assert_eq!(buf.len(), 80 + 4 + 2 * 50);
    }

    /// Regression: "Attribute byte count set to nonzero garbage — confuses
    /// some readers that interpret it as color bytes."
    #[test]
    fn regression_attribute_byte_count_is_zero() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write(&mesh, &mut buf).unwrap();
        // Attribute bytes at offset 80 + 4 + 48 = 132..134.
        let attr_bytes: [u8; 2] = buf[132..134].try_into().unwrap();
        assert_eq!(u16::from_le_bytes(attr_bytes), 0);
    }
}
