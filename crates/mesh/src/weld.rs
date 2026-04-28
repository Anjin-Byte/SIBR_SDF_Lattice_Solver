//! Vertex welding — deduplicates coincident positions in a [`Mesh`].
//!
//! Marching Cubes emits each triangle with 3 fresh vertices (see
//! [`super`]'s module-level note on welding). Shared grid edges between
//! neighboring voxels produce bit-exact identical positions at distinct
//! indices. Welding collapses these duplicates so downstream consumers
//! (STL exporters, `PreForm`, GPU paths) see a properly indexed mesh.
//!
//! The weld is **position-based**: two vertices merge if their positions
//! quantize to the same integer triple under the caller-provided
//! `tolerance`. Choice of tolerance is load-bearing — see
//! [`weld_by_position`]'s doc comment for guidance.

use glam::Vec3;
use std::collections::HashMap;

use super::Mesh;

/// Welds vertices within `tolerance` into a single representative and
/// drops triangles that become degenerate as a result. Operates in place.
///
/// # Algorithm
///
/// Single pass over `mesh.vertices`, building a hashmap keyed by
/// quantized position → canonical index. Each old index maps to either a
/// newly-assigned index (for the first representative in its bucket) or
/// the canonical index of an earlier vertex in the same bucket. Indices
/// are rewritten; triangles whose three remapped indices are not pairwise
/// distinct are dropped (they would be zero-area after the weld). Orphan
/// vertices (not referenced by any surviving triangle) are compacted out.
///
/// # Tolerance
///
/// Must be positive. Choose below the smallest legitimate vertex
/// separation in the mesh (otherwise distinct-but-close vertices
/// erroneously merge) and above f32 rounding noise at the mesh's
/// position magnitudes (otherwise bit-identical vertices land in
/// different buckets — unlikely at realistic scales).
///
/// Empirical guidance for Marching-Cubes output: `grid.cell_size() * 1e-4`.
/// This is the quantum at which the cubic-lattice mesh's topology
/// converges to 0 non-manifold edges by index count.
///
/// # Panics
///
/// `debug_assert!`s that `tolerance > 0.0`. Non-positive tolerance has no
/// meaningful geometric interpretation.
///
/// # Determinism
///
/// Output is byte-identical for identical inputs: vertices are visited
/// in original insertion order, and the first representative per bucket
/// wins the canonical-index slot.
#[allow(clippy::cast_possible_truncation, clippy::expect_used)]
pub fn weld_by_position(mesh: &mut Mesh, tolerance: f32) {
    debug_assert!(
        tolerance > 0.0,
        "weld tolerance must be positive, got {tolerance}"
    );

    if mesh.vertices.is_empty() {
        // Also drop any triangles (there shouldn't be any, but be defensive).
        mesh.indices.clear();
        return;
    }

    // 1. Quantize + dedup: build canonical vertices and old→new remap.
    let inv_tol = 1.0 / tolerance;
    let quantize = |v: Vec3| -> (i64, i64, i64) {
        (
            (v.x * inv_tol).round() as i64,
            (v.y * inv_tol).round() as i64,
            (v.z * inv_tol).round() as i64,
        )
    };

    let mut canonical: HashMap<(i64, i64, i64), u32> = HashMap::with_capacity(mesh.vertices.len());
    let mut new_vertices: Vec<Vec3> = Vec::with_capacity(mesh.vertices.len());
    let mut old_to_new: Vec<u32> = Vec::with_capacity(mesh.vertices.len());

    for v in &mesh.vertices {
        let key = quantize(*v);
        if let Some(&new_idx) = canonical.get(&key) {
            old_to_new.push(new_idx);
        } else {
            let new_idx =
                u32::try_from(new_vertices.len()).expect("vertex count exceeds u32 during weld");
            canonical.insert(key, new_idx);
            new_vertices.push(*v);
            old_to_new.push(new_idx);
        }
    }

    // 2. Rewrite indices; retain only non-degenerate triangles.
    let mut new_indices: Vec<[u32; 3]> = Vec::with_capacity(mesh.indices.len());
    for tri in &mesh.indices {
        let a = old_to_new[tri[0] as usize];
        let b = old_to_new[tri[1] as usize];
        let c = old_to_new[tri[2] as usize];
        if a != b && b != c && a != c {
            new_indices.push([a, b, c]);
        }
    }

    // 3. Compact orphans: keep only vertices referenced by some surviving
    //    triangle. Rebuild a second remap in one pass.
    let mut referenced = vec![false; new_vertices.len()];
    for tri in &new_indices {
        referenced[tri[0] as usize] = true;
        referenced[tri[1] as usize] = true;
        referenced[tri[2] as usize] = true;
    }
    let compaction_needed = referenced.iter().any(|&r| !r);
    if compaction_needed {
        let mut compact_remap = vec![u32::MAX; new_vertices.len()];
        let mut compact_vertices: Vec<Vec3> = Vec::with_capacity(new_vertices.len());
        for (old_idx, (&keep, v)) in referenced.iter().zip(new_vertices.iter()).enumerate() {
            if keep {
                let new_idx = u32::try_from(compact_vertices.len())
                    .expect("vertex count exceeds u32 during compaction");
                compact_remap[old_idx] = new_idx;
                compact_vertices.push(*v);
            }
        }
        for tri in &mut new_indices {
            tri[0] = compact_remap[tri[0] as usize];
            tri[1] = compact_remap[tri[1] as usize];
            tri[2] = compact_remap[tri[2] as usize];
        }
        mesh.vertices = compact_vertices;
    } else {
        mesh.vertices = new_vertices;
    }
    mesh.indices = new_indices;
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    missing_docs
)]
mod tests {
    use super::*;

    // --------------------------------------------------------------
    // Ordinary behavior tests.
    // --------------------------------------------------------------

    #[test]
    fn weld_deduplicates_bit_identical_vertices() {
        // Three triangles, each referencing the same geometric vertex at
        // (0, 0, 0) via three different old indices.
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::ZERO,               // 0
                Vec3::new(1.0, 0.0, 0.0), // 1
                Vec3::new(0.0, 1.0, 0.0), // 2
                Vec3::ZERO,               // 3 — dup of 0
                Vec3::new(2.0, 0.0, 0.0), // 4
                Vec3::new(0.0, 2.0, 0.0), // 5
                Vec3::ZERO,               // 6 — dup of 0
                Vec3::new(3.0, 0.0, 0.0), // 7
                Vec3::new(0.0, 3.0, 0.0), // 8
            ],
            indices: vec![[0, 1, 2], [3, 4, 5], [6, 7, 8]],
        };
        weld_by_position(&mut mesh, 1e-6);
        // Expected: (0,0,0) dedups to one vertex, others remain distinct.
        assert_eq!(mesh.vertices.len(), 7, "expected 1 deduped + 6 unique");
        assert_eq!(mesh.triangle_count(), 3, "no degenerates");
        // All triangles should share their first index (the origin).
        assert_eq!(mesh.indices[0][0], mesh.indices[1][0]);
        assert_eq!(mesh.indices[1][0], mesh.indices[2][0]);
    }

    #[test]
    fn weld_preserves_non_coincident_vertices() {
        let vertices = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
        ];
        let indices = vec![[0, 1, 2], [3, 4, 5]];
        let mut mesh = Mesh {
            vertices: vertices.clone(),
            indices: indices.clone(),
        };
        weld_by_position(&mut mesh, 1e-6);
        assert_eq!(mesh.vertices, vertices);
        assert_eq!(mesh.indices, indices);
    }

    #[test]
    fn weld_drops_degenerate_triangle_after_merge() {
        // Triangle whose vertices 0 and 1 are close enough to merge under
        // the tolerance — becomes degenerate after the weld and is dropped.
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),  // 0
                Vec3::new(1e-9, 0.0, 0.0), // 1 — merges with 0 under tol=1e-6
                Vec3::new(0.0, 1.0, 0.0),  // 2
                Vec3::new(5.0, 0.0, 0.0),  // 3
                Vec3::new(0.0, 5.0, 0.0),  // 4
                Vec3::new(5.0, 5.0, 0.0),  // 5
            ],
            indices: vec![
                [0, 1, 2], // will become degenerate
                [3, 4, 5], // survives
            ],
        };
        weld_by_position(&mut mesh, 1e-6);
        assert_eq!(mesh.triangle_count(), 1, "degenerate triangle dropped");
        // Only the surviving triangle's 3 vertices should remain
        // (0 was merged but orphan-compacted).
        assert_eq!(mesh.vertices.len(), 3);
    }

    #[test]
    fn weld_compacts_orphan_vertices() {
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(9.0, 9.0, 9.0), // orphan — no triangle references it
            ],
            indices: vec![[0, 1, 2]],
        };
        weld_by_position(&mut mesh, 1e-6);
        assert_eq!(mesh.vertices.len(), 3, "orphan dropped");
        assert_eq!(mesh.triangle_count(), 1);
    }

    #[test]
    fn weld_is_identity_on_already_welded_input() {
        let vertices = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let indices = vec![[0, 1, 2]];
        let mut mesh = Mesh {
            vertices: vertices.clone(),
            indices: indices.clone(),
        };
        weld_by_position(&mut mesh, 1e-6);
        assert_eq!(mesh.vertices, vertices);
        assert_eq!(mesh.indices, indices);
    }

    #[test]
    fn weld_is_deterministic() {
        let mesh_src = Mesh {
            vertices: vec![
                Vec3::ZERO,
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::ZERO,
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(0.0, 2.0, 0.0),
            ],
            indices: vec![[0, 1, 2], [3, 4, 5]],
        };
        let mut a = mesh_src.clone();
        let mut b = mesh_src.clone();
        weld_by_position(&mut a, 1e-6);
        weld_by_position(&mut b, 1e-6);
        assert_eq!(a.vertices, b.vertices);
        assert_eq!(a.indices, b.indices);
    }

    // --------------------------------------------------------------
    // Edge / boundary tests.
    // --------------------------------------------------------------

    #[test]
    fn weld_handles_empty_mesh() {
        let mut mesh = Mesh::default();
        weld_by_position(&mut mesh, 1e-6);
        assert!(mesh.vertices.is_empty());
        assert!(mesh.indices.is_empty());
    }

    #[test]
    fn weld_handles_mesh_with_vertices_but_no_triangles() {
        let mut mesh = Mesh {
            vertices: vec![Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)],
            indices: vec![],
        };
        weld_by_position(&mut mesh, 1e-6);
        // All vertices are orphans → dropped during compaction.
        assert!(mesh.vertices.is_empty());
        assert!(mesh.indices.is_empty());
    }

    #[test]
    fn weld_with_tolerance_larger_than_mesh_extent_collapses_everything() {
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            indices: vec![[0, 1, 2]],
        };
        // Tolerance 100 quantizes all three to (0, 0, 0) → single vertex
        // → triangle becomes degenerate → dropped.
        weld_by_position(&mut mesh, 100.0);
        assert!(mesh.vertices.is_empty());
        assert!(mesh.indices.is_empty());
    }

    #[test]
    #[should_panic(expected = "weld tolerance must be positive")]
    fn weld_debug_rejects_non_positive_tolerance() {
        let mut mesh = Mesh::default();
        weld_by_position(&mut mesh, 0.0);
    }

    // --------------------------------------------------------------
    // Adversarial: welding must preserve manifoldness of a valid
    // hand-built mesh.
    // --------------------------------------------------------------

    /// A square split into two triangles sharing an edge. After welding,
    /// by-index counting finds 5 total edges: 4 boundary (the square's
    /// perimeter) and 1 shared (the diagonal, count 2).
    #[test]
    fn weld_allows_index_based_edge_counting_to_find_shared_edge() {
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                // Duplicates at positions 0 and 2 to force welding work.
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
            ],
            // Second triangle uses duplicates of (0, 0, 0) and (1, 1, 0).
            indices: vec![[0, 1, 2], [4, 5, 3]],
        };
        weld_by_position(&mut mesh, 1e-6);
        assert_eq!(mesh.vertices.len(), 4, "4 unique corners of a square");
        assert_eq!(mesh.triangle_count(), 2);

        // Count edges by index pair (min, max). Shared edge shows count=2.
        let mut edge_counts: HashMap<(u32, u32), u32> = HashMap::new();
        for tri in &mesh.indices {
            for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_counts.entry(key).or_insert(0) += 1;
            }
        }
        let shared = edge_counts.values().filter(|&&c| c == 2).count();
        let boundary = edge_counts.values().filter(|&&c| c == 1).count();
        assert_eq!(shared, 1, "two triangles share exactly one edge");
        assert_eq!(boundary, 4, "square perimeter has 4 boundary edges");
    }

    // --------------------------------------------------------------
    // Property tests (proptest).
    // --------------------------------------------------------------

    use proptest::prelude::*;

    /// Generate a small mesh with up to 10 random vertices and up to 5
    /// random triangles referencing those vertices. All positions are
    /// finite. Triangle indices always reference a valid vertex.
    fn arb_mesh() -> impl Strategy<Value = Mesh> {
        let vertex_strat = prop::collection::vec(
            (-10.0_f32..10.0, -10.0_f32..10.0, -10.0_f32..10.0)
                .prop_map(|(x, y, z)| Vec3::new(x, y, z)),
            1..=10,
        );
        vertex_strat.prop_flat_map(|verts| {
            let n = verts.len();
            let tri_strat = prop::collection::vec(
                (0..n, 0..n, 0..n).prop_map(|(a, b, c)| [a as u32, b as u32, c as u32]),
                0..=5,
            );
            tri_strat.prop_map(move |indices| Mesh {
                vertices: verts.clone(),
                indices,
            })
        })
    }

    proptest! {
        #[test]
        fn prop_weld_preserves_index_validity(mesh in arb_mesh()) {
            let mut m = mesh;
            weld_by_position(&mut m, 1e-6);
            let n = m.vertices.len() as u32;
            for tri in &m.indices {
                for &idx in tri {
                    prop_assert!(idx < n, "index {idx} out of range {n}");
                }
                prop_assert!(tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);
            }
        }

        #[test]
        fn prop_weld_never_grows_vertex_count(mesh in arb_mesh()) {
            let original_len = mesh.vertices.len();
            let mut m = mesh;
            weld_by_position(&mut m, 1e-6);
            prop_assert!(m.vertices.len() <= original_len);
        }

        #[test]
        fn prop_weld_never_grows_triangle_count(mesh in arb_mesh()) {
            let original_tris = mesh.indices.len();
            let mut m = mesh;
            weld_by_position(&mut m, 1e-6);
            prop_assert!(m.indices.len() <= original_tris);
        }

        #[test]
        fn prop_weld_is_idempotent(mesh in arb_mesh()) {
            let mut once = mesh.clone();
            weld_by_position(&mut once, 1e-6);
            let mut twice = once.clone();
            weld_by_position(&mut twice, 1e-6);
            prop_assert_eq!(once.vertices, twice.vertices);
            prop_assert_eq!(once.indices, twice.indices);
        }
    }
}
