//! Butterfly surface subdivision (Dyn, Levin, Gregory 1990).
//!
//! Butterfly is an **interpolating** subdivision scheme for triangle
//! meshes: each iteration splits every triangle into four and inserts a
//! new vertex at the midpoint of every edge, but **original vertex
//! positions are preserved exactly**. The new midpoint vertex is placed
//! using the classic 8-point stencil — a weighted average of the two
//! edge endpoints, the two "wings" (opposite vertices in the two
//! triangles sharing the edge), and four "tails" (vertices completing
//! the triangles adjacent to each wing).
//!
//! Unlike Loop or Catmull-Clark (both *approximating* schemes that
//! contract the mesh toward a smoothed limit surface), Butterfly keeps
//! every vertex the caller started with, so it's appropriate when the
//! input already represents the desired surface and subdivision is
//! purely for adding surface density.
//!
//! # Scope
//!
//! This is a **Tier-1** implementation: classic 8-point Butterfly, no
//! valence-dependent stencil correction. Mild artifacts are visible at
//! extraordinary vertices (valence ≠ 6 for interior). Addressing those
//! via Zorin 1996's "Modified Butterfly" is deferred to a follow-on
//! session. A linear-midpoint fallback handles boundary edges and any
//! case where the 8-point stencil can't be fully evaluated.
//!
//! # Prerequisite: welded mesh
//!
//! Subdivision operates on the vertex-index graph, so coincident
//! vertices must share an index. Welding ([`super::weld`]) provides this
//! guarantee. On an unwelded mesh, every triangle is topologically
//! isolated — edges have no shared faces and every stencil falls back
//! to linear midpoint.
//!
//! # Reference
//!
//! N. Dyn, D. Levin, J. A. Gregory, "A butterfly subdivision scheme
//! for surface interpolation with tension control," ACM TOG 1990.

use glam::Vec3;
use std::collections::HashMap;

use super::Mesh;

/// Tension parameter `w` from Dyn-Levin-Gregory. Standard choice; the
/// value used throughout the literature.
const TENSION: f32 = 1.0 / 16.0;

/// Sentinel representing "no second face" for a boundary edge.
const NO_FACE: u32 = u32::MAX;

/// Parameters for [`butterfly`] subdivision. Construct via
/// [`ButterflyParams::new`]. Zero iterations is a no-op.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ButterflyParams {
    iterations: u32,
}

impl ButterflyParams {
    /// Constructs a parameter set. No validation — any `u32` is
    /// acceptable as an iteration count, including zero (no-op) and
    /// values large enough to exhaust memory. Memory use grows ~4× per
    /// iteration; callers are responsible for sizing.
    pub fn new(iterations: u32) -> Self {
        Self { iterations }
    }

    /// Number of Butterfly iterations to apply.
    pub fn iterations(&self) -> u32 {
        self.iterations
    }
}

/// Applies `params.iterations` iterations of classic 8-point Butterfly
/// subdivision to `mesh`. Original vertex positions are preserved
/// exactly; each iteration quadruples triangle count and adds
/// edge-count new vertices.
///
/// # Performance
///
/// Per iteration: `O(T)` time for edge-to-face map construction +
/// `O(E)` stencil evaluations + `O(T)` triangle splits, where `T` is
/// triangle count and `E` is edge count (`≈ 1.5 T` for a closed
/// manifold). Memory: output mesh is ~4× the input size.
pub fn butterfly(mesh: &mut Mesh, params: ButterflyParams) {
    for _ in 0..params.iterations {
        butterfly_once(mesh);
    }
}

/// Applies a single Butterfly iteration in place.
#[allow(clippy::cast_possible_truncation, clippy::expect_used)]
fn butterfly_once(mesh: &mut Mesh) {
    if mesh.indices.is_empty() {
        return;
    }

    let edge_to_face = build_edge_to_face(mesh);

    // Preserve every original vertex at its existing index; new vertices
    // are appended. This is what makes the scheme interpolating.
    let mut new_vertices: Vec<Vec3> = mesh.vertices.clone();
    let mut edge_to_new_idx: HashMap<(u32, u32), u32> = HashMap::with_capacity(edge_to_face.len());

    // Iterate edges in sorted order so vertex-index assignment (and
    // therefore the output mesh) is deterministic.
    let mut edges: Vec<(u32, u32)> = edge_to_face.keys().copied().collect();
    edges.sort_unstable();

    for edge in edges {
        let new_pos = compute_midpoint(mesh, &edge_to_face, edge.0, edge.1);
        let new_idx =
            u32::try_from(new_vertices.len()).expect("subdivide: vertex count exceeds u32");
        new_vertices.push(new_pos);
        edge_to_new_idx.insert(edge, new_idx);
    }

    // Split every triangle into four, preserving winding.
    let mut new_indices: Vec<[u32; 3]> = Vec::with_capacity(mesh.indices.len() * 4);
    for tri in &mesh.indices {
        let [a, b, c] = *tri;
        let m_ab = edge_to_new_idx[&canonical_edge(a, b)];
        let m_bc = edge_to_new_idx[&canonical_edge(b, c)];
        let m_ca = edge_to_new_idx[&canonical_edge(c, a)];
        new_indices.push([a, m_ab, m_ca]);
        new_indices.push([m_ab, b, m_bc]);
        new_indices.push([m_ca, m_bc, c]);
        new_indices.push([m_ab, m_bc, m_ca]);
    }

    mesh.vertices = new_vertices;
    mesh.indices = new_indices;
}

/// Returns `(min(a, b), max(a, b))` so edges hash identically regardless
/// of traversal direction.
fn canonical_edge(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// Builds a map from each edge to the (up to two) triangle indices that
/// contain it. Boundary edges have `[face_idx, NO_FACE]`.
///
/// The face index stored is the position in `mesh.indices`, castable to
/// `u32` for meshes with fewer than 4 billion triangles (covers all
/// realistic workloads).
#[allow(clippy::cast_possible_truncation)]
fn build_edge_to_face(mesh: &Mesh) -> HashMap<(u32, u32), [u32; 2]> {
    let mut map: HashMap<(u32, u32), [u32; 2]> = HashMap::with_capacity(mesh.indices.len() * 3 / 2);
    for (face_idx, tri) in mesh.indices.iter().enumerate() {
        let face_idx = face_idx as u32;
        for &(i, j) in &[(0, 1), (1, 2), (2, 0)] {
            let edge = canonical_edge(tri[i], tri[j]);
            let entry = map.entry(edge).or_insert([NO_FACE, NO_FACE]);
            if entry[0] == NO_FACE {
                entry[0] = face_idx;
            } else if entry[1] == NO_FACE {
                entry[1] = face_idx;
            }
            // Three+ faces sharing an edge: non-manifold input. We
            // record only the first two; the stencil for this edge
            // will still evaluate (incorrectly) but won't panic.
        }
    }
    map
}

/// Returns the vertex of `tri` that is neither `v1` nor `v2`.
///
/// Returns `None` for a degenerate triangle (two or three identical
/// indices) or one that doesn't contain both `v1` and `v2` — both
/// shouldn't happen for a well-formed input, but yield a safe
/// linear-fallback path rather than a panic.
fn third_vertex(tri: &[u32; 3], v1: u32, v2: u32) -> Option<u32> {
    tri.iter().copied().find(|&v| v != v1 && v != v2)
}

/// Finds the "tail" vertex: the third vertex of the triangle that shares
/// edge `(wing, corner)` with the face `excluded_face_idx`, picking the
/// OTHER face (not the one we're computing the midpoint from).
fn find_tail(
    edge_to_face: &HashMap<(u32, u32), [u32; 2]>,
    mesh_indices: &[[u32; 3]],
    wing: u32,
    corner: u32,
    excluded_face_idx: u32,
) -> Option<u32> {
    let edge = canonical_edge(wing, corner);
    let faces = edge_to_face.get(&edge)?;
    let tail_face_idx = if faces[0] == excluded_face_idx {
        faces[1]
    } else if faces[1] == excluded_face_idx {
        faces[0]
    } else {
        // excluded_face_idx should always be one of the two faces on
        // this edge; if it isn't, the edge_to_face map has an
        // inconsistency — fall back.
        return None;
    };
    if tail_face_idx == NO_FACE {
        return None;
    }
    let tail_face = &mesh_indices[tail_face_idx as usize];
    third_vertex(tail_face, wing, corner)
}

/// Computes the new midpoint vertex position for edge `(v1, v2)` using
/// the 8-point Butterfly stencil, falling back to the linear midpoint
/// if any stencil point is unavailable (e.g., boundary).
#[allow(clippy::many_single_char_names)]
fn compute_midpoint(
    mesh: &Mesh,
    edge_to_face: &HashMap<(u32, u32), [u32; 2]>,
    v1: u32,
    v2: u32,
) -> Vec3 {
    let p1 = mesh.vertices[v1 as usize];
    let p2 = mesh.vertices[v2 as usize];
    let linear = (p1 + p2) * 0.5;

    let Some(faces) = edge_to_face.get(&canonical_edge(v1, v2)) else {
        return linear;
    };
    if faces[0] == NO_FACE || faces[1] == NO_FACE {
        return linear;
    }

    let face0_idx = faces[0];
    let face1_idx = faces[1];
    let Some(w1) = third_vertex(&mesh.indices[face0_idx as usize], v1, v2) else {
        return linear;
    };
    let Some(w2) = third_vertex(&mesh.indices[face1_idx as usize], v1, v2) else {
        return linear;
    };

    // Tails around wing1 (on the two edges of triangle face0 that don't
    // include the edge being subdivided), and likewise for wing2.
    let Some(t1) = find_tail(edge_to_face, &mesh.indices, w1, v1, face0_idx) else {
        return linear;
    };
    let Some(t2) = find_tail(edge_to_face, &mesh.indices, w1, v2, face0_idx) else {
        return linear;
    };
    let Some(t3) = find_tail(edge_to_face, &mesh.indices, w2, v1, face1_idx) else {
        return linear;
    };
    let Some(t4) = find_tail(edge_to_face, &mesh.indices, w2, v2, face1_idx) else {
        return linear;
    };

    let pw1 = mesh.vertices[w1 as usize];
    let pw2 = mesh.vertices[w2 as usize];
    let pt1 = mesh.vertices[t1 as usize];
    let pt2 = mesh.vertices[t2 as usize];
    let pt3 = mesh.vertices[t3 as usize];
    let pt4 = mesh.vertices[t4 as usize];

    // Stencil: (1/2)(v1+v2) + 2w·(w1+w2) − w·(t1+t2+t3+t4). Weights sum
    // to 1 (affine invariance): 1/2 + 1/2 + 4w − 4w = 1.
    0.5 * (p1 + p2) + (2.0 * TENSION) * (pw1 + pw2) - TENSION * (pt1 + pt2 + pt3 + pt4)
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

    /// A regular tetrahedron: 4 vertices, 4 triangular faces, every
    /// edge shared by exactly 2 faces. All vertices have valence 3.
    fn tetrahedron() -> Mesh {
        Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
            indices: vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
        }
    }

    // --------------------------------------------------------------
    // ButterflyParams construction.
    // --------------------------------------------------------------

    #[test]
    fn params_stores_iteration_count() {
        assert_eq!(ButterflyParams::new(0).iterations(), 0);
        assert_eq!(ButterflyParams::new(7).iterations(), 7);
    }

    // --------------------------------------------------------------
    // Ordinary tests: core invariants (plan §Invariants 1-6, 8).
    // --------------------------------------------------------------

    #[test]
    fn butterfly_zero_iterations_is_identity() {
        let before = tetrahedron();
        let mut mesh = before.clone();
        butterfly(&mut mesh, ButterflyParams::new(0));
        assert_eq!(mesh.vertices, before.vertices);
        assert_eq!(mesh.indices, before.indices);
    }

    #[test]
    fn butterfly_handles_empty_mesh() {
        let mut mesh = Mesh::default();
        butterfly(&mut mesh, ButterflyParams::new(5));
        assert!(mesh.vertices.is_empty());
        assert!(mesh.indices.is_empty());
    }

    #[test]
    fn butterfly_quadruples_triangle_count_per_iteration() {
        // Tetrahedron has 4 triangles. After 1 iteration: 16. After 2: 64.
        let mut mesh = tetrahedron();
        butterfly(&mut mesh, ButterflyParams::new(1));
        assert_eq!(mesh.triangle_count(), 16);
        butterfly(&mut mesh, ButterflyParams::new(1));
        assert_eq!(mesh.triangle_count(), 64);
    }

    #[test]
    fn butterfly_vertex_count_grows_by_edge_count() {
        // Tetrahedron: V=4, T=4, E=6. After 1 iter: V = 4 + 6 = 10.
        let mut mesh = tetrahedron();
        butterfly(&mut mesh, ButterflyParams::new(1));
        assert_eq!(mesh.vertex_count(), 10);
    }

    #[test]
    fn butterfly_preserves_original_vertex_positions() {
        // The defining property of interpolating subdivision.
        let before = tetrahedron();
        let mut mesh = before.clone();
        butterfly(&mut mesh, ButterflyParams::new(2));
        for (i, original) in before.vertices.iter().enumerate() {
            assert_eq!(
                mesh.vertices[i], *original,
                "original vertex {i} moved from {original:?} to {:?}",
                mesh.vertices[i]
            );
        }
    }

    #[test]
    fn butterfly_indices_are_in_range_after_subdivision() {
        let mut mesh = tetrahedron();
        butterfly(&mut mesh, ButterflyParams::new(2));
        let n = u32::try_from(mesh.vertices.len()).unwrap();
        for tri in &mesh.indices {
            for &idx in tri {
                assert!(idx < n, "index {idx} out of range {n}");
            }
        }
    }

    #[test]
    fn butterfly_is_deterministic() {
        let src = tetrahedron();
        let mut a = src.clone();
        let mut b = src.clone();
        butterfly(&mut a, ButterflyParams::new(2));
        butterfly(&mut b, ButterflyParams::new(2));
        assert_eq!(a.vertices, b.vertices);
        assert_eq!(a.indices, b.indices);
    }

    // --------------------------------------------------------------
    // Edge tests: boundary handling, degenerate inputs.
    // --------------------------------------------------------------

    #[test]
    fn butterfly_single_triangle_uses_linear_fallback() {
        // A single triangle has 3 boundary edges → all midpoints use
        // linear fallback → each midpoint is exactly the average of
        // its endpoints.
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            indices: vec![[0, 1, 2]],
        };
        butterfly(&mut mesh, ButterflyParams::new(1));
        // Output: 4 triangles, 6 vertices (3 original + 3 midpoints).
        assert_eq!(mesh.triangle_count(), 4);
        assert_eq!(mesh.vertex_count(), 6);
        // Originals unchanged.
        assert_eq!(mesh.vertices[0], Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(mesh.vertices[1], Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh.vertices[2], Vec3::new(0.0, 1.0, 0.0));
        // The three midpoints are at pairwise averages (order by
        // canonical edge sort: (0,1), (0,2), (1,2)).
        let midpoints: Vec<Vec3> = mesh.vertices[3..].to_vec();
        let expected = [
            Vec3::new(0.5, 0.0, 0.0), // edge 0-1
            Vec3::new(0.0, 0.5, 0.0), // edge 0-2
            Vec3::new(0.5, 0.5, 0.0), // edge 1-2
        ];
        for (got, want) in midpoints.iter().zip(expected.iter()) {
            assert!(
                (*got - *want).length() < 1e-6,
                "midpoint {got:?} != expected {want:?}"
            );
        }
    }

    #[test]
    fn butterfly_isolated_orphan_vertex_stays_in_place() {
        // A triangle plus an unreferenced vertex. The orphan's position
        // is preserved at its original index (interpolation property).
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(99.0, 99.0, 99.0), // orphan
            ],
            indices: vec![[0, 1, 2]],
        };
        butterfly(&mut mesh, ButterflyParams::new(1));
        assert_eq!(mesh.vertices[3], Vec3::new(99.0, 99.0, 99.0));
    }

    // --------------------------------------------------------------
    // Adversarial: geometric invariants (plan §Adversarial).
    // --------------------------------------------------------------

    #[test]
    fn butterfly_preserves_flat_triangulated_plane() {
        // A 3x3 grid of vertices, triangulated flat. All z=0. After
        // subdivision, every new vertex must also have z=0 (weighted
        // sum of coplanar points is coplanar).
        let mut verts = Vec::new();
        for y in 0..3 {
            for x in 0..3 {
                verts.push(Vec3::new(x as f32, y as f32, 0.0));
            }
        }
        let idx = |x: u32, y: u32| y * 3 + x;
        let mut indices = Vec::new();
        for y in 0..2 {
            for x in 0..2 {
                indices.push([idx(x, y), idx(x + 1, y), idx(x, y + 1)]);
                indices.push([idx(x + 1, y), idx(x + 1, y + 1), idx(x, y + 1)]);
            }
        }
        let mut mesh = Mesh {
            vertices: verts,
            indices,
        };
        butterfly(&mut mesh, ButterflyParams::new(2));
        for v in &mesh.vertices {
            assert!(
                v.z.abs() < 1e-5,
                "flat plane lost planarity: vertex {v:?} has z = {}",
                v.z
            );
        }
    }

    // --------------------------------------------------------------
    // Property tests (plan §Property).
    // --------------------------------------------------------------

    use proptest::prelude::*;

    /// Small random meshes with valid index ranges.
    fn arb_small_mesh() -> impl Strategy<Value = Mesh> {
        let vert_strat = prop::collection::vec(
            (-5.0_f32..5.0, -5.0_f32..5.0, -5.0_f32..5.0).prop_map(|(x, y, z)| Vec3::new(x, y, z)),
            3..=10,
        );
        vert_strat.prop_flat_map(|verts| {
            let n = verts.len();
            let tri_strat = prop::collection::vec(
                (0..n, 0..n, 0..n).prop_map(|(a, b, c)| [a as u32, b as u32, c as u32]),
                1..=6,
            );
            tri_strat.prop_map(move |indices| Mesh {
                vertices: verts.clone(),
                indices,
            })
        })
    }

    proptest! {
        #[test]
        fn prop_butterfly_preserves_original_vertex_positions(
            mesh in arb_small_mesh(),
            iters in 0_u32..=2,
        ) {
            let before = mesh.clone();
            let mut m = mesh;
            butterfly(&mut m, ButterflyParams::new(iters));
            // First V_in positions must match input exactly.
            for (i, original) in before.vertices.iter().enumerate() {
                prop_assert_eq!(m.vertices[i], *original);
            }
        }

        #[test]
        fn prop_butterfly_indices_are_in_range(
            mesh in arb_small_mesh(),
            iters in 0_u32..=2,
        ) {
            let mut m = mesh;
            butterfly(&mut m, ButterflyParams::new(iters));
            let n = u32::try_from(m.vertices.len()).unwrap();
            for tri in &m.indices {
                for &idx in tri {
                    prop_assert!(idx < n);
                }
            }
        }

        #[test]
        fn prop_butterfly_quadruples_triangle_count(
            mesh in arb_small_mesh(),
            iters in 0_u32..=2,
        ) {
            let tri_before = mesh.triangle_count();
            let mut m = mesh;
            butterfly(&mut m, ButterflyParams::new(iters));
            prop_assert_eq!(m.triangle_count(), tri_before * 4_usize.pow(iters));
        }

        #[test]
        fn prop_butterfly_is_deterministic(
            mesh in arb_small_mesh(),
            iters in 0_u32..=2,
        ) {
            let src = mesh;
            let mut a = src.clone();
            let mut b = src.clone();
            butterfly(&mut a, ButterflyParams::new(iters));
            butterfly(&mut b, ButterflyParams::new(iters));
            prop_assert_eq!(a.vertices, b.vertices);
            prop_assert_eq!(a.indices, b.indices);
        }
    }
}
