//! Taubin low-pass mesh smoothing (Taubin 1995).
//!
//! Classical Laplacian smoothing moves each vertex toward the centroid of
//! its graph neighbors. Applied iteratively, it low-pass-filters the mesh
//! but also causes global shrinkage (a sphere gets smaller). Taubin's
//! contribution is a two-pass iteration that cancels the shrinkage while
//! retaining the low-pass property:
//!
//! ```text
//! Pass 1:  v ← v + λ · Δ(v)       (λ > 0 — smoothing, shrinks)
//! Pass 2:  v ← v + μ · Δ(v)       (μ < 0, |μ| > λ — anti-shrinks)
//! ```
//!
//! where `Δ(v_i) = mean(v_j : j ∈ neighbors(i)) - v_i` is the umbrella
//! operator (discrete Laplacian). `Δ` is recomputed between passes on the
//! updated positions.
//!
//! # Prerequisite: welded mesh
//!
//! Smoothing operates via **graph adjacency**, which requires vertices that
//! are geometrically coincident to share an index. The upstream welding
//! stage ([`super::weld`]) provides that guarantee. Calling [`taubin`] on
//! an unwelded mesh does nothing useful — each triangle's three vertices
//! have no graph neighbors outside the triangle itself, so the umbrella
//! operator is near-zero everywhere.
//!
//! # Parameter choice
//!
//! The `(λ, μ)` pair controls the filter's cutoff frequency. Taubin 1995
//! recommends `λ = 0.5`, `μ ≈ -0.5263`, giving a pass-band frequency
//! `k_PB ≈ 0.1` (frequencies below 10% of the Nyquist limit pass with
//! near-unit gain; frequencies above are attenuated). Our default
//! constructor [`TaubinParams::default_with_iterations`] uses these
//! values; callers who want different behavior can use
//! [`TaubinParams::new`].
//!
//! # Reference
//!
//! G. Taubin, "Curve and Surface Smoothing without Shrinkage",
//! ICCV 1995.

use glam::Vec3;
use thiserror::Error;

use super::Mesh;
use crate::progress::Progress;

/// Taubin defaults from the 1995 paper: pass-band frequency `k_PB ≈ 0.1`.
const DEFAULT_LAMBDA: f32 = 0.5;
const DEFAULT_MU: f32 = -0.5263;

/// Construction failure for [`TaubinParams`]. All variants originate from
/// the validating [`TaubinParams::new`] constructor.
#[derive(Debug, Clone, Copy, Error, PartialEq)]
pub enum SmoothError {
    /// A float parameter was NaN or ±∞.
    #[error("smooth: parameter {field} must be finite, got {value}")]
    NonFinite {
        /// Name of the offending parameter.
        field: &'static str,
        /// The received value.
        value: f32,
    },
    /// λ must be strictly positive (Taubin's smoothing-pass scalar).
    #[error("smooth: lambda must be > 0, got {value}")]
    LambdaNonPositive {
        /// The received value.
        value: f32,
    },
    /// μ must be strictly negative (Taubin's anti-shrink pass scalar).
    #[error("smooth: mu must be < 0, got {value}")]
    MuNonNegative {
        /// The received value.
        value: f32,
    },
    /// Taubin's anti-shrinkage criterion requires `|μ| > λ`.
    #[error("smooth: |mu| must be > lambda (|{mu}| vs {lambda})")]
    InsufficientAntiShrink {
        /// The λ parameter.
        lambda: f32,
        /// The μ parameter.
        mu: f32,
    },
}

/// Parameters for Taubin smoothing. Construct via [`TaubinParams::new`]
/// so validation runs at construction time — the algorithm assumes a
/// valid parameter set and does not re-check at runtime.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TaubinParams {
    lambda: f32,
    mu: f32,
    iterations: u32,
}

impl TaubinParams {
    /// Constructs a validated Taubin parameter set.
    ///
    /// # Errors
    ///
    /// - [`SmoothError::NonFinite`] if `lambda` or `mu` is NaN or ±∞.
    /// - [`SmoothError::LambdaNonPositive`] if `lambda <= 0`.
    /// - [`SmoothError::MuNonNegative`] if `mu >= 0`.
    /// - [`SmoothError::InsufficientAntiShrink`] if `|mu| <= lambda`.
    pub fn new(lambda: f32, mu: f32, iterations: u32) -> Result<Self, SmoothError> {
        if !lambda.is_finite() {
            return Err(SmoothError::NonFinite {
                field: "lambda",
                value: lambda,
            });
        }
        if !mu.is_finite() {
            return Err(SmoothError::NonFinite {
                field: "mu",
                value: mu,
            });
        }
        if lambda <= 0.0 {
            return Err(SmoothError::LambdaNonPositive { value: lambda });
        }
        if mu >= 0.0 {
            return Err(SmoothError::MuNonNegative { value: mu });
        }
        if mu.abs() <= lambda {
            return Err(SmoothError::InsufficientAntiShrink { lambda, mu });
        }
        Ok(Self {
            lambda,
            mu,
            iterations,
        })
    }

    /// Taubin 1995's recommended defaults (`λ = 0.5`, `μ = -0.5263`,
    /// pass-band `k_PB ≈ 0.1`) with the requested iteration count.
    ///
    /// # Errors
    ///
    /// Cannot fail for any finite `iterations`. Returns
    /// `Result<Self, SmoothError>` for signature consistency with
    /// [`Self::new`].
    pub fn default_with_iterations(iterations: u32) -> Result<Self, SmoothError> {
        Self::new(DEFAULT_LAMBDA, DEFAULT_MU, iterations)
    }

    /// Smoothing-pass scalar (λ > 0).
    pub fn lambda(&self) -> f32 {
        self.lambda
    }
    /// Anti-shrink-pass scalar (μ < 0, |μ| > λ).
    pub fn mu(&self) -> f32 {
        self.mu
    }
    /// Number of Taubin iterations (each iteration is one λ pass + one
    /// μ pass).
    pub fn iterations(&self) -> u32 {
        self.iterations
    }
}

/// Applies `params.iterations` iterations of Taubin low-pass smoothing
/// to `mesh`. Topology (`mesh.indices`) is unchanged; only vertex
/// positions are mutated.
///
/// See the module-level documentation for algorithm details.
///
/// # Performance
///
/// One-time adjacency build: `O(T)` where T is triangle count.
/// Per-iteration cost: `O(V + E)` where V is vertex count and E is the
/// total adjacency edges. On a ~300k-vertex welded mesh, 10 iterations
/// complete in well under 100 ms on a modern CPU.
pub fn taubin(mesh: &mut Mesh, params: TaubinParams) {
    taubin_with_progress(mesh, params, &mut ());
}

/// Like [`taubin`], but reports progress to `progress`. Ticks once per
/// λ+μ iteration; declared length is `params.iterations()`.
pub fn taubin_with_progress(mesh: &mut Mesh, params: TaubinParams, progress: &mut impl Progress) {
    progress.set_len(u64::from(params.iterations));
    if params.iterations == 0 || mesh.vertices.is_empty() {
        progress.finish();
        return;
    }

    let adjacency = build_adjacency(mesh);

    // Scratch buffer reused across iterations. Allocated once.
    let mut deltas: Vec<Vec3> = vec![Vec3::ZERO; mesh.vertices.len()];

    for _ in 0..params.iterations {
        // Pass 1: λ (smoothing, shrinks).
        compute_laplacian(mesh, &adjacency, &mut deltas);
        apply_scaled(mesh, &deltas, params.lambda);
        // Pass 2: μ (anti-shrink). Recomputed Δ on updated positions.
        compute_laplacian(mesh, &adjacency, &mut deltas);
        apply_scaled(mesh, &deltas, params.mu);
        progress.inc(1);
    }
    progress.finish();
}

/// Builds per-vertex neighbor lists from the mesh's triangle indices.
/// Neighbors are sorted and deduplicated so subsequent centroid
/// computations touch each neighbor exactly once.
fn build_adjacency(mesh: &Mesh) -> Vec<Vec<u32>> {
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); mesh.vertices.len()];
    for tri in &mesh.indices {
        let [a, b, c] = *tri;
        adj[a as usize].push(b);
        adj[a as usize].push(c);
        adj[b as usize].push(a);
        adj[b as usize].push(c);
        adj[c as usize].push(a);
        adj[c as usize].push(b);
    }
    for list in &mut adj {
        list.sort_unstable();
        list.dedup();
    }
    adj
}

/// Writes `deltas[i] = centroid(neighbors of i) - mesh.vertices[i]`
/// for every vertex. Vertices with no neighbors get `Δ = Vec3::ZERO`.
#[allow(clippy::cast_precision_loss)]
fn compute_laplacian(mesh: &Mesh, adjacency: &[Vec<u32>], deltas: &mut [Vec3]) {
    for (i, neighbors) in adjacency.iter().enumerate() {
        if neighbors.is_empty() {
            deltas[i] = Vec3::ZERO;
            continue;
        }
        let mut sum = Vec3::ZERO;
        for &j in neighbors {
            sum += mesh.vertices[j as usize];
        }
        let centroid = sum / neighbors.len() as f32;
        deltas[i] = centroid - mesh.vertices[i];
    }
}

/// In-place `mesh.vertices[i] += scale * deltas[i]` for all i.
fn apply_scaled(mesh: &mut Mesh, deltas: &[Vec3], scale: f32) {
    for (v, d) in mesh.vertices.iter_mut().zip(deltas.iter()) {
        *v += scale * *d;
    }
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

    fn tetrahedron() -> Mesh {
        // 4 vertices, 4 triangular faces (outward winding).
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
    // TaubinParams validation tests.
    // --------------------------------------------------------------

    #[test]
    fn params_default_values_match_taubin_1995() {
        let p = TaubinParams::default_with_iterations(1).unwrap();
        assert_eq!(p.lambda(), DEFAULT_LAMBDA);
        assert_eq!(p.mu(), DEFAULT_MU);
        // Anti-shrink criterion satisfied.
        assert!(p.mu().abs() > p.lambda());
    }

    #[test]
    fn params_new_rejects_nan_lambda() {
        assert!(matches!(
            TaubinParams::new(f32::NAN, -0.5, 1),
            Err(SmoothError::NonFinite {
                field: "lambda",
                ..
            })
        ));
    }

    #[test]
    fn params_new_rejects_infinite_mu() {
        assert!(matches!(
            TaubinParams::new(0.5, f32::NEG_INFINITY, 1),
            Err(SmoothError::NonFinite { field: "mu", .. })
        ));
    }

    #[test]
    fn params_new_rejects_lambda_zero_or_negative() {
        assert!(matches!(
            TaubinParams::new(0.0, -0.6, 1),
            Err(SmoothError::LambdaNonPositive { .. })
        ));
        assert!(matches!(
            TaubinParams::new(-0.1, -0.6, 1),
            Err(SmoothError::LambdaNonPositive { .. })
        ));
    }

    #[test]
    fn params_new_rejects_mu_zero_or_positive() {
        assert!(matches!(
            TaubinParams::new(0.5, 0.0, 1),
            Err(SmoothError::MuNonNegative { .. })
        ));
        assert!(matches!(
            TaubinParams::new(0.5, 0.1, 1),
            Err(SmoothError::MuNonNegative { .. })
        ));
    }

    #[test]
    fn params_new_rejects_abs_mu_leq_lambda() {
        // λ=0.5, μ=-0.4 → |μ|=0.4 < λ=0.5 → fails.
        assert!(matches!(
            TaubinParams::new(0.5, -0.4, 1),
            Err(SmoothError::InsufficientAntiShrink { .. })
        ));
        // λ=0.5, μ=-0.5 → |μ|=λ=0.5 → fails (strict inequality).
        assert!(matches!(
            TaubinParams::new(0.5, -0.5, 1),
            Err(SmoothError::InsufficientAntiShrink { .. })
        ));
    }

    #[test]
    fn params_new_accepts_valid_inputs() {
        let p = TaubinParams::new(0.5, -0.5263, 10).unwrap();
        assert_eq!(p.iterations(), 10);
    }

    // --------------------------------------------------------------
    // Topology preservation and degenerate-input behavior.
    // --------------------------------------------------------------

    #[test]
    fn taubin_preserves_topology() {
        let before = tetrahedron();
        let mut mesh = before.clone();
        let params = TaubinParams::default_with_iterations(5).unwrap();
        taubin(&mut mesh, params);
        assert_eq!(mesh.indices, before.indices);
        assert_eq!(mesh.vertices.len(), before.vertices.len());
    }

    #[test]
    fn taubin_identity_on_zero_iterations() {
        let before = tetrahedron();
        let mut mesh = before.clone();
        let params = TaubinParams::default_with_iterations(0).unwrap();
        taubin(&mut mesh, params);
        assert_eq!(mesh.vertices, before.vertices);
        assert_eq!(mesh.indices, before.indices);
    }

    #[test]
    fn taubin_handles_empty_mesh() {
        let mut mesh = Mesh::default();
        let params = TaubinParams::default_with_iterations(5).unwrap();
        taubin(&mut mesh, params);
        assert!(mesh.vertices.is_empty());
        assert!(mesh.indices.is_empty());
    }

    #[test]
    fn taubin_leaves_isolated_vertex_in_place() {
        // A triangle (0,1,2) + an orphan vertex at index 3.
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(99.0, 99.0, 99.0), // orphan
            ],
            indices: vec![[0, 1, 2]],
        };
        let orphan_before = mesh.vertices[3];
        let params = TaubinParams::default_with_iterations(10).unwrap();
        taubin(&mut mesh, params);
        // The orphan has no neighbors → Δ = 0 → position unchanged.
        assert_eq!(mesh.vertices[3], orphan_before);
    }

    // --------------------------------------------------------------
    // Smoothing behavior: high-frequency noise attenuated.
    // --------------------------------------------------------------

    #[test]
    fn taubin_moves_noisy_vertex_toward_neighbor_plane() {
        // Six vertices forming a regular pattern in the z=0 plane, with
        // one center vertex lifted to z = 1 (high-frequency "spike").
        // Taubin should pull the spike down toward its neighbors' plane.
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(0.0, 0.0, 1.0),  // 0 — the spike
                Vec3::new(1.0, 0.0, 0.0),  // 1
                Vec3::new(-1.0, 0.0, 0.0), // 2
                Vec3::new(0.0, 1.0, 0.0),  // 3
                Vec3::new(0.0, -1.0, 0.0), // 4
            ],
            indices: vec![[0, 1, 3], [0, 3, 2], [0, 2, 4], [0, 4, 1]],
        };
        let z_before = mesh.vertices[0].z;
        let params = TaubinParams::default_with_iterations(5).unwrap();
        taubin(&mut mesh, params);
        let z_after = mesh.vertices[0].z;
        // The spike started at z=1; neighbors all at z=0. Taubin should
        // pull it strictly toward zero (without overshooting past zero
        // at these modest parameters).
        assert!(
            z_after < z_before,
            "spike should move toward neighbors' plane: before={z_before}, after={z_after}"
        );
        assert!(
            z_after > 0.0,
            "5 iterations shouldn't overshoot past the target plane: after={z_after}"
        );
    }

    // --------------------------------------------------------------
    // Determinism.
    // --------------------------------------------------------------

    // --------------------------------------------------------------
    // Progress plumbing (Level 1).
    // --------------------------------------------------------------

    #[test]
    fn taubin_with_progress_ticks_per_iteration() {
        use crate::progress::Spy;
        let mut mesh = tetrahedron();
        let params = TaubinParams::default_with_iterations(5).unwrap();
        let mut spy = Spy::default();
        taubin_with_progress(&mut mesh, params, &mut spy);
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, 5);
        assert_eq!(spy.inc_sum, 5);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn taubin_with_progress_zero_iterations_still_finishes() {
        use crate::progress::Spy;
        let mut mesh = tetrahedron();
        let params = TaubinParams::default_with_iterations(0).unwrap();
        let mut spy = Spy::default();
        taubin_with_progress(&mut mesh, params, &mut spy);
        assert_eq!(spy.set_len_calls, 1);
        assert_eq!(spy.total, 0);
        assert_eq!(spy.inc_sum, 0);
        assert_eq!(spy.finish_calls, 1);
    }

    #[test]
    fn taubin_with_progress_empty_mesh_still_finishes() {
        use crate::progress::Spy;
        let mut mesh = Mesh::default();
        let params = TaubinParams::default_with_iterations(5).unwrap();
        let mut spy = Spy::default();
        taubin_with_progress(&mut mesh, params, &mut spy);
        assert_eq!(spy.finish_calls, 1);
        assert_eq!(spy.inc_sum, 0);
    }

    #[test]
    fn taubin_is_deterministic_on_identical_input() {
        let src = tetrahedron();
        let params = TaubinParams::default_with_iterations(3).unwrap();
        let mut a = src.clone();
        let mut b = src.clone();
        taubin(&mut a, params);
        taubin(&mut b, params);
        assert_eq!(a.vertices, b.vertices);
        assert_eq!(a.indices, b.indices);
    }

    // --------------------------------------------------------------
    // Property tests (proptest).
    // --------------------------------------------------------------

    use proptest::prelude::*;

    /// Generate small welded meshes: 4–10 vertices + up to 8 triangles,
    /// triangle indices always in range. Not every generated mesh is
    /// manifold; Taubin's topology-preservation invariant holds
    /// regardless.
    fn arb_small_mesh() -> impl Strategy<Value = Mesh> {
        let vert_strat = prop::collection::vec(
            (-5.0_f32..5.0, -5.0_f32..5.0, -5.0_f32..5.0).prop_map(|(x, y, z)| Vec3::new(x, y, z)),
            4..=10,
        );
        vert_strat.prop_flat_map(|verts| {
            let n = verts.len();
            let tri_strat = prop::collection::vec(
                (0..n, 0..n, 0..n).prop_map(|(a, b, c)| [a as u32, b as u32, c as u32]),
                1..=8,
            );
            tri_strat.prop_map(move |indices| Mesh {
                vertices: verts.clone(),
                indices,
            })
        })
    }

    fn arb_valid_params() -> impl Strategy<Value = TaubinParams> {
        (0.1_f32..0.9_f32, 1_u32..=5).prop_map(|(lambda, iters)| {
            // Choose μ so |μ| > λ by a small margin (Taubin-style).
            let mu = -(lambda + 0.05_f32.max(lambda * 0.1));
            TaubinParams::new(lambda, mu, iters).unwrap()
        })
    }

    proptest! {
        /// Topology (index list) must survive any valid smoothing run.
        #[test]
        fn prop_taubin_preserves_indices(mesh in arb_small_mesh(), params in arb_valid_params()) {
            let indices_before = mesh.indices.clone();
            let verts_before_len = mesh.vertices.len();
            let mut m = mesh;
            taubin(&mut m, params);
            prop_assert_eq!(m.indices, indices_before);
            prop_assert_eq!(m.vertices.len(), verts_before_len);
        }

        /// No vertex should move more than `iterations * (λ + |μ|) *
        /// max_original_edge_length` in total — a loose-but-useful
        /// sanity bound on displacement.
        #[test]
        fn prop_taubin_vertex_motion_bounded(mesh in arb_small_mesh(), params in arb_valid_params()) {
            let before = mesh.clone();
            let mut m = mesh;
            // Longest edge in the *original* mesh bounds the Laplacian
            // magnitude at the first iteration.
            let max_edge = before.indices.iter().flat_map(|tri| {
                let v0 = before.vertices[tri[0] as usize];
                let v1 = before.vertices[tri[1] as usize];
                let v2 = before.vertices[tri[2] as usize];
                [(v1 - v0).length(), (v2 - v1).length(), (v0 - v2).length()]
            }).fold(0.0_f32, f32::max);
            // Displacement bound (generous to accommodate growing Δ
            // under repeated iteration — uses 3× slack).
            let bound = params.iterations() as f32
                * (params.lambda() + params.mu().abs())
                * max_edge
                * 3.0;
            taubin(&mut m, params);
            for (a, b) in before.vertices.iter().zip(m.vertices.iter()) {
                let dist = (*b - *a).length();
                prop_assert!(
                    dist <= bound + 1e-4,
                    "vertex moved {dist} > bound {bound}"
                );
            }
        }

        /// Running Taubin twice on the same input must produce identical
        /// output.
        #[test]
        fn prop_taubin_is_deterministic(mesh in arb_small_mesh(), params in arb_valid_params()) {
            let src = mesh;
            let mut a = src.clone();
            let mut b = src.clone();
            taubin(&mut a, params);
            taubin(&mut b, params);
            prop_assert_eq!(a.vertices, b.vertices);
            prop_assert_eq!(a.indices, b.indices);
        }
    }
}
