---
title: GPU Compute Pipeline
tags: [domain, gpu, performance, wgpu]
---

# GPU Compute Pipeline

The end-to-end structure of SDF evaluation → isosurface extraction → mesh export on GPU. This note explains **why the approach is performant** and, more importantly, **the three places performance is lost if the implementation is sloppy**.

## The two-pass structure

The production path is **two wgpu compute passes**, plus pre- and post-processing stages:

```
   [CLI input: LatticeJob]
           │
           ▼
   [Pass 0 (CPU): Job compilation]
   Lower LatticeJob into a GPU-friendly SDF graph
           │
           ▼
   [Pass 1 (GPU compute): Field evaluation]
   For each active voxel, evaluate composed SDF
   Output: sparse bricked scalar field
           │
           ▼
   [Pass 2 (GPU compute): Isosurface extraction (MC33)]
   For each sign-change voxel, emit triangles
   Output: vertex buffer + index buffer (with duplicates)
           │
           ▼
   [Pass 3 (GPU compute or CPU): Vertex welding]
   Scan + hash dedup
   Output: compact indexed mesh
           │
           ▼
   [Pass 4 (CPU): Mesh decimation]
   QEM to printer target resolution
           │
           ▼
   [Pass 5 (CPU): Export]
   STL / 3MF write
```

Passes 1 and 2 are the heavy compute. Passes 3 and 4 are mesh-domain operations — covered in [[Mesh Decimation]].

## Why this is performant

### Pass 1: SDF field evaluation

Per voxel:
- Evaluate `O(k)` primitive SDFs (struts near the query point, found via `opRepetition`)
- Combine via SDF Booleans — `min`/`max`/`smin` (see [[SDF vs Mesh Booleans]])
- Intersect against primitive boundary SDF

**Hardware fit:**
- Per-voxel uniform work — no lane divergence within a SIMD group.
- No inter-voxel dependencies — `O(G_s)` work distributed over tens of thousands of GPU threads.
- Memory access pattern is write-mostly (output field), with small read of constants (strut parameters, cell layout).

**Rough numbers on a 10 TFLOPS GPU:**
- ~10⁸ voxels evaluated per second for `k ≈ 3`, cylindrical strut SDFs.
- For the [[Meshing Complexity Analysis|Woodward Fig. 4A example]] (~10⁷ active voxels): ~100 ms.

### Pass 2: MC33 extraction

Per voxel:
- Read 8 corner field values
- Compute 8-bit sign index
- Table lookup: edges with sign changes + triangle configuration
- Interpolate vertex positions on sign-change edges
- Emit up to 5 triangles

**Hardware fit:**
- Per-voxel independent — even more parallel than field evaluation.
- Memory-bound: field reads + triangle writes dominate, not ALU.
- Table lookups are constant-time indexed reads.

**Rough numbers:**
- ~10⁸–10⁹ voxels extracted per second; often faster than Pass 1 on the same data.
- For 10⁷ active voxels: 10–100 ms.

### Together

End-to-end compute time for the worked example: **~100–300 ms**. Output mesh transfer to host (for export) adds tens of ms more. Total **well under a second** for Woodward-scale parts on consumer hardware.

## The three places performance is lost

If you see 10× slower numbers than the above, the cause is almost always one of these:

### 1. Not using `opRepetition`

**Symptom:** per-voxel cost scales with total strut count `N_struts`, not struts per cell `k`.

For a Woodward Fig. 4A cylinder with 21,000 struts, a naive implementation evaluates all 21,000 per voxel. With `opRepetition` + intersection with primitive boundary, each voxel evaluates only the ~3 struts in its local cell neighborhood. **Factor: ~7,000×.**

The full expression for a repeated lattice body:
```glsl
float sdLatticeBody(vec3 p) {
    // Shift point into fundamental cell of infinite repetition
    vec3 cell = round(p / L_c);
    vec3 q = p - L_c * cell;
    // Evaluate k struts in this cell
    float d = /* min over cell's k strut SDFs evaluated at q */;
    // Intersect against primitive boundary
    return max(d, sdPrimitive(p));
}
```

No broadphase, no data structures, no per-cell iteration — the `round()` puts you in the right cell, period.

### 2. Dense grid over mostly-empty volume

**Symptom:** memory usage and eval time scale with primitive volume `V`, not surface area `S`.

For a lattice with porosity 0.9, >99% of voxels are either fully inside void or fully inside the primitive — no sign change at the corners. Evaluating the SDF at all of them is waste.

**Fix:** bricked sparse activation. Coarsely pre-evaluate the SDF on an 8× or 16× coarser grid; activate only bricks where `|f_coarse| < threshold × h_coarse`. Then run the fine pass only on activated bricks. Typical savings: **10–50×** depending on porosity.

This also means the field memory scales with `S/h²`, not `V/h³` — the difference between "fits in 8 GB" and "needs 200 GB." See [[Meshing Complexity Analysis]] Method B.

### 3. Naive vertex output

**Symptom:** the output mesh has 3× more vertices than needed; later mesh processing (decimation, export) is disproportionately slow.

Raw MC emits 3 vertices per triangle, independently per voxel. An edge shared by adjacent voxels is interpolated twice, producing duplicate vertices at the shared edge. For a typical mesh, this inflates vertex count ~3× and leaves the topology implicit (the mesh has no shared vertices, only coincident duplicates).

**Fix:** a parallel hash-welding pass. Hash each vertex position (quantized to ~h/10 to handle floating-point jitter); deduplicate; rewrite index buffer. This is Pass 3 in the pipeline above.

Alternative: per-edge vertex allocation — each grid edge with a sign change owns exactly one vertex, shared across all triangles that touch it. Cleaner but needs a scan pass to assign vertex indices.

Either way: **do not emit vertices per triangle; emit per sign-change edge.**

## Handoff to CPU

After Pass 3, the mesh is a compact indexed structure in GPU memory. It must move to host for:
- Decimation (QEM requires a priority queue; GPU-friendly but not needed at the scale involved)
- Export (file I/O)
- Optional preview rendering (stays on GPU)

Transfer cost is `O(|V| + |F|) × bytes` — typically 50–500 MB for Woodward-scale parts. On PCIe 4.0 x16 (~32 GB/s peak, ~20 GB/s effective), this is 25–250 ms. Non-negligible but bounded.

## wgpu-specific considerations

- **Compute dispatch size** limited to 65,535 workgroups per dimension per dispatch. For >65k-brick workloads, split into multiple dispatches or use indirect dispatch. Not a hard blocker at Woodward scales.
- **Storage buffer size limits** vary by backend (Vulkan: typically 2 GB+ per buffer; D3D12: 4 GB per buffer). Bricked storage avoids single-buffer ceilings anyway.
- **Atomic operations** on u32 are universally supported; on f32 more patchy. Design the welding pass around u32 atomics (hash tables).
- **Subgroup operations** (wave-level primitives like ballot, shuffle) are useful for welding and broadphase but not universally available. Keep them behind capability checks.

All of this is implementation hygiene, not design-level concern — but worth flagging so the gpu crate's capability-detection layer is built in early.

## Related

- [[Meshing Complexity Analysis]]
- [[Isosurface Extraction Methods]]
- [[SDF vs Mesh Booleans]]
- [[Mesh Decimation]]
- [[SDF Primitives Catalog]]
- [[gpu crate]]
- [[wgpu]]
