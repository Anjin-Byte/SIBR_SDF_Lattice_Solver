#!/usr/bin/env python3
"""Manifoldness analyzer for binary STL files.

Reports edge / vertex topology after welding vertices to a quantization
grid. Distinguishes the four edge classes Blender lumps under "non-manifold":

  - boundary    : exactly one incident triangle (open mesh edge)
  - manifold    : exactly two incident triangles (the goal)
  - non-manifold: three or more incident triangles (genuine topology bug)
  - orient-err  : two triangles share an edge in the same direction
                  (one of them is wound backwards)

Usage:
    python3 check_manifold.py PATH.stl [PATH.stl ...]

Welds at three quantization scales (1e-6, 1e-5, 1e-4 absolute) so you can
see whether tightening / loosening the grid changes the topology
classification — useful for diagnosing sub-ULP welding issues.

Pure stdlib. Designed for wafer / small lattice samples (sub-million
triangles). For multi-million-triangle meshes it'll work but will take
multiple minutes — install numpy and use a numpy-based analyzer for that
scale.
"""

import collections
import struct
import sys
from pathlib import Path


def read_stl_binary(path):
    """Yields (v0, v1, v2) tuples for every triangle in the binary STL."""
    with open(path, "rb") as f:
        f.read(80)  # header
        n = struct.unpack("<I", f.read(4))[0]
        # Read everything in one shot for speed.
        raw = f.read(50 * n)
    # Each tri is 50 bytes: 12 normal + 36 verts (9 floats) + 2 attribute.
    # We only care about the vertex floats — at offset 12, length 36.
    unpack_verts = struct.Struct("<9f").unpack_from
    for i in range(n):
        off = i * 50 + 12
        v = unpack_verts(raw, off)
        yield (v[0], v[1], v[2]), (v[3], v[4], v[5]), (v[6], v[7], v[8])


def _find(parent, x):
    """Union-Find with path compression."""
    root = x
    while parent[root] != root:
        root = parent[root]
    while parent[x] != root:
        parent[x], x = root, parent[x]
    return root


def analyze(path, eps):
    """Read STL at `path`, weld at quantization `eps`, classify edges
    and count connected components (shells)."""
    inv = 1.0 / eps
    keymap = {}
    n_tris = 0
    edge_counts = collections.Counter()  # frozenset({a,b}) -> incident triangle count
    edge_dirs = collections.Counter()  # (a,b) directed -> count
    n_degenerate = 0
    parent = []  # Union-Find structure for connected components

    for v0, v1, v2 in read_stl_binary(path):
        n_tris += 1
        # Quantize each vertex via integer triple as dict key.
        idx = []
        for vx, vy, vz in (v0, v1, v2):
            key = (int(vx * inv + 0.5 if vx >= 0 else vx * inv - 0.5),
                   int(vy * inv + 0.5 if vy >= 0 else vy * inv - 0.5),
                   int(vz * inv + 0.5 if vz >= 0 else vz * inv - 0.5))
            i = keymap.get(key)
            if i is None:
                i = len(keymap)
                keymap[key] = i
                parent.append(i)
            idx.append(i)
        a, b, c = idx
        if a == b or b == c or a == c:
            n_degenerate += 1
            continue
        for u, v in ((a, b), (b, c), (c, a)):
            edge_dirs[(u, v)] += 1
            edge_counts[frozenset((u, v))] += 1
        # Union-Find: merge the three vertex sets via this triangle.
        ra, rb, rc = _find(parent, a), _find(parent, b), _find(parent, c)
        if ra != rb:
            parent[rb] = ra
        if rc != _find(parent, a):
            parent[rc] = ra

    boundary = manifold = nonmanifold = 0
    for n in edge_counts.values():
        if n == 1:
            boundary += 1
        elif n == 2:
            manifold += 1
        else:
            nonmanifold += 1

    # Orientation errors: for each manifold edge, the two triangles using
    # it should traverse it in opposite directions. Count edges where
    # forward count != reverse count.
    orient_errors = 0
    for e, n in edge_counts.items():
        if n != 2:
            continue
        u, v = tuple(e)
        f = edge_dirs.get((u, v), 0)
        r = edge_dirs.get((v, u), 0)
        if f != r:
            orient_errors += 1

    # Count distinct components: only consider vertices that appear in
    # at least one non-degenerate triangle (orphans don't count).
    referenced = set()
    for e in edge_counts:
        referenced.update(e)
    roots = set(_find(parent, v) for v in referenced)
    shells = len(roots)

    return {
        "verts": len(keymap),
        "tris": n_tris,
        "degenerate": n_degenerate,
        "edges_total": len(edge_counts),
        "boundary": boundary,
        "manifold": manifold,
        "nonmanifold": nonmanifold,
        "orient_errors": orient_errors,
        "shells": shells,
    }


def fmt(r):
    return (
        f"verts={r['verts']:>9,}  tris={r['tris']:>9,}  degen={r['degenerate']:>5,}  "
        f"shells={r['shells']:>5,}  boundary={r['boundary']:>7,}  "
        f"manifold={r['manifold']:>9,}  non-manifold={r['nonmanifold']:>5,}  "
        f"orient-err={r['orient_errors']:>5,}"
    )


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__ or "")
        return 1
    for path in argv[1:]:
        p = Path(path)
        if not p.exists():
            sys.stderr.write(f"warn: {path} not found, skipping\n")
            continue
        print(f"=== {path} ({p.stat().st_size / 1e6:.1f} MB) ===")
        for eps in (1e-6, 1e-5, 1e-4):
            r = analyze(p, eps)
            print(f"  weld eps={eps:.0e}  {fmt(r)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
