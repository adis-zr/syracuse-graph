#!/usr/bin/env python3
"""
bijection.py — Verify the history bijection (Theorem 5.3 of the paper).

For each depth d and distance k, the history map

    v  ->  (eps(v), Sigma(v))

is a bijection from the depth-d, distance-k layer to {+,-} x S_{d,k},
where

    S_{d,k} = {Sigma subset [d] : d in Sigma, |Sigma| = 2k+2}.

Consequently the layer has size 2 * C(d-1, 2k+1).

This script:
  1. Builds G_D and computes distances to C_D.
  2. For every vertex v outside C_D, computes (eps(v), Sigma(v)).
  3. Checks that the map is a bijection onto {+,-} x S_{d,k}.
  4. Checks layer sizes against the binomial formula.

Usage:
    python bijection.py [D]          # default D = 16
"""

from __future__ import annotations

import sys
from math import comb
from itertools import combinations
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

from graph import (
    Vertex, build, tarjan_scc, bfs_distance_to, boundary_set, v2
)


# ---------------------------------------------------------------------------
# Seed and orientation
# ---------------------------------------------------------------------------

def seed_and_orientation(v0: Vertex, g) -> Tuple[FrozenSet[int], str]:
    """Return (Seed(v0), eps(v0)) for a vertex v0 in C_D (§5.1).

    Seed(v0) = {a, d} where a is the tau_d-orbit label (1 <= a < d).
    eps(v0)  = '-' if r < 2^{d-1}, '+' if r > 2^{d-1}.

    The d-1 tau_d-orbits at depth d are (Proposition 5.1):
      - For a = 1..d-2: {unresolved type a, direct-return sibling} -> label {a, d}
      - Remaining orbit: {unresolved type d-1, unresolved type d} -> label {d-1, d}

    So:
      - Unresolved type t (1 <= t <= d-1): seed = {t, d}
      - Unresolved type d (maximal): seed = {d-1, d}
      - Direct return: parent was unresolved type a at depth d-1, seed = {a, d}
        where a = v2(r+1) (the t-value of the direct return's forward edge)
    """
    r, d = v0
    eps = '-' if r < (1 << (d - 1)) else '+'

    t_type = g.exc_type.get(v0)
    if t_type is not None:
        # Unresolved vertex: type t
        # If t == d (maximal), it shares orbit {d-1, d} with type d-1
        a = t_type if t_type < d else d - 1
        return frozenset({a, d}), eps

    # Direct return: resolved vertex whose forward edge lands at (1,1)
    # t = v2(r+1) is the type of its unresolved sibling = orbit label a
    t = v2(r + 1)
    a = t
    return frozenset({a, d}), eps


# ---------------------------------------------------------------------------
# History map
# ---------------------------------------------------------------------------

def history(v: Vertex, g, C: set) -> Optional[Tuple[str, FrozenSet[int]]]:
    """Compute (eps(v), Sigma(v)) for a vertex v outside C_D.

    Follows the valuation chain v = v_k -> ... -> v_0 in C_D and assembles:
        Sigma(v) = Seed(v_0) union {d_{m-1}+t_m, d_m}_{m=1..k}

    Returns None if v is in C_D or if the chain cannot be traced.
    """
    if v in C:
        return None

    chain: List[Tuple[Vertex, int, int]] = []  # (source_vertex, t, rho)
    cur = v
    seen: set = set()
    while cur not in C:
        if cur in seen:
            return None  # cycle outside C (shouldn't happen for finite-dist vertices)
        seen.add(cur)
        nbrs = g.adj.get(cur, [])
        if not nbrs:
            return None
        nxt = nbrs[0]
        edge = g.edge_data.get((cur, nxt))
        if edge is None or edge.kind != "forward":
            return None
        chain.append((cur, edge.t, edge.rho))
        cur = nxt

    v0 = cur
    seed, eps = seed_and_orientation(v0, g)

    sigma_elements = set(seed)
    prev_depth = v0[1]  # d_0
    for (src, t, rho) in reversed(chain):
        d_m = src[1]
        i = prev_depth + t    # d_{m-1} + t_m
        j = d_m               # d_m
        sigma_elements.add(i)
        sigma_elements.add(j)
        prev_depth = d_m

    return eps, frozenset(sigma_elements)


# ---------------------------------------------------------------------------
# Admissible subsets S_{d,k}
# ---------------------------------------------------------------------------

def admissible_subsets(d: int, k: int) -> Set[FrozenSet[int]]:
    """Enumerate S_{d,k}: all (2k+2)-element subsets of [d] containing d."""
    return {frozenset(rest) | {d}
            for rest in combinations(range(1, d), 2 * k + 1)}


# ---------------------------------------------------------------------------
# Main: verify bijection theorem
# ---------------------------------------------------------------------------

def verify_bijection(D: int = 16):
    """
    History bijection (§6): for each (d, k), the history map is a
    bijection from depth-d distance-k vertices to {+,-} x S_{d,k},
    and the layer size equals 2 * C(d-1, 2k+1).
    """
    g = build(D)
    comps = tarjan_scc(g.vertices, g.adj)
    C = set(next(c for c in comps if (1, 1) in c))
    dist = bfs_distance_to(g.vertices, g.adj, C)

    print(f"Verifying History Bijection (D={D})")
    print(f"{'d':>3} {'k':>3}  {'|layer|':>8}  {'2*C(d-1,2k+1)':>14}  "
          f"{'injective':>10}  {'surjective':>11}  {'result':>6}")
    print("-" * 62)

    all_pass = True

    # Group vertices outside C_D by (depth, distance)
    layers: Dict[Tuple[int, int], List[Vertex]] = {}
    for v in g.vertices:
        if v in C:
            continue
        k = dist.get(v)
        if k is None:
            continue  # boundary-trapped
        d = v[1]
        layers.setdefault((d, k), []).append(v)

    for (d, k), verts in sorted(layers.items()):
        target = admissible_subsets(d, k)
        expected_size = comb(d - 1, 2 * k + 1)

        image: Dict[Tuple[str, FrozenSet[int]], Vertex] = {}
        injective = True
        bad_image = False

        for v in verts:
            coord = history(v, g, C)
            if coord is None:
                bad_image = True
                continue
            eps, sigma = coord
            key = (eps, sigma)
            if key in image:
                injective = False
            image[key] = v
            # Check sigma is in the right family
            if d not in sigma or len(sigma) != 2 * k + 2:
                bad_image = True

        # Check image equals {+,-} x S_{d,k}
        image_sigmas = {sigma for _, sigma in image}
        surjective = (image_sigmas == target) and not bad_image
        size_ok = (len(verts) == 2 * expected_size)

        ok = injective and surjective and size_ok and not bad_image
        all_pass = all_pass and ok

        print(f"  d={d:2d}, k={k}:  {len(verts):>8}  {2*expected_size:>14}  "
              f"{'yes' if injective else 'NO':>10}  "
              f"{'yes' if surjective else 'NO':>11}  "
              f"{'PASS' if ok else 'FAIL':>6}")

    print()
    print("All layers passed." if all_pass else "SOME CHECKS FAILED.")


if __name__ == "__main__":
    D = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    verify_bijection(D)
