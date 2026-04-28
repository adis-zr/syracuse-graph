#!/usr/bin/env python3
"""
weights.py — Verify the weight factorization theorem (Theorem 7.1).

For a vertex v with distance k and depth d, the total chain weight

    Omega(v) = sum_{m=1}^{k} (t_m * alpha - (t_m + rho_m))

equals the history weight

    Omega(Sigma) = Lambda(Sigma) * alpha - (d - sigma_2)

where:
    alpha   = log2(3)
    Sigma   = Sigma(v) = {sigma_1 < ... < sigma_{2k+2} = d}
    Lambda  = sum of odd-indexed gaps: (sigma_3-sigma_2) + (sigma_5-sigma_4) + ...
    sigma_2 = second-smallest element of Sigma (seed depth)

This means the chain weight depends only on the history subset, not on
the actual residue values traversed.

Usage:
    python weights.py [D]          # default D = 16
"""

from __future__ import annotations

import sys
from math import log2
from typing import Dict, List, Optional, Tuple, FrozenSet

from graph import Vertex, build, tarjan_scc, bfs_distance_to, v2
from bijection import history, seed_and_orientation

ALPHA = log2(3)


# ---------------------------------------------------------------------------
# Weight from chain steps
# ---------------------------------------------------------------------------

def chain_weight(v: Vertex, g, C: set) -> Optional[float]:
    """Compute Omega(v) = sum_m (t_m * alpha - (t_m + rho_m)) directly."""
    if v in C:
        return None
    cur = v
    seen: set = set()
    omega = 0.0
    while cur not in C:
        if cur in seen:
            return None
        seen.add(cur)
        nbrs = g.adj.get(cur, [])
        if not nbrs:
            return None
        nxt = nbrs[0]
        edge = g.edge_data.get((cur, nxt))
        if edge is None or edge.kind != "forward":
            return None
        t, rho = edge.t, edge.rho
        omega += t * ALPHA - (t + rho)
        cur = nxt
    return omega


# ---------------------------------------------------------------------------
# Weight from history subset
# ---------------------------------------------------------------------------

def subset_weight(sigma: FrozenSet[int]) -> float:
    """Compute Omega(Sigma) = Lambda(Sigma) * alpha - (d - sigma_2).

    sigma_1 < sigma_2 < ... < sigma_{2k+2} = d
    Lambda = (sigma_3 - sigma_2) + (sigma_5 - sigma_4) + ...
           = sum of odd-indexed gaps after the seed depth sigma_2
    c = d - sigma_2
    """
    elems = sorted(sigma)
    k = len(elems) // 2 - 1  # |sigma| = 2k+2
    lam = sum(elems[2 * m] - elems[2 * m - 1] for m in range(1, k + 1))
    c = elems[-1] - elems[1]   # d - sigma_2
    return lam * ALPHA - c


# ---------------------------------------------------------------------------
# Main: verify weight factorization
# ---------------------------------------------------------------------------

def verify_weights(D: int = 16):
    """
    Weight factorization (§7): Omega(v) = Omega(Sigma(v)) for all
    positive-distance vertices v in G_D.
    """
    g = build(D)
    comps = tarjan_scc(g.vertices, g.adj)
    C = set(next(c for c in comps if (1, 1) in c))
    dist = bfs_distance_to(g.vertices, g.adj, C)

    print(f"Verifying Weight Factorization Theorem (D={D})")
    print(f"  Omega(v) from chain steps  ==  Omega(Sigma(v)) from history subset")
    print()

    mismatches = 0
    checked = 0
    max_err = 0.0

    for v in g.vertices:
        if v in C or v not in dist:
            continue

        w_chain = chain_weight(v, g, C)
        coord = history(v, g, C)
        if w_chain is None or coord is None:
            continue

        _, sigma = coord
        w_subset = subset_weight(sigma)

        err = abs(w_chain - w_subset)
        max_err = max(max_err, err)
        checked += 1
        if err > 1e-10:
            mismatches += 1
            print(f"  MISMATCH at v={v}: chain={w_chain:.10f}  subset={w_subset:.10f}  err={err:.2e}")

    print(f"Checked {checked} vertices.")
    print(f"Max numerical error: {max_err:.2e}")
    if mismatches == 0:
        print("All weights match. PASS.")
    else:
        print(f"{mismatches} mismatches. FAIL.")

    # Also verify top-bit independence: tau_d partner has same Omega
    print()
    print("Checking top-bit independence: Omega(v) == Omega(tau_d(v))")
    tau_mismatches = 0
    tau_checked = 0

    for v in g.vertices:
        r, d = v
        if v in C or v not in dist or d < 2:
            continue
        # tau_d partner
        half = 1 << (d - 1)
        r2 = r + half if r < half else r - half
        v2v = (r2, d)
        if v2v not in dist:
            continue

        w1 = chain_weight(v, g, C)
        w2 = chain_weight(v2v, g, C)
        if w1 is None or w2 is None:
            continue

        err = abs(w1 - w2)
        tau_checked += 1
        if err > 1e-10:
            tau_mismatches += 1
            print(f"  TAU MISMATCH at v={v} vs tau={v2v}: {w1:.10f} vs {w2:.10f}")

    print(f"Checked {tau_checked} tau-pairs.")
    if tau_mismatches == 0:
        print("All tau-pairs match. PASS.")
    else:
        print(f"{tau_mismatches} tau-pair mismatches. FAIL.")


if __name__ == "__main__":
    D = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    verify_weights(D)
