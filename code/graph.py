#!/usr/bin/env python3
"""
graph.py — G_D construction and core analysis.

Builds the Collatz residue-refinement graph G_D (§2 of the paper) and
verifies the SCC theorem (§4): G_D has a unique nontrivial strongly
connected component C_D of size D(D-2), and every non-boundary vertex
reaches it.

Usage:
    python graph.py [D_max]          # default D_max = 18
"""

from __future__ import annotations

import sys
from collections import deque, defaultdict
from dataclasses import dataclass
from math import comb
from typing import Dict, List, Optional, Tuple

Vertex = Tuple[int, int]  # (r, d)


# ---------------------------------------------------------------------------
# Arithmetic primitives
# ---------------------------------------------------------------------------

def v2(n: int) -> int:
    """2-adic valuation of n (n > 0)."""
    a = 0
    while n % 2 == 0:
        n //= 2
        a += 1
    return a


def exceptional_residue(d: int, t: int) -> int:
    """Unique unresolved residue of type t at depth d (1 <= t < d).

    r_{d,t} = 2^t * (3^t)^{-1} - 1  mod 2^d.
    """
    mod = 1 << (d - t)
    inv = pow(pow(3, t), -1, mod)
    return ((1 << t) * inv - 1) % (1 << d)


def exceptional_residues(d: int) -> Dict[int, int]:
    """Return {type: residue} for all d unresolved residues at depth d.

    Types 1..d-1 are non-maximal; type d is the maximal residue 2^d - 1.
    """
    out = {t: exceptional_residue(d, t) for t in range(1, d)}
    out[d] = (1 << d) - 1
    return out


def forward_edge(r: int, d: int) -> Tuple[int, int, int, int]:
    """Compute the valuation edge from a resolved vertex (r, d).

    Returns (s, e, t, rho) where:
        t   = v2(r + 1)
        rho = v2(3^t * q - 1),  q = (r+1) / 2^t
        e   = d - t - rho
        s   = ((3^t q - 1) / 2^rho) mod 2^e
    """
    t = v2(r + 1)
    q = (r + 1) >> t
    y = pow(3, t) * q - 1
    rho = v2(y)
    e = d - t - rho
    s = (y >> rho) % (1 << e)
    return s, e, t, rho


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

@dataclass
class Edge:
    src: Vertex
    dst: Vertex
    kind: str            # "refine" or "forward"
    t: Optional[int] = None
    rho: Optional[int] = None


@dataclass
class Graph:
    D: int
    vertices: List[Vertex]
    edges: List[Edge]
    adj: Dict[Vertex, List[Vertex]]
    edge_data: Dict[Tuple[Vertex, Vertex], Edge]
    exc_type: Dict[Vertex, int]           # vertex -> exceptional type (or absent)
    exc_by_type: Dict[int, Dict[int, int]] # depth -> {type: residue}


def build(D: int) -> Graph:
    """Construct G_D."""
    exc_by_type = {d: exceptional_residues(d) for d in range(1, D + 1)}
    exc_type: Dict[Vertex, int] = {}
    for d, tr in exc_by_type.items():
        for t, r in tr.items():
            exc_type[(r, d)] = t

    vertices: List[Vertex] = []
    edges: List[Edge] = []
    adj: Dict[Vertex, List[Vertex]] = {}
    edge_data: Dict[Tuple[Vertex, Vertex], Edge] = {}

    for d in range(1, D + 1):
        for r in range(1, 1 << d, 2):
            u = (r, d)
            vertices.append(u)
            t_type = exc_type.get(u)

            if t_type is not None:
                # Exceptional: refine to depth d+1 (unless at boundary)
                if d < D:
                    nbrs = []
                    for bit in (0, 1):
                        v = (r + bit * (1 << d), d + 1)
                        e = Edge(src=u, dst=v, kind="refine")
                        edges.append(e)
                        edge_data[(u, v)] = e
                        nbrs.append(v)
                    adj[u] = nbrs
                else:
                    adj[u] = []
            else:
                s, e, t, rho = forward_edge(r, d)
                v = (s, e)
                edge = Edge(src=u, dst=v, kind="forward", t=t, rho=rho)
                edges.append(edge)
                edge_data[(u, v)] = edge
                adj[u] = [v]

    return Graph(D=D, vertices=vertices, edges=edges, adj=adj,
                 edge_data=edge_data, exc_type=exc_type, exc_by_type=exc_by_type)


# ---------------------------------------------------------------------------
# SCC (iterative Tarjan)
# ---------------------------------------------------------------------------

def tarjan_scc(vertices: List[Vertex], adj: Dict[Vertex, List[Vertex]]) -> List[List[Vertex]]:
    index = 0
    stack: List[Vertex] = []
    on_stack: set = set()
    idx: Dict[Vertex, int] = {}
    low: Dict[Vertex, int] = {}
    comps: List[List[Vertex]] = []

    for root in vertices:
        if root in idx:
            continue
        work = [(root, iter(adj[root]))]
        idx[root] = low[root] = index
        index += 1
        stack.append(root)
        on_stack.add(root)

        while work:
            v, nbrs = work[-1]
            advanced = False
            for w in nbrs:
                if w not in idx:
                    idx[w] = low[w] = index
                    index += 1
                    stack.append(w)
                    on_stack.add(w)
                    work.append((w, iter(adj[w])))
                    advanced = True
                    break
                elif w in on_stack:
                    low[v] = min(low[v], idx[w])
            if not advanced:
                work.pop()
                if work:
                    low[work[-1][0]] = min(low[work[-1][0]], low[v])
                if low[v] == idx[v]:
                    comp: List[Vertex] = []
                    while True:
                        w = stack.pop()
                        on_stack.remove(w)
                        comp.append(w)
                        if w == v:
                            break
                    comps.append(comp)
    return comps


def bfs_distance_to(vertices: List[Vertex], adj: Dict[Vertex, List[Vertex]],
                    target: set) -> Dict[Vertex, int]:
    """Shortest distance from each vertex to the target set (reverse BFS)."""
    radj: Dict[Vertex, List[Vertex]] = {v: [] for v in vertices}
    for u in vertices:
        for v in adj[u]:
            radj[v].append(u)

    dist = {v: 0 for v in target}
    q: deque[Vertex] = deque(target)
    while q:
        v = q.popleft()
        for u in radj[v]:
            if u not in dist:
                dist[u] = dist[v] + 1
                q.append(u)
    return dist


def boundary_set(g: Graph) -> set:
    """The D+1 boundary-trapped vertices (§4)."""
    D = g.D
    out = {(r, D) for r in g.exc_by_type[D].values()}
    if D >= 2:
        out.add(((1 << (D - 1)) - 1, D - 1))
    return out


# ---------------------------------------------------------------------------
# Main: verify the SCC theorem
# ---------------------------------------------------------------------------

def verify_scc_theorem(D_max: int = 18):
    """
    SCC Theorem (§4): For D >= 3, G_D has a unique nontrivial SCC C_D
    of size D(D-2).  A vertex reaches C_D iff it is not in the boundary.
    """
    print("Verifying SCC Theorem: |C_D| = D(D-2), unique nontrivial SCC")
    print(f"{'D':>3}  {'|C_D|':>8}  {'D(D-2)':>8}  {'nontrivial SCCs':>16}  "
          f"{'boundary':>10}  {'all reach C_D':>14}  {'result':>6}")
    print("-" * 80)

    all_pass = True
    for D in range(3, D_max + 1):
        g = build(D)
        comps = tarjan_scc(g.vertices, g.adj)
        core_comp = next(c for c in comps if (1, 1) in c)
        C = set(core_comp)
        dist = bfs_distance_to(g.vertices, g.adj, C)
        bd = boundary_set(g)

        nontrivial = [c for c in comps if len(c) > 1]
        predicted = D * (D - 2)

        # Every non-boundary, non-core vertex should have finite distance to C
        all_reach = all(v in dist for v in g.vertices if v not in bd)

        ok = (len(C) == predicted and len(nontrivial) == 1 and all_reach)
        all_pass = all_pass and ok

        print(f"{D:>3}  {len(C):>8}  {predicted:>8}  {len(nontrivial):>16}  "
              f"{len(bd):>10}  {str(all_reach):>14}  {'PASS' if ok else 'FAIL':>6}")

    print()
    print("All D passed." if all_pass else "SOME CHECKS FAILED.")


if __name__ == "__main__":
    D_max = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    verify_scc_theorem(D_max)
