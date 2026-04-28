# A finite graph from the 2-adic valuation of the 3n+1 map

**Author:** Aditya Sriram

## Abstract

This paper builds a finite directed graph $G_D$ from the 2-adic valuation data of the 3n+1 map.

The main result is an explicit bijection for the finite-distance layers of this adaptive graph. After quotienting by the top-bit symmetry, the depth-$d$, distance-$k$ layer is in canonical bijection with

$$S_{d,k} = \{\Sigma \subseteq [d] : d \in \Sigma,\ |\Sigma| = 2k+2\}.$$

The same history subset also carries the logarithmic drift weight of the corresponding valuation chain, so the bijection turns the graph into a usable coordinate system rather than only a counting result.

The paper does not prove a new convergence result for the 3n+1 problem. Its contribution is a finite graph-theoretic and enumerative classification of the adaptive 2-adic valuation structure that the odd-to-odd map exposes.

## Files

- `paper/collatz_graph_paper.tex` — LaTeX source (current)
- `paper/history/` — earlier versions v1–v9
- `code/` — Python scripts that verify the paper's main results

## Code

Three scripts reproduce the paper's core theorems. Each is self-contained and requires only the Python standard library (Python 3.8+).

### `code/graph.py` — SCC Theorem (§4)

Builds $G_D$ and verifies that it has a unique nontrivial strongly connected component $C_D$ of size $D(D-2)$, and that every non-boundary vertex reaches $C_D$.

```
python3 code/graph.py [D_max]        # default D_max = 18
```

### `code/bijection.py` — History Bijection (§6)

For every depth $d$ and distance $k$, verifies that the history map $v \mapsto (\varepsilon(v), \Sigma(v))$ is a bijection from the depth-$d$, distance-$k$ layer to $\{+,-\} \times S_{d,k}$, and that the layer size equals $2\binom{d-1}{2k+1}$.

```
python3 code/bijection.py [D]        # default D = 16
```

### `code/weights.py` — Weight Factorization (§7)

Verifies that $\Omega(v) = \Omega(\Sigma(v))$ for every positive-distance vertex: the total log-drift of the chain equals the alternating-gap functional on the history subset, independent of the actual residue values traversed.

```
python3 code/weights.py [D]          # default D = 16
```
