"""Unit test for common.find_bridges (iterative Tarjan on multigraphs) against brute force.

Brute force: an edge is a bridge iff deleting it increases the number of connected components.
Random multigraphs include parallel edges and self-loops (never bridges). Seed fixed.
Usage: python discovery/test_bridges.py
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from common import find_bridges


def n_comp(n, u, v):
    if len(u) == 0:
        return n
    A = coo_matrix((np.ones(len(u)), (u, v)), shape=(n, n))
    return connected_components(A, directed=False)[0]


def main():
    rng = np.random.default_rng(20261012)
    for trial in range(300):
        n = int(rng.integers(2, 60))
        m = int(rng.integers(1, 90))
        u = rng.integers(0, n, m)
        v = rng.integers(0, n, m)
        k = int(rng.integers(0, 6))  # parallel copies of some edges
        u, v = np.r_[u, u[:k]], np.r_[v, v[:k]]
        got = find_bridges(n, u, v)
        base = n_comp(n, u, v)
        want = np.zeros(len(u), bool)
        for i in range(len(u)):
            if u[i] == v[i]:
                continue
            keep = np.arange(len(u)) != i
            want[i] = n_comp(n, u[keep], v[keep]) > base
        assert (got == want).all(), f"mismatch in trial {trial}"
    print("find_bridges: 300 random multigraphs match brute force")


if __name__ == "__main__":
    main()
