#!/usr/bin/env python
"""Routines on PD codes for band moves, and the Smith normal form; a module imported by the band searches and checks.

smith(A) and invariants(A) give the Smith normal form of an integer matrix and the invariant factors of its cokernel.
from_pd(pd) turns a PD code into the half-edge structure described below; faces lists the faces of a diagram,
components and graph_components count the components of the link and of the diagram, band(nbr, under, h1, h3, t) is
the band move between the half-edges h1 and h3 of a face with t crossings between the long sides of the band (for
t > 0 the strand from lower left to upper right is over, as in Figure 1 of the computational paper), pd_code turns
the result back into a PD code, and band_pd(pd, face, edges, t) is the band move on a PD code.
Usage: imported; no command line.   Output: none
"""
from collections import defaultdict


def smith(A):
    """Smith normal form over Z: (d, U, V) with U A V = diag(d), U and V unimodular, d_i | d_(i+1), d_i >= 0."""
    m = len(A); n = len(A[0]) if m else 0
    M = [list(r) for r in A]
    U = [[int(i == j) for j in range(m)] for i in range(m)]
    V = [[int(i == j) for j in range(n)] for i in range(n)]

    def rowop(i, j, q):                   # row_i -= q row_j
        M[i] = [a - q * b for a, b in zip(M[i], M[j])]; U[i] = [a - q * b for a, b in zip(U[i], U[j])]

    def colop(i, j, q):                   # col_i -= q col_j
        for R in M + V:
            R[i] -= q * R[j]

    def rowswap(i, j):
        M[i], M[j] = M[j], M[i]; U[i], U[j] = U[j], U[i]

    def colswap(i, j):
        for R in M + V:
            R[i], R[j] = R[j], R[i]

    for t in range(min(m, n)):
        nz = [(abs(M[i][j]), i, j) for i in range(t, m) for j in range(t, n) if M[i][j]]
        if not nz:
            break
        _, i, j = min(nz)
        rowswap(t, i); colswap(t, j)
        while True:
            p = M[t][t]
            r = next((i for i in range(t + 1, m) if M[i][t] % p), None)
            if r is not None:
                rowop(r, t, M[r][t] // p); rowswap(t, r); continue
            for i in range(t + 1, m):
                rowop(i, t, M[i][t] // p)
            c = next((j for j in range(t + 1, n) if M[t][j] % p), None)
            if c is not None:
                colop(c, t, M[t][c] // p); colswap(t, c); continue
            for j in range(t + 1, n):
                colop(j, t, M[t][j] // p)
            bad = next(((i, j) for i in range(t + 1, m) for j in range(t + 1, n) if M[i][j] % p), None)
            if bad is None:
                break
            rowop(t, bad[0], -1)
        if M[t][t] < 0:
            M[t] = [-a for a in M[t]]; U[t] = [-a for a in U[t]]
    d = [M[i][i] for i in range(min(m, n))] + [0] * (n - min(m, n))
    UA = [[sum(U[i][a] * A[a][b] for a in range(m)) for b in range(n)] for i in range(m)]
    assert all(sum(UA[i][b] * V[b][j] for b in range(n)) == (d[i] if i == j else 0) for i in range(m) for j in range(n))
    return d, U, V


def invariants(A):
    """Invariant factors > 1 of coker A (0 for a free summand)."""
    return tuple(x for x in smith(A)[0] if x != 1)


def from_pd(pd):
    pos = defaultdict(list)
    for c, X in enumerate(pd):
        for j, lab in enumerate(X):
            pos[lab].append((c, j))
    nbr, label = {}, {}
    for lab, (h, k) in pos.items():
        nbr[h] = k; nbr[k] = h; label[h] = label[k] = lab
    return nbr, [0] * len(pd), label


def faces(nbr):
    """Faces as lists of departing half-edges, face on the left: after arriving at (c, i), leave by (c, i - 1)."""
    seen, out = set(), []
    for h0 in sorted(nbr):
        if h0 in seen:
            continue
        face, h = [], h0
        while h not in seen:
            seen.add(h); face.append(h)
            c, i = nbr[h]
            h = (c, (i - 1) % 4)
        assert h == h0
        out.append(face)
    return out


def components(nbr):
    seen, comps = set(), []
    for h0 in sorted(nbr):
        if h0 in seen:
            continue
        comp, h = [], h0
        while h not in seen:
            a = nbr[h]; seen.add(h); seen.add(a); comp.append(h)
            h = (a[0], (a[1] + 2) % 4)
        comps.append(comp)
    return comps


def graph_components(nbr, ncross):
    parent = list(range(ncross))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (c, _), (d, _) in nbr.items():
        parent[find(c)] = find(d)
    return len({find(c) for c in range(ncross)})


def band(nbr, under, h1, h3, t):
    """Band move along an arc in the face containing the darts h1 -> nbr[h1] and h3 -> nbr[h3] (see the docstring)."""
    h2, h4 = nbr[h1], nbr[h3]
    assert len({h1, h2, h3, h4}) == 4
    N, U = dict(nbr), list(under)
    for h in (h1, h2, h3, h4):
        del N[h]

    def link(a, b):
        N[a] = b; N[b] = a
    if t == 0:
        link(h1, h4); link(h2, h3)
    else:                                 # new crossings: slots 0 NE, 1 NW, 2 SW, 3 SE
        cs = list(range(len(U), len(U) + abs(t)))
        U += [1 if t > 0 else 0] * abs(t)
        link(h1, (cs[0], 1)); link(h2, (cs[0], 2))
        for a, b in zip(cs, cs[1:]):
            link((a, 0), (b, 1)); link((a, 3), (b, 2))
        link((cs[-1], 0), h4); link((cs[-1], 3), h3)
    V, E, F = len(U), len(N) // 2, len(faces(N))
    assert V - E + F == 2 * graph_components(N, V)      # every component of the rotation system is spherical
    return N, U


def pd_code(nbr, under):
    """KnotTheory PD code (incoming under-strand first, counterclockwise), edges numbered along the components."""
    lab, arriving, k = {}, set(), 0
    for comp in components(nbr):
        for h in comp:
            k += 1; lab[h] = lab[nbr[h]] = k; arriving.add(nbr[h])
    pd = []
    for c in range(len(under)):
        s0 = under[c] if (c, under[c]) in arriving else under[c] + 2
        assert (c, s0) in arriving and (c, (s0 + 2) % 4) not in arriving
        pd.append(tuple(lab[(c, (s0 + j) % 4)] for j in range(4)))
    return pd


def band_pd(pd, face, edges, t):
    nbr, under, label = from_pd(pd)
    f = faces(nbr)[face]
    h1, h3 = (next(h for h in f if label[h] == e) for e in edges)
    return pd_code(*band(nbr, under, h1, h3, t))
