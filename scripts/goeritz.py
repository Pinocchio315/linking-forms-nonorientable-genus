#!/usr/bin/env python
"""Goeritz matrices with the incidence numbers of Gordon and Litherland; a module imported by small_examples.py.

pd_faces(pd) traces the faces of a PD code, the face on the left of each edge; tait(pd, white_colour) colours them
(face 0 gets colour 0) and returns the faces, the colouring, the white faces and, for each crossing, the white faces
at two opposite corners with eta: +1 if these are the corners 0 and 2 (corner k lies between the slots k and k + 1),
-1 otherwise.  This is the eta of Gordon and Litherland (their Fig. 1): with the white regions to the left and right
of the crossing, eta = +1 exactly when the over strand runs from bottom left to top right.
laplacian(vertices, edges, delete) is the Goeritz matrix of these edges (entries -sum eta, row sums 0) with the vertex
delete removed, gl_goeritz(pd) that of the checkerboard colouring with the fewer white regions with no region deleted,
delete(M, r) removes the row and column r, and region_order(M, F) finds an order of the regions in which M is F.
Usage: imported; no command line.   Output: none
"""
from collections import defaultdict
from itertools import permutations


def slots_of(pd):
    s = defaultdict(list)
    for x, cr in enumerate(pd):
        for i, e in enumerate(cr):
            s[e].append((x, i))
    bad = [e for e, v in s.items() if len(v) != 2]
    assert not bad, f'edges not appearing twice: {bad[:5]}'
    return s


def other_end(S, e, x, i):
    a, b = S[e]
    return b if a == (x, i) else a


def pd_faces(pd):
    """Faces as lists of directed edge traversals (e, (X, i) -> (Y, j)), face on the left: leave Y via slot j - 1."""
    S = slots_of(pd)
    seen, out = set(), []
    for x in range(len(pd)):
        for i in range(4):
            if (x, i) in seen:
                continue
            f, cur = [], (x, i)
            while cur not in seen:
                seen.add(cur)
                X, I = cur
                e = pd[X][I]
                Y, J = other_end(S, e, X, I)
                f.append((e, (X, I), (Y, J)))
                cur = (Y, (J - 1) % 4)
            out.append(f)
    return out


def tait(pd, white_colour=0):
    """Returns (faces, colour, white list, edges) with edges = list of (crossing, white face a, white face b, eta).
    Corner k of crossing x lies between the slots k and k + 1; eta = +1 if the white corners are {0, 2}, else -1."""
    F = pd_faces(pd)
    face_of_corner = {}
    for fi, f in enumerate(F):
        for e, (X, I), (Y, J) in f:
            face_of_corner[(Y, (J - 1) % 4)] = fi
    adj = defaultdict(set); side = defaultdict(list)
    for fi, f in enumerate(F):
        for e, a, b in f:
            side[e].append(fi)
    for e, fs in side.items():
        adj[fs[0]].add(fs[1]); adj[fs[1]].add(fs[0])
    col = {0: 0}; q = [0]
    while q:
        u = q.pop()
        for w in adj[u]:
            if w not in col:
                col[w] = 1 - col[u]; q.append(w)
            else:
                assert col[w] != col[u]
    white = sorted(fi for fi in col if col[fi] == white_colour)
    edges = []
    for x in range(len(pd)):
        cf = {k: face_of_corner[(x, k)] for k in range(4)}
        wc = [k for k in range(4) if col[cf[k]] == white_colour]
        assert len(wc) == 2 and wc[1] - wc[0] == 2
        eta = +1 if wc[0] == 0 else -1     # Gordon-Litherland: white corners {0, 2}
        edges.append((x, cf[wc[0]], cf[wc[1]], eta))
    return F, col, white, edges


def laplacian(vertices, edges, delete):
    idx = {v: k for k, v in enumerate(v for v in vertices if v != delete)}
    n = len(idx)
    G = [[0] * n for _ in range(n)]
    for x, a, b, eta in edges:
        if a == b:
            continue
        for u, v in ((a, b), (b, a)):
            if u != delete:
                G[idx[u]][idx[u]] += eta
                if v != delete:
                    G[idx[u]][idx[v]] -= eta
    return G, idx


def delete(M, r):
    return [[M[i][j] for j in range(len(M)) if j != r] for i in range(len(M)) if i != r]


def gl_goeritz(pd):
    """The Goeritz matrix, with the incidence numbers of Gordon and Litherland, of the checkerboard colouring with the
    fewer white regions, with no region deleted."""
    colourings = [tait([list(x) for x in pd], w) for w in (0, 1)]
    _, _, white, edges = min(colourings, key=lambda c: len(c[2]))
    return laplacian(white, edges, None)[0]


def region_order(M, F):
    """An order of the regions of M in which M is F, or None."""
    n = len(F)
    for perm in permutations(range(n)):
        if all(M[perm[i]][perm[j]] == F[i][j] for i in range(n) for j in range(n)):
            return perm
    return None
