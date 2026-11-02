#!/usr/bin/env python
"""The knots of Examples 3.5 and 3.6 of the theory paper, checked from PD codes (computational paper, Section 2.1).

J = P(3, -3, 1), drawn as in the theory paper's Figure 2(a) with the white regions X_0, X_1, X_2, is the numerator
closure of spherogram's rational tangles -1/3, 1/3 and -1/1: spherogram's tangle 1/n is the mirror image of a twist
region of the theory paper with |n| crossings of incidence number eta = sign n.  There b_n
is the band at the top of the third twist region that joins its two upper strands, with a box labelled n: |n|
crossings between its long sides, all with incidence number eta = sign(n) of Gordon and Litherland.  The band is
inserted with bands.band, whose t counts the crossings with the sign of Figure 1 of the computational
paper; b_n is the band with t = n.  Goeritz matrices are those of the white regions with the incidence numbers of
Gordon and Litherland (goeritz.gl_goeritz) and are compared with those of the theory paper up to the order of the
regions only, so that the sign fixes the chirality.  The script checks:
- J has the complement of 6_1, determinant 9 and the Goeritz matrix [[0, 3], [3, -2]] of the theory paper (basis X_1,
  X_2), and the band move along b_0 gives a diagram that Reidemeister simplification reduces to two unlinked unknots;
- the band move along b_-3 gives a knot K whose Goeritz matrix, with all five white regions, is that of the
  incidences of the theory paper, with H_K = Z/27 and the complement of 8_10; the band move along b_3 gives the
  complement of 10_143 (a knot complement determines the knot up to mirror image), and the one along b_-2 a
  two-component link L;
- in the diagram of L, the band b' in the bigon between the two crossings of b_-2, with two crossings of incidence
  number 1 between its long sides (t = -2), gives a knot K' whose Goeritz matrix is the matrix G' of the proof of
  Example 3.6, [[0, 3, 0, 0], [3, -3, 1, 0], [0, 1, 1, -2], [0, 0, -2, 1]] in the basis X_1, X_2, Y_1, Y_2, with Smith
  normal form (1, 1, 1, 27), and whose complement is that of 8_11.
For each matrix the script also records that minus the matrix of the paper does not occur (the mirror image).  It
records the PD codes of K and K' quoted in the computational paper, in the convention of KnotInfo (edges numbered from 1
along the orientation, each crossing listed counterclockwise from the incoming under-strand), and repeats the checks of
the Goeritz matrix and of the complement on them.
Usage: small_examples.py   Output: data/small_examples.txt
"""
import os, sys, warnings
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from bands import from_pd, faces, band, pd_code, components, smith
from smooth_top_gap import det_link
from goeritz import gl_goeritz, region_order, delete

G_J = [[0, 3], [3, -2]]
G_PRIME = [[0, 3, 0, 0], [3, -3, 1, 0], [0, 1, 1, -2], [0, 0, -2, 1]]


def goeritz_from_incidences(pairs, n):
    """The Goeritz matrix of the theory paper, before a region is deleted: entries -sum(eta) over the crossings between
    two white regions, row sums 0; pairs maps (i, j) to the sum of eta over the crossings between the regions i, j."""
    M = [[0] * n for _ in range(n)]
    for (i, j), eta in pairs.items():
        M[i][j] -= eta; M[j][i] -= eta
    for i in range(n):
        M[i][i] = -sum(M[i][j] for j in range(n) if j != i)
    return M


def restore(G):
    """Add back the deleted region of a Goeritz matrix, with row sums 0."""
    row0 = [-sum(r) for r in G]
    return [[-sum(row0)] + row0] + [[row0[i]] + G[i] for i in range(len(G))]


def compare(pd, F):
    """The Goeritz matrix of the white regions of pd in the order of the regions of F (or as computed), whether it is
    F, and whether it is -F."""
    M = gl_goeritz(pd)
    order = region_order(M, F)
    mirror = region_order(M, [[-x for x in r] for r in F]) is not None
    return ([[M[i][j] for j in order] for i in order] if order else M), order is not None, mirror


def knotinfo_pd(code):
    """The PD code of the same diagram in the convention of KnotInfo: spherogram's, with the edges numbered from 1."""
    import spherogram
    return [tuple(e + 1 for e in X) for X in spherogram.Link([list(x) for x in code]).PD_code()]


def complement_is(L, name):
    import snappy
    M = L.exterior()
    try:
        return M.is_isometric_to(snappy.Manifold(name))
    except RuntimeError:
        return any(name == str(x).split('(')[0] for x in M.identify())


if __name__ == '__main__':
    import snappy, spherogram
    from spherogram import RationalTangle
    lines = [f'Examples 3.5 and 3.6 of the theory paper; snappy {snappy.__version__}, spherogram {spherogram.__version__}']
    J = (RationalTangle(1, -3) + RationalTangle(1, 3) + RationalTangle(1, -1)).numerator_closure()   # P(3, -3, 1)
    pd = [list(x) for x in J.PD_code()]
    third = len(pd) - 1                                  # the crossing of the third twist region
    nbr, under, label = from_pd([tuple(x) for x in pd])
    where = {}
    for c, X in enumerate(pd):
        for e in X:
            where.setdefault(e, set()).add(c)
    # the face above the third twist region: its two edges at that crossing are the two upper strands
    arc = None
    for f in faces(nbr):
        ix = [k for k, h in enumerate(f) if third in where[label[h]]]
        if len(ix) == 2:
            N, U = band(nbr, under, f[ix[0]], f[ix[1]], 0)
            S = spherogram.Link([list(x) for x in pd_code(N, U)]); S.simplify('global')
            if len(components(N)) == 2 and len(S.crossings) == 0:
                arc = (f, ix[0], ix[1])
                break
    assert arc is not None
    f, i, j = arc
    MJ, sameJ, mirrorJ = compare(pd, restore(G_J))
    GJ = delete(MJ, 0)
    c61 = complement_is(J, '6_1')
    ok_J = c61 and det_link(pd) == 9 and sameJ and not mirrorJ
    lines.append(f'J = P(3,-3,1): {len(pd)} crossings, determinant {det_link(pd)}, complement of 6_1: '
                 f'{c61}, Goeritz matrix {GJ} (that of the paper: {sameJ}; minus it: {mirrorJ}); '
                 f'b_0 gives the two-component unlink: True')

    def b(n):
        """PD code of the band move along b_n (t = n) and its number of components."""
        N, U = band(nbr, under, f[i], f[j], n)
        return N, U, pd_code(N, U)

    # K = J + b_-3: regions X_0, X_1, X_2 and the two bigons Y_1, Y_2 of the band, eta(b_-3) = -1
    N3, U3, P3 = b(-3)
    full_K = goeritz_from_incidences({(0, 1): 3, (1, 2): -3, (2, 0): 1, (2, 3): -1, (3, 4): -1, (4, 0): -1}, 5)
    MK, sameK, mirrorK = compare(P3, full_K)
    snf_K = smith(delete(MK, 0))[0]
    K = spherogram.Link([list(x) for x in P3])
    c810 = complement_is(K, '8_10')
    ok_K = (len(components(N3)) == 1 and sameK and not mirrorK and snf_K[-1] == 27 and all(x == 1 for x in snf_K[:-1])
            and c810)
    lines.append(f'b_-3 (t = -3): {len(components(N3))} component(s), Goeritz matrix of the five white regions that '
                 f'of the incidences of the paper: {sameK} (minus it: {mirrorK}), Smith normal form {snf_K}, '
                 f'complement of 8_10: {c810}')
    N3m, U3m, P3m = b(3)
    c10143 = complement_is(spherogram.Link([list(x) for x in P3m]), '10_143')
    N2, U2, P2 = b(-2)
    two = len(components(N2)) == 2
    lines.append(f'b_3 (t = 3): determinant {det_link(P3m)}, complement of 10_143: {c10143}; '
                 f'b_-2 (t = -2): two components: {two}')
    # b' in the bigon between the two crossings of b_-2
    n2, u2, l2 = from_pd([tuple(x) for x in P2])
    new = set(range(len(pd), len(P2)))
    where2 = {}
    for c, X in enumerate(P2):
        for e in X:
            where2.setdefault(e, set()).add(c)
    bigons = [F for F in faces(n2) if len(F) == 2 and all(where2[l2[h]] == new for h in F)]
    assert len(bigons) == 1
    F = bigons[0]
    NN, UU = band(n2, u2, F[0], F[1], -2)
    PK = pd_code(NN, UU)
    Kp = spherogram.Link([list(x) for x in PK])
    c811 = complement_is(Kp, '8_11')
    MKp, same, mirrorKp = compare(PK, restore(G_PRIME))
    G = delete(MKp, 0)
    snf = smith(G)[0]
    lines.append(f"b' in the bigon of b_-2 with 2 crossings (t = -2): {len(components(NN))} component(s), determinant "
                 f"{det_link(PK)}, Goeritz matrix in the basis X_1, X_2, Y_1, Y_2 {G}, the matrix G' of the paper: "
                 f"{same} (minus it: {mirrorKp}); Smith normal form {snf}; complement of 8_11: {c811}")
    # the PD codes quoted in the computational paper, with the checks repeated on them
    quoted = []
    for name, code, F, ref in (('K', P3, full_K, '8_10'), ("K'", PK, restore(G_PRIME), '8_11')):
        pk = knotinfo_pd(code)
        _, s, m = compare([[e - 1 for e in X] for X in pk], F)
        c = complement_is(spherogram.Link([list(X) for X in pk]), ref)
        quoted.append(s and not m and c)
        lines.append(f'PD code of {name} (KnotInfo convention): {pk}; Goeritz matrix that of the paper: {s} '
                     f'(minus it: {m}); complement of {ref}: {c}')
    ok = (ok_J and ok_K and c10143 and two and c811 and same and not mirrorKp and snf == [1, 1, 1, 27]
          and len(components(NN)) == 1 and all(quoted))
    lines.append(f'all checks pass: {ok}')
    print('\n'.join(lines))
    with open(os.path.join(ROOT, 'data', 'small_examples.txt'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
