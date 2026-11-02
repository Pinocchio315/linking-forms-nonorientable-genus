#!/usr/bin/env python
"""The linking forms of rank1_linking.py recomputed from the KnotInfo PD codes ([C] Section 3.5, Linking forms).

Each checkerboard (Goeritz) matrix G of a diagram presents H_1(Sigma_K) and its linking form up to a global sign, as
V + V^T does.  For every knot of data/rank1_targets.json and data/rank1_controls.json the script compares the isometry
classes of the forms of both Goeritz matrices of the KnotInfo diagram with that of V + V^T, allowing one sign for all
primes; the class records, for each odd p and level v, the rank and the square class of the discriminant of the Jordan
block (Jordan data of rank1_linking.invariants).  For the obstructed knots whose Jordan blocks all have level 1 it also
records which standard basis vectors generate coker G; there are such only when det K is square-free, the case decided
by Algorithm 1 of Lee and Sabloff.
goeritz_pair() traces the faces of a PD code through its rotation system and allows crossings of both signs; the other
scripts take their Goeritz matrices from it.
Usage: python scripts/rank1_goeritz.py [workers]
Output: data/rank1_goeritz.json
"""
import json, os, sys
from collections import defaultdict
from fractions import Fraction
from math import lcm
from multiprocessing import Pool
from rank1_linking import invariants, legendre
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')


def goeritz_pair(pd):
    """The two checkerboard matrices of a PD code (KnotTheory convention), with one face deleted from each."""
    darts = defaultdict(list)
    for c, row in enumerate(pd):
        for j, label in enumerate(row):
            darts[label].append((c, j))
    assert all(len(x) == 2 for x in darts.values())
    twin = {a: b for ab in darts.values() for a, b in (ab, ab[::-1])}
    remaining, face_of, nf = set(twin), {}, 0
    while remaining:                      # the face of (c, j) contains the corner between darts j and j+1 of c
        here = start = min(remaining)
        while True:
            remaining.remove(here); face_of[here] = nf
            c, j = here; here = twin[c, (j + 1) % 4]
            if here == start:
                break
        nf += 1
    assert nf == len(pd) + 2
    edges = [(face_of[c, j], face_of[c, (j + 2) % 4], sign) for c in range(len(pd)) for j, sign in ((0, 1), (1, -1))]
    colour = {0: 0}
    todo = [0]
    while todo:                           # faces opposite at a crossing get the same colour
        a = todo.pop()
        for x, y, _ in edges:
            for u, v in ((x, y), (y, x)):
                if u == a and v not in colour:
                    colour[v] = 0; todo.append(v)
    for f in range(nf):
        colour.setdefault(f, 1)
    out = []
    for col in (0, 1):
        comp = [f for f in range(nf) if colour[f] == col]
        ix = {v: i for i, v in enumerate(comp)}
        L = [[0] * len(comp) for _ in comp]
        for a, b, sign in edges:
            if a in ix and a != b:
                i, j = ix[a], ix[b]
                L[i][j] += sign; L[j][i] += sign; L[i][i] -= sign; L[j][j] -= sign
        out.append([row[1:] for row in L[1:]])
    return out


def det(A):
    n, M, d = len(A), [[Fraction(x) for x in row] for row in A], Fraction(1)
    for k in range(n):
        p = next((r for r in range(k, n) if M[r][k] != 0), None)
        if p is None:
            return 0
        if p != k:
            M[k], M[p] = M[p], M[k]; d = -d
        d *= M[k][k]
        for r in range(k + 1, n):
            f = M[r][k] / M[k][k]
            for c in range(k, n):
                M[r][c] -= f * M[k][c]
    return int(d)


def inverse(A):
    n = len(A)
    M = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for k in range(n):
        p = next(r for r in range(k, n) if M[r][k] != 0)
        M[k], M[p] = M[p], M[k]
        piv = M[k][k]
        M[k] = [x / piv for x in M[k]]
        for r in range(n):
            if r != k and M[r][k] != 0:
                f = M[r][k]
                M[r] = [x - f * y for x, y in zip(M[r], M[k])]
    return [row[n:] for row in M]


def iso_class(inv, eta):
    """For each p and level v: (v, rank, Legendre symbol of the discriminant of eta times the block)."""
    out = {}
    for p, J in inv.items():
        levels = defaultdict(list)
        for v, u in J:
            levels[v].append(u)
        cls = []
        for v in sorted(levels):
            d = 1
            for u in levels[v]:
                d = d * u % p
            cls.append((v, len(levels[v]), legendre(d * eta ** len(levels[v]), p)))
        out[p] = cls
    return out


def one(item):
    name, pd, V, rec = item
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    D = rec['det']
    ref = iso_class(invariants(A, D), 1)
    checks = []
    for G in goeritz_pair(pd):
        dG = abs(det(G))
        if dG == D:
            invG = invariants(G, D)
            eta = [e for e in (1, -1) if iso_class(invG, e) == ref]
        else:
            eta = []
        c = dict(size=len(G), det=dG, eta=eta)
        if rec['obstructed'] and all(v == 1 for J in rec['jordan'].values() for v, _ in J):
            Gi = inverse(G)                # order of e_i in coker G = lcm of the denominators of column i of G^-1
            c['generating_basis_vectors'] = [i for i in range(len(G)) if lcm(*(Gi[r][i].denominator for r in range(len(G)))) == D]
        checks.append(c)
    return name, dict(det=D, obstructed=rec['obstructed'], goeritz=checks,
                      agree=all(c['det'] == D and c['eta'] for c in checks))


if __name__ == '__main__':
    from collections import Counter
    from database_knotinfo import link_list
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    BY = {r['name']: r for r in link_list()[1:]}
    recs = {}
    for mode in ('targets', 'controls'):
        recs.update(json.load(open(os.path.join(ROOT, 'data', f'rank1_{mode}.json'))))
    items = [(n, json.loads(BY[n]['pd_notation']), json.loads(BY[n]['seifert_matrix']), r) for n, r in recs.items()]
    with Pool(workers) as p:
        out = dict(p.imap_unordered(one, items, chunksize=10))
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', 'rank1_goeritz.json'), 'w'), indent=0)
    print('knots', len(out), '| both Goeritz matrices isometric to +-(V + V^T):', sum(r['agree'] for r in out.values()),
          '| failures', [n for n, r in out.items() if not r['agree']][:10])
    print('sign pattern (eta of the two colourings):', Counter(tuple(tuple(c['eta']) for c in r['goeritz']) for r in out.values()))
    for n, r in sorted(out.items()):
        if any('generating_basis_vectors' in c for c in r['goeritz']):
            print(n, r['det'], [(c['size'], c['generating_basis_vectors']) for c in r['goeritz']])
