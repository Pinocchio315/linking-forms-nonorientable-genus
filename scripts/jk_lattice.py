#!/usr/bin/env python
"""The lattice criterion of Jabuka and Kelly (arXiv:1708.03000, Thms. 2.10 and 2.12) for alternating knots, compared
with the linking form and the correction terms of rank1_dinv.py ([C] Sections 2.3 and 3.5, Donaldson's theorem).

Thm. 2.10 has three hypotheses: sigma + 4 Arf = 2 eps (mod 8), G is a -eps-definite Goeritz form, and G + [-l] embeds
in (Z^{rank G + 1}, eps Id) for no divisor l of det K with det K / l a square; then gamma_4 >= 2.  In this statement
the diagonal form has the sign eps, in which no (-eps)-definite form embeds; the argument gives the sign of [C]
Proposition 2.7, with a positive definite lattice L, chosen by the sign rules below, in place of G and the standard
lattice Z^{rank L + 1}.  Thm. 2.12 (sigma + 4 Arf = 0) is the same with G+ + [l] and G- + [-l] both failing.
If the definite lattice L embeds in Z^{n+1} at all, let u be a primitive vector spanning its orthogonal complement and m
the index of L in its saturation.  Then |m u|^2 = det K, so l = det K works.  The test is therefore whether L embeds in
Z^{rank L + 1}.
Signs, with G_neg and G_pos the negative and positive definite forms of rank1_goeritz.goeritz_pair:
- s = sigma + 4 Arf = 2 needs L = G_pos, s = 6 needs L = -G_neg, and s = 0 needs one of the two (obstructed_signed);
- the unsigned version accepts either lattice in every case (obstructed).
This is the convention of rank1_dinv.py, where the boundary of X_{G_neg} is -Sigma_K and W is positive definite when
s = 2.  It is checked on the controls.
Embedding search (embeds): rows are placed in a most-constrained order.  Each new row is generated coordinate by
coordinate with Cauchy-Schwarz pruning.  It is canonical under the stabiliser of the rows already placed (permutations
of equal columns, sign changes of zero columns), so the search is complete.
Modes: controls, the alternating knots with KnotInfo value 1; targets, the alternating knots whose KnotInfo range
starts at 1, including those that the linking form obstructs; ls14, the alternating knots with 14 crossings that the
table of Lee and Sabloff leaves open (no Arf obstruction, no band, linking form not obstructed), unsigned test only.
Usage: python scripts/jk_lattice.py controls | targets [workers]
       python scripts/jk_lattice.py ls14 path/to/knot_3-14data.csv [workers]
Output: data/jk_<mode>.json
"""
import json, os, sys
from math import isqrt
from multiprocessing import Pool
from rank1_goeritz import goeritz_pair
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, os.path.join(ROOT, 'vendor', 'code', 'lower_bounds', 'owens'))
from owens_obstruction import positive_definite


def order_rows(A):
    """Most-constrained order: start at the largest norm, then the row with most placed neighbours."""
    n = len(A)
    out = [max(range(n), key=lambda i: A[i][i])]
    while len(out) < n:
        rest = [i for i in range(n) if i not in out]
        out.append(max(rest, key=lambda i: (sum(A[i][j] != 0 for j in out), A[i][i])))
    return out


def embeds(A, N=None):
    """Is there an integer matrix V (n x N) with V V^T = A?  N defaults to n + 1."""
    n = len(A); N = n + 1 if N is None else N
    order = order_rows(A)
    rows = []

    def extend(k):
        if k == n:
            return True
        i = order[k]
        a = A[i][i]
        targets = [A[i][order[j]] for j in range(k)]
        cols = sorted(range(N), key=lambda c: tuple(r[c] for r in rows))
        key = [tuple(r[c] for r in rows) for c in cols]
        suffix = [[0] * (N + 1) for _ in range(k)]          # remaining norm of each placed row on cols[p:]
        for j in range(k):
            for p in range(N - 1, -1, -1):
                suffix[j][p] = suffix[j][p + 1] + rows[j][cols[p]] ** 2
        v = [0] * N

        def coords(p, rem, part):
            if p == N:
                if rem == 0 and part == targets:
                    rows.append(list(v))
                    if extend(k + 1):
                        return True
                    rows.pop()
                return False
            c = cols[p]
            hi = isqrt(rem)
            if p and key[p] == key[p - 1]:
                hi = min(hi, v[cols[p - 1]])                   # equal columns: non-increasing entries
            lo = 0 if not any(key[p]) else -isqrt(rem)          # zero columns: nonnegative entries
            for x in range(hi, lo - 1, -1):
                r2 = rem - x * x
                q = [t + x * rows[j][c] for j, t in enumerate(part)]
                if all((targets[j] - q[j]) ** 2 <= suffix[j][p + 1] * r2 for j in range(k)):
                    v[c] = x
                    if coords(p + 1, r2, q):
                        return True
            v[c] = 0
            return False

        return coords(0, a, [0] * k)

    return extend(0)


def one(item):
    name, pd, s = item
    Gs = goeritz_pair(pd)
    neg = [G for G in Gs if positive_definite([[-v for v in row] for row in G])]
    pos = [G for G in Gs if positive_definite(G)]
    if len(neg) != 1 or len(pos) != 1:
        return name, dict(definite=False)
    L1 = [[-v for v in row] for row in neg[0]]
    L2 = pos[0]
    e1, e2 = embeds(L1), embeds(L2)
    need = {2: (e2,), 6: (e1,)}.get(s, (e1, e2))
    return name, dict(definite=True, s=s, rank=(len(L1), len(L2)), embeds=(e1, e2),
                      obstructed=not (e1 or e2), obstructed_signed=not any(need))


if __name__ == '__main__':
    from collections import Counter
    from database_knotinfo import link_list
    mode = sys.argv[1]
    if mode == 'ls14':                    # 14-crossing alternating knots left open by Lee-Sabloff (no s: unsigned only)
        import csv, re, snappy
        path, workers = sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 8
        L14 = json.load(open(os.path.join(ROOT, 'data', 'rank1_ls14.json')))
        R = {n: dict(obstructed=r['obstructed']) for n, r in L14.items()}
        names = [r['Name'] for r in csv.DictReader(open(path, newline='')) if re.fullmatch(r'K14a\d+', r['Name'])
                 and r['Arf obstructed'] == 'False' and not r['band'] and r['Linking Form'] != 'Obstructed']
        items = [(n, [list(c) for c in snappy.Link(n).PD_code()], None) for n in names]
    else:
        workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
        BY = {r['name']: r for r in link_list()[1:]}
        G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))
        R = {}
        for f in ('rank1_targets.json', 'rank1_controls.json'):
            R.update(json.load(open(os.path.join(ROOT, 'data', f))))
        if mode == 'controls':
            names = [n for n, v in G4.items() if v == '1' and BY[n]['alternating'] == 'Y']
        else:                             # every alternating open knot, including those obstructed by the linking form
            names = [n for n, v in G4.items() if v.startswith('[1,') and BY[n]['alternating'] == 'Y']
        items = [(n, json.loads(BY[n]['pd_notation']), R[n]['s']) for n in names]
    with Pool(workers) as p:
        out = dict(p.imap_unordered(one, items, chunksize=4))
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', f'jk_{mode}.json'), 'w'), indent=0)
    D = json.load(open(os.path.join(ROOT, 'data', f'rank1_dinv_{mode}.json')))
    print(mode, 'knots', len(out), '| definite', sum(r['definite'] for r in out.values()),
          '| JK obstructed', sum(r.get('obstructed', False) for r in out.values()),
          '| JK obstructed with sign', sum(r.get('obstructed_signed', False) for r in out.values()))
    cmp = Counter((('linking' if R[n]['obstructed'] else 'd' if D.get(n, {}).get('obstructed') else 'passes'),
                   'JK signed' if r['obstructed_signed'] else 'JK unsigned' if r['obstructed'] else 'JK passes')
                  for n, r in out.items() if r['definite'])
    for k in sorted(cmp):
        print('  (ours, JK)', k, cmp[k])
