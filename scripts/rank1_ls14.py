#!/usr/bin/env python
"""The linking-form tests for the knots with 14 crossings of the table of Lee and Sabloff ([C] Theorem D
and Section 3.7).

For every 14-crossing row of knot_3-14data.csv (github.com/LeeJoPing/Mobius-Band-Search), the knot is built from its
Hoste-Thistlethwaite name with SnapPy.  The tests are run on a Goeritz matrix of that diagram; both checkerboard
matrices must have the same determinant and isometric forms up to one global sign (field consistent).
Two verdicts are stored.  obstructed is the Witt-level test, h_lambda >= 2: condition (L_1) of [T] Theorem 4.4 fails
for both signs ((R_1) always holds), the unsigned test that matches rank1_linking.py.  obstructed_summand is the
orthogonal-sum test of rank1_linking.survivors.  They can differ only at an exceptional prime ([T] Proposition 4.9),
which needs p^9 | det K and does not occur here.  h_witt is h_lambda (linking_genus.h_witt).
The verdicts are compared with the table: knots with a band must pass, knots obstructed by Lee-Sabloff's linking form
must fail, and the open rows (no band, no Arf obstruction, linking form not obstructed) are the new cases.
Usage: python scripts/rank1_ls14.py path/to/knot_3-14data.csv [workers]
Output: data/rank1_ls14.json
"""
import csv, json, os, re, sys
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)          # the CSV is untrusted data: with python -I the local modules are still found
from rank1_linking import invariants, survivors, legendre
from linking_genus import h_witt, excess, anisotropic
from rank1_goeritz import goeritz_pair, det, iso_class


def witt_rank1(inv, D):
    """Unsigned Witt-level rank-one test: a_p <= 1 at every p, and (disc A_p / p) = (tau D'_p / p) wherever a_p = 1, for
    one tau in {1, -1}."""
    A = {p: anisotropic(J, p) for p, J in inv.items()}
    if any(a > 1 for a, _ in A.values()):
        return False
    return any(all(a == 0 or legendre(d, p) == legendre(tau * (D // p ** sum(v for v, _ in inv[p])), p)
                   for p, (a, d) in A.items()) for tau in (1, -1))


def one(item):
    name, ls = item
    import snappy
    L = snappy.Link(name)
    assert len(L.link_components) == 1
    G1, G2 = goeritz_pair([list(c) for c in L.PD_code()])
    D = abs(det(G1))
    rec = dict(ls, det=D, det2=abs(det(G2)))
    if D > 1 and rec['det2'] == D:
        inv = invariants(G1, D)
        c1, c2 = iso_class(inv, 1), iso_class(invariants(G2, D), 1)
        rec['consistent'] = c2 in (c1, iso_class(inv, -1))
        surv = survivors(inv, D)
        witt = witt_rank1(inv, D)
        hw = h_witt(inv, D, (-sum(excess(J, p) for p, J in inv.items())) % 8)
    else:
        rec['consistent'] = rec['det2'] == D
        surv, witt, hw = [(1, 1), (1, -1)], True, (0 if D == 1 else None)
    rec.update(survivors=surv, obstructed=not witt, obstructed_summand=not surv, h_witt=hw)
    return name, rec


if __name__ == '__main__':
    from collections import Counter
    path = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    items = []
    for r in csv.DictReader(open(path, newline='')):
        if re.fullmatch(r'K14[an]\d+', r['Name']):
            items.append((r['Name'], dict(arf=r['Arf obstructed'] == 'True', linking=r['Linking Form'], band=bool(r['band']))))
    with Pool(workers) as p:
        out = dict(p.imap_unordered(one, items, chunksize=50))
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', 'rank1_ls14.json'), 'w'), indent=0)
    print('14-crossing knots', len(out), '| inconsistent Goeritz pairs', [n for n, r in out.items() if not r['consistent']][:10])
    print('Witt-level vs orthogonal-sum verdict (obstructed, obstructed_summand):',
          sorted(Counter((r['obstructed'], r['obstructed_summand']) for r in out.values()).items()))
    print('h_witt:', sorted(Counter(r['h_witt'] for r in out.values()).items(), key=str))
    tab = Counter((('band' if r['band'] else 'Arf' if r['arf'] else 'open'), r['linking'], r['obstructed']) for r in out.values())
    for k in sorted(tab):
        print('  (LS status, LS linking form, ours obstructed)', k, tab[k])
