#!/usr/bin/env python
"""Check of h_lambda >= 3 for the knots with 14 crossings of the table of Lee and Sabloff by binary forms ([C]
Theorem D).

For every knot with h_witt >= 2 in data/rank1_ls14.json, the Jordan data of a Goeritz matrix of the SnapPy diagram are
tested by the binary-form enumeration of rank2.linking_signs, with the signature rule of rank2.py (s = sigma + 4 Arf,
sigma mod 8 from the excesses of lambda by [T] Lemma 2.5, Arf from det mod 8).  The knots with no odd binary Q
(field obstructed_odd) are compared with those with h_witt = 3 ([T] Corollary 4.5(c)).
Usage: python scripts/rank2_ls14.py [workers]
Output: data/rank2_ls14.json
"""
import json, os, sys
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
from rank1_goeritz import goeritz_pair, det
from rank1_linking import invariants
from linking_genus import excess
from rank2 import linking_signs, ALLOWED


def one(name):
    import snappy
    G1, _ = goeritz_pair([list(c) for c in snappy.Link(name).PD_code()])
    D = abs(det(G1))
    inv = invariants(G1, D)
    s = ((-sum(excess(J, p) for p, J in inv.items())) % 8 + 4 * (0 if D % 8 in (1, 7) else 1)) % 8
    ok = linking_signs(inv, D)
    return name, dict(det=D, s=s, obstructed_odd=not any(par == 'odd' for S, par in ok),
                      obstructed_odd_signed=not any(par == 'odd' and S in ALLOWED[s] for S, par in ok))


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    L14 = json.load(open(os.path.join(ROOT, 'data', 'rank1_ls14.json')))
    names = sorted(k for k, r in L14.items() if r['h_witt'] is not None and r['h_witt'] >= 2)
    with Pool(workers) as p:
        out = dict(p.imap_unordered(one, names, chunksize=20))
    for k in out:
        out[k]['h_witt'] = L14[k]['h_witt']
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', 'rank2_ls14.json'), 'w'), indent=0)
    print('knots with h_witt >= 2:', len(out), '| (h_witt >= 3, no odd Q, no odd Q with allowed sign):',
          sorted(Counter((r['h_witt'] >= 3, r['obstructed_odd'], r['obstructed_odd_signed']) for r in out.values()).items()))
    h3 = [k for k, r in out.items() if r['h_witt'] == 3]
    print('h_witt = 3:', len(h3), '| alternating', sum(k[3] == 'a' for k in h3), '| nonalternating', sum(k[3] == 'n' for k in h3))
