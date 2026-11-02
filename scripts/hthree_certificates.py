#!/usr/bin/env python
"""Certificates of h_lambda = 3 for the 49 knots with at most 13 crossings ([C] Table 3, by [C] Corollary 2.3(c)).

For each knot with h_witt = 3 in data/linking_genus.json, from the Seifert matrix V of the KnotInfo package
(database_knotinfo 2026.8.1) and A = V + V^T: the determinant, the group H_K, the signature, s = sigma + 4 Arf mod 8
(Arf from the determinant mod 8, by Levine and Murasugi), the primes p with a_p = 2 with their signs epsilon_p, and the
conditions (i)-(iii) of [C] Corollary 2.3(c) ([T] Corollary 4.5(c)) that hold; and the value of gamma_4 in KnotInfo
2026.10.5.  [C] Table 3 is made from this file alone, so the table does not depend on the installed package.
Usage: python scripts/hthree_certificates.py
Output: data/hthree_certificates.json
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
from database_knotinfo import link_list
from rank1_linking import invariants, signature, legendre
from linking_genus import anisotropic


def load(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


if __name__ == '__main__':
    lg, g4 = load('linking_genus.json'), load('gamma4_knotinfo_2026.10.5.json')
    OPEN = load('open_gamma4_knotinfo_2026.10.5.json')
    h3 = sorted(k for k, v in lg.items() if v['h_witt'] == 3)
    by = {r['name']: r for r in link_list()[1:] if r['name'] in h3}
    out = {}
    for k in h3:
        V = json.loads(by[k]['seifert_matrix'])
        n = len(V)
        A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
        D, sig = lg[k]['det'], signature(A)
        assert sig == lg[k]['sigma']
        inv = invariants(A, D)
        arf = 0 if D % 8 in (1, 7) else 1
        if k in OPEN:
            assert int(OPEN[k]['arf_invariant']) == arf
        s4 = (sig + 4 * arf) % 8
        eps = []
        for p in sorted(inv):
            a, d = anisotropic(inv[p], p)
            if a == 2:
                eps.append([p, legendre(d * (D // p ** sum(v for v, _ in inv[p])), p)])
        crit = []
        if any(p % 4 == 1 and e == -1 for p, e in eps):
            crit.append('(i)')
        if len({e for p, e in eps if p % 4 == 3}) == 2:
            crit.append('(ii)')
        if s4 == 4 and any(p % 4 == 3 and e == -1 for p, e in eps):
            crit.append('(iii)')
        assert crit, k
        out[k] = dict(det=D, group=[[p, v] for p in sorted(inv) for v, _ in sorted(inv[p])], sigma=sig, arf=arf, s=s4,
                      eps=eps, criteria=crit, knotinfo=str(g4[k]).replace(' ', ''))
    meta = dict(seifert_matrices='database_knotinfo 2026.8.1', gamma4_values='KnotInfo data file 2026.10.5')
    json.dump(dict(meta=meta, knots=out), open(os.path.join(ROOT, 'data', 'hthree_certificates.json'), 'w'), indent=1)
    print(len(out), 'knots; criteria', sorted({c for v in out.values() for c in v['criteria']}))
