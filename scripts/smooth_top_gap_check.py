#!/usr/bin/env python
"""Recheck of the bands of [C] Theorem B found by scripts/smooth_top_gap.py ([C] Section 4, independent checks).

For every knot K that data/smooth_top_gap.json certifies (the twelve knots of [C] Proposition 3.1(c) in Theorem B)
one hit is rechecked, the first one in the KnotInfo diagram if there is one, otherwise the first hit: (1) the stored
diagram is K (SnapPy is_isometric_to the exterior of the KnotInfo PD code); (2) the band rebuilt by
bands.band_pd from the stored face, edges and twist t gives a knot J isometric to the stored J; (3)
spherogram's own band code (add_one_band of Dunfield-Gong) gives a knot isometric to J; (4) J has one component,
Delta_J is a unit (SymPy determinant of V - t V^T for a Seifert matrix V of a simplified diagram) and the Euler
characteristic of knot Floer homology is 1; (5) tau(J) from knot Floer homology.  It also tabulates tau over all hits:
a smoothly slice J would contradict the smooth obstruction, and tau(J) != 0 shows that J is not smoothly slice.
Usage: python scripts/smooth_top_gap_check.py
Output: data/smooth_top_gap_check.json
"""
import json, os, sys, warnings
from collections import Counter
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
import snappy, spherogram, sympy
from spherogram.links.bands.core import add_one_band
from database_knotinfo import link_list
from bands import band_pd, from_pd, faces


def iso(A, B, tries=6):
    for _ in range(tries):
        try:
            if A.is_isometric_to(B):
                return True
        except Exception:
            pass
        A.randomize(); B.randomize()
    return False


def check(K, r, h, BY):
    diag = r['diagrams'][h['diagram']]
    pd = diag['pd']
    EK = snappy.Link(json.loads(BY[K]['pd_notation'])).exterior()
    J = band_pd([tuple(x) for x in pd], h['face'], h['edges'], h['t'])
    LJ = spherogram.Link([list(x) for x in J])
    nbr, under, label = from_pd([tuple(x) for x in pd])
    f = faces(nbr)[h['face']]
    B = add_one_band(spherogram.Link([list(x) for x in pd]), ([f[h['darts'][0]], nbr[f[h['darts'][1]]]], [], -h['t']))
    S = LJ.copy(); S.simplify('global')
    V = sympy.Matrix(S.seifert_matrix()); t = sympy.symbols('t')
    poly = sympy.Poly(sympy.expand((V - t * V.T).det(method='bareiss')), t) if len(V) else sympy.Poly(1, t)
    H = S.knot_floer_homology()
    chi = Counter()
    for (a, m), rk in H['ranks'].items():
        chi[a] += (-1) ** m * rk
    return dict(diagram_source=diag['source'], face_edges=h['face_edges'], band_edges=h['edges'], t=h['t'],
                J_label=h['J_label'], diagram_is_K=iso(spherogram.Link([list(x) for x in pd]).exterior(), EK),
                rebuilt_J_is_stored_J=iso(LJ.exterior(), spherogram.Link([list(x) for x in h['J_pd']]).exterior()),
                add_one_band_is_J=len(B.link_components) == 1 and iso(B.exterior(), LJ.exterior()),
                J_components=len(LJ.link_components), J_crossings_simplified=len(S.crossings),
                alexander_seifert_unit=len(poly.terms()) == 1 and abs(poly.terms()[0][1]) == 1,
                alexander_hfk=sorted((a, int(c)) for a, c in chi.items() if c), tau=H['tau'])


if __name__ == '__main__':
    d = json.load(open(os.path.join(ROOT, 'data', 'smooth_top_gap.json')))
    BY = {r['name']: r for r in link_list()[1:]}
    out = {}
    for K in d['summary']['certified']:
        r = d['candidates'][K]
        ki = [h for h in r['hits'] if r['diagrams'][h['diagram']]['source'] == 'knotinfo']
        rec = check(K, r, (ki or r['hits'])[0], BY)
        rec['tau_over_all_hits'] = {f'{lab} tau={tau}': n for (lab, tau), n in
                                    Counter((h['J_label'], h['hfk']['tau'] if h.get('hfk') else None) for h in r['hits']).items()}
        out[K] = rec
        print(K, rec, flush=True)
    ok = all(v['diagram_is_K'] and v['rebuilt_J_is_stored_J'] and v['add_one_band_is_J'] and v['J_components'] == 1
             and v['alexander_seifert_unit'] and v['alexander_hfk'] == [(0, 1)] and v['tau'] != 0 for v in out.values())
    zero_tau = sum(n for v in out.values() for key, n in v['tau_over_all_hits'].items() if key.endswith('tau=0'))
    json.dump(dict(all_checks_pass=ok, hits_with_tau_zero=zero_tau, knots=out),
              open(os.path.join(ROOT, 'data', 'smooth_top_gap_check.json'), 'w'), indent=1)
    print('all checks pass:', ok, '| hits with tau = 0:', zero_tau)
