#!/usr/bin/env sage -python
"""[C] Proposition 3.2(b), with the bounds of [C] Section 2.2: the test of [C] Proposition 2.6 on the surgery
descriptions of scripts/surgery_dinv.py for every V allowed by the bounds of Hom-Wu and Rasmussen, so that no value of
V_j(C) is needed.

Setting.  scripts/surgery_dinv.py identifies, for a knot K with u(K) = 1, a knot C (kappa_hat there) with
Sigma_K = S^3_{D/2}(C) up to orientation (D = det K), and records the knot Floer homology of C from Szabo's program and
a diagram of C.  The correction terms of Y = S^3_{p/q}(C), p, q > 0, are
    d(Y, i) = d(L(p, q), i) - 2 max{V_{floor(i/q)}(C), V_{floor((p+q-1-i)/q)}(C)}
(Ni-Wu, Prop. 1.6 and Remark 2.10), with the labelling checked in surgery_dinv.lens_check.  Each V_j is a nonnegative
integer, and
- V_j >= 1 for 0 <= j < tau(C): nu^+(C) = min{j : V_j = 0} (Hom-Wu, Def. 2.1 and the sentence after it) and
  nu^+(C) >= tau(C) (Hom-Wu, Prop. 2.3(3));
- V_j <= g_4(C) - j (Rasmussen, thesis, Cor. 7.4, where h_j is V_j), and g_4(C) <= g_alg, the genus of the Seifert
  surface that Seifert's algorithm builds from the stored diagram of C.
So V lies in the box V_min <= V <= V_max (componentwise), V_min_j = 1 for j < tau(C) and 0 otherwise, V_max_j =
max(g_alg - j, 0).  Each configuration (sign, l, M, g) of the test of Proposition 2.6 (scripts/rank1_dinv.py, test)
requires inequalities 4 sign d(Z, t) >= b, and d is monotone in each V_j; so a configuration with Z_0 = +Z passes for
some V in the box if and only if it passes for the corner that makes d(Z) largest, and one with Z_0 = -Z if and only if
it passes for the corner that makes d(Z) smallest.  Here Z = -e Y' as in surgery_dinv.py (Y' = S^3_{D/2}(C), e the
sign of the slope).  The knot is "box-obstructed" if neither orientation passes at its corner; then every V allowed by
the bounds is obstructed, and K bounds no smooth Mobius band.
The lower bound V_min needs the chirality of C, which surgery_dinv.py fixes by an orientation-preserving isometry of the
exterior of the stored diagram with E(C), or by Regina's names of the two candidate surgeries when E(C) is not hyperbolic
(field kappa.diagram.orientation); for a route whose chirality is not fixed, V_min = 0 is used.
For each knot every route of surgery_dinv.py that identified C (not a duplicate) and has the slope +-D/2 is tested: the
routes of the unknotting crossings, and the drilling routes of that slope.  The box is that of S^3_{D/2}(C), so a
drilling route of another slope p/q is not used here (surgery_drill.py treats it with the exact V).  This is done for
the three parts of
data/surgery_dinv.json: targets (the 50 knots of Proposition 3.2(b)), calibration, and pilot (the record of 11n_176,
[C] Section 3.8).  Checks: a box-obstructed knot is also obstructed with the exact V where surgery_dinv.py determined V,
the exact V lies in the box, and no calibration knot with KnotInfo gamma_4 = 1 (version 2026.10.5) is box-obstructed.
Usage: sage -python scripts/surgery_box.py
Output: data/surgery_box.json
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'vendor', 'code', 'lower_bounds', 'owens'))
from rank1_dinv import test
import surgery_dinv as sd


def routes(rec):
    """The routes of surgery_dinv.py that identified C, without duplicates, with the slope +-D/2 (module docstring)."""
    return [r for r in rec.get('routes', []) + rec.get('drill_routes', [])
            if r.get('kappa_hat') and not r.get('duplicate_of') and slope_q(r, rec['det']) == 2]


def slope_q(r, D):
    """The denominator q of the slope p/q of a route; the routes of unknotting crossings have the slope +-D/2."""
    if not r.get('slope_p_q'):
        return 2
    p, q = r['slope_p_q']
    assert abs(p) == D, (p, D)
    return q


def g_alg(route):
    """The genus of the Seifert surface of Seifert's algorithm on the stored diagram of kappa (that of kappa_hat is its
    mirror image or the same diagram, with the same genus)."""
    import spherogram
    L = spherogram.Link([list(c) for c in route['kappa']['hfk_diagram']['pd']])
    V = L.seifert_matrix()
    assert V.nrows() == V.ncols() and V.nrows() % 2 == 0
    return V.nrows() // 2


def box(D, e, tau, g):
    """The test at the corners of the box V_min <= V <= V_max (module docstring)."""
    lo = [1 if j < max(tau, 0) else 0 for j in range(g + 1)]
    hi = [max(g - j, 0) for j in range(g + 1)]
    Z = lambda V: (sd.negate(sd.surgery_form(D, 2, V)) if e == 1 else sd.surgery_form(D, 2, V))
    plus_corner, minus_corner = (hi, lo) if e == 1 else (lo, hi)    # d(Z) largest, smallest
    plus, minus = test(Z(plus_corner), (1,)), test(Z(minus_corner), (-1,))
    return dict(V_min=lo, V_max=hi, plus_corner='V_max' if e == 1 else 'V_min', passing_plus=plus, passing_minus=minus,
                obstructed=plus is None and minus is None)


def one(rec):
    """The box test on every route of one knot."""
    out = []
    for r in routes(rec):
        hh = r['kappa_hat']['hfk']
        g = g_alg(r)
        assert g >= hh['seifert_genus'] >= max(hh['tau'], 0)
        ori = r['kappa'].get('diagram', {}).get('orientation')
        b = box(rec['det'], r['e'], hh['tau'] if ori is not None else 0, g)
        V = r.get('V')
        if V is not None:                                      # the exact V lies in the box
            Vp = V + [0] * (g + 1 - len(V))
            b['exact_V_in_box'] = all(lo <= v <= hi for v, lo, hi in zip(Vp, b['V_min'], b['V_max']))
        out.append(dict(crossing=r.get('crossing'), diagram_index=r.get('diagram_index'), e=r['e'], tau=hh['tau'],
                        chirality=ori, g_hfk=hh['seifert_genus'], g_alg=g, kappa=r['kappa'].get('identify'),
                        which=r['kappa_hat']['which'], **b))
    return dict(det=rec['det'], verdict=rec['verdict'], routes=out, box_obstructed=any(x['obstructed'] for x in out))


if __name__ == '__main__':
    S = json.load(open(os.path.join(ROOT, 'data', 'surgery_dinv.json')))
    G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))   # gamma_4 in KnotInfo 2026.10.5
    res = {part: {n: one(rec) for n, rec in S[part].items() if isinstance(rec, dict) and 'det' in rec}
           for part in ('targets', 'calibration', 'pilot')}
    T, C = res['targets'], res['calibration']
    exact = {n for n, r in T.items() if r['verdict'] == 'obstructed'}
    boxed = sorted(n for n, r in T.items() if r['box_obstructed'])
    checks = dict(
        box_implies_exact=[n for n in boxed if S['targets'][n]['verdict'] not in ('obstructed', 'obstructed_by_bounds')],
        exact_V_outside_box=[n for p in res.values() for n, r in p.items() for x in r['routes']
                             if x.get('exact_V_in_box') is False],
        calibration_gamma4_one_box_obstructed=[n for n, r in C.items() if r['box_obstructed'] and G4.get(n) == '1'])
    assert not any(checks.values()), checks
    summary = dict(targets=len(T), box_obstructed=boxed, exact_obstructed=sorted(exact),
                   exact_not_box=sorted(exact - set(boxed)),
                   box_not_exact=sorted(set(boxed) - exact),
                   calibration=dict(knots=len(C), box_obstructed=sorted(n for n, r in C.items() if r['box_obstructed'])),
                   pilot={n: r['box_obstructed'] for n, r in res['pilot'].items()}, checks=checks)
    json.dump(dict(summary=summary, **res), open(os.path.join(ROOT, 'data', 'surgery_box.json'), 'w'), indent=1)
    print(json.dumps(summary, indent=1))
