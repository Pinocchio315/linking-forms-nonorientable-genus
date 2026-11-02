#!/usr/bin/env python
"""The Fox-Milnor condition for the knots with Mobius bands that neither h_lambda nor the signature shows to be not
topologically slice ([C] Section 3.4).

These are the knots with a Moebius band of the first search ([C] Proposition 3.3(a);
data/klein_bottle_search_mobius.json) that have h_lambda = 0 and signature 0 (data/linking_genus.json) and, with the
argument wide, those of the second search ([C] Section 3.3; data/mobius_open_search.json.gz).
If K is topologically slice, then Delta_K(t) = f(t) f(1/t) up to a unit, with f in Z[t] (Fox and Milnor, Osaka J.
Math. 3 (1966), Thm 2, for piecewise linear discs; Friedl, Nagel, Orson and Powell, NYJM Monographs 6 (2025), Thm 15.2,
for locally flat discs).  Since Delta_K(1/t) = Delta_K(t) up to a unit, this holds if and only if every irreducible
factor g of Delta_K over Z with g(1/t) = g(t) up to a unit (a reciprocal factor) has even multiplicity.
For each knot the script computes Delta_K by Fox calculus from the KnotInfo PD code (smooth_top_gap.alexander_fox),
compares it with the Alexander polynomial of KnotInfo, factors it over Z in PARI and in SymPy, and records whether the
condition holds.  The smooth and topological 4-genus listed by KnotInfo are recorded for comparison.
Usage: python scripts/fox_milnor.py [wide]
Output: data/fox_milnor.json (wide: data/fox_milnor_wide.json)
"""
import json, os, sys, warnings
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from smooth_top_gap import alexander_fox, pari
from klein_bottle_search import data, knot_order, knotinfo_alexander, fox_milnor


def knots():
    """The knots with a Moebius band (mode mobius of klein_bottle_search.py) with h_lambda = 0 and signature 0."""
    km, lg = data('klein_bottle_search_mobius.json')['knots'], data('linking_genus.json')
    mb = [k for k, v in km.items() if v['group'] == 'mobius' and v['moebius_band_of_K']]
    return sorted((k for k in mb if lg[k]['h_witt'] == 0 and lg[k]['sigma'] == 0), key=knot_order)


def knots_wide():
    """The knots with a Moebius band of the second construction (scripts/mobius_open_search.py and its extension to every
    knot J, scripts/mobius_J_search.py) with h_lambda = 0 and signature 0."""
    import gzip
    m = json.load(gzip.open(os.path.join(ROOT, 'data', 'mobius_open_search.json.gz'), 'rt'))
    j = json.load(gzip.open(os.path.join(ROOT, 'data', 'mobius_J_search.json.gz'), 'rt'))
    lg = data('linking_genus.json')
    mb = [k for k, v in m['knots'].items() if v['group'] == 'target' and v['moebius_band']]
    mb += [k for k, v in j['knots'].items() if v['group'] == 'target' and v['certificates']]
    return sorted((k for k in mb if lg[k]['h_witt'] == 0 and lg[k]['sigma'] == 0), key=knot_order)


def factors_pari(A):
    """The irreducible factors of positive degree over Z of the polynomial with coefficients A (lowest degree first), as
    (factor, degree, multiplicity, reciprocal)."""
    P = pari(); t = P('t')
    F = P.factor(sum(c * t ** e for e, c in enumerate(A)))
    out = []
    for i in range(int(P.matsize(F)[0])):
        g, m = F[i, 0], int(F[i, 1])
        if P.poldegree(g) > 0:
            out.append((str(g), int(P.poldegree(g)), m, bool(g == P.polrecip(g) or g == -P.polrecip(g))))
    return out


def factors_sympy(A):
    """The same as (degree, multiplicity, reciprocal), sorted, from SymPy's factorisation over Z."""
    import sympy
    t = sympy.Symbol('t')
    out = []
    for g, m in sympy.factor_list(sum(c * t ** e for e, c in enumerate(A)), t)[1]:
        p = sympy.Poly(g, t)
        r = sympy.Poly(p.all_coeffs()[::-1], t)
        out.append((p.degree(), m, p == r or p == -r))
    return sorted(out)


if __name__ == '__main__':
    import cypari, sympy, database_knotinfo
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}
    wide = sys.argv[1:2] == ['wide']
    out = {}
    for k in (knots_wide() if wide else knots()):
        r = BY[k]
        A = alexander_fox(json.loads(r['pd_notation']))
        F = factors_pari(A)
        holds = not any(m % 2 and rec for _, _, m, rec in F)
        assert holds == fox_milnor(A) and sorted((d, m, rec) for _, d, m, rec in F) == factors_sympy(A)
        assert A == A[::-1] and sum(A) == 1                     # symmetric, Delta(1) = 1
        out[k] = dict(alexander=A, agrees_with_knotinfo=A == knotinfo_alexander(r),
                      factors=[dict(factor=g, degree=d, multiplicity=m, reciprocal=rec) for g, d, m, rec in F],
                      fox_milnor=holds, knotinfo_smooth_four_genus=r['smooth_four_genus'].strip(),
                      knotinfo_topological_four_genus=r['topological_four_genus'].strip())
    summary = dict(knots=len(out), agree_with_knotinfo=sum(v['agrees_with_knotinfo'] for v in out.values()),
                   fox_milnor_fails=[k for k, v in out.items() if not v['fox_milnor']],
                   fox_milnor_holds=[k for k, v in out.items() if v['fox_milnor']],
                   knotinfo_topological_four_genus=dict(sorted(
                       (val, [k for k, v in out.items() if v['knotinfo_topological_four_genus'] == val])
                       for val in {v['knotinfo_topological_four_genus'] for v in out.values()})),
                   versions=dict(database_knotinfo=database_knotinfo.version(), cypari=cypari.__version__,
                                 sympy=sympy.__version__))
    with open(os.path.join(ROOT, 'data', 'fox_milnor_wide.json' if wide else 'fox_milnor.json'), 'w') as f:
        json.dump(dict(summary=summary, knots=out), f, indent=1)
    print(json.dumps(summary))
