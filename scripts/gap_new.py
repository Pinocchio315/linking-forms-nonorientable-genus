#!/usr/bin/env python
"""[C] Theorem B for the knots of [C] Proposition 3.2(a), (b): the knot Floer homology of the knot J of the band
that gives each of them a locally flat Mobius band.

The knots are those that Proposition 3.2(a), (b) exclude: data/plumbing_dinv.json (summary field obstructed),
data/surgery_box.json (summary field box_obstructed) and data/surgery_cfk.json (summary field excluded).  For each of
them the band search of scripts/smooth_top_gap.py (one band move from K to a knot J with Delta_J = 1) was run by mode
gap of plumbing_dinv.py, surgery_dinv.py or surgery_cfk.py, and its result is read from data/plumbing_dinv.json,
data/surgery_dinv.json or data/surgery_cfk_gap.json; every knot must have been searched.  For the representative band
of each knot with a certificate, the script checks that J (a stored PD code) has one component, removes the
Reidemeister I and II reductions of the diagram, and computes the knot Floer homology of J with Szabo's program: tau,
epsilon, the Seifert genus and the Alexander polynomial (the graded Euler characteristic), which must be 1.  A J with
tau(J) != 0 is not smoothly slice (Ozsvath-Szabo); for a J with tau(J) = 0, Theorem B itself shows that J is not
smoothly slice, since a smooth slice disc of J and the band would form a smooth Mobius band for K.
Usage: python scripts/gap_new.py
Output: data/gap_new.json
"""
import json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')


def load(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


def hfk_of(pd):
    """Knot Floer homology of the knot J of a PD code (Szabo's program): tau, epsilon, the total rank, the Seifert
    genus, and whether Delta_J = 1."""
    import spherogram
    L = spherogram.Link([list(c) for c in pd])
    assert len(L.link_components) == 1
    L.simplify('basic')                   # Szabo's program refuses diagrams with a reducing Reidemeister I move
    h = L.knot_floer_homology()
    euler = Counter()
    for (a, m), r in h['ranks'].items():
        euler[a] += (-1) ** m * r
    alex = {a: c for a, c in euler.items() if c}
    return dict(tau=h['tau'], epsilon=h['epsilon'], total_rank=h['total_rank'], seifert_genus=h['seifert_genus'],
                alexander_is_one=alex == {0: 1}, crossings=len(L.crossings))


if __name__ == '__main__':
    P, S, B = load('plumbing_dinv.json'), load('surgery_dinv.json'), load('surgery_box.json')
    pob = set(P['summary']['obstructed'])
    sob = set(B['summary']['box_obstructed'])
    out = {}
    for k in sorted(pob):
        g = P['gap'][k]
        assert g['contradictions'] == 0 and g['hits_failing_a_recheck'] == 0
        if g['certificate']:
            out[k] = dict(source='plumbing', representative=g['representative'], J_found=g['J_found'])
    for k in sorted(sob - pob):
        g = S['gap'][k]
        assert g['contradictions'] == 0 and g['hits_failing_a_recheck'] == 0
        if g['certificate']:
            out[k] = dict(source='surgery', representative=g['representative'], J_found=g['J_found'])
    C, CG = load('surgery_cfk.json'), load('surgery_cfk_gap.json')
    cob = set(C['summary']['excluded']) - pob - sob               # excluded with the exact V_j(C) only
    assert cob == set(CG['summary']['knots']) and CG['summary']['contradictions'] == 0
    for k in sorted(cob):
        g = CG['summary']['knots_summary'][k]
        assert g['contradictions'] == 0 and g['hits_failing_a_recheck'] == 0
        if g['certificate']:
            out[k] = dict(source='surgery', representative=g['representative'], J_found=g['J_found'])
    for k, r in out.items():
        r['J'] = hfk_of(r['representative']['J_pd'])
        assert r['J']['alexander_is_one'], k
        print(k, r['source'], r['representative']['J_label'], 'tau', r['J']['tau'], 'epsilon', r['J']['epsilon'], flush=True)
    searched = sorted(set(P['gap']) & pob) + sorted(set(S['gap']) & (sob - pob)) + sorted(cob)
    summary = dict(candidates=len(pob | sob | cob), searched=len(searched), certified=sorted(out),
                   by_source=dict(Counter(r['source'] for r in out.values())),
                   representative_from_knotinfo_diagram=sorted(k for k, r in out.items()
                                                               if r['representative']['diagram_source'] == 'knotinfo'),
                   tau_of_representative_J=dict(sorted(Counter(str(r['J']['tau']) for r in out.values()).items())),
                   tau_zero=sorted(k for k, r in out.items() if r['J']['tau'] == 0),
                   tau_zero_J=sorted({r['representative']['J_label'] for r in out.values() if r['J']['tau'] == 0}),
                   # knots J told apart by an invariant: the volume label (a numerical value, used only to group) and the
                   # Seifert genus from knot Floer homology; J of 13n_1281 and 13n_1307 share the volume, not the genus
                   tau_zero_J_distinct=len({(r['representative']['J_label'], r['J']['seifert_genus'])
                                            for r in out.values() if r['J']['tau'] == 0}))
    assert summary['searched'] == summary['candidates']
    json.dump(dict(summary=summary, knots=out), open(os.path.join(ROOT, 'data', 'gap_new.json'), 'w'), indent=1)
    print(json.dumps(summary, indent=1))
