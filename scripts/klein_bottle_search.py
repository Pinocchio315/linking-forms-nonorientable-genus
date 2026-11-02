#!/usr/bin/env python
"""Surfaces built from band moves: the first search of [C] Section 3.3, for [C] Proposition 3.3 and Section 3.8.

Modes and knots.  The ranges are those of KnotInfo 2026.10.5 (data/open_gamma4_knotinfo_2026.10.5.json), h_lambda is
taken from data/linking_genus.json, and the tests of [C] Proposition 3.1(b), (c) from data/rank1_targets.json (linking
form) and data/rank1_dinv_targets.json, data/rank1_dinv_mtargets.json (correction terms), field obstructed.
  default  Punctured Klein bottles for the 283 knots of range [1,3] that Proposition 3.1(b) or (c) obstructs and that
           have h_lambda <= 2 (Proposition 3.3(b)), and a surface with b_1 = 3 for 13a_2890 (range [2,4], h_lambda = 3;
           Proposition 3.3(c)).  The other 4 knots of range [1,3] that Proposition 3.1(b) or (c) obstructs have
           h_lambda = 3; they and the other 45 knots with h_lambda = 3 are controls, which must give no band to a knot
           with KnotInfo gamma_4 = 1 ([C] Section 3.5).  Field proposition_5_1b: the knot is one of these 287.
  ranges   Punctured Klein bottles for the 1715 knots of range [2,3] (1710) or [2,4] (5) with h_lambda = 2 (Proposition
           3.3(b)); with h_lambda = 2 <= gamma_4^top <= gamma_4 they give gamma_4 = gamma_4^top = 2.  The 41 knots of
           range [2,3] with h_lambda = 3 are controls of the default mode.  Group upper: the 4 knots of range [1,4] or
           [1,5], for which a Mobius band of K by our own bands gives gamma_4 <= 1 (Proposition 3.3(a)); without one,
           the surfaces of the targets would give gamma_4 <= 2 or <= 3.  11n_133 (KnotInfo gamma_4 = 2, but
           h_lambda = 3) is treated as 13a_2890: a surface with b_1 = 3 gives gamma_4 = 3 ([C] Section 3.8).
  mobius   Mobius bands for the 1280 knots of range [1,2] (794) or [1,3] (486) that neither test obstructs (group
           mobius; Proposition 3.3(a)); a Mobius band gives gamma_4 <= 1.  For the knots of group mobius with range
           [1,3] and no Mobius band, punctured Klein bottles as for the targets of mode ranges, which give the range
           [1,2] (the last group of Proposition 3.3(b)).  Controls: the 484 knots of these ranges that a test obstructs
           (197 and 287); any Mobius band of a control is a contradiction ([C] Section 3.5).
  replay   The stored diagrams of a mode scanned again (see Reproducibility).

Surfaces.  A band move from a knot K to a knot J is nonorientable ([C] Section 2.4): pushed into a collar of
S^3 = boundary of B^4, the band together with a surface F' of J in B^4 is a nonorientable surface of K with
b_1 = b_1(F') + 1.  The surface of J is built from our own bands: a face band of a diagram of J to the unknot (band and
disc; unknotting_band), or a face band J -> J' and a face band from J' to the two-component unlink (J' is a ribbon knot;
band, band and two discs; ribbon_moebius_band), a Mobius band of J (moebius_band); then K bounds a punctured Klein
bottle and gamma_4(K) <= 2.  For J with KnotInfo gamma_4 = 2 (13a_2890 and 11n_133, and group upper), two bands
J -> J' -> unknot give a punctured Klein bottle of J (klein_bottle_of), and K bounds a surface with b_1 = 3.  A Mobius
band of K (group upper, mode mobius) is a band of K to the unknot or a band K -> J' and a ribbon band of J'.  The
unknot is recognised by Reidemeister simplification to no crossings or by knot Floer homology of rank one
(unknot_check), the two-component unlink by Reidemeister simplification to no crossings (is_unlink2).
Verification of a certificate (LEVELS): 'self-contained' when the surface of J lies in the diagram of J made by the
band, so that every band lies in a diagram stored here (K -> J -> unknot, K -> J -> J' -> unlink or
K -> J -> J' -> unknot) and the surface depends neither on an identification of J nor on any value of KnotInfo;
'self-contained up to Reidemeister moves of J' when it lies in a diagram that spherogram's Reidemeister moves reach from
the diagram of J made by the band (stage C; the diagram is stored, and scripts/klein_bottle_rebuild.py reaches it
again); 'surface of J in another diagram of J' when it lies in the KnotInfo diagram of J or, as a last resort, in
another stored diagram of J (tied to J by the isometry check); otherwise 'KnotInfo value of gamma_4(J)'.
Reference knots J: the 6479 knots with KnotInfo gamma_4 '1' (the script checks that the tests of [C] Proposition 3.1
obstruct none of them), and, for 13a_2890, 11n_133 and group upper, also the knots with KnotInfo gamma_4 '2' and
h_lambda <= 2 (9_49 and 11n_133 are listed as 2 by KnotInfo, but have h_lambda = 3).  KnotInfo gives slice knots the
value 1 too.  For a target a band to the unknot or to a knot that KnotInfo lists as slice (smooth 4-genus 0) would give
a Mobius band of K, against gamma_4(K) >= 2, and for a control (h_lambda = 3) a band to a knot with KnotInfo gamma_4 = 1
would give gamma_4(K) <= 2; such bands are flagged as contradictions (not for the groups upper and mobius).  KnotInfo's
value of gamma_4(J) is checked independently by the surface of J: KnotInfo is not free of errors (its upper bound 2 for
9_49 and 11n_133 is below h_lambda = 3).

Search.  Face bands of scripts/bands.py, as in scripts/smooth_top_gap.py: for a face of a diagram of K and
two distinct edges on it, J_t is the band move along the arc in that face whose long sides cross t times ([C] Figure 1),
|t| <= T_MAX.  Stage A scans the KnotInfo diagram.  Stage B, for the controls and for the knots without a self-contained
certificate after stage A (group upper: without a Mobius band of K), adds the diagrams of smooth_top_gap.diagrams_of
(spherogram many_diagrams, then nonminimal diagrams with at most c + 6 crossings by random Reidemeister II and III
moves), each checked against the KnotInfo exterior (SnapPy is_isometric_to); the stage-B diagrams stored in the outputs
given as earlier_output.json are scanned there again (field carried).
A knot J_t passes to identification if det J_t (Goeritz matrix, exact) is the determinant of a reference knot and Delta
(Fox calculus, PARI) is the Alexander polynomial of a reference knot; Delta = 1 always passes, so that an unknot is
seen.  Identification: SnapPy identify of the exterior of J, KnotInfo name from the DT name; if SnapPy finds nothing,
the signed Gauss code of the globally simplified diagram is compared with those of the KnotInfo diagrams of the
reference knots with that polynomial and of their mirrors.  A hit: J identified as a reference knot whose KnotInfo
polynomial is Delta_J.
Rechecks of the chosen hits (up to REPS per knot with distinct J, more while none is self-contained, at most
MAX_RECHECK; then, while none is self-contained, further hits whose diagram of J has an unknotting band, SC_TRIES hits
tried, or a band to a ribbon knot, RB_HITS hits tried): one component (traversal and spherogram); Delta_J by Fox
calculus and from a Seifert matrix equal to the KnotInfo polynomial; knot Floer homology of J equal to that of the
KnotInfo diagram of the name or of its mirror; the exterior of J isometric to the exterior of the KnotInfo diagram (by
Gordon-Luecke J is the named knot or its mirror, and gamma_4 is mirror-invariant), or for a non-hyperbolic J the
identification by the diagram; the dual band (cocore) of J gives back K; spherogram's own band code
(links/bands/core.add_one_band of Dunfield-Gong, num_twist = -t) rebuilds J; the surface of J.
Stage C (mode ranges, the targets still without a self-contained certificate after stage B): first an unknotting band of
the diagram of J made by the band, for every hit not tried yet (J with KnotInfo gamma_4 '1'); then, for the hits (the
knots J in turn, at most C_HITS), a Mobius band of J in diagrams reached from the diagram of J made by the band by
spherogram's Reidemeister moves (simplification, and random moves followed by simplification; C_VARIANTS attempts; the
diagram found is stored).  In mode ranges the first band of klein_bottle_of also records the edges of its face
(face_edges), as the other bands do, for scripts/klein_bottle_rebuild.py.
Consistency checks: det J_t = |alpha + n t| on every arc (n = det K, Montesinos trick); |Delta(-1)| = det for every
Alexander polynomial computed; the Fox polynomial of the KnotInfo diagram of every reference knot is its KnotInfo
polynomial.
Mode mobius.  A Mobius band of K is looked for in the KnotInfo diagram (stage A), then in the diagrams of stage B until
one is found (controls: in all of them), and the diagrams searched are stored.  moebius_band_fast does the search of
moebius_band with the same candidates in the same order, and then tries the band to the unknot with the det-1 value of s
outside |s| <= S_MAX (as in scripts/smooth_top_gap.py).  It is faster: by det = |alpha + n s| on each arc (arc_alphas)
it builds only the bands whose det can be 1, a nonzero square or (bands of J' to the unlink) 0, and it passes a link to
is_unlink2 only if two conditions that hold for the two-component unlink do (a minor of the Fox matrix vanishes, every
component has det 1).  For the punctured Klein bottles the diagrams of stage B searched for a Mobius band are carried
into stage B, and a band to the unknot or to a slice knot is no contradiction.

Reproducibility.  spherogram's random moves are not reproducible, so the output stores every diagram scanned.  Mode
replay scans the stored diagrams again and compares the counts and the hits identified by SnapPy diagram by diagram
(identification by the simplified diagram, used for non-hyperbolic J such as 3_1, 8_19 and 10_124, depends on
spherogram's random simplification, so those hits are only counted).  It rebuilds the J of every stored certificate
(same PD code, exterior isometric to the KnotInfo exterior of its name), the surface of J and every Mobius band of K
(the unknot or the two-component unlink recognised again).  In mode mobius it compares the Mobius band search diagram by
diagram, repeats the punctured Klein bottle searches, and also runs moebius_band in the KnotInfo diagram of every knot,
which must find the band of moebius_band_fast there (none if that has |s| > S_MAX).  Replay writes nothing;
data/klein_bottle_search_mobius_replay.txt is the printed output of mode mobius replay.
Mode mobius also writes the records, as they arrive, to the output path plus .partial; a knot whose computation raises
gets the record {error}, and then no output is written; with KLEIN_RESUME=1 a new run takes the other records of the
.partial file and runs only the remaining knots.  Its workers may grow their PARI stack up to PARI_SIZEMAX: the default
maximum, 8 MB, is too small for some knots.
Usage: python scripts/klein_bottle_search.py [ranges | mobius] [workers] [earlier_output.json ...] | [ranges | mobius] replay [workers]
Environment, for tests: KLEIN_ONLY=knot,knot,... restricts the knots (not in replay); KLEIN_OUT=path replaces the output
file, and in replay the file read; KLEIN_RESUME=1 in mode mobius, see above.
Output: data/klein_bottle_search.json, in mode ranges data/klein_bottle_search_ranges.json, in mode mobius
data/klein_bottle_search_mobius.json
"""
import json, os, re, sys, time, warnings, zlib
from collections import Counter, defaultdict
from itertools import combinations
from math import isqrt
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from bands import from_pd, faces, band, pd_code, components, graph_components
from smooth_top_gap import (det_link, alexander_fox, alexander_seifert, normalize, poly_coeffs, pari, gauss_key,
                            isometric, diagrams_of, verify_diagrams, orientation)

T_MAX = 6                                 # |t| <= T_MAX for the band K -> J
S_MAX = 6                                 # |s| <= S_MAX for the band J -> unknot
REPS, MAX_RECHECK = 3, 12                 # rechecked hits per knot, with distinct J
SC_TRIES = 200                            # further hits whose J is tested for an unknotting band (see certify)
RB_HITS = 30                              # further hits whose J is tested for a band to a ribbon knot (see certify)
RB_TRIES = 20                             # knots J' tried for a ribbon band J' -> two-component unlink
KB_TRIES = 40                             # knots J' tried for a band J' -> unknot when gamma_4(J) = 2
J_DIAGRAMS = (10, 20)                     # minimal and nonminimal diagrams of J tried as a last resort
C_HITS, C_VARIANTS = 40, 2000             # stage C (mode ranges): hits tried, attempts to reach diagrams of J
PARI_STACK, PARI_SIZEMAX = 2 ** 25, 2 ** 31  # mode mobius: PARI stack of the workers, and its maximum
REACHED = 'diagram of J reached from J made by the band'


def data(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


def knot_order(name):
    a, b = name.split('_')
    c = re.match(r'\d+', a).group()
    return (int(c), a[len(c):], int(b))


def knots():
    """Targets and controls (see the docstring), ordered by crossing number, type and index (knot_order)."""
    op = data('open_gamma4_knotinfo_2026.10.5.json')
    rng = {k: v['smooth_4d_crosscap_number'].replace(' ', '') for k, v in op.items()}
    r1, da, dm = data('rank1_targets.json'), data('rank1_dinv_targets.json'), data('rank1_dinv_mtargets.json')
    LG, G = data('linking_genus.json'), data('gamma4_knotinfo_2026.10.5.json')
    lam = {k for k, v in r1.items() if v['obstructed']}
    dd = {k for k, v in da.items() if v['obstructed']} | {k for k, v in dm.items() if v['obstructed'] and k in rng}
    three = [k for k in lam | dd if rng.get(k) == '[1,3]']
    assert len(three) == 287 and not lam & dd
    out = {}
    for k in three:
        h = LG[k]['h_witt']
        out[k] = dict(group='target' if h <= 2 else 'control', proposition_5_1b=True, knotinfo_gamma4=G[k], h_witt=h,
                      test='linking form' if k in lam else 'd-invariant', det=LG[k]['det'], J_gamma4=['1'])
    for k, v in LG.items():
        if v['h_witt'] == 3 and k not in out:
            out[k] = dict(group='control', proposition_5_1b=False, knotinfo_gamma4=G[k], h_witt=3, test='h_lambda = 3',
                          det=v['det'], J_gamma4=['1'])
    assert sum(v['group'] == 'target' for v in out.values()) == 283 and sum(v['group'] == 'control' for v in out.values()) == 49
    assert [k for k, v in out.items() if v['knotinfo_gamma4'] == '[2,4]' and v['h_witt'] == 3] == ['13a_2890']
    out['13a_2890']['J_gamma4'] = ['1', '2']
    return dict(sorted(out.items(), key=lambda kv: knot_order(kv[0])))


def ranges_knots():
    """Mode ranges: targets of range [2,3] or [2,4] with h_lambda = 2, 11n_133, and the group upper (see the docstring)."""
    op = data('open_gamma4_knotinfo_2026.10.5.json')
    rng = {k: v['smooth_4d_crosscap_number'].replace(' ', '') for k, v in op.items()}
    LG, G = data('linking_genus.json'), data('gamma4_knotinfo_2026.10.5.json')
    out = {}
    for k, r in rng.items():
        h = LG[k]['h_witt']
        if r in ('[2,3]', '[2,4]') and h == 2:
            out[k] = dict(group='target', range=r, knotinfo_gamma4=G[k], h_witt=h, det=LG[k]['det'], J_gamma4=['1'],
                          stage_C=True)
        elif r in ('[1,4]', '[1,5]'):
            out[k] = dict(group='upper', range=r, knotinfo_gamma4=G[k], h_witt=h, det=LG[k]['det'], J_gamma4=['1', '2'])
    groups = Counter((v['group'], v['range']) for v in out.values())
    assert groups == {('target', '[2,3]'): 1710, ('target', '[2,4]'): 5, ('upper', '[1,4]'): 3, ('upper', '[1,5]'): 1}
    assert all(v['h_witt'] <= 1 for v in out.values() if v['group'] == 'upper')
    assert sum(r == '[2,3]' and LG[k]['h_witt'] == 3 for k, r in rng.items()) == 41
    assert [k for k, r in rng.items() if r == '[2,4]' and LG[k]['h_witt'] == 3] == ['13a_2890']
    assert '11n_133' not in rng and G['11n_133'] == '2' and LG['11n_133']['h_witt'] == 3
    out['11n_133'] = dict(group='control', range=None, knotinfo_gamma4=G['11n_133'], h_witt=3, det=LG['11n_133']['det'],
                          J_gamma4=['1', '2'])
    return dict(sorted(out.items(), key=lambda kv: knot_order(kv[0])))


def mobius_knots():
    """Mode mobius: the knots of range [1,2] or [1,3], group mobius if no test of the paper obstructs them, else control
    (see the docstring)."""
    op = data('open_gamma4_knotinfo_2026.10.5.json')
    rng = {k: v['smooth_4d_crosscap_number'].replace(' ', '') for k, v in op.items()}
    r1, da, dm = data('rank1_targets.json'), data('rank1_dinv_targets.json'), data('rank1_dinv_mtargets.json')
    LG, G = data('linking_genus.json'), data('gamma4_knotinfo_2026.10.5.json')
    lam = {k for k, v in r1.items() if v['obstructed']}
    dd = {k for k, v in da.items() if v['obstructed']} | {k for k, v in dm.items() if v['obstructed'] and k in rng}
    assert not lam & dd
    out = {}
    for k, r in rng.items():
        if r in ('[1,2]', '[1,3]'):
            test = 'linking form' if k in lam else 'd-invariant' if k in dd else None
            out[k] = dict(group='control' if test else 'mobius', range=r, test=test, knotinfo_gamma4=G[k],
                          h_witt=LG[k]['h_witt'], det=LG[k]['det'], J_gamma4=['1'])
    groups = Counter((v['group'], v['range']) for v in out.values())
    assert groups == {('mobius', '[1,2]'): 794, ('mobius', '[1,3]'): 486, ('control', '[1,2]'): 197,
                      ('control', '[1,3]'): 287}
    assert Counter((v['range'], v['test']) for v in out.values() if v['test']) == {
        ('[1,2]', 'linking form'): 116, ('[1,2]', 'd-invariant'): 81, ('[1,3]', 'linking form'): 135,
        ('[1,3]', 'd-invariant'): 152}
    assert all(v['knotinfo_gamma4'] == v['range'] for v in out.values())
    assert all(v['h_witt'] <= 1 and v['det'] >= 3 for v in out.values() if v['group'] == 'mobius')
    return dict(sorted(out.items(), key=lambda kv: knot_order(kv[0])))


def knotinfo_alexander(r):
    if r['alexander_polynomial_vector'].strip():
        v = json.loads(r['alexander_polynomial_vector'])
        return normalize(dict(zip(range(v[0], v[1] + 1), v[2:])))
    return normalize(poly_coeffs(pari()(r['alexander_polynomial'])))


def references(BY, G, LG):
    """Reference knots J: KnotInfo gamma_4 '1', and '2' with h_lambda <= 2 (used for 13a_2890 only)."""
    r1, da, dm = data('rank1_targets.json'), data('rank1_dinv_targets.json'), data('rank1_dinv_mtargets.json')
    obstructed = ({k for k, v in r1.items() if v['obstructed']} | {k for k, v in da.items() if v['obstructed']}
                  | {k for k, v in dm.items() if v['obstructed']})
    ref = {}
    for k, v in G.items():
        if v == '1':
            assert LG[k]['h_witt'] <= 1 and k not in obstructed      # the tests of [C] Prop. 3.1 allow a Mobius band
        elif not (v == '2' and LG[k]['h_witt'] <= 2):
            continue
        r = BY[k]
        ref[k] = dict(g4=v, det=int(r['determinant']), alex=knotinfo_alexander(r), pd=json.loads(r['pd_notation']),
                      crossings=int(r['crossing_number']), smooth_four_genus=r['smooth_four_genus'])
    return ref


# ---------- worker state

REF, DT, G4 = {}, {}, {}                  # reference knots, SnapPy census name -> KnotInfo name, KnotInfo gamma_4
BYPOLY, DETS, GKEYS = defaultdict(list), {}, {}


def init_worker(ref, dt, g4):
    import snappy, spherogram             # spherogram's exterior() needs snappy imported
    REF.update(ref); DT.update(dt); G4.update(g4)
    for k, v in ref.items():
        BYPOLY[tuple(v['alex'])].append(k)
    for allowed in (('1',), ('1', '2')):
        DETS[allowed] = {v['det'] for v in ref.values() if v['g4'] in allowed} | {1}


def gauss_keys(name):
    """Signed Gauss codes of the KnotInfo diagram of a reference knot and of its mirror."""
    if name not in GKEYS:
        import spherogram
        L = spherogram.Link(REF[name]['pd'])
        GKEYS[name] = {gauss_key([tuple(x) for x in L.PD_code()]), gauss_key([tuple(x) for x in L.mirror().PD_code()])}
    return GKEYS[name]


def census_names(E):
    for _ in range(4):
        try:
            return [str(m).split('(')[0] for m in E.identify()]
        except Exception:
            E.randomize()
    return []


def identify(P, alex):
    """KnotInfo name of the knot with PD code P and Alexander polynomial alex (see the docstring)."""
    import spherogram
    L = spherogram.Link([list(x) for x in P])
    ids = census_names(L.exterior())
    names = sorted({DT[x] for x in ids if x in DT})
    if names:
        return dict(how='snappy', names=names)
    if ids:
        return dict(how='snappy', names=[], census=ids[0])
    S = L.copy(); S.simplify('global')
    if not S.crossings:
        return dict(how='unknot', names=[])
    key = gauss_key([tuple(x) for x in S.PD_code()])
    names = [k for k in BYPOLY.get(tuple(alex), ()) if key in gauss_keys(k)]
    return dict(how='diagram' if names else 'unidentified', names=names, crossings=len(S.crossings))


# ---------- the search

def scan(pd, n, allowed):
    """Face bands of one diagram of K, |t| <= T_MAX: counts, fit failures and the knots J that pass the det and
    Alexander filters for the reference knots with KnotInfo gamma_4 in allowed."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    st, cands, fit_fail = Counter(), [], []
    for fi, f in enumerate(faces(nbr)):
        for i, j in combinations(range(len(f)), 2):
            h1, h3 = f[i], f[j]
            if label[h1] == label[h3]:
                continue
            st['arcs'] += 1
            fam = {}
            for t in range(-T_MAX, T_MAX + 1):
                N, U = band(nbr, under, h1, h3, t)
                st['bands'] += 1
                if graph_components(N, len(U)) > 1:
                    fam[t] = (2, 0); st['split'] += 1; continue
                P = pd_code(N, U)
                fam[t] = (len(components(N)), det_link(P))
                if fam[t][0] != 1:
                    continue
                st['knots'] += 1
                if fam[t][1] not in DETS[allowed]:
                    continue
                st['det_pass'] += 1
                A = alexander_fox(P)
                st['fox_det_checks'] += 1
                st['fox_det_failures'] += abs(sum(x * (-1) ** e for e, x in enumerate(A))) != fam[t][1]
                if A == [1] or any(REF[k]['g4'] in allowed for k in BYPOLY.get(tuple(A), ())):
                    st['alexander_pass'] += 1
                    cands.append(dict(face=fi, darts=[i, j], edges=[label[h1], label[h3]], t=t, alexander=A, pd=P))
            if not any(all(abs(a + n * t) == d for t, (_, d) in fam.items()) for a in {fam[0][1], -fam[0][1]}):
                fit_fail.append([fi, i, j]); st['fit_failures'] += 1
    return st, cands, fit_fail


def classify(c, allowed, group):
    """Identify the knot J of a candidate band and decide whether it is a hit (and whether it is a contradiction)."""
    idf = identify(c['pd'], c['alexander'])
    names = [k for k in idf['names'] if k in REF and REF[k]['alex'] == c['alexander'] and REF[k]['g4'] in allowed]
    out = dict(how=idf['how'], names=idf['names'], census=idf.get('census'), J=names[0] if len(names) == 1 else None,
               unknot=idf['how'] == 'unknot')
    out['alexander_mismatch'] = any(k in REF and REF[k]['alex'] != c['alexander'] for k in idf['names'])
    if out['J']:
        r = REF[out['J']]
        out['contradiction'] = group not in ('upper', 'mobius') and (r['smooth_four_genus'] == '0' or
                                                                     (group == 'control' and r['g4'] == '1'))
    else:
        out['contradiction'] = out['unknot'] and group not in ('upper', 'mobius')
    return out


def unknot_check(P):
    import spherogram
    S = spherogram.Link([list(x) for x in P]); S.simplify('global')
    if not S.crossings:
        return 'simplifies to no crossings'
    if S.knot_floer_homology()['total_rank'] == 1:
        return 'knot Floer homology of rank one'
    return None


def unknotting_band(pd):
    """A face band of the diagram pd of a knot to the unknot, the least |s| <= S_MAX first; None if there is none."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = [(fi, f, i, j) for fi, f in enumerate(faces(nbr)) for i, j in combinations(range(len(f)), 2)
            if label[f[i]] != label[f[j]]]
    for s in sorted(range(-S_MAX, S_MAX + 1), key=lambda x: (abs(x), x)):
        for fi, f, i, j in arcs:
            N, U = band(nbr, under, f[i], f[j], s)
            if graph_components(N, len(U)) > 1 or len(components(N)) != 1:
                continue
            P = pd_code(N, U)
            if det_link(P) != 1 or alexander_fox(P) != [1]:
                continue
            how = unknot_check(P)
            if how:
                return dict(face=fi, darts=[i, j], face_edges=[label[h] for h in f], edges=[label[f[i]], label[f[j]]],
                            s=s, recognised_by=how, crossings=len(P))
    return None


def fox_milnor(A):
    """Delta = f(t) f(1/t) up to units: every reciprocal irreducible factor has even multiplicity (true for slice knots)."""
    P = pari(); t = P('t')
    F = P.factor(sum(c * t ** e for e, c in enumerate(A)))
    for i in range(int(P.matsize(F)[0])):
        g, m = F[i, 0], int(F[i, 1])
        if P.poldegree(g) > 0 and m % 2 and (g == P.polrecip(g) or g == -P.polrecip(g)):
            return False
    return True


def is_unlink2(P):
    import spherogram
    L = spherogram.Link([list(x) for x in P]); L.simplify('global')
    return len(L.crossings) == 0 and L.unlinked_unknot_components == 2


def ribbon_band(pd):
    """A face band of the diagram pd of a knot to the two-component unlink (recognised by Reidemeister simplification
    to no crossings): with two discs it is a ribbon disc of the knot.  The least |r| <= S_MAX first; None if none."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = [(fi, f, i, j) for fi, f in enumerate(faces(nbr)) for i, j in combinations(range(len(f)), 2)
            if label[f[i]] != label[f[j]]]
    for r in sorted(range(-S_MAX, S_MAX + 1), key=lambda x: (abs(x), x)):
        for fi, f, i, j in arcs:
            N, U = band(nbr, under, f[i], f[j], r)
            if len(components(N)) != 2:
                continue
            P = pd_code(N, U)
            if graph_components(N, len(U)) == 1 and det_link(P) != 0:
                continue
            if is_unlink2(P):
                return dict(face=fi, darts=[i, j], face_edges=[label[h] for h in f], edges=[label[f[i]], label[f[j]]],
                            r=r, crossings=len(P))
    return None


def ribbon_moebius_band(pd):
    """A Mobius band of the knot J with diagram pd: a face band J -> J' and a ribbon disc of J' (ribbon_band).  J' runs
    over the knots of the face bands of pd (least |s| <= S_MAX first) whose det is a square and whose Alexander
    polynomial has the Fox-Milnor form, at most RB_TRIES of them."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = [(fi, f, i, j) for fi, f in enumerate(faces(nbr)) for i, j in combinations(range(len(f)), 2)
            if label[f[i]] != label[f[j]]]
    tried = 0
    for s in sorted(range(-S_MAX, S_MAX + 1), key=lambda x: (abs(x), x)):
        for fi, f, i, j in arcs:
            N, U = band(nbr, under, f[i], f[j], s)
            if graph_components(N, len(U)) > 1 or len(components(N)) != 1:
                continue
            P = pd_code(N, U)
            d = det_link(P)
            if isqrt(d) ** 2 != d:
                continue
            A = alexander_fox(P)
            if not fox_milnor(A):
                continue
            rb = ribbon_band(P)
            if rb:
                return dict(face=fi, darts=[i, j], face_edges=[label[h] for h in f], edges=[label[f[i]], label[f[j]]],
                            s=s, J1_pd=[list(x) for x in P], J1_det=d, J1_alexander=A, ribbon_band_of_J1=rb)
            tried += 1
            if tried >= RB_TRIES:
                return None
    return None


def moebius_band(pd):
    """A Mobius band of the knot with diagram pd by our own bands: an unknotting band, else a band to a ribbon knot."""
    mb = unknotting_band(pd)
    if mb:
        return dict(mb, kind='band to the unknot')
    rb = ribbon_moebius_band(pd)
    return dict(rb, kind='band to a ribbon knot') if rb else None


# ---------- mode mobius: the search of moebius_band with the dets of the bands predicted

def arc_alphas(nbr, under, h1, h3, n):
    """The alpha with det = |alpha + n s| for the bands with s = -1, 0, 1 half twists on the arc of the darts h1, h3
    (split diagrams: det 0).  By the Montesinos trick det = |alpha + n s| for every s (n = det of the knot; checked on
    every arc by scan and by scripts/smooth_top_gap.py).  [] if no alpha fits."""
    dets = {}
    for s in (-1, 0, 1):
        N, U = band(nbr, under, h1, h3, s)
        dets[s] = 0 if graph_components(N, len(U)) > 1 else det_link(pd_code(N, U))
    return sorted(a for a in {dets[0], -dets[0]} if all(abs(a + n * s) == d for s, d in dets.items()))


def fox_minor_vanishes(P):
    """Whether the minor without the last row and column of the Fox matrix of the Wirtinger presentation of the link
    diagram P (n crossings, rows as in smooth_top_gap.alexander_fox, components oriented by the traversal) vanishes;
    True also when it is not formed (a component without undercrossing).  For the two-component unlink it vanishes: it
    lies in the elementary ideal E_1, an invariant of the group (Fox), and E_1 = 0 for the free group of rank two."""
    n = len(P)
    entering, comps = orientation(P)
    if n < 2 or not all((c, 0) in entering and ((c, 1) in entering) != ((c, 3) in entering) for c in range(n)):
        return True
    labels = sorted({x for X in P for x in X})
    parent = {x: x for x in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for X in P:
        parent[find(X[1])] = find(X[3])
    roots = sorted({find(x) for x in labels})
    if len(roots) != n:
        return True
    arc = {x: roots.index(find(x)) for x in labels}
    p = pari(); t = p('t')
    M = [[p(0)] * n for _ in range(n)]
    for c, X in enumerate(P):
        k, i, j = arc[X[1]], arc[X[0]], arc[X[2]]
        M[c][k] += 1 - t
        if (c, 3) in entering:
            M[c][i] += t; M[c][j] -= 1
        else:
            M[c][i] -= 1; M[c][j] += t
    return p.matdet(p.matrix(n - 1, n - 1, [M[r][s] for r in range(n - 1) for s in range(n - 1)])) == 0


def components_det_one(P):
    """Whether every component of the link diagram P, drawn alone (its crossings with the other components removed),
    has det 1, as every component of an unlink has.  A component with at most one crossing of its own counts as det 1."""
    parent = {}

    def find(q, x):
        q.setdefault(x, x)
        while q[x] != x:
            q[x] = q[q[x]]; x = q[x]
        return x
    for a, b, c, d in P:
        parent[find(parent, a)] = find(parent, c); parent[find(parent, b)] = find(parent, d)
    comp = {x: find(parent, x) for X in P for x in X}
    for C in set(comp.values()):
        q, own = {}, []
        for a, b, c, d in P:
            if comp[a] == C and comp[b] == C:
                own.append((a, b, c, d))
            elif comp[a] == C:                    # C passes under another component: one edge of C
                q[find(q, a)] = find(q, c)
            elif comp[b] == C:
                q[find(q, b)] = find(q, d)
        if len(own) > 1 and det_link([tuple(find(q, x) for x in X) for X in own]) != 1:
            return False
    return True


def ribbon_band_fast(pd, n1, st):
    """ribbon_band for a knot J' with det n1: on an arc only the band whose det can be 0 (r = -alpha / n1) is built, every
    r if no alpha fits, in the order of ribbon_band; a link goes to is_unlink2 only if fox_minor_vanishes and
    components_det_one, both true for the two-component unlink."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = [(fi, f, i, j) for fi, f in enumerate(faces(nbr)) for i, j in combinations(range(len(f)), 2)
            if label[f[i]] != label[f[j]]]
    cands = []
    for k, (fi, f, i, j) in enumerate(arcs):
        al = arc_alphas(nbr, under, f[i], f[j], n1)
        st['J1_arcs'] += 1
        if not al:
            st['J1_fit_failures'] += 1
            cands += [(r, k, False) for r in range(-S_MAX, S_MAX + 1)]
        cands += [(-a // n1, k, True) for a in al if a % n1 == 0 and abs(a // n1) <= S_MAX]
    for r, k, predicted in sorted(cands, key=lambda x: (abs(x[0]), x[0], x[1])):
        fi, f, i, j = arcs[k]
        N, U = band(nbr, under, f[i], f[j], r)
        if len(components(N)) != 2:
            st['det_mismatch'] += predicted
            continue
        P = pd_code(N, U)
        if graph_components(N, len(U)) == 1 and det_link(P) != 0:
            st['det_mismatch'] += predicted
            continue
        if not (fox_minor_vanishes(P) and components_det_one(P)):
            continue
        st['unlink_checks'] += 1
        if is_unlink2(P):
            return dict(face=fi, darts=[i, j], face_edges=[label[h] for h in f], edges=[label[f[i]], label[f[j]]],
                        r=r, crossings=len(P))
    return None


def moebius_band_fast(pd, n):
    """Mode mobius: the search of moebius_band, then a band to the unknot with |s| > S_MAX (the det-1 value of s, as in
    scripts/smooth_top_gap.py), with the same candidates in the same order but faster.  On every arc det = |alpha + n s|
    (arc_alphas), so only the bands whose det can be 1 (bands to the unknot) or a nonzero square (bands to a ribbon knot
    J') are built, and ribbon_band_fast does the same for J'; an arc where no alpha fits is searched in full, and the det
    of every band built is compared with the prediction (det_mismatch).  Returns the Mobius band or None, counts, and
    the arcs of pd where no alpha fits."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = [(fi, f, i, j) for fi, f in enumerate(faces(nbr)) for i, j in combinations(range(len(f)), 2)
            if label[f[i]] != label[f[j]]]
    st, fit_fail, alphas = Counter(arcs=len(arcs)), [], []
    for fi, f, i, j in arcs:
        alphas.append(arc_alphas(nbr, under, f[i], f[j], n))
        if not alphas[-1]:
            fit_fail.append([fi, i, j])
    dets = lambda k, s: {abs(a + n * s) for a in alphas[k]}
    order = [(s, k) for s in sorted(range(-S_MAX, S_MAX + 1), key=lambda x: (abs(x), x)) for k in range(len(arcs))]

    def knot(k, s):
        """The PD code and det of the knot of the band, or None if it is not a knot; det_mismatch when det is not
        |alpha + n s|."""
        fi, f, i, j = arcs[k]
        N, U = band(nbr, under, f[i], f[j], s)
        if graph_components(N, len(U)) > 1:
            st['det_mismatch'] += bool(alphas[k]) and 0 not in dets(k, s)
            return None, None
        P = pd_code(N, U)
        d = det_link(P)
        st['det_mismatch'] += bool(alphas[k]) and d not in dets(k, s)
        return (P, d) if len(components(N)) == 1 else (None, None)

    def found(k, s, **kw):
        fi, f, i, j = arcs[k]
        return dict(face=fi, darts=[i, j], face_edges=[label[h] for h in f], edges=[label[f[i]], label[f[j]]], s=s, **kw)
    for s, k in order:                            # unknotting_band
        if alphas[k] and 1 not in dets(k, s):
            continue
        P, d = knot(k, s)
        if d == 1 and alexander_fox(P) == [1]:
            how = unknot_check(P)
            if how:
                return found(k, s, recognised_by=how, crossings=len(P), kind='band to the unknot'), st, fit_fail
    tried = 0
    for s, k in order:                            # ribbon_moebius_band
        if alphas[k] and not any(x and isqrt(x) ** 2 == x for x in dets(k, s)):
            continue
        P, d = knot(k, s)
        if d is None or isqrt(d) ** 2 != d:
            continue
        A = alexander_fox(P)
        if not fox_milnor(A):
            continue
        rb = ribbon_band_fast(P, d, st)
        if rb:
            return found(k, s, J1_pd=[list(x) for x in P], J1_det=d, J1_alexander=A, ribbon_band_of_J1=rb,
                         kind='band to a ribbon knot'), st, fit_fail
        tried += 1
        if tried >= RB_TRIES:
            break
    far = [((e - a) // n, k) for k in range(len(arcs)) for a in alphas[k] for e in (1, -1)
           if (e - a) % n == 0 and abs((e - a) // n) > S_MAX]
    for s, k in sorted(far, key=lambda x: (abs(x[0]), x)):     # the det-1 value of s outside |s| <= S_MAX
        st['far_bands'] += 1
        P, d = knot(k, s)
        if d == 1 and alexander_fox(P) == [1]:
            how = unknot_check(P)
            if how:
                return found(k, s, recognised_by=how, crossings=len(P), kind='band to the unknot'), st, fit_fail
    return None, st, fit_fail


def klein_bottle_of(pd):
    """For a knot J with gamma_4(J) = 2: bands J -> J' -> unknot, the first in the diagram pd of J and the second in
    the diagram of J' it makes, so that J bounds a punctured Klein bottle.  J' runs over the knots of scan (det and
    Alexander polynomial of a reference knot with gamma_4 = 1), the least |t| first, at most KB_TRIES of them."""
    _, cands, _ = scan(pd, det_link(pd), ('1',))
    for c in sorted(cands, key=lambda c: (abs(c['t']), len(c['pd']), c['face'], c['darts']))[:KB_TRIES]:
        mb = unknotting_band(c['pd'])
        if mb:
            return dict(face=c['face'], darts=c['darts'], edges=c['edges'], t=c['t'], J1_pd=[list(x) for x in c['pd']],
                        J1_alexander=c['alexander'], unknotting_band_of_J1=mb)
    return None


def other_diagrams(pd, seed):
    """Diagrams of a knot other than pd: many_diagrams (minimal), then random Reidemeister II and III moves on those."""
    import random, spherogram
    random.seed(seed)
    minimal = spherogram.Link([list(x) for x in pd]).many_diagrams(target=J_DIAGRAMS[0], tries=60, method='backtrack')
    out = list(minimal)
    c = min(len(D.crossings) for D in minimal)
    for i in range(J_DIAGRAMS[1]):
        random.seed(seed + 1 + i)
        D = minimal[i % len(minimal)].copy()
        D.backtrack(steps=1 + i % 4, prob_type_1=0, prob_type_2=0.5)
        if c < len(D.crossings) <= c + 4:
            out.append(D)
    return out


def moebius_band_elsewhere(name):
    """Last resort for a reference knot J without a Mobius band in the two diagrams of recheck: a Mobius band in
    another diagram of J (stored, with its isometry check against the KnotInfo exterior)."""
    import spherogram
    E0 = spherogram.Link(REF[name]['pd']).exterior()
    for D in other_diagrams(REF[name]['pd'], zlib.crc32(name.encode())):
        P = [tuple(x) for x in D.PD_code()]
        mb = moebius_band(P)
        if mb:
            return dict(mb, diagram='another diagram of J', diagram_pd=[list(x) for x in P],
                        diagram_isometric_to_knotinfo=isometric(D.exterior(), E0))
    return None


def moebius_band_reached(pd, seed):
    """Stage C: a Mobius band of the knot with diagram pd in a diagram reached from pd by spherogram's Reidemeister
    moves: simplify('global'), then random moves that do not decrease the crossing number (backtrack, 1 to 12 steps)
    followed by simplify('global'); C_VARIANTS attempts, each new diagram searched once (the diagram is stored, the
    moves are not, and scripts/klein_bottle_rebuild.py finds the diagram from pd again)."""
    import random, spherogram
    keys = set()
    for i in range(C_VARIANTS):
        random.seed(seed + i)
        L = spherogram.Link([list(x) for x in pd])
        if i:
            L.backtrack(steps=1 + i % 12, prob_type_1=0.2, prob_type_2=0.4)
        L.simplify('global')
        V = [tuple(x) for x in L.PD_code()]
        if not V or gauss_key(V) in keys:
            continue
        keys.add(gauss_key(V))
        mb = moebius_band(V)
        if mb:
            return dict(mb, diagram=REACHED, diagram_pd=[list(x) for x in V], variant=i)
    return None


def is_hyperbolic(E):
    """True if SnapPy finds a solution with all tetrahedra positively oriented (retriangulating up to five times)."""
    for _ in range(6):
        if E.solution_type() == 'all tetrahedra positively oriented':
            return True
        E.randomize()
    return False


def recheck(pd, hit, E_K, surface=None):
    """Independent checks of a hit (see the docstring).  surface: a Mobius band of J found by stage C, used when the
    diagram of J made by the band has none."""
    import spherogram
    from spherogram.links.bands.core import add_one_band
    nbr, under, label = from_pd([tuple(x) for x in pd])
    f = faces(nbr)[hit['face']]
    h1, h3 = f[hit['darts'][0]], f[hit['darts'][1]]
    t, ref = hit['t'], REF[hit['J']]
    N, U = band(nbr, under, h1, h3, t)
    P = pd_code(N, U)
    L = spherogram.Link([list(x) for x in P])
    rec = dict(face_edges=[label[h] for h in f], J_pd=[list(x) for x in P], J_crossings=len(P),
               components_traversal=len(components(N)), components_spherogram=len(L.link_components),
               alexander_fox=alexander_fox(P), alexander_seifert=alexander_seifert(L.seifert_matrix()),
               alexander_knotinfo=ref['alex'], J_knotinfo_crossings=ref['crossings'])
    S = L.copy(); S.simplify('global')
    rec['J_crossings_simplified'] = len(S.crossings)
    H = (S.knot_floer_homology() if S.crossings else
         dict(ranks={(0, 0): 1}, tau=0, seifert_genus=0, total_rank=1, fibered=True))
    Href = spherogram.Link(ref['pd']).knot_floer_homology()['ranks']
    rec['hfk'] = {k: H[k] for k in ('tau', 'seifert_genus', 'total_rank', 'fibered')}
    rec['hfk_equals_knotinfo'] = ('same' if H['ranks'] == Href else
                                  'mirror' if H['ranks'] == {(-a, -m): x for (a, m), x in Href.items()} else False)
    E_J = L.exterior()
    hyperbolic = is_hyperbolic(E_J)
    rec['J_solution_type'] = E_J.solution_type()
    rec['identified_by'] = hit['how']
    rec['isometric_to_knotinfo'] = isometric(E_J, spherogram.Link(ref['pd']).exterior()) if hyperbolic else None
    # the cocore of the band: a face band with t = 0 between the dart h2 = nbr[h1] and the dart running back to h1
    h2, d = nbr[h1], ((len(under), 1) if t else nbr[h3])
    assert N[d] == h1 and any(h2 in g and d in g for g in faces(N))
    N2, U2 = band(N, U, h2, d, 0)
    rec['dual_band_gives_K'] = isometric(spherogram.Link([list(x) for x in pd_code(N2, U2)]).exterior(), E_K)
    B = add_one_band(spherogram.Link([list(x) for x in pd]), ([h1, nbr[h3]], [], -t))
    sb = dict(num_twist=-t, components=len(B.link_components))
    if sb['components'] == 1:
        sb['alexander_seifert'] = alexander_seifert(B.seifert_matrix())
        sb['isometric_to_J'] = isometric(B.exterior(), E_J) if hyperbolic else None
    rec['spherogram_add_one_band'] = sb
    rec['gamma4_knotinfo'], rec['smooth_four_genus_knotinfo'] = ref['g4'], ref['smooth_four_genus']
    # the surface of J: a Mobius band (gamma_4(J) = 1) or a punctured Klein bottle (gamma_4(J) = 2), found by us
    field, find = (('moebius_band_of_J', moebius_band) if ref['g4'] == '1' else ('klein_bottle_of_J', klein_bottle_of))
    surf = None
    for where, Q in (('J made by the band', P), ('KnotInfo diagram of J', ref['pd'])):
        if surface is not None and where != 'J made by the band':
            surf = dict(surface, diagram_isometric_to_J=isometric(spherogram.Link(surface['diagram_pd']).exterior(), E_J)
                        if hyperbolic else None)
            break
        surf = find(Q)
        if surf:
            surf['diagram'] = where
            break
    rec[field] = surf
    rec['self_contained'] = bool(surf) and surf['diagram'] == 'J made by the band'
    rec['contradiction'] = hit['contradiction']
    ident = rec['isometric_to_knotinfo'] is True or (not hyperbolic and hit['how'] == 'diagram')
    rebuilt = sb.get('isometric_to_J') is True or (not hyperbolic and sb.get('alexander_seifert') == ref['alex'])
    rec['passes'] = (rec['components_traversal'] == rec['components_spherogram'] == 1
                     and rec['alexander_fox'] == rec['alexander_seifert'] == ref['alex'] == hit['alexander']
                     and bool(rec['hfk_equals_knotinfo']) and ident and rec['dual_band_gives_K'] is True and rebuilt)
    return rec


LEVELS = ('self-contained', 'self-contained up to Reidemeister moves of J', 'surface of J in another diagram of J',
          'KnotInfo value of gamma_4(J)')


def verification(c):
    """How the surface bounded by J is known: by our bands in the diagram of J made by the band (K -> J -> unknot,
    self-contained), by our bands in a diagram reached from it by Reidemeister moves (stage C), by our bands in
    another diagram of J (linked to J by the isometry check), or only by KnotInfo."""
    s = c.get('moebius_band_of_J') or c.get('klein_bottle_of_J')
    return (LEVELS[0] if c['self_contained'] else LEVELS[1] if s and s['diagram'] == REACHED else
            LEVELS[2] if s else LEVELS[3])


def band_pd_darts(pd, face, darts, t):
    nbr, under, label = from_pd([tuple(x) for x in pd])
    f = faces(nbr)[face]
    return pd_code(*band(nbr, under, f[darts[0]], f[darts[1]], t))


def preference(hit, diags):
    """KnotInfo diagram first, then small |t|, small J, small diagram."""
    return (hit['diagram'] > 0, abs(hit['t']), REF[hit['J']]['crossings'], len(diags[hit['diagram']]['pd']) + abs(hit['t']),
            hit['diagram'], hit['face'], hit['darts'])


def scan_diagram(d, n, allowed, group):
    st, cands, fit_fail = scan(d['pd'], n, allowed)
    ident, hits, g4_seen, contra = Counter(), [], Counter(), []
    for c in cands:
        cl = classify(c, allowed, group)
        ident[cl['how'] if cl['names'] or cl['how'] in ('unknot', 'unidentified') else 'outside_knotinfo'] += 1
        ident['alexander_mismatch'] += cl['alexander_mismatch']
        for k in cl['names']:
            g4_seen[G4.get(k)] += 1
        if cl['J'] or (cl['unknot'] and group not in ('upper', 'mobius')):   # groups upper, mobius: see moebius_band
            h = dict(face=c['face'], darts=c['darts'], edges=c['edges'], t=c['t'], J=cl['J'] or 'unknot', how=cl['how'],
                     alexander=c['alexander'], contradiction=cl['contradiction'])
            (contra if cl['contradiction'] else hits).append(h)
    d.update({k: st[k] for k in ('arcs', 'bands', 'knots', 'det_pass', 'alexander_pass')}, hits=len(hits),
             contradictions=len(contra))
    return st, hits, contra, fit_fail, ident, g4_seen


def one(item):
    """One knot: stage A, stage B if needed, rechecks.  In replay (stored diagrams given) the scan only."""
    name, info, pd0, stored, carried = item
    import spherogram
    t0 = time.time()
    n, allowed, group = info['det'], tuple(info['J_gamma4']), info['group']
    E_K = spherogram.Link([list(x) for x in pd0]).exterior()
    diags = ([dict(source='knotinfo', pd=[list(x) for x in pd0])] if stored is None
             else [dict(source=d['source'], pd=d['pd'], **({'carried': True} if d.get('carried') else {})) for d in stored])
    verify_diagrams(diags, pd0)
    tot, ident, g4_seen = Counter(), Counter(), Counter()
    hits, contra, fit_fail = [], [], []

    def run(first):
        for di in range(first, len(diags)):
            d = diags[di]
            assert det_link(d['pd']) == n
            st, h, c, ff, idc, g4 = scan_diagram(d, n, allowed, group)
            tot.update(st); ident.update(idc); g4_seen.update(g4)
            for x in h + c:
                x['diagram'] = di
            hits.extend(h); contra.extend(c); fit_fail.extend([di] + x for x in ff)
    run(0)
    if stored is not None:
        return name, dict(diagrams=diags, hits=hits, contradiction_hits=contra, fit_failures=fit_fail,
                          identification=dict(ident))
    certs = []
    good = lambda: [c for c in certs if c['passes'] and not c['contradiction']]
    self_contained = lambda: any(c['self_contained'] for c in good())
    # group upper: J with gamma_4 '1' (b_1 = 2) before J with gamma_4 '2' (b_1 = 3); the other groups keep their order
    b1_three = lambda h: group == 'upper' and REF[h['J']]['g4'] != '1'
    mobius = []                                   # group upper: a Mobius band of K by our own bands (gamma_4 = 1)
    tried_unknotting = set()                      # hits whose diagram of J was searched for an unknotting band

    def mobius_search(first):
        for di in range(first, len(diags)):
            if mobius:
                return
            mb = moebius_band(diags[di]['pd'])
            if mb:
                mobius.append(dict(mb, diagram=di, diagram_source=diags[di]['source'], diagram_pd=diags[di]['pd']))

    def add(h, surface=None):
        rec = dict(h, diagram_source=diags[h['diagram']]['source'], diagram_pd=diags[h['diagram']]['pd'])
        rec.update(recheck(diags[h['diagram']]['pd'], h, E_K, surface))
        certs.append(rec)
        tried_unknotting.add((h['diagram'], h['face'], tuple(h['darts']), h['t']))

    def certify_reached():
        """Stage C (mode ranges, targets): first an unknotting band of the diagram of J made by the band for every hit
        not tried yet (J with gamma_4 '1'); then hits whose J bounds a Mobius band in a diagram reached from the
        diagram of J made by the band by Reidemeister moves (moebius_band_reached), the knots J in turn, at most
        C_HITS hits, until one passes the rechecks."""
        for h in sorted(hits, key=lambda h: preference(h, diags)):
            if self_contained():
                return
            key = (h['diagram'], h['face'], tuple(h['darts']), h['t'])
            if REF[h['J']]['g4'] == '1' and key not in tried_unknotting and \
                    unknotting_band(band_pd_darts(diags[h['diagram']]['pd'], h['face'], h['darts'], h['t'])):
                add(h)
        if self_contained():
            return
        rank, order = Counter(), []
        for h in sorted(hits, key=lambda h: preference(h, diags)):
            if REF[h['J']]['g4'] == '1':
                order.append(((rank[h['J']], preference(h, diags)), h)); rank[h['J']] += 1
        for i, (_, h) in enumerate(sorted(order, key=lambda x: x[0])[:C_HITS]):
            P = band_pd_darts(diags[h['diagram']]['pd'], h['face'], h['darts'], h['t'])
            surf = moebius_band_reached(P, zlib.crc32(name.encode()) + 1000 * i)
            if surf:
                add(h, surf)
                if any(c.get('moebius_band_of_J') and c['moebius_band_of_J']['diagram'] in ('J made by the band', REACHED)
                       for c in good()):
                    return

    def certify():
        """Hits with distinct J first: REPS of them, more while none is self-contained, at most MAX_RECHECK in all.
        Then, while none is self-contained, further hits in the same order whose diagram of J (made by the band) bounds
        a surface by our own bands: an unknotting band (SC_TRIES hits; klein_bottle_of when gamma_4(J) = 2, REPS
        hits), then a band to a ribbon knot (RB_HITS hits).  Those found are rechecked."""
        ordered = sorted(hits, key=lambda h: (b1_three(h), preference(h, diags)))
        done = {c['J'] for c in certs}
        for h in ordered:
            if len(certs) >= MAX_RECHECK or (len(good()) >= REPS and self_contained()):
                break
            if h['J'] not in done:
                done.add(h['J']); add(h)
        sc_budget = SC_TRIES if allowed == ('1',) or group == 'upper' else REPS
        for find, budget in ((unknotting_band, sc_budget), (ribbon_moebius_band, RB_HITS)):
            seen = {(c['diagram'], c['face'], tuple(c['darts']), c['t']) for c in certs}
            tried = 0
            for h in ordered:
                if self_contained() or tried >= budget:
                    break
                if (h['diagram'], h['face'], tuple(h['darts']), h['t']) in seen:
                    continue
                if find is ribbon_moebius_band and REF[h['J']]['g4'] != '1':
                    continue
                tried += 1
                if find is unknotting_band:
                    tried_unknotting.add((h['diagram'], h['face'], tuple(h['darts']), h['t']))
                P = band_pd_darts(diags[h['diagram']]['pd'], h['face'], h['darts'], h['t'])
                if (find if REF[h['J']]['g4'] == '1' else klein_bottle_of)(P):
                    add(h)
    stage = 'A'
    certify()
    if group == 'upper':
        mobius_search(0)
    if group == 'control' or (not mobius if group == 'upper' else not self_contained()):  # controls: always stage B
        stage = 'B'
        more = diagrams_of(name, pd0, zlib.crc32(name.encode()), carried)[1:]
        verify_diagrams(more, pd0)
        k = len(diags)
        diags.extend(more)
        run(k)
        certify()
        if group == 'upper':
            mobius_search(k)
    if info.get('stage_C') and not self_contained():
        stage = 'C'
        certify_reached()
    if good() and not any(c.get('moebius_band_of_J') or c.get('klein_bottle_of_J') for c in good()):
        for c in sorted(good(), key=lambda c: preference(c, diags))[:REPS]:   # last resort: other diagrams of J
            if REF[c['J']]['g4'] == '1':
                c['moebius_band_of_J'] = moebius_band_elsewhere(c['J'])
                if c['moebius_band_of_J']:
                    break
    for h in contra[:MAX_RECHECK]:                # a contradiction is rechecked like a hit
        if h['J'] != 'unknot':
            add(h)
    for c in certs:
        c['verification'] = verification(c)
    rep = min(good(), key=lambda c: (LEVELS.index(c['verification']), b1_three(c), preference(c, diags)), default=None)
    J = Counter(h['J'] for h in hits)
    extra = dict(moebius_band_of_K=mobius[0] if mobius else None) if group == 'upper' else {}
    return name, dict(info, **extra, stage=stage, diagrams=diags, arcs=tot['arcs'], bands_tried=tot['bands'],
                      bands_split=tot['split'], knot_results=tot['knots'], det_pass=tot['det_pass'],
                      alexander_pass=tot['alexander_pass'], identification=dict(ident), J_gamma4_seen=dict(g4_seen),
                      fox_det_checks=[tot['fox_det_checks'], tot['fox_det_checks'] - tot['fox_det_failures']],
                      fit_failures=fit_fail, hits=hits, contradiction_hits=contra, J_found=dict(J.most_common()),
                      certificates=certs, certified=rep is not None, self_contained=self_contained(),
                      verification=None if rep is None else rep['verification'],
                      representative=None if rep is None else {k: rep[k] for k in (
                          'diagram', 'diagram_source', 'face', 'face_edges', 'edges', 't', 'J', 'J_pd',
                          'J_crossings_simplified', 'self_contained', 'verification', 'moebius_band_of_J',
                          'klein_bottle_of_J') if k in rep},
                      contradictions=len(contra), seconds=round(time.time() - t0, 1))


def mobius_init_worker(ref, dt, g4):
    """init_worker, and a PARI stack that may grow to PARI_SIZEMAX (mode mobius: the default maximum, 8 MB, is too
    small for some knots)."""
    init_worker(ref, dt, g4)
    pari().allocatemem(PARI_STACK, PARI_SIZEMAX, silent=True)


def mobius_task(item):
    """mobius_one in a worker; an exception becomes the record {error: traceback} instead of stopping the run (the main
    process then writes no output; KLEIN_RESUME=1 runs those knots again)."""
    try:
        return mobius_one(item)
    except Exception:
        import traceback
        return item[0], dict(error=traceback.format_exc())


def klein_task(item):
    """one() in a worker of mobius_replay, with the errors of mobius_task."""
    try:
        return one(item)
    except Exception:
        import traceback
        return item[0], dict(error=traceback.format_exc())


def mobius_one(item):
    """Mode mobius, one knot (see the docstring): Mobius bands of K in the KnotInfo diagram, then in the diagrams of
    stage B (group mobius until one is found, controls all); for a knot of group mobius with range [1,3] and none, the
    punctured Klein bottle search of one().  In replay (stored diagrams given) the search of the stored diagrams only."""
    name, info, pd0, stored, carried = item
    t0 = time.time()
    n, control = info['det'], info['group'] == 'control'
    diags = ([dict(source='knotinfo', pd=[list(x) for x in pd0])] if stored is None
             else [dict(source=d['source'], pd=d['pd'], **({'carried': True} if d.get('carried') else {})) for d in stored])
    verify_diagrams(diags, pd0)
    bands, tot, fit_fail = [], Counter(), []

    def search(first):
        for di in range(first, len(diags)):
            if bands and not control and stored is None:
                del diags[di:]                    # not searched, not stored
                return
            d = diags[di]
            assert det_link(d['pd']) == n
            mb, st, ff = moebius_band_fast(d['pd'], n)
            tot.update(st); fit_fail.extend([di] + x for x in ff)
            d['moebius_band'] = mb['kind'] if mb else None
            if mb:
                bands.append(dict(mb, diagram=di, diagram_source=d['source'], diagram_pd=d['pd']))
    search(0)
    if stored is not None:
        return name, dict(diagrams=diags, moebius_bands=bands, fit_failures=fit_fail)
    stage = 'A'
    if control or not bands:
        stage = 'B'
        more = diagrams_of(name, pd0, zlib.crc32(name.encode()), carried)[1:]
        verify_diagrams(more, pd0)
        k = len(diags)
        diags.extend(more)
        search(k)
    kb = None
    if info['group'] == 'mobius' and info['range'] == '[1,3]' and not bands:
        kb = one((name, dict(info, stage_C=True), pd0, None, diags[1:]))[1]
        kb = {k: v for k, v in kb.items() if k not in info and k != 'stage_C'}
    return name, dict(info, stage=stage, diagrams=diags, moebius_band_of_K=bands[0] if bands else None,
                      contradiction_bands=bands if control else [], contradictions=len(bands) if control else 0,
                      search_counts=dict(tot), fit_failures=fit_fail,
                      certificates=kb.pop('certificates') if kb else [],
                      representative=kb.pop('representative') if kb else None,
                      certified=bool(kb and kb['certified']), self_contained=bool(kb and kb['self_contained']),
                      verification=kb and kb['verification'], klein_bottle=kb, seconds=round(time.time() - t0, 1))


def validate(item):
    """Fox polynomial of the KnotInfo diagram of a reference knot against its KnotInfo polynomial."""
    name, pd, alex, det = item
    A = alexander_fox([tuple(x) for x in pd])
    return name, A == alex and abs(sum(x * (-1) ** e for e, x in enumerate(A))) == det


def snappy_hits(hs):
    """The hits identified by SnapPy, as comparable keys.  An identification by the simplified diagram (non-hyperbolic
    J, the unknot) depends on spherogram's random global simplification, so those hits are only counted."""
    return sorted((h['diagram'], h['face'], tuple(h['darts']), h['t'], h['J']) for h in hs if h['how'] == 'snappy')


def compare_scan(o, new):
    """A scan of one() replayed, against its stored record o: whether the counts of every diagram, the hits identified
    by SnapPy and the fit failures agree, and the numbers of the other hits, [stored, now]."""
    same = (len(o['diagrams']) == len(new['diagrams'])
            and all(all(a[k] == b[k] for k in ('arcs', 'bands', 'knots', 'det_pass', 'alexander_pass',
                                               'isometric_to_knotinfo'))
                    for a, b in zip(o['diagrams'], new['diagrams']))
            and snappy_hits(o['hits']) == snappy_hits(new['hits'])
            and snappy_hits(o['contradiction_hits']) == snappy_hits(new['contradiction_hits'])
            and o['fit_failures'] == new['fit_failures'])
    other = lambda hs: sum(h['how'] != 'snappy' for h in hs)
    return same, [other(o['hits'] + o['contradiction_hits']), other(new['hits'] + new['contradiction_hits'])]


def rebuild_certificates(old, BY):
    """The J of every stored certificate rebuilt in its diagram of K, with the isometry check against the KnotInfo
    exterior of its name, and its surface of J (surface_check): tuples (knot, J, same PD code as stored, same isometry
    check as stored, surface of J rebuilt or None if there is none)."""
    import snappy, spherogram
    rebuilt = []
    for name, o in old.items():
        for c in o['certificates']:
            P = [list(x) for x in band_pd_darts(c['diagram_pd'], c['face'], c['darts'], c['t'])]
            E = spherogram.Link(P).exterior()
            iso = (isometric(E, spherogram.Link(json.loads(BY[c['J']]['pd_notation'])).exterior())
                   if is_hyperbolic(E) else None)
            rebuilt.append((name, c['J'], P == c['J_pd'], iso == c['isometric_to_knotinfo'], surface_check(c, BY)))
    return rebuilt


def print_rebuilt(rebuilt):
    """The counts of rebuild_certificates, and the certificates that fail."""
    print('certificates rebuilt:', len(rebuilt), 'with a different J or isometry check:',
          [x[:4] for x in rebuilt if not (x[2] and x[3])])
    print('surfaces of J rebuilt:', sum(x[4] is not None for x in rebuilt), 'failing:',
          [x[:2] for x in rebuilt if x[4] is False])


def replay(C, BY, workers, init, path):
    """Default mode and mode ranges: the stored diagrams scanned again and compared with the stored output, and the J
    and the surface of J of every stored certificate rebuilt (see the docstring); writes nothing."""
    old = json.load(open(path))['knots']
    items = [(k, v, json.loads(BY[k]['pd_notation']), old[k]['diagrams'], ()) for k, v in C.items() if k in old]
    bad, varying = [], []
    with Pool(workers, initializer=init_worker, initargs=init) as p:
        for name, new in p.imap_unordered(one, items):
            same, counts = compare_scan(old[name], new)
            print(name, 'diagrams', len(new['diagrams']), 'hits', len(new['hits']), 'same as stored:', same,
                  'hits identified by diagram (stored, now):', counts, flush=True)
            if not same:
                bad.append(name)
            if counts[0] != counts[1]:
                varying.append([name] + counts)
    rebuilt = rebuild_certificates(old, BY)
    print('replay: knots', len(items), 'differing from the stored output:', bad)
    print('knots whose number of hits identified by the simplified diagram changed:', varying)
    print_rebuilt(rebuilt)
    mob = {name: mobius_check(o['moebius_band_of_K']) for name, o in old.items() if o.get('moebius_band_of_K')}
    if any(o['group'] == 'upper' for o in old.values()):
        print('Mobius bands of K rebuilt:', len(mob), 'failing:', [k for k, v in mob.items() if not v])


def original_search(item):
    """moebius_band in the KnotInfo diagram of a knot, against which mobius_replay checks moebius_band_fast."""
    name, pd = item
    return name, moebius_band(pd)


def mobius_replay(C, BY, workers, init, path):
    """Mode mobius: the Mobius band search of the stored diagrams again (per diagram the kind of band found, or none, and
    the first band), the punctured Klein bottle searches as in replay, the J and the surface of J of every stored
    certificate, and every stored Mobius band of K rebuilt (mobius_check); writes nothing.  Also moebius_band in the
    KnotInfo diagram of every knot, which must find the band of moebius_band_fast there (none if that has |s| > S_MAX)."""
    old = json.load(open(path))['knots']
    pd = lambda k: json.loads(BY[k]['pd_notation'])
    items = [(k, v, pd(k), old[k]['diagrams'], ()) for k, v in C.items() if k in old]
    kitems = [(k, dict(v, stage_C=True), pd(k), old[k]['klein_bottle']['diagrams'], ()) for k, v in C.items()
              if k in old and old[k]['klein_bottle']]
    band_key = lambda m: m and (m['diagram'], m['face'], tuple(m['darts']), m['s'], m['kind'],
                                m.get('ribbon_band_of_J1') and (m['ribbon_band_of_J1']['face'],
                                                                tuple(m['ribbon_band_of_J1']['darts']),
                                                                m['ribbon_band_of_J1']['r']))
    bad, kbad, varying, fast0 = [], [], [], {}
    with Pool(workers, initializer=mobius_init_worker, initargs=init) as p:
        for name, new in p.imap_unordered(mobius_task, items):
            o = old[name]
            if 'error' in new:
                print(name, 'ERROR', new['error'], flush=True); bad.append(name); fast0[name] = None
                continue
            fast0[name] = next((m for m in new['moebius_bands'] if m['diagram'] == 0), None)
            same = (len(o['diagrams']) == len(new['diagrams'])
                    and all(a[k] == b[k] for a, b in zip(o['diagrams'], new['diagrams'])
                            for k in ('moebius_band', 'isometric_to_knotinfo'))
                    and band_key(o['moebius_band_of_K']) == band_key(new['moebius_bands'][0] if new['moebius_bands']
                                                                     else None)
                    and len(o['contradiction_bands']) == (len(new['moebius_bands']) if o['group'] == 'control' else 0)
                    and o['fit_failures'] == new['fit_failures'])
            print(name, o['group'], o['range'], 'diagrams', len(new['diagrams']), 'Mobius bands', len(new['moebius_bands']),
                  'same as stored:', same, flush=True)
            if not same:
                bad.append(name)
        for name, new in p.imap_unordered(klein_task, kitems):
            if 'error' in new:
                print(name, 'ERROR', new['error'], flush=True); kbad.append(name)
                continue
            same, counts = compare_scan(old[name]['klein_bottle'], new)
            print(name, 'punctured Klein bottle search: diagrams', len(new['diagrams']), 'hits', len(new['hits']),
                  'same as stored:', same, 'hits identified by diagram (stored, now):', counts, flush=True)
            if not same:
                kbad.append(name)
            if counts[0] != counts[1]:
                varying.append([name] + counts)
        orig = dict(p.imap_unordered(original_search, [(k, pd(k)) for k in fast0], chunksize=4))
    near = lambda m: m if m and abs(m['s']) <= S_MAX else None
    differ = [k for k in fast0 if band_key(dict(orig[k], diagram=0) if orig[k] else None) != band_key(near(fast0[k]))]
    rebuilt = rebuild_certificates(old, BY)
    mob = {name: mobius_check(o['moebius_band_of_K']) for name, o in old.items() if o['moebius_band_of_K']}
    print('replay: knots', len(items), 'Mobius band search differing from the stored output:', bad)
    print('punctured Klein bottle searches', len(kitems), 'differing from the stored output:', kbad)
    print('knots whose number of hits identified by the simplified diagram changed:', varying)
    print_rebuilt(rebuilt)
    print('Mobius bands of K rebuilt:', len(mob), 'failing:', [k for k, v in mob.items() if not v])
    print('moebius_band in the KnotInfo diagrams:', len(orig), 'bands', sum(bool(v) for v in orig.values()),
          'differing from moebius_band_fast:', differ)


def mobius_check(m):
    """Rebuild a stored Mobius band of K (group upper) in its diagram of K: the unknot, or J' and the two-component
    unlink, recognised again."""
    if m['kind'] == 'band to the unknot':
        return bool(unknot_check(band_pd_darts(m['diagram_pd'], m['face'], m['darts'], m['s'])))
    P1 = band_pd_darts(m['diagram_pd'], m['face'], m['darts'], m['s'])
    r = m['ribbon_band_of_J1']
    return [list(x) for x in P1] == m['J1_pd'] and is_unlink2(band_pd_darts(P1, r['face'], r['darts'], r['r']))


def surface_check(c, BY):
    """Rebuild the stored surface of J (unknotting band; band to a ribbon knot and its ribbon band; or bands
    J -> J' -> unknot) and recognise the unknot or the two-component unlink again."""
    s = c.get('moebius_band_of_J') or c.get('klein_bottle_of_J')
    if not s:
        return None
    pd = {'J made by the band': c['J_pd'], 'KnotInfo diagram of J': json.loads(BY[c['J']]['pd_notation']),
          'another diagram of J': s.get('diagram_pd'), REACHED: s.get('diagram_pd')}[s['diagram']]
    if 'unknotting_band_of_J1' in s:
        P1 = band_pd_darts(pd, s['face'], s['darts'], s['t'])
        u = s['unknotting_band_of_J1']
        return [list(x) for x in P1] == s['J1_pd'] and bool(unknot_check(band_pd_darts(P1, u['face'], u['darts'], u['s'])))
    if 'ribbon_band_of_J1' in s:
        P1 = band_pd_darts(pd, s['face'], s['darts'], s['s'])
        r = s['ribbon_band_of_J1']
        return [list(x) for x in P1] == s['J1_pd'] and is_unlink2(band_pd_darts(P1, r['face'], r['darts'], r['r']))
    return bool(unknot_check(band_pd_darts(pd, s['face'], s['darts'], s['s'])))


def add_face_edges(out, BY):
    """Mode ranges: the edges of the face of the first band of klein_bottle_of (J -> J'), as for the other bands."""
    for r in out.values():
        for c in r['certificates'] + ([r['representative']] if r['representative'] else []):
            s = c.get('klein_bottle_of_J')
            if s and 'face_edges' not in s:
                pd = c['J_pd'] if s['diagram'] == 'J made by the band' else json.loads(BY[c['J']]['pd_notation'])
                nbr, under, label = from_pd([tuple(x) for x in pd])
                s['face_edges'] = [label[h] for h in faces(nbr)[s['face']]]


def default_summary(out, vres):
    """Default mode: the summary of the targets, the controls and the checks."""
    T = {k: r for k, r in out.items() if r['group'] == 'target'}
    rep = {k: r['representative'] for k, r in T.items() if r['certified']}
    return dict(
        targets=len(T), certified=len(rep), self_contained=sum(r['self_contained'] for r in T.values()),
        not_certified=[k for k, r in T.items() if not r['certified']],
        verification=dict(Counter(r['verification'] for r in T.values() if r['certified'])),
        not_self_contained={k: r['verification'] for k, r in T.items() if r['certified'] and not r['self_contained']},
        certified_by_test={t: [sum(r['certified'] for r in T.values() if r['test'] == t),
                               sum(r['test'] == t for r in T.values())] for t in ('linking form', 'd-invariant')},
        certified_by_crossings={c: [sum(r['certified'] for k, r in T.items() if knot_order(k)[0] == c),
                                    sum(knot_order(k)[0] == c for k in T)] for c in (11, 12, 13)},
        stage_B=[k for k, r in T.items() if r['stage'] == 'B'],
        representative_t=dict(sorted(Counter(abs(r['t']) for r in rep.values()).items())),
        representative_from_knotinfo_diagram=sum(r['diagram_source'] == 'knotinfo' for r in rep.values()),
        representative_J=dict(Counter(r['J'] for r in rep.values()).most_common()),
        controls={k: dict(knotinfo_gamma4=r['knotinfo_gamma4'], proposition_5_1b=r['proposition_5_1b'],
                          diagrams=len(r['diagrams']), knot_results=r['knot_results'], hits=len(r['hits']),
                          contradictions=r['contradictions']) for k, r in out.items() if r['group'] == 'control'},
        knot_13a_2890=dict(certified=out['13a_2890']['certified'] if '13a_2890' in out else None,
                           representative=out['13a_2890']['representative'] if '13a_2890' in out else None),
        contradictions=sum(r['contradictions'] for r in out.values()),
        fit_failures=sum(len(r['fit_failures']) for r in out.values()),
        fox_det_checks=[sum(r['fox_det_checks'][0] for r in out.values()), sum(r['fox_det_checks'][1] for r in out.values())],
        alexander_mismatches=sum(r['identification'].get('alexander_mismatch', 0) for r in out.values()),
        diagrams=sum(len(r['diagrams']) for r in out.values()),
        diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for r in out.values() for d in r['diagrams']),
        bands_tried=sum(r['bands_tried'] for r in out.values()),
        knot_results=sum(r['knot_results'] for r in out.values()),
        candidates=sum(r['alexander_pass'] for r in out.values()),
        hits=sum(len(r['hits']) for r in out.values()),
        certificates_rechecked=sum(len(r['certificates']) for r in out.values()),
        certificates_failing_a_recheck=sum(not c['passes'] for r in out.values() for c in r['certificates']),
        reference_fox_validation=[sum(vres.values()), len(vres)])


def ranges_summary(out, vres, ref):
    """Mode ranges: the summary of the targets, of 11n_133, of the group upper and of the checks."""
    T = {k: r for k, r in out.items() if r['group'] == 'target'}
    rep = {k: r['representative'] for k, r in T.items() if r['certified']}
    sc = {k for k, r in T.items() if r['self_contained']}
    kind = lambda r: (r['moebius_band_of_J']['kind'] if r.get('moebius_band_of_J') else
                      'band J -> J\' -> unknot' if r.get('klein_bottle_of_J') else None)
    b1 = lambda r: None if not r['self_contained'] else 2 if ref[r['representative']['J']]['g4'] == '1' else 3
    x = out.get('11n_133')
    return dict(
        targets=len(T), certified=len(rep), self_contained=len(sc),
        self_contained_by_range={g: [sum(k in sc for k, r in T.items() if r['range'] == g),
                                     sum(r['range'] == g for r in T.values())] for g in ('[2,3]', '[2,4]')},
        not_self_contained=[k for k in T if k not in sc],
        not_certified=[k for k, r in T.items() if not r['certified']],
        verification=dict(Counter(r['verification'] for r in T.values() if r['certified'])),
        self_contained_by_crossings={c: [sum(k in sc for k in T if knot_order(k)[0] == c),
                                         sum(knot_order(k)[0] == c for k in T)] for c in (11, 12, 13)},
        stage_B=[k for k, r in T.items() if r['stage'] in ('B', 'C')],
        stage_B_self_contained=[k for k, r in T.items() if r['stage'] == 'B' and r['self_contained']],
        stage_C={k: r['verification'] for k, r in T.items() if r['stage'] == 'C'},
        representative_t=dict(sorted(Counter(abs(r['t']) for r in rep.values()).items())),
        representative_from_knotinfo_diagram=sum(r['diagram_source'] == 'knotinfo' for r in rep.values()),
        representative_surface_of_J=dict(Counter(kind(r) for r in rep.values())),
        representative_J=dict(Counter(r['J'] for r in rep.values()).most_common()),
        knot_11n_133=None if x is None else dict(
            certified=x['certified'], self_contained=x['self_contained'], b1=b1(x), verification=x['verification'],
            representative=x['representative'], diagrams=len(x['diagrams']), hits=len(x['hits']), J_found=x['J_found'],
            J_gamma4_seen=x['J_gamma4_seen'], contradictions=x['contradictions']),
        upper={k: dict(range=r['range'], h_witt=r['h_witt'], stage=r['stage'], diagrams=len(r['diagrams']),
                       moebius_band_of_K=r['moebius_band_of_K'] and dict(
                           kind=r['moebius_band_of_K']['kind'], diagram=r['moebius_band_of_K']['diagram'],
                           diagram_source=r['moebius_band_of_K']['diagram_source']),
                       certified=r['certified'], self_contained=r['self_contained'], b1=b1(r),
                       J=r['representative']['J'] if r['representative'] else None,
                       upper_bound=1 if r['moebius_band_of_K'] else b1(r), hits=len(r['hits']),
                       J_gamma4_seen=r['J_gamma4_seen'])
               for k, r in out.items() if r['group'] == 'upper'},
        contradictions=sum(r['contradictions'] for r in out.values()),
        fit_failures=sum(len(r['fit_failures']) for r in out.values()),
        fox_det_checks=[sum(r['fox_det_checks'][0] for r in out.values()), sum(r['fox_det_checks'][1] for r in out.values())],
        alexander_mismatches=sum(r['identification'].get('alexander_mismatch', 0) for r in out.values()),
        diagrams=sum(len(r['diagrams']) for r in out.values()),
        diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for r in out.values() for d in r['diagrams']),
        bands_tried=sum(r['bands_tried'] for r in out.values()),
        knot_results=sum(r['knot_results'] for r in out.values()),
        candidates=sum(r['alexander_pass'] for r in out.values()),
        hits=sum(len(r['hits']) for r in out.values()),
        certificates_rechecked=sum(len(r['certificates']) for r in out.values()),
        certificates_failing_a_recheck=sum(not c['passes'] for r in out.values() for c in r['certificates']),
        reference_fox_validation=[sum(vres.values()), len(vres)])


def mobius_summary(out, vres):
    """Mode mobius: the summary of the Mobius bands, of the punctured Klein bottles, of the controls and of the checks."""
    M = {k: r for k, r in out.items() if r['group'] == 'mobius'}
    mob = {k: r['moebius_band_of_K'] for k, r in M.items() if r['moebius_band_of_K']}
    kb = {k: r for k, r in M.items() if r['klein_bottle']}
    by_range = lambda keys: {g: [sum(M[k]['range'] == g for k in keys), sum(r['range'] == g for r in M.values())]
                             for g in ('[1,2]', '[1,3]')}
    ctl = {g: {k: r for k, r in out.items() if r['group'] == 'control' and r['range'] == g} for g in ('[1,2]', '[1,3]')}
    rep = {k: r['representative'] for k, r in kb.items() if r['certified']}
    return dict(
        knots=len(out), group_mobius=len(M), moebius_band_of_K=by_range(mob),
        moebius_band_of_K_stage={g: dict(Counter(M[k]['stage'] for k in mob if M[k]['range'] == g))
                                 for g in ('[1,2]', '[1,3]')},
        moebius_band_kind=dict(Counter(m['kind'] for m in mob.values())),
        moebius_band_diagram_source=dict(Counter(m['diagram_source'] for m in mob.values())),
        moebius_band_s=dict(sorted(Counter(abs(m['s']) for m in mob.values()).items())),
        moebius_band_r=dict(sorted(Counter(abs(m['ribbon_band_of_J1']['r']) for m in mob.values()
                                           if 'ribbon_band_of_J1' in m).items())),
        moebius_band_J1_det=dict(sorted(Counter(m['J1_det'] for m in mob.values() if 'J1_det' in m).items())),
        moebius_band_far=[k for k, m in mob.items() if abs(m['s']) > S_MAX],
        moebius_band_stage_B=[k for k, r in M.items() if k in mob and r['stage'] == 'B'],
        no_moebius_band={g: [k for k, r in M.items() if r['range'] == g and k not in mob] for g in ('[1,2]', '[1,3]')},
        klein_bottle=dict(
            knots=len(kb), certified=len(rep), self_contained=[sum(r['self_contained'] for r in kb.values()), len(kb)],
            not_self_contained=[k for k, r in kb.items() if not r['self_contained']],
            stage=dict(Counter(r['klein_bottle']['stage'] for r in kb.values())),
            verification=dict(Counter(r['verification'] for r in kb.values() if r['certified'])),
            representative_t=dict(sorted(Counter(abs(r['t']) for r in rep.values()).items())),
            representative_from_knotinfo_diagram=sum(r['diagram_source'] == 'knotinfo' for r in rep.values()),
            representative_J=dict(Counter(r['J'] for r in rep.values()).most_common()),
            certificates_rechecked=sum(len(r['certificates']) for r in kb.values()),
            certificates_failing_a_recheck=sum(not c['passes'] for r in kb.values() for c in r['certificates']),
            contradictions=sum(r['klein_bottle']['contradictions'] for r in kb.values()),
            fit_failures=sum(len(r['klein_bottle']['fit_failures']) for r in kb.values()),
            fox_det_checks=[sum(r['klein_bottle']['fox_det_checks'][0] for r in kb.values()),
                            sum(r['klein_bottle']['fox_det_checks'][1] for r in kb.values())],
            alexander_mismatches=sum(r['klein_bottle']['identification'].get('alexander_mismatch', 0)
                                     for r in kb.values()),
            diagrams=sum(len(r['klein_bottle']['diagrams']) for r in kb.values()),
            diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for r in kb.values()
                                       for d in r['klein_bottle']['diagrams'])),
        range_one_three_to_one_two=[k for k, r in kb.items() if r['self_contained']],
        controls={g: dict(knots=len(C), diagrams=sum(len(r['diagrams']) for r in C.values()),
                          moebius_bands=sum(r['contradictions'] for r in C.values())) for g, C in ctl.items()},
        contradictions=sum(r['contradictions'] for r in out.values()),
        contradiction_knots=[k for k, r in out.items() if r['contradictions']],
        diagrams=sum(len(r['diagrams']) for r in out.values()),
        diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for r in out.values() for d in r['diagrams']),
        search_counts=dict(sum((Counter(r['search_counts']) for r in out.values()), Counter())),
        fit_failures=sum(len(r['fit_failures']) for r in out.values()),
        reference_fox_validation=[sum(vres.values()), len(vres)])


def versions():
    """The versions of the programs, for the field meta of the output."""
    import snappy, spherogram, cypari, database_knotinfo
    return dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                database_knotinfo=database_knotinfo.version(), cypari=cypari.__version__,
                pari='.'.join(map(str, pari().version())), python=sys.version.split()[0])


if __name__ == '__main__':
    from database_knotinfo import link_list
    args = sys.argv[1:]
    mode_ranges, mode_mobius = args[:1] == ['ranges'], args[:1] == ['mobius']
    args = args[1:] if mode_ranges or mode_mobius else args
    mode_replay = args[:1] == ['replay']
    args = args[1:] if mode_replay else args
    workers = int(args[0]) if args else 10
    path = os.environ.get('KLEIN_OUT', os.path.join(ROOT, 'data', 'klein_bottle_search_ranges.json' if mode_ranges
                                                    else 'klein_bottle_search_mobius.json' if mode_mobius
                                                    else 'klein_bottle_search.json'))
    t_start = time.time()
    LL = link_list()[1:]
    BY = {r['name']: r for r in LL}
    G, LG = data('gamma4_knotinfo_2026.10.5.json'), data('linking_genus.json')
    C = ranges_knots() if mode_ranges else mobius_knots() if mode_mobius else knots()
    ref = references(BY, G, LG)
    dt = {'K' + r['dt_name'].replace('_', ''): r['name'] for r in LL if r.get('dt_name')}
    init = (ref, dt, G)
    if mode_replay:
        (mobius_replay if mode_mobius else replay)(C, BY, workers, init, path)
        sys.exit()
    only = os.environ.get('KLEIN_ONLY')          # comma-separated knots, for tests
    if only:
        C = {k: C[k] for k in only.split(',')}
    earlier = [json.load(open(f))['knots'] for f in args[1:]]
    items = [(k, v, json.loads(BY[k]['pd_notation']), None,
              [d for E in earlier if k in E for d in E[k]['diagrams'][1:]]) for k, v in C.items()]
    items.sort(key=lambda it: it[1]['group'] in ('target', 'mobius'))  # the controls (and group upper) run longest: first
    resumed, partial = {}, path + '.partial'      # mode mobius: [name, record] per line, written as the records arrive
    if mode_mobius and os.environ.get('KLEIN_RESUME') and os.path.exists(partial):
        for line in open(partial):
            name, rec = json.loads(line)
            if 'error' not in rec:
                resumed[name] = rec
        items = [it for it in items if it[0] not in resumed]
        print('records taken from', partial, len(resumed), 'knots to run', len(items), flush=True)
    elif mode_mobius:
        open(partial, 'w').close()
    out = {}
    with Pool(workers, initializer=mobius_init_worker if mode_mobius else init_worker, initargs=init) as p:
        vres = dict(p.imap_unordered(validate, [(k, v['pd'], v['alex'], v['det']) for k, v in ref.items()],
                                     chunksize=50))
        print('Fox polynomials of the reference knots:', sum(vres.values()), 'of', len(vres), 'agree', flush=True)
        for name, rec in p.imap_unordered(mobius_task if mode_mobius else one, items):
            out[name] = rec
            r, m = rec.get('representative'), rec.get('moebius_band_of_K')
            if mode_mobius:
                with open(partial, 'a') as f:
                    f.write(json.dumps([name, rec], separators=(',', ':'), default=str) + '\n')
                print(f'{name}: ERROR {rec["error"]}' if 'error' in rec else
                      f'{name}: {rec["group"]} {rec["range"]} h {rec["h_witt"]} det {rec["det"]} stage {rec["stage"]} '
                      f'diagrams {len(rec["diagrams"])} '
                      + (f'Mobius band of K: {m["kind"]}, diagram {m["diagram"]} ({m["diagram_source"]}), s {m["s"]} '
                         if m else 'no Mobius band of K ')
                      + (f'punctured Klein bottle: stage {rec["klein_bottle"]["stage"]} certified {rec["certified"]} '
                         f'self-contained {rec["self_contained"]} ' if rec['klein_bottle'] else '')
                      + (f'rep {r["diagram_source"]} t {r["t"]} J {r["J"]} [{r["verification"]}] ' if r else '')
                      + f'contradictions {rec["contradictions"]} ({rec["seconds"]}s)', flush=True)
                continue
            print(f'{name}: {rec["group"]} h {rec["h_witt"]} det {rec["det"]} stage {rec["stage"]} diagrams '
                  f'{len(rec["diagrams"])} knots {rec["knot_results"]} candidates {rec["alexander_pass"]} hits '
                  f'{len(rec["hits"])} J {len(rec["J_found"])} certified {rec["certified"]} self-contained '
                  f'{rec["self_contained"]} contradictions {rec["contradictions"]} '
                  + (f'rep {r["diagram_source"]} t {r["t"]} J {r["J"]} [{r["verification"]}] ' if r else '')
                  + (f'Mobius band of K: {rec["moebius_band_of_K"]["kind"]} ' if rec.get('moebius_band_of_K') else '')
                  + f'({rec["seconds"]}s)', flush=True)
    out = dict(sorted(dict(out, **resumed).items(), key=lambda kv: knot_order(kv[0])))
    if mode_mobius:
        errors = [k for k, r in out.items() if 'error' in r]
        if errors:                                # the records so far stay in the .partial file
            print('knots with an error (run again with KLEIN_RESUME=1):', errors)
            sys.exit(1)
        summary = mobius_summary(out, vres)
        meta = dict(mode='mobius', T_max=T_MAX, S_max=S_MAX, reps=REPS, max_recheck=MAX_RECHECK, sc_tries=SC_TRIES,
                    rb_hits=RB_HITS, rb_tries=RB_TRIES, kb_tries=KB_TRIES, j_diagrams=J_DIAGRAMS, c_hits=C_HITS,
                    c_variants=C_VARIANTS, pari_stack=[PARI_STACK, PARI_SIZEMAX], resumed_records=len(resumed),
                    reference_knots=len(ref), reference_gamma4_1=sum(v['g4'] == '1' for v in ref.values()),
                    earlier_outputs=len(args[1:]), only=only, versions=versions(), seconds=round(time.time() - t_start))
        json.dump(dict(meta=meta, summary=summary, knots=out), open(path, 'w'), separators=(',', ':'), default=str)
        os.remove(partial)
        print('summary:', {k: v for k, v in summary.items() if k not in ('no_moebius_band', 'klein_bottle')})
        print('punctured Klein bottles:', summary['klein_bottle'])
        print('time %.0fs' % (time.time() - t_start))
        sys.exit()
    if mode_ranges:
        add_face_edges(out, BY)
        summary = ranges_summary(out, vres, ref)
        meta = dict(mode='ranges', T_max=T_MAX, S_max=S_MAX, reps=REPS, max_recheck=MAX_RECHECK, sc_tries=SC_TRIES,
                    rb_hits=RB_HITS, rb_tries=RB_TRIES, kb_tries=KB_TRIES, j_diagrams=J_DIAGRAMS, reference_knots=len(ref),
                    reference_gamma4_1=sum(v['g4'] == '1' for v in ref.values()), earlier_outputs=len(args[1:]),
                    only=only, versions=versions(), seconds=round(time.time() - t_start))
        json.dump(dict(meta=meta, summary=summary, knots=out), open(path, 'w'), separators=(',', ':'), default=str)
        print('summary:', {k: v for k, v in summary.items() if k not in ('representative_J', 'upper', 'knot_11n_133')})
        print('11n_133:', summary['knot_11n_133'])
        print('upper:', summary['upper'])
        print('time %.0fs' % (time.time() - t_start))
        sys.exit()
    summary = default_summary(out, vres)
    meta = dict(T_max=T_MAX, S_max=S_MAX, reps=REPS, max_recheck=MAX_RECHECK, sc_tries=SC_TRIES, rb_hits=RB_HITS,
                rb_tries=RB_TRIES, kb_tries=KB_TRIES, j_diagrams=J_DIAGRAMS, reference_knots=len(ref),
                reference_gamma4_1=sum(v['g4'] == '1' for v in ref.values()), earlier_outputs=len(args[1:]),
                only=only, versions=versions(), seconds=round(time.time() - t_start))
    json.dump(dict(meta=meta, summary=summary, knots=out), open(path, 'w'), indent=0, default=str)
    print('summary:', {k: v for k, v in summary.items() if k not in ('controls',)})
    print('time %.0fs' % (time.time() - t_start))
