#!/usr/bin/env python
"""Locally flat Mobius bands for the knots of range [1,2] for which only the correction terms of [C] Proposition 3.1(c)
exclude a smooth one: the search of [C] Section 4, which gives the twelve knots of Proposition 3.1(c) in Theorem B.

Candidates: KnotInfo gamma_4 '[1,2]' (data/gamma4_knotinfo_2026.10.5.json), not obstructed by the linking-form test
(field obstructed False in data/rank1_targets.json), obstructed by the correction terms of [C] Proposition 3.1(c)
(field obstructed True in data/rank1_dinv_targets.json or data/rank1_dinv_mtargets.json).  These knots have
gamma_4 = 2; there are 81 of them (72 alternating, 9 Montesinos), all with h_witt = 1 (data/linking_genus.json), and
the output records these counts (summary, field candidates).  The knots of [C] Proposition 3.2(a), (b) in Theorem B
are searched with the function one of this script by plumbing_dinv.py, surgery_dinv.py and surgery_cfk.py (mode gap).

Certificate.  If one band move turns K into a knot J with Delta_J = 1, then J bounds a locally flat disc D in B^4
(Freedman), the band is nonorientable (both ends are knots), and the band in a collar together with D is a locally flat
Mobius band.  So gamma_4^top(K) = 1 < gamma_4(K) = 2.

Search.  Face bands of scripts/bands.py: for a face of a diagram of K and two distinct edges on it, J_t is
the band move along the arc in that face whose long sides cross t times ([C] Figure 1).  Diagrams of K: the KnotInfo
PD code, the diagrams with distinct dual graphs found by spherogram many_diagrams (random Reidemeister moves, then
simplify('global'); almost all minimal), and nonminimal diagrams with at most c + EXTRA crossings (random Reidemeister
II and III moves from those).  Each diagram is checked against the KnotInfo exterior (SnapPy is_isometric_to; mirror
images are not distinguished, and gamma_4, gamma_4^top and the certificate are mirror-invariant).
For every arc, det J_t (|t| <= T_MAX) is |det| of a Goeritz matrix (exact; split diagrams have det 0).  By the
Montesinos trick Sigma_2(J_t) is integral surgery on the lift of the arc, with slope moving by one per half twist, so
det J_t = |alpha + n t|, n = det K; this is checked for every arc.  Delta_J = 1 forces det J = 1, which happens for at
most one t per arc (n >= 3).  That t may lie outside |t| <= T_MAX; it is then built and tested too.  For every knot
with det 1, Delta is computed exactly by Fox calculus on the Wirtinger presentation (PARI matdet over Z[t]).
Linking-form screen, used only as a check: det J = 1 makes Sigma_2(J) a ZHS, so Sigma_2(K) is (+-n)-surgery on a knot
in it and lambda_K = <+-1/n>, i.e. [n, +-1] is among the survivors of rank1_linking.py.  Candidates without it must
give no knot of det 1.

Rechecks of every hit (band with Delta_J = 1): Delta_J from a Seifert matrix of J (spherogram seifert_matrix, PARI det
of V - t V^T) and from the Euler characteristic of knot Floer homology; one component (traversal and spherogram);
H_1(Sigma_2 J) = 0 from both Goeritz matrices; SnapPy identify of J and the KnotInfo record of the name found; the dual
band of J (its cocore, a face band with t = 0) gives back K (is_isometric_to the KnotInfo exterior); and J is rebuilt
from the diagram of K by spherogram's own band code (links/bands/core.add_one_band of Dunfield-Gong, band [h1, h4] with
num_twist = -t, its sign convention), which must be a knot with Delta = 1 (Seifert matrix) isometric to J.
J is never smoothly slice when the correction terms are right (a smooth slice disc would give a smooth Mobius band), so
a hit whose J is known to be smoothly slice is flagged as a contradiction.
Validation of the three Alexander routines: all KnotInfo knots up to 13 crossings (Fox and Seifert), a random sample
for knot Floer homology.  Also, for sampled knots J_t, |Delta(-1)| from Fox calculus equals the Goeritz det.

Reproducibility.  spherogram's random moves are not reproducible across runs even with a fixed seed (it iterates over
sets of objects hashed by address), so the output stores every diagram used.  Mode replay scans the stored diagrams
again and compares the counts and the bands with det 1 diagram by diagram; it writes nothing.  The diagrams stored in
the outputs given as earlier_output.json are scanned again (field carried), before the new random diagrams: rare hits
occur in only a few random diagrams, so one run can miss them.
Usage: python scripts/smooth_top_gap.py [workers] [earlier_output.json ...] | replay [workers]
Output: data/smooth_top_gap.json
"""
import json, os, random, re, sys, time, warnings, zlib
from collections import Counter, defaultdict
from itertools import combinations
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from bands import from_pd, faces, band, pd_code, components, graph_components, invariants
from rank1_goeritz import goeritz_pair

T_MAX = 6                                 # |t| <= T_MAX for every arc, plus the det-1 value of t outside that range
EXTRA = 6                                 # nonminimal diagrams have at most c + EXTRA crossings
MINIMAL_TARGET, NONMINIMAL_TRIES = 30, 300   # diagrams wanted from many_diagrams, attempts at nonminimal diagrams
KI_FIELDS = ('alexander_polynomial', 'smooth_four_genus', 'topological_four_genus', 'ozsvath_szabo_tau_invariant',
             'rasmussen_invariant', 'signature', 'three_genus', 'montesinos_notation', 'geometric_type')


def data(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


def candidates():
    G, R1 = data('gamma4_knotinfo_2026.10.5.json'), data('rank1_targets.json')
    D, DM, LG = data('rank1_dinv_targets.json'), data('rank1_dinv_mtargets.json'), data('linking_genus.json')
    out = {}
    for k, v in G.items():
        if v != '[1,2]' or k not in R1 or R1[k]['obstructed'] is not False:
            continue
        dinv = D[k] if k in D and D[k]['obstructed'] is True else DM[k] if k in DM and DM[k]['obstructed'] is True else None
        if dinv is None:
            continue
        r = R1[k]
        out[k] = dict(det=r['det'], test='alternating' if k in D else 'montesinos', h_witt=LG[k]['h_witt'],
                      survivors=r['survivors'], s=r['s'], jordan=r['jordan'], d_spin=dinv['d_spin'],
                      rank_one_full=any(x[0] == r['det'] for x in r['survivors']))
    return dict(sorted(out.items(), key=lambda kv: (int(re.match(r'\d+', kv[0]).group()), kv[0][2], int(kv[0].split('_')[1]))))


# ---------- exact invariants

def bareiss_det(A):
    """Determinant of an integer matrix by fraction-free (Bareiss) elimination."""
    n = len(A)
    if n == 0:
        return 1
    M, sign, prev = [list(r) for r in A], 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            p = next((r for r in range(k + 1, n) if M[r][k] != 0), None)
            if p is None:
                return 0
            M[k], M[p] = M[p], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def det_link(P):
    """Determinant of a link diagram: |det| of the smaller of its two Goeritz matrices."""
    return abs(bareiss_det(min(goeritz_pair([list(x) for x in P]), key=len)))


def pari():
    from cypari import pari as p
    return p


def normalize(coeffs):
    """Laurent polynomial {exponent: coefficient} -> coefficient list from the lowest degree, sign with Delta(1) > 0."""
    c = {e: v for e, v in coeffs.items() if v}
    if not c:
        return [0]
    L = [c.get(e, 0) for e in range(min(c), max(c) + 1)]
    return [-x for x in L] if sum(L) < 0 else L


def poly_coeffs(p):
    P = pari(); t = P('t')
    return {} if p == 0 else {e: int(P.polcoef(p, e, t)) for e in range(int(P.poldegree(p, t)) + 1)}


def orientation(pd):
    """Slots where the link enters each crossing (slot 0 = incoming under-strand), and the number of components."""
    occ = defaultdict(list)
    for c, X in enumerate(pd):
        for s, lab in enumerate(X):
            occ[lab].append((c, s))
    assert all(len(v) == 2 for v in occ.values())
    entering, seen, comps = set(), set(), 0
    for start in [(c, 0) for c in range(len(pd))]:
        if start in seen:
            continue
        comps += 1
        c, s = start
        while (c, s) not in seen:
            seen.add((c, s)); entering.add((c, s)); seen.add((c, (s + 2) % 4))
            a, b = occ[pd[c][(s + 2) % 4]]
            c, s = b if a == (c, (s + 2) % 4) else a
    return entering, comps


def alexander_fox(pd):
    """Delta of a knot from the Wirtinger presentation: row (1-t, t, -1) at (over, in, out) arcs of a positive crossing,
    (1-t, -1, t) for a negative one; det of the minor without the last row and column."""
    n = len(pd)
    if n <= 1:
        return [1]
    entering, comps = orientation(pd)
    assert comps == 1 and all((c, 0) in entering and ((c, 1) in entering) != ((c, 3) in entering) for c in range(n))
    labels = sorted({x for X in pd for x in X})
    parent = {x: x for x in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for X in pd:
        parent[find(X[1])] = find(X[3])
    roots = sorted({find(x) for x in labels})
    arc = {x: roots.index(find(x)) for x in labels}
    P = pari(); t = P('t')
    M = [[P(0)] * n for _ in range(n)]
    for c, X in enumerate(pd):
        k, i, j = arc[X[1]], arc[X[0]], arc[X[2]]
        M[c][k] += 1 - t
        if (c, 3) in entering:            # over-strand west to east, under-strand south to north: positive
            M[c][i] += t; M[c][j] -= 1
        else:
            M[c][i] -= 1; M[c][j] += t
    return normalize(poly_coeffs(P.matdet(P.matrix(n - 1, n - 1, [M[r][s] for r in range(n - 1) for s in range(n - 1)]))))


def alexander_seifert(V):
    n = len(V)
    if n == 0:
        return [1]
    P = pari(); t = P('t')
    return normalize(poly_coeffs(P.matdet(P.matrix(n, n, [V[i][j] - t * V[j][i] for i in range(n) for j in range(n)]))))


def alexander_hfk(ranks):
    c = defaultdict(int)
    for (a, m), r in ranks.items():
        c[a] += r if m % 2 == 0 else -r
    return normalize(c)


def gauss_key(pd):
    """Signed Gauss code up to starting point and direction (for removing repeated diagrams)."""
    entering, comps = orientation(pd)
    assert comps == 1
    occ = defaultdict(list)
    for c, X in enumerate(pd):
        for s, lab in enumerate(X):
            occ[lab].append((c, s))
    seq, (c, s) = [], (0, 0)
    for _ in range(2 * len(pd)):
        seq.append((c, 'U' if s % 2 == 0 else 'O'))
        a, b = occ[pd[c][(s + 2) % 4]]
        c, s = b if a == (c, (s + 2) % 4) else a
    sign = {c: 1 if (c, 3) in entering else -1 for c in range(len(pd))}
    best = None
    for word in (seq, seq[::-1]):
        for r in range(len(word)):
            w, ren = word[r:] + word[:r], {}
            code = tuple((ren.setdefault(c, len(ren)), ou, sign[c]) for c, ou in w)
            best = code if best is None or code < best else best
    return best


# ---------- diagrams of K

def isometric(A, B, tries=8):
    """Whether the exteriors A and B are isometric (SnapPy), retriangulating after a failure; None if every try fails."""
    for _ in range(tries):
        try:
            return bool(A.is_isometric_to(B))
        except Exception:
            A.randomize(); B.randomize()
    return None


def diagrams_of(name, pd0, seed, carried=()):
    """Diagrams of K without repetition (signed Gauss code): the KnotInfo diagram pd0, the carried diagrams, those of
    many_diagrams and the nonminimal diagrams made from them by random Reidemeister II and III moves."""
    import spherogram
    out, keys = [], set()

    def add(src, pd, **kw):
        key = gauss_key(pd)
        if key not in keys:
            keys.add(key); out.append(dict(source=src, pd=[list(x) for x in pd], **kw))
    add('knotinfo', [tuple(x) for x in pd0])
    for d in carried:
        add(d['source'], [tuple(x) for x in d['pd']], carried=True)
    random.seed(seed)
    minimal = spherogram.Link([list(x) for x in pd0]).many_diagrams(target=MINIMAL_TARGET, tries=150, method='backtrack')
    for D in minimal:
        add('many_diagrams', [tuple(x) for x in D.PD_code()])
    c = min(len(D.crossings) for D in minimal)
    for i in range(NONMINIMAL_TRIES):
        random.seed(seed + 1 + i)
        D = minimal[i % len(minimal)].copy()
        D.backtrack(steps=1 + i % 4, prob_type_1=0, prob_type_2=0.5)
        if c < len(D.crossings) <= c + EXTRA:
            add('reidemeister_II_III', [tuple(x) for x in D.PD_code()])
    return out


def verify_diagrams(diags, pd0):
    """Record the crossing number of each diagram and whether its exterior is isometric to that of pd0."""
    import snappy, spherogram             # spherogram's exterior() needs snappy imported
    E0 = spherogram.Link([list(x) for x in pd0]).exterior()
    for d in diags:
        d['crossings'] = len(d['pd'])
        d['isometric_to_knotinfo'] = isometric(spherogram.Link(d['pd']).exterior(), E0)


# ---------- face bands

def arc_darts(pd, face, i, j):
    nbr, under, label = from_pd([tuple(x) for x in pd])
    f = faces(nbr)[face]
    return nbr, under, label, f[i], f[j]


def scan(pd, n, rng):
    """All face bands of one diagram: counts, fit failures, the det-1 knots with Delta (Fox), and sampled checks."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    st, det1, fit_fail, checks = Counter(), [], [], []
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
                if fam[t][0] == 1:
                    st['knots'] += 1
                    if rng.random() < 0.004:  # Fox |Delta(-1)| against the Goeritz det
                        a = alexander_fox(P)
                        checks.append(abs(sum(x * (-1) ** e for e, x in enumerate(a))) == fam[t][1])
            alphas = sorted(a for a in {fam[0][1], -fam[0][1]} if all(abs(a + n * t) == d for t, (_, d) in fam.items()))
            if not alphas:
                fit_fail.append([fi, i, j]); st['fit_failures'] += 1
            ts = {t for t, (c, d) in fam.items() if c == 1 and d == 1}
            for a in alphas:
                for s in (1, -1):
                    if (s - a) % n == 0:
                        t = (s - a) // n
                        if abs(t) <= T_MAX:
                            assert fam[t] == (1, 1), (fi, i, j, t, fam[t])
                        else:
                            ts.add(t); st['bands_outside_range'] += 1
            for t in sorted(ts):
                N, U = band(nbr, under, h1, h3, t)
                P = pd_code(N, U)
                if abs(t) > T_MAX:
                    assert graph_components(N, len(U)) == 1 and len(components(N)) == 1 and det_link(P) == 1
                    st['knots'] += 1
                st['det1_knots'] += 1
                det1.append(dict(face=fi, darts=[i, j], edges=[label[h1], label[h3]], t=t, alexander_fox=alexander_fox(P)))
    return st, det1, fit_fail, checks


def recheck(pd, hit, E_K, BY, DT):
    """Independent checks of a band with Delta_J = 1 (see the docstring)."""
    import spherogram
    nbr, under, label, h1, h3 = arc_darts(pd, hit['face'], *hit['darts'])
    t = hit['t']
    N, U = band(nbr, under, h1, h3, t)
    P = pd_code(N, U)
    L = spherogram.Link([list(x) for x in P])
    rec = dict(face_edges=[label[h] for h in faces(nbr)[hit['face']]], J_pd=[list(x) for x in P], J_crossings=len(P),
               components_traversal=len(components(N)),
               components_spherogram=len(L.link_components),
               h1_goeritz=[list(invariants(G)) if G else [] for G in goeritz_pair([list(x) for x in P])],
               alexander_seifert=alexander_seifert(L.seifert_matrix()))
    S = L.copy(); S.simplify('global')
    rec['J_crossings_simplified'] = len(S.crossings)
    if len(S.crossings) == 0:
        rec.update(alexander_hfk=[1], unknot=True)
    else:
        H = S.knot_floer_homology()
        rec['alexander_hfk'] = alexander_hfk(H['ranks'])
        rec['hfk'] = {k: H[k] for k in ('tau', 'epsilon', 'seifert_genus', 'total_rank', 'fibered', 'L_space_knot')}
        rec['unknot'] = H['total_rank'] == 1
        E = S.exterior()
        try:
            ids = [str(m).split('(')[0] for m in E.identify()]
        except Exception:
            ids = []
        rec['snappy_identify'] = ids
        rec['snappy_solution_type'], rec['snappy_volume'] = E.solution_type(), round(float(E.volume()), 6)
        names = sorted({DT[x] for x in ids if x in DT})
        rec['knotinfo_name'] = names[0] if len(names) == 1 else None
        if rec['knotinfo_name']:
            rec['knotinfo'] = {k: BY[names[0]][k] for k in KI_FIELDS}
        hyp = rec['snappy_solution_type'] == 'all tetrahedra positively oriented'
        rec['J_label'] = (rec['knotinfo_name'] or (ids[0] if ids else None)
                          or (f'hyperbolic, volume {rec["snappy_volume"]:.5f}' if hyp else 'unidentified'))
    # the cocore of the band: a face band with t = 0 between the dart h2 = nbr[h1] and the dart running back to h1
    h2, d = nbr[h1], ((len(under), 1) if t else nbr[h3])
    assert N[d] == h1 and any(h2 in f and d in f for f in faces(N))
    N2, U2 = band(N, U, h2, d, 0)
    rec['dual_band_gives_K'] = isometric(spherogram.Link([list(x) for x in pd_code(N2, U2)]).exterior(), E_K)
    from spherogram.links.bands.core import add_one_band
    B = add_one_band(spherogram.Link([list(x) for x in pd]), ([h1, nbr[h3]], [], -t))
    sb = dict(num_twist=-t, components=len(B.link_components))
    if sb['components'] == 1:
        sb.update(alexander_seifert=alexander_seifert(B.seifert_matrix()), isometric_to_J=isometric(B.exterior(), L.exterior()))
    rec['spherogram_add_one_band'] = sb
    ki = rec.get('knotinfo', {})
    rec['contradiction'] = bool(rec['unknot'] or ki.get('smooth_four_genus') == '0')
    rec['passes'] = (rec['components_traversal'] == rec['components_spherogram'] == 1 and hit['alexander_fox'] == [1]
                     and rec['alexander_seifert'] == [1] and rec['alexander_hfk'] == [1]
                     and all(h == [] for h in rec['h1_goeritz']) and rec['dual_band_gives_K'] is True
                     and sb.get('alexander_seifert') == [1] and sb.get('isometric_to_J') is True)
    return rec


KI, DT = {}, {}                           # KnotInfo fields by name, SnapPy census name -> KnotInfo name (set per worker)


def init_worker(ki, dt):
    KI.update(ki); DT.update(dt)


def one(item):
    """One knot (name, info, KnotInfo PD code, stored diagrams or None, carried diagrams): its diagrams, the scan of
    each, and the rechecks of the hits; given stored diagrams (replay), the scan of those only."""
    name, info, pd0, stored, carried = item
    import spherogram
    t0 = time.time()
    seed = zlib.crc32(name.encode())
    rng = random.Random(seed)
    n = info['det']
    diags = (diagrams_of(name, pd0, seed, carried) if stored is None
             else [dict(source=d['source'], pd=d['pd'], **({'carried': True} if d.get('carried') else {})) for d in stored])
    verify_diagrams(diags, pd0)
    E_K = spherogram.Link([list(x) for x in pd0]).exterior()
    tot, hits, det1_other, fit_fail, checks = Counter(), [], [], [], []
    for di, d in enumerate(diags):
        assert det_link(d['pd']) == n
        st, det1, ff, ch = scan(d['pd'], n, rng)
        tot.update(st); checks += ch
        d.update({k: st[k] for k in ('arcs', 'bands', 'knots', 'det1_knots')}, hits=sum(h['alexander_fox'] == [1] for h in det1))
        fit_fail += [[di] + x for x in ff]
        for h in det1:
            h['diagram'] = di
            if h['alexander_fox'] == [1]:
                if stored is None:
                    h.update(recheck(d['pd'], h, E_K, KI, DT))
                hits.append(h)
            else:
                det1_other.append(h)
    if stored is not None:                # replay: the scan only
        return name, dict(diagrams=diags, hits=hits, det1_not_hits=det1_other, fit_failures=fit_fail,
                          fox_det_checks=[len(checks), sum(checks)])
    J = Counter(h['J_label'] for h in hits)
    good = [h for h in hits if h['passes'] and not h['contradiction']]
    rep = min(good, key=lambda h: (h['diagram'] > 0, h['J_crossings_simplified'], abs(h['t']), h['J_crossings'], h['diagram']),
              default=None)
    if rep is not None:
        rep = dict(diagram=rep['diagram'], diagram_source=diags[rep['diagram']]['source'],
                   **{k: rep[k] for k in ('face', 'face_edges', 'edges', 't', 'J_label', 'J_crossings', 'J_pd')})
    return name, dict(info, diagrams=diags, arcs=tot['arcs'], bands_tried=tot['bands'] + tot['bands_outside_range'],
                      bands_split=tot['split'], bands_outside_range=tot['bands_outside_range'], knot_results=tot['knots'],
                      det1_knots=tot['det1_knots'], fit_failures=fit_fail,
                      fox_det_checks=[len(checks), sum(checks)], hits=hits, det1_not_hits=det1_other,
                      J_found=dict(J.most_common()), certificate=bool(good), representative=rep,
                      contradictions=sum(h['contradiction'] for h in hits), seconds=round(time.time() - t0, 1))


def validate(item):
    """Fox, Seifert (and for sampled knots HFK) Delta against KnotInfo."""
    name, pd, ref_str, ref_vec, do_hfk = item
    import spherogram
    P = pari()
    if ref_vec.strip():
        v = json.loads(ref_vec)
        ref = normalize(dict(zip(range(v[0], v[1] + 1), v[2:])))
    else:
        ref = normalize(poly_coeffs(P(ref_str)))
    L = spherogram.Link(pd)
    res = [alexander_fox([tuple(x) for x in pd]) == ref, alexander_seifert(L.seifert_matrix()) == ref]
    if do_hfk:
        res.append(alexander_hfk(L.knot_floer_homology()['ranks']) == ref)
    return name, res


def replay(C, BY, workers):
    """Mode replay (see the docstring): the stored diagrams of every candidate scanned again and compared."""
    old = json.load(open(os.path.join(ROOT, 'data', 'smooth_top_gap.json')))['candidates']
    items = [(k, v, json.loads(BY[k]['pd_notation']), old[k]['diagrams'], ()) for k, v in C.items()]
    key = lambda h: (h['diagram'], h['face'], tuple(h['darts']), h['t'], tuple(h['alexander_fox']))
    bad = []
    with Pool(workers) as p:
        for name, new in p.imap_unordered(one, items):
            o = old[name]
            same = (all(all(a[k] == b[k] for k in ('arcs', 'bands', 'knots', 'det1_knots', 'hits', 'isometric_to_knotinfo'))
                        for a, b in zip(o['diagrams'], new['diagrams'])) and len(o['diagrams']) == len(new['diagrams'])
                    and sorted(map(key, o['hits'])) == sorted(map(key, new['hits']))
                    and sorted(map(key, o['det1_not_hits'])) == sorted(map(key, new['det1_not_hits']))
                    and o['fit_failures'] == new['fit_failures'] and o['fox_det_checks'] == new['fox_det_checks'])
            print(name, 'diagrams', len(new['diagrams']), 'hits', len(new['hits']), 'same as stored:', same, flush=True)
            if not same:
                bad.append(name)
    print('replay: candidates', len(items), 'differing from the stored output:', bad)


if __name__ == '__main__':
    import snappy, spherogram, cypari, database_knotinfo
    from database_knotinfo import link_list
    mode_replay = len(sys.argv) > 1 and sys.argv[1] == 'replay'
    args = sys.argv[2:] if mode_replay else sys.argv[1:]
    workers = int(args[0]) if args else 8
    t_start = time.time()
    LL = link_list()[1:]
    BY = {r['name']: r for r in LL}
    ki = {r['name']: {k: r[k] for k in KI_FIELDS} for r in LL}
    dt = {'K' + r['dt_name'].replace('_', ''): r['name'] for r in LL if r.get('dt_name')}
    C = candidates()
    if mode_replay:
        replay(C, BY, workers)
        sys.exit()
    check = dict(count=len(C), alternating=sum(v['test'] == 'alternating' for v in C.values()),
                 montesinos=sum(v['test'] == 'montesinos' for v in C.values()),
                 h_witt_all_one=all(v['h_witt'] == 1 for v in C.values()),
                 rank_one_full=sum(v['rank_one_full'] for v in C.values()))
    print('candidates:', check, flush=True)
    earlier = [json.load(open(f))['candidates'] for f in args[1:]]
    items = [(k, v, json.loads(BY[k]['pd_notation']), None, [d for E in earlier for d in E[k]['diagrams']])
             for k, v in C.items()]
    rnd = random.Random(20261018)
    vitems = [(r['name'], json.loads(r['pd_notation']), r['alexander_polynomial'], r['alexander_polynomial_vector'],
               rnd.random() < 0.05) for r in LL if r['crossing_number'] and int(r['crossing_number']) >= 3]
    out = {}
    with Pool(workers, initializer=init_worker, initargs=(ki, dt)) as p:
        for name, rec in p.imap_unordered(one, items):
            out[name] = rec
            print(f'{name}: det {rec["det"]} rank-one {rec["rank_one_full"]} diagrams {len(rec["diagrams"])} '
                  f'bands {rec["bands_tried"]} knots {rec["knot_results"]} det1 {rec["det1_knots"]} hits {len(rec["hits"])} '
                  f'J {rec["J_found"]} cert {rec["certificate"]} fit failures {len(rec["fit_failures"])} '
                  f'contradictions {rec["contradictions"]} ({rec["seconds"]}s)', flush=True)
        vres = dict(p.imap_unordered(validate, vitems, chunksize=50))
    val = dict(knots=len(vres), fox_ok=sum(r[0] for r in vres.values()), seifert_ok=sum(r[1] for r in vres.values()),
               hfk_checked=sum(len(r) == 3 for r in vres.values()), hfk_ok=sum(len(r) == 3 and r[2] for r in vres.values()),
               failures=sorted(k for k, r in vres.items() if not all(r)))
    print('validation of the Alexander routines:', val)
    out = dict(sorted(out.items(), key=lambda kv: list(C).index(kv[0])))
    certified = {k: r['J_found'] for k, r in out.items() if r['certificate']}
    J_all = sum((Counter(r['J_found']) for r in out.values()), Counter())
    summary = dict(candidates=check, certified=certified, not_certified=[k for k in out if k not in certified],
                   not_certified_without_rank_one_form=[k for k, r in out.items() if not r['rank_one_full']],
                   J_found=dict(J_all.most_common()), contradictions=sum(r['contradictions'] for r in out.values()),
                   fit_failures=sum(len(r['fit_failures']) for r in out.values()),
                   det1_knots_without_rank_one_form=sum(r['det1_knots'] for r in out.values() if not r['rank_one_full']),
                   diagrams=sum(len(r['diagrams']) for r in out.values()),
                   diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for r in out.values() for d in r['diagrams']),
                   bands_tried=sum(r['bands_tried'] for r in out.values()),
                   knot_results=sum(r['knot_results'] for r in out.values()),
                   det1_knots=sum(r['det1_knots'] for r in out.values()),
                   hits=sum(len(r['hits']) for r in out.values()),
                   hits_failing_a_recheck=sum(not h['passes'] for r in out.values() for h in r['hits']),
                   fox_det_checks=[sum(r['fox_det_checks'][0] for r in out.values()),
                                   sum(r['fox_det_checks'][1] for r in out.values())])
    meta = dict(T_max=T_MAX, extra_crossings=EXTRA, many_diagrams_target=MINIMAL_TARGET, nonminimal_tries=NONMINIMAL_TRIES,
                earlier_outputs=len(args[1:]),
                carried_diagrams=sum(bool(d.get('carried')) for r in out.values() for d in r['diagrams']),
                versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                              database_knotinfo=database_knotinfo.version(), cypari=cypari.__version__,
                              pari='.'.join(map(str, pari().version())), python=sys.version.split()[0]),
                seconds=round(time.time() - t_start))
    json.dump(dict(meta=meta, summary=summary, validation=val, candidates=out),
              open(os.path.join(ROOT, 'data', 'smooth_top_gap.json'), 'w'), indent=0, default=str)
    print('certified:', len(certified), certified)
    print('summary:', {k: v for k, v in summary.items() if k not in ('certified', 'not_certified')})
    print('time %.0fs' % (time.time() - t_start))
