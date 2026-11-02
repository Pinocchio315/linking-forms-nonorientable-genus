#!/usr/bin/env python
"""The second search for Mobius bands of [C] Section 3.3 (Proposition 3.3(a)), with its controls ([C] Section 3.5).

Knots.  Targets: the 305 knots that the first search (mode mobius of scripts/klein_bottle_search.py) and [C]
Propositions 3.1 and 3.2(a), (b) leave open (the set left0 - newd of scripts/paper_counts.py, taken from
scripts/open_sets.py): 58 alternating, the Montesinos knot 12n_219 and 246 that are neither.  The first search found no
Mobius band for them in about 200 diagrams each (data/klein_bottle_search_mobius.json).  Controls: CONTROLS knots drawn
by random.Random(0) from the 484 controls of mode mobius, for which [C] Proposition 3.1(b) or (c) gives gamma_4 >= 2,
and the two knots of SETTLED, which the exact correction terms of [C] Proposition 3.2(b) exclude; a Mobius band of a
control, or a band of a control to a knot that KnotInfo lists as slice, is a contradiction.
Mobius bands.  A face band turns K into a knot J (knot to knot, so the band is nonorientable), and J bounds a disc made
of our bands: one band to the unknot and a disc; or one band to the two-component unlink and two discs; or two bands
J -> L -> three-component unlink and three discs (a ribbon disc: chi = 1, one boundary component).  The band and that
disc form a Mobius band (chi = 0, one boundary component, nonorientable), so gamma_4(K) <= 1.  Verification
'self-contained': every band lies in the diagram made by the band before it, starting from a stored diagram of K, and
the unknot or unlink at the end is recognised by Reidemeister simplification to no crossings; no identification of
knots and no value of KnotInfo is used.  Verification 'self-contained up to Reidemeister moves of J': the disc of J lies
in a stored diagram of J that spherogram's Reidemeister moves relate to the diagram of J made by the band
(surface_of_J); scripts/mobius_open_rebuild.py relates them again.
KnotInfo-slice J (reported apart, never counted as self-contained).  A knot J of the search (square det, Delta_J of
Fox-Milnor form, signature 0, no disc found) that SnapPy identifies with a knot X that KnotInfo lists as smoothly slice
(smooth four-genus 0 in database_knotinfo 2026.8.1, the same as in the KnotInfo 2026.10.5 CSV for every knot), with
KnotInfo's Alexander polynomial: then K bounds a smooth Mobius band if KnotInfo is right.  Rechecks (slice_recheck):
Delta_J by Fox calculus and from a Seifert matrix; the exterior of J isometric to that of the KnotInfo diagram of X,
with the determinant of the cusp map of each isometry that extends to the knots (+1: J = X, -1: J = mirror image of
X); the Jones polynomial of J against that of the KnotInfo PD code of X and its mirror image (jones: own Kauffman
bracket, which gives t + t^3 - t^4 for the positive trefoil and KnotInfo's listed polynomials on the KnotInfo PD codes,
so spherogram and KnotInfo read those PD codes with the same chirality; HFK tau of spherogram agrees with KnotInfo's,
and only spherogram's Seifert matrix has the opposite sign); the dual band gives back K; spherogram's band code
rebuilds J.  For at most SLICE_EXTRA such X per target the search then looks for a disc of J by our bands in diagrams
of J related by Reidemeister moves (surface_of_J); one found makes a Mobius band of the second verification above.
Search in a diagram of K (search_diagram).  Face bands of scripts/bands.py with |s| <= S_MAX half twists and
the det-1 value of s beyond; det J = |alpha + n s| is predicted on every arc (n = det K, Montesinos trick,
klein_bottle_search.arc_alphas), only bands that can give det 1 or a nonzero square are built, and every det built is
compared with the prediction.  In turn: a band to the unknot; for the J with square det and Delta_J of Fox-Milnor form,
SnapPy's identification, and a J that KnotInfo lists as not slice (smooth four-genus at least 1) is not searched
further (this only orders the search: the certificates do not use it); the other J with signature 0 (Seifert matrix,
exact) get a ribbon disc by one band (ribbon_one), all of them (the first search tried at most 20 J per diagram, and
only discs of one band); bands to the unknot with |s| > S_MAX; a ribbon disc by two bands (ribbon_two) for the J with
|s| <= S2_MAX and the KnotInfo-slice ones.  A band of a ribbon disc splits a component; its number r of half twists is
fixed by linking number 0 of the two new components (split_twist), |r| <= R_MAX.  Before the Reidemeister
simplification an unlink must satisfy unlink_conditions: components of det 1 and Alexander polynomial 1, two-component
sublinks split or of det 0 with a vanishing minor of the Fox matrix, and the Kauffman bracket of the unlink.
Diagrams of K.  First the diagrams of mode mobius (field stored: the KnotInfo diagram, spherogram many_diagrams, and
random Reidemeister moves with at most c + 6 crossings; their isometry check was done there), scanned again with this
search; then rounds of CHUNK new diagrams made by longer sequences of spherogram's random Reidemeister moves from base
diagrams (those of least crossing number): kind 'walk' (10 to 60 moves that do not decrease the crossing number, of
types II and III, on average 1 to 5 of type II), 'walk_level' (20 to 100 moves of types I, II, III, then
simplify('level')), 'walk_global' (20 to 200 moves, then simplify('global')); kept with c to c + X_MAX crossings and a
signed Gauss code new for the knot.  Every new diagram is checked against the KnotInfo exterior (SnapPy
is_isometric_to) and stored with its PD code and the index of the diagram it was made from (field parent):
spherogram's random moves cannot be reproduced from a seed.  A target gets NEW_ROUNDS rounds of new diagrams, a control
one round of CHUNK stored diagrams and one round of new diagrams.
Rounds and checkpoints.  Each knot is scanned in rounds of CHUNK diagrams (stored diagrams first), the knots in turn
round by round; a target stops at its first Mobius band, a control does not.  The main process appends the result of
each round to the .partial file (its first line holds the parameters) and rewrites the .progress.txt file.  Run again
with the same parameters, the script resumes from the .partial file: finished rounds are kept, a round that was running
is done again (with other new diagrams).  The amount of work is set by numbers of diagrams, not by time, so a sleep of
the computer only pauses the run.  The workers enlarge the PARI stack (allocatemem), turn exceptions into the record
{error} and end when the main process is gone; a round is tried at most three times, a pool whose worker died is
replaced, and SIGTERM stops the workers and the main process.
Usage: python scripts/mobius_open_search.py [fresh] [workers] | status [workers] | finalize
  fresh: move an existing .partial file aside and start again; status: print the progress recorded in the .partial
  file; finalize: write the output from the .partial file of a run that was stopped before all its rounds were done,
  each knot with the rounds done (data/mobius_open_search.json.gz is the output of such a run; meta.stopped records
  the rounds done).
Environment, for tests: MOS_ONLY=knot,knot,... restricts the knots, MOS_OUT=path.json.gz replaces the output, MOS_CHUNK
and MOS_NEW_ROUNDS replace CHUNK and NEW_ROUNDS.
Output: data/mobius_open_search.json.gz; while running data/mobius_open_search.partial and .progress.txt
"""
import gzip, heapq, json, os, random, signal, sys, time, warnings, zlib
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from bands import from_pd, faces, band, pd_code, components, graph_components
from smooth_top_gap import det_link, alexander_fox, alexander_seifert, normalize, poly_coeffs, pari, gauss_key, \
    isometric, verify_diagrams, orientation
from klein_bottle_search import arc_alphas, fox_milnor, fox_minor_vanishes, unknot_check, \
    census_names, is_hyperbolic, knot_order
from rank1_linking import signature

S_MAX = 6                                 # |s| <= S_MAX for the band K -> J (and the det-1 value of s beyond)
R_MAX = 12                                # |r| <= R_MAX for the bands of a ribbon disc (r is fixed by linking number 0)
S2_MAX = 0                                # ribbon discs of two bands for the J with |s| <= S2_MAX (and KnotInfo-slice J)
UNLINK_TRIES = 3                          # attempts of simplify('global') that recognise an unlink
CHUNK = int(os.environ.get('MOS_CHUNK', 40))            # diagrams per round
NEW_ROUNDS = int(os.environ.get('MOS_NEW_ROUNDS', 6))    # rounds of new diagrams per target
X_MAX = 10                                # new nonminimal diagrams have at most c + X_MAX crossings
BASES = 40                                # base diagrams of the random walks (the minimal ones first)
ATTEMPTS = 8                              # attempts per new diagram wanted
SLICE_EXTRA = 2                           # KnotInfo-slice knots X per target given the search for a surface of J
J_VARIANTS = 100                          # attempts at other diagrams of such a J
X_DIAGRAMS, X_FOUND = 300, 3              # attempts at diagrams of X with a ribbon disc, and how many are kept
J_REACH = 1000                            # attempts to relate the diagram of J to such a diagram of X
SLICE_PER_J = 2                           # KnotInfo-slice hits rechecked and stored per round and knot X
CONTROLS = 30                             # control knots (gamma_4 >= 2 by [C] Proposition 3.1(b) or (c))
SETTLED = ('13n_3844', '13n_3855')        # controls: gamma_4 = 2 by the exact correction terms of [C] Prop. 3.2(b)
PARI_STACK, PARI_SIZEMAX = 2 ** 25, 2 ** 31
REACHED = 'diagram of J reached from J made by the band'
KI_X = 'diagram of X with a ribbon disc'
KI_X_REFLECTED = 'reflection of a diagram of X with a ribbon disc'
OUT = os.environ.get('MOS_OUT', os.path.join(ROOT, 'data', 'mobius_open_search.json.gz'))
PARTIAL, PROGRESS = OUT.replace('.json.gz', '') + '.partial', OUT.replace('.json.gz', '') + '.progress.txt'


def data(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


# ---------- knots

def targets():
    """The knots that the first search and [C] Propositions 3.1 and 3.2(a), (b) leave open (left0 - newd of
    scripts/paper_counts.py), with their kind: alternating, Montesinos, or neither."""
    from open_sets import upto
    g = upto('newd')
    oleft, alt, mont, nd = g['left0'] - g['newd'], g['alt'], g['mont'], g['nd']
    kind = {k: 'alternating' if alt(k) else 'Montesinos' if mont(k) else 'not alternating, not Montesinos' for k in oleft}
    assert {k for k in oleft if kind[k] == 'Montesinos'} == {'12n_219'} and Counter(kind.values())['alternating'] == 58
    assert {k for k in oleft if kind[k].startswith('not')} <= nd and not oleft & set(SETTLED)
    return {k: kind[k] for k in sorted(oleft, key=knot_order)}


def controls(old):
    """CONTROLS knots of mode mobius whose gamma_4 >= 2 by [C] Proposition 3.1(b) or (c) (group control there), and the
    knots of SETTLED."""
    C = sorted((k for k, v in old.items() if v['group'] == 'control'), key=knot_order)
    assert len(C) == 484 and all(old[k]['group'] == 'mobius' and not old[k]['moebius_band_of_K'] for k in SETTLED)
    return sorted(random.Random(0).sample(C, CONTROLS) + list(SETTLED), key=knot_order)


# ---------- worker state

G = {}                                    # KnotInfo data in the workers: slice, alexander, pd, jones, dt


def worker_init(init):
    """KnotInfo data, a larger PARI stack, and a thread that ends the worker when the main process is gone."""
    import snappy, spherogram             # spherogram's exterior() needs snappy imported
    import threading
    G.update(init)
    pari().allocatemem(PARI_STACK, PARI_SIZEMAX, silent=True)
    parent = os.getppid()

    def watch():
        while True:
            time.sleep(30)
            if os.getppid() != parent:
                os._exit(1)
    threading.Thread(target=watch, daemon=True).start()


def knotinfo_data():
    """Smooth 4-genus, Alexander polynomial, PD code and Jones polynomial of the KnotInfo knots, and the KnotInfo name
    of SnapPy's census names (DT names)."""
    from database_knotinfo import link_list
    LL = link_list()[1:]
    sl, alex, pd, jo = {}, {}, {}, {}
    for r in LL:
        sl[r['name']] = r['smooth_four_genus'].strip()
        if r['smooth_four_genus'].strip() == '0' and r['crossing_number'] != '0':
            v = json.loads(r['alexander_polynomial_vector']) if r['alexander_polynomial_vector'].strip() else None
            alex[r['name']] = (normalize(dict(zip(range(v[0], v[1] + 1), v[2:]))) if isinstance(v, list) else
                               normalize(poly_coeffs(pari()(r['alexander_polynomial']))))
            pd[r['name']] = json.loads(r['pd_notation'])
            jo[r['name']] = r['jones_polynomial']
    dt = {'K' + r['dt_name'].replace('_', ''): r['name'] for r in LL if r.get('dt_name')}
    return dict(slice=sl, alex=alex, pd=pd, jones=jo, dt=dt)


# ---------- Jones polynomial

def bracket(pd):
    """Kauffman bracket <D> of a link diagram (PD code, slots counterclockwise from the incoming under-strand) as
    {exponent of A: coefficient}, <unknot> = 1: summed over the crossings in turn (each next crossing the one sharing
    most edges with those done), the partial states being pairings of the open edges; the A-smoothing of X[a, b, c, d]
    joins a to b and c to d."""
    n = len(pd)
    if n == 0:
        return {0: 1}
    order, count = [], Counter()
    while len(order) < n:
        c = max((c for c in range(n) if c not in order), key=lambda c: (sum(count[x] == 1 for x in pd[c]), -c))
        order.append(c); count.update(pd[c])
    states = {frozenset(): {0: 1}}
    for c in order:
        a, b, cc, d = pd[c]
        new = defaultdict(Counter)
        for key, poly in states.items():
            for arcs, w in ((((a, b), (cc, d)), 1), (((a, d), (b, cc)), -1)):
                m, loops = {}, 0
                for u, v in key:
                    m[u] = v; m[v] = u
                for u, v in arcs:
                    if u == v:                    # both ends of edge u at this crossing
                        loops += 1
                    elif u in m and v in m:
                        if m[u] == v:
                            del m[u], m[v]; loops += 1
                        else:
                            pu, pv = m.pop(u), m.pop(v); m[pu] = pv; m[pv] = pu
                    elif u in m or v in m:
                        u, v = (u, v) if u in m else (v, u)
                        pu = m.pop(u); m[pu] = v; m[v] = pu
                    else:
                        m[u] = v; m[v] = u
                q = Counter({e + w: x for e, x in poly.items()})
                for _ in range(loops):            # d = -A^2 - A^-2
                    q2 = Counter()
                    for e, x in q.items():
                        q2[e + 2] -= x; q2[e - 2] -= x
                    q = q2
                new[frozenset(tuple(sorted(p)) for p in m.items())].update(q)
        states = {k: {e: x for e, x in p.items() if x} for k, p in new.items()}
    assert set(states) == {frozenset()}
    return divide_by_d(states[frozenset()], 1)


def divide_by_d(p, k):
    """p / d^k for d = -A^2 - A^-2 = -A^-2 (A^4 + 1) (Laurent polynomials in A as {exponent: coefficient}), or None if
    d^k does not divide p."""
    for _ in range(k):
        if not p:
            return {}
        lo, hi = min(p), max(p)
        r = [p.get(lo + i, 0) for i in range(hi - lo + 1)]
        q = [0] * max(len(r) - 4, 0)
        for i in range(len(r) - 1, 3, -1):            # long division by x^4 + 1
            if r[i]:
                q[i - 4] = r[i]; r[i - 4] -= r[i]; r[i] = 0
        if any(r[:4]):
            return None
        p = {lo + 2 + i: -x for i, x in enumerate(q) if x}
    return p


def jones(pd):
    """Jones polynomial of a knot diagram as {exponent of t: coefficient}: V = (-A^3)^(-w) <D> with A = t^(-1/4)
    (bracket; w the writhe, the orientation by traversal).  Gives t + t^3 - t^4 for the positive trefoil, and on the
    KnotInfo PD codes the polynomials that KnotInfo lists."""
    br = bracket(pd)
    entering, comps = orientation([tuple(x) for x in pd])
    assert comps == 1
    w = sum(1 if (c, 3) in entering else -1 for c in range(len(pd)))
    out = Counter()
    for e, x in br.items():
        assert (e - 3 * w) % 4 == 0
        out[-(e - 3 * w) // 4] += (-1) ** (w % 2) * x
    return {e: x for e, x in sorted(out.items()) if x}


def bracket_of_unlink(pd, ncomp):
    """Whether the Kauffman bracket of the diagram is that of the ncomp-component unlink up to a unit: +-A^k d^(ncomp-1)."""
    q = divide_by_d(bracket(pd), ncomp - 1)
    return q is not None and len(q) == 1 and abs(next(iter(q.values()))) == 1


def mirror_jones(V):
    return {-e: x for e, x in sorted(V.items(), key=lambda kv: -kv[0])}


def knotinfo_jones(s):
    """KnotInfo's Jones polynomial string (in t, with t^(-k) or 1/t^k) as {exponent: coefficient}."""
    import sympy
    t = sympy.Symbol('t')
    e = sympy.expand(sympy.sympify(s.replace('^', '**'), locals={'t': t}))
    k = 60
    p = sympy.Poly(sympy.expand(e * t ** k), t)
    return {m[0] - k: int(c) for m, c in zip(p.monoms(), p.coeffs()) if c}


# ---------- links: linking numbers, unlinks

def face_arcs(nbr, label):
    return [(fi, f, i, j) for fi, f in enumerate(faces(nbr)) for i in range(len(f)) for j in range(i + 1, len(f))
            if label[f[i]] != label[f[j]]]


def band_record(arc, label, **kw):
    fi, f, i, j = arc
    return dict(face=fi, darts=[i, j], face_edges=[label[h] for h in f], edges=[label[f[i]], label[f[j]]], **kw)


def arrivals(N, h0):
    """The ends at which the component through the dart h0 arrives at crossings, traversed starting along h0."""
    arr, h = set(), h0
    while True:
        a = N[h]; arr.add(a); h = (a[0], (a[1] + 2) % 4)
        if h == h0:
            return arr


def linking(U, A, B):
    """Linking number of the components with arriving ends A and B: half the sum of the signs of their crossings, a
    crossing positive when the over-strand arrives one slot before the under-strand (slots counterclockwise, as in
    smooth_top_gap.alexander_fox)."""
    at, total = {c: s for c, s in A}, 0
    for c, s in B:
        if c in at:
            ui, oi = (at[c], s) if (at[c] - U[c]) % 2 == 0 else (s, at[c])
            total += 1 if oi == (ui - 1) % 4 else -1
    assert total % 2 == 0
    return total // 2


def split_twist(nbr, under, h1, h3, forward, third=False):
    """For the face band on the darts h1, h3 of one component: the band splits the component when it is coherent,
    that is, for even r when the component runs the same way around the face along both edges (forward: the darts left
    by the traversal of components), and for odd r otherwise; the r of that parity with linking number 0 between the
    two new components (traversed from h1 and from h2 = nbr[h1], along the band): two more half twists change it by
    +-1.  third (a two-component link): also linking number 0 of both with the other component and det 0 of the result
    (det L_r = |alpha + r det L| = |alpha| when det L = 0).  Returns r, or None."""
    h2 = nbr[h1]
    par = 0 if (h1 in forward) == (h3 in forward) else 1
    N, U = band(nbr, under, h1, h3, par)
    A, B = arrivals(N, h1), arrivals(N, h2)
    if third:
        used = A | B | {(c, (s + 2) % 4) for c, s in A | B}
        Y = arrivals(N, next(h for h in sorted(N) if h not in used))
        if linking(U, A, Y) or linking(U, B, Y):
            return None
        if graph_components(N, len(U)) == 1 and det_link(pd_code(N, U)) != 0:
            return None
    l0 = linking(U, A, B)
    if l0 == 0:
        return par
    N2, U2 = band(nbr, under, h1, h3, par + 2)
    l2 = linking(U2, arrivals(N2, h1), arrivals(N2, h2))
    if abs(l2 - l0) != 1:
        return None
    return par - 2 * l0 * (l2 - l0)


def label_components(P):
    """Edge label -> component, components numbered by their least label (the order of pd_code)."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b, c, d in P:
        parent[find(a)] = find(c); parent[find(b)] = find(d)
    least = {}
    for X in P:
        for x in X:
            r = find(x)
            least[r] = min(least.get(r, x), x)
    rank = {r: i for i, r in enumerate(sorted(least, key=least.get))}
    return {x: rank[find(x)] for X in P for x in X}


def sublink(P, keep, comp):
    """The diagram of the sublink of the components in keep (comp = label_components(P)), the other components erased,
    and whether it is split (its crossings fall into several pieces, or a component of it has no crossing with the
    others)."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    own = []
    for a, b, c, d in P:
        ka, kb = comp[a] in keep, comp[b] in keep
        if ka and kb:
            own.append((a, b, c, d))
        elif ka:
            parent[find(a)] = find(c)
        elif kb:
            parent[find(b)] = find(d)
    Q = [tuple(find(x) for x in X) for X in own]
    cross = list(range(len(Q)))

    def root(x):
        while cross[x] != x:
            cross[x] = cross[cross[x]]; x = cross[x]
        return x
    where = defaultdict(list)
    for i, X in enumerate(Q):
        for x in X:
            where[x].append(i)
    for v in where.values():
        for i in v[1:]:
            cross[root(i)] = root(v[0])
    split = len({root(i) for i in range(len(Q))}) > 1 or len({comp[x] for X in own for x in X}) < len(keep)
    return Q, split


def component_alexander(P, comp):
    """Alexander polynomial (Fox calculus) of each component of the link diagram P drawn alone."""
    return [alexander_fox(sublink(P, {k}, comp)[0]) for k in range(max(comp.values()) + 1)]


def unlink_conditions(P, ncomp):
    """Necessary conditions for the ncomp-component unlink, checked before the Reidemeister moves, the cheap ones first:
    every component alone has det 1, every two-component sublink is split or has det 0; every component has Alexander
    polynomial 1, every non-split two-component sublink has a vanishing minor of its Fox matrix
    (klein_bottle_search.fox_minor_vanishes: for the unlink the first elementary ideal is 0); the Kauffman bracket is
    that of the unlink up to a unit (bracket_of_unlink)."""
    comp = label_components(P)
    one = [sublink(P, {k}, comp)[0] for k in range(ncomp)]
    two = [sublink(P, {x, y}, comp) for x in range(ncomp) for y in range(x + 1, ncomp)]
    if any(len(Q) > 1 and det_link(Q) != 1 for Q in one) or any(not split and det_link(Q) != 0 for Q, split in two):
        return False
    if any(alexander_fox(Q) != [1] for Q in one) or any(not split and not fox_minor_vanishes(Q) for Q, split in two):
        return False
    return bracket_of_unlink(P, ncomp)


def is_unlink(P, n):
    import spherogram
    for _ in range(UNLINK_TRIES):
        L = spherogram.Link([list(x) for x in P]); L.simplify('global')
        if not L.crossings and L.unlinked_unknot_components == n:
            return True
    return False


def ribbon_one(pd, st):
    """A face band from the knot diagram pd to the two-component unlink (r fixed by linking number 0), the least |r|
    first; with two discs a ribbon disc.  Checked before the Reidemeister moves: det 0 and unlink_conditions."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs, forward = face_arcs(nbr, label), {h for comp in components(nbr) for h in comp}
    cands = []
    for k, (fi, f, i, j) in enumerate(arcs):
        r = split_twist(nbr, under, f[i], f[j], forward)
        if r is not None and abs(r) <= R_MAX:
            cands.append((abs(r), r, k))
    st['ribbon_one_arcs'] += len(arcs)
    for _, r, k in sorted(cands):
        fi, f, i, j = arcs[k]
        N, U = band(nbr, under, f[i], f[j], r)
        P = pd_code(N, U)
        if graph_components(N, len(U)) == 1 and det_link(P) != 0:
            continue
        if not unlink_conditions(P, 2):
            continue
        st['unlink_checks'] += 1
        if is_unlink(P, 2):
            return band_record(arcs[k], label, r=r, crossings=len(P), components=2)
    return None


def ribbon_two(pd, st):
    """Two face bands from the knot diagram pd to the three-component unlink, J -> L -> unlink, the second in the
    diagram of L made by the first; with three discs a ribbon disc (chi = 1, one boundary component).  Each band splits
    a component, r fixed by linking number 0.  Checked first: det L = 0, a component of L with Alexander polynomial 1
    and the other of Fox-Milnor form; for the second band (on the other component), linking number 0 with the first
    component and det 0 (det L_r = |alpha| for every r, as det L = 0); then unlink_conditions."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs, forward = face_arcs(nbr, label), {h for comp in components(nbr) for h in comp}
    firsts = []
    for k, (fi, f, i, j) in enumerate(arcs):
        r = split_twist(nbr, under, f[i], f[j], forward)
        if r is None or abs(r) > R_MAX:
            continue
        N, U = band(nbr, under, f[i], f[j], r)
        Q = pd_code(N, U)
        if graph_components(N, len(U)) == 1 and det_link(Q) != 0:
            continue
        ca = component_alexander(Q, label_components(Q))
        cd = [A == [1] for A in ca]                 # a component that may be the unknot
        if any(cd[a] and fox_milnor(ca[1 - a]) for a in (0, 1)):
            firsts.append((abs(r), r, k, Q, cd))
    st['ribbon_two_links'] += len(firsts)
    for _, r, k, Q, cd in sorted(firsts, key=lambda x: x[:3]):
        nb2, un2, lab2 = from_pd([tuple(x) for x in Q])
        comp = label_components(Q)
        arcs2, fwd2 = face_arcs(nb2, lab2), {h for c in components(nb2) for h in c}
        cands = []
        for k2, (gi, g, a, b) in enumerate(arcs2):
            X = comp[lab2[g[a]]]
            if comp[lab2[g[b]]] != X or not cd[1 - X]:
                continue
            r2 = split_twist(nb2, un2, g[a], g[b], fwd2, third=True)
            if r2 is not None and abs(r2) <= R_MAX:
                cands.append((abs(r2), r2, k2))
        st['ribbon_two_arcs'] += len(arcs2)
        for _, r2, k2 in sorted(cands):
            gi, g, a, b = arcs2[k2]
            R = pd_code(*band(nb2, un2, g[a], g[b], r2))
            if not unlink_conditions(R, 3):
                continue
            st['unlink_checks'] += 1
            if is_unlink(R, 3):
                return [band_record(arcs[k], label, r=r, after=[list(x) for x in Q], components=2),
                        band_record(arcs2[k2], lab2, r=r2, crossings=len(R), components=3)]
    return None


def ribbon(pd, st):
    """A ribbon disc of the knot diagram pd by our bands: ribbon_one, else ribbon_two (stored as the list of bands)."""
    r = ribbon_one(pd, st)
    return ([r] if r else None) or ribbon_two(pd, st)


# ---------- the search in one diagram of K

def sigma(P):
    import spherogram
    V = spherogram.Link([list(x) for x in P]).seifert_matrix()
    return signature([[V[i][j] + V[j][i] for j in range(len(V))] for i in range(len(V))]) if len(V) else 0


def identify(P):
    """SnapPy's census names of the exterior of the knot P, and the KnotInfo names among them (DT names)."""
    import spherogram
    ids = census_names(spherogram.Link([list(x) for x in P]).exterior())
    return ids, sorted({G['dt'][x] for x in ids if x in G['dt']})


def search_diagram(pd, n, st, slice_hits):
    """Face bands K -> J in the diagram pd of K (det n), |s| <= S_MAX and the det-1 value of s beyond; det J =
    |alpha + n s| is predicted on every arc (arc_alphas) and only bands that can give det 1 or a nonzero square are
    built.  In turn: J the unknot (Alexander polynomial 1, then Reidemeister simplification or knot Floer homology);
    for the J with square det and Delta_J of Fox-Milnor form, SnapPy's identification: a J that KnotInfo lists as not
    slice (smooth four-genus at least 1) is not searched further; the others (KnotInfo-slice with KnotInfo's Alexander
    polynomial, or not identified, or of open sliceness) with signature 0 get a ribbon disc by one band (ribbon_one);
    then bands to the unknot with |s| > S_MAX; then a ribbon disc by two bands (ribbon_two) for those with |s| <= S2_MAX
    and the KnotInfo-slice ones.  The first found is
    returned (a Mobius band of K).  Otherwise the KnotInfo-slice J are appended to slice_hits.  Counts go to st."""
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = face_arcs(nbr, label)
    alphas = [arc_alphas(nbr, under, f[i], f[j], n) for fi, f, i, j in arcs]
    st['arcs'] += len(arcs); st['fit_failures'] += sum(not a for a in alphas)
    dets = lambda k, s: {abs(a + n * s) for a in alphas[k]}
    square = lambda x: x > 0 and isqrt(x) ** 2 == x

    def knot(k, s):
        fi, f, i, j = arcs[k]
        N, U = band(nbr, under, f[i], f[j], s)
        st['bands_built'] += 1
        if graph_components(N, len(U)) > 1:
            st['det_mismatch'] += bool(alphas[k]) and 0 not in dets(k, s)
            return None, None
        P = pd_code(N, U)
        d = det_link(P)
        st['det_mismatch'] += bool(alphas[k]) and d not in dets(k, s)
        return (P, d) if len(components(N)) == 1 else (None, None)

    def found(k, s, P, **kw):
        return band_record(arcs[k], label, s=s, J_pd=[list(x) for x in P], **kw)
    Js = []
    for s in sorted(range(-S_MAX, S_MAX + 1), key=lambda x: (abs(x), x)):
        for k in range(len(arcs)):
            if alphas[k] and not any(square(x) for x in dets(k, s)):
                continue
            P, d = knot(k, s)
            if P is None or not square(d):
                continue
            st['J_square_det'] += 1
            A = alexander_fox(P)
            if A == [1]:
                how = unknot_check(P)
                if how:
                    return found(k, s, P, kind='band to the unknot', recognised_by=how)
                st['J_alexander_one'] += 1
            if not fox_milnor(A):
                continue
            st['J_fox_milnor'] += 1
            ids, names = identify(P)
            g4 = {G['slice'].get(x) for x in names}
            sl = [x for x in names if G['slice'].get(x) == '0' and G['alex'][x] == A]
            if names and not sl and not g4 & {'0', '[0,1]', ''}:
                st['J_knotinfo_not_slice'] += 1
                continue
            st['J_alexander_mismatch'] += bool(names) and '0' in g4 and not sl
            if sigma(P) != 0:
                st['J_signature_nonzero'] += 1
                continue
            st['J_searched'] += 1
            st['J_knotinfo_slice'] += bool(sl)
            Js.append((s, k, P, d, A, sl, names, ids))
    for s, k, P, d, A, sl, names, ids in Js:
        st['ribbon_one'] += 1
        rb = ribbon_one(P, st)
        if rb:
            return found(k, s, P, kind='band to a ribbon knot', J_det=d, J_alexander=A, ribbon_bands=[rb])
    far = [((e - a) // n, k) for k in range(len(arcs)) for a in alphas[k] for e in (1, -1)
           if (e - a) % n == 0 and abs((e - a) // n) > S_MAX]
    for s, k in sorted(set(far), key=lambda x: (abs(x[0]), x)):
        st['far_bands'] += 1
        P, d = knot(k, s)
        if d == 1 and alexander_fox(P) == [1]:
            how = unknot_check(P)
            if how:
                return found(k, s, P, kind='band to the unknot', recognised_by=how)
    for s, k, P, d, A, sl, names, ids in Js:
        if abs(s) > S2_MAX and not sl:
            continue
        st['ribbon_two'] += 1
        rb = ribbon_two(P, st)
        if rb:
            return found(k, s, P, kind='band to a ribbon knot', J_det=d, J_alexander=A, ribbon_bands=rb)
    for s, k, P, d, A, sl, names, ids in Js:
        if sl:
            slice_hits.append(found(k, s, P, J=sl[0], J_names=names, census=ids[:3], J_det=d, J_alexander=A))
    return None


# ---------- the surface of a KnotInfo-slice J by our bands, in diagrams related by Reidemeister moves

def reflection(P):
    """The reflected diagram (counterclockwise order reversed): a diagram of the mirror image."""
    return [[a, d, c, b] for a, b, c, d in P]


def variant(P, i):
    """Attempt i of turning the diagram P into others by spherogram's Reidemeister moves (as reaches in
    scripts/klein_bottle_rebuild.py): simplify('level'), simplify('global'), then random moves (100 as in
    many_diagrams, or 1 to 12 that do not decrease the crossing number) followed by simplify('global')."""
    import spherogram
    L = spherogram.Link([list(x) for x in P])
    if i >= 2 and i % 2:
        L.backtrack(100)
    elif i >= 2:
        L.backtrack(steps=1 + i % 12, prob_type_1=0.2, prob_type_2=0.4)
    L.simplify('level' if i == 0 else 'global')
    return [tuple(x) for x in L.PD_code()]


X_CACHE = {}                              # worker: KnotInfo-slice knot X -> diagrams of X with a ribbon disc


def ribbon_diagrams_of_X(X, seed, st):
    """Diagrams of the KnotInfo-slice knot X with a ribbon disc by our bands (ribbon): the KnotInfo diagram, then
    diagrams made from it by spherogram's random moves (1 to 8 moves of types II and III that do not decrease the
    crossing number, every third time after 100 moves and simplify('global')), X_DIAGRAMS attempts, at most X_FOUND of
    them; kept in the worker."""
    import spherogram
    if X not in X_CACHE:
        D0, out, keys = [tuple(x) for x in G['pd'][X]], [], set()
        for i in range(X_DIAGRAMS):
            random.seed(seed + i)
            D = D0
            if i:
                L = spherogram.Link([list(x) for x in D0])
                if i % 3 == 0:
                    L.backtrack(100); L.simplify('global')
                L.backtrack(steps=1 + i % 8, prob_type_1=0, prob_type_2=0.5)
                D = [tuple(x) for x in L.PD_code()]
            if not D or gauss_key(D) in keys:
                continue
            keys.add(gauss_key(D))
            st['X_diagrams'] += 1
            rb = ribbon(D, st)
            if rb:
                out.append(D)
                if len(out) >= X_FOUND:
                    break
        X_CACHE[X] = out
    return X_CACHE[X]


def surface_of_J(P, X, seed, st):
    """For the diagram P of J made by the band, J identified with the KnotInfo-slice knot X: a ribbon disc of J by our
    bands (ribbon) in (a) a diagram V reached from P by spherogram's Reidemeister moves (variant), or W made from such a
    V by 1 to 8 random moves of types II and III that do not decrease the crossing number, J_VARIANTS attempts; or (b) a
    diagram D of X with such a disc (ribbon_diagrams_of_X), reflected when the Jones polynomials say that J is the mirror
    image of X, once variants of P and of D meet (a common signed Gauss code, at most J_REACH attempts each).  The
    diagrams are stored; scripts/mobius_open_rebuild.py relates them to P by Reidemeister moves again.  None if neither
    is found."""
    import spherogram
    keys = set()
    for i in range(J_VARIANTS):
        random.seed(seed + i)
        V = variant(P, i)
        W, steps = V, 0
        if V and i % 2:
            L = spherogram.Link([list(x) for x in V])
            steps = 1 + (i // 2) % 8
            L.backtrack(steps=steps, prob_type_1=0, prob_type_2=0.5)
            W = [tuple(x) for x in L.PD_code()]
        if not W or gauss_key(W) in keys:
            continue
        keys.add(gauss_key(W))
        st['J_variants'] += 1
        rb = ribbon(W, st)
        if rb:
            return dict(diagram=REACHED, diagram_pd=[list(x) for x in W], variant=i, walk_steps=steps,
                        simplified_pd=[list(x) for x in V] if steps else None, ribbon_bands=rb)
    VJ, VX = jones(P), jones([tuple(x) for x in G['pd'][X]])
    sides = [s for s, ok in ((False, VJ == VX), (True, VJ == mirror_jones(VX))) if ok]
    for D0 in ribbon_diagrams_of_X(X, zlib.crc32(X.encode()), st):
        for reflected in sides:
            D = reflection(D0) if reflected else D0
            rb = ribbon([tuple(x) for x in D], st)
            if not rb:
                continue
            side, other = {gauss_key(P)}, {gauss_key([tuple(x) for x in D])}
            for i in range(J_REACH):
                random.seed(seed + 10 ** 6 + i)
                VP, VD = variant(P, i), variant(D, i)
                side.add(VP and gauss_key(VP)); other.add(VD and gauss_key(VD))
                if (side & other) - {None, ()}:
                    return dict(diagram=KI_X_REFLECTED if reflected else KI_X, diagram_pd=[list(x) for x in D],
                                met_after=i + 1, ribbon_bands=rb)
            st['X_diagram_not_met'] += 1
    return None


# ---------- rechecks

def isometry_signs(A, B, tries=8):
    """Determinants of the cusp maps of the isometries between two knot exteriors that extend to the knots (+1: the
    knots are equal, -1: mirror images); [] if not isometric, None if SnapPy fails."""
    for _ in range(tries):
        try:
            return sorted({int(m.cusp_maps()[0].det()) for m in A.is_isometric_to(B, return_isometries=True)
                           if m.extends_to_link()})
        except Exception:
            A.randomize(); B.randomize()
    return None


def chirality(VJ, VX):
    if VJ == VX == mirror_jones(VX):
        return 'equal, Jones polynomial symmetric'
    return 'same' if VJ == VX else 'mirror' if VJ == mirror_jones(VX) else 'different'


def cocore_and_rebuild(pd, rec, t, E_K, E_J, hyperbolic):
    """The dual band (cocore, a face band with t = 0) of the band gives back K; spherogram's band code
    (links/bands/core.add_one_band of Dunfield-Gong, num_twist = -t) rebuilds J."""
    import spherogram
    from spherogram.links.bands.core import add_one_band
    nbr, under, label = from_pd([tuple(x) for x in pd])
    f = faces(nbr)[rec['face']]
    h1, h3 = f[rec['darts'][0]], f[rec['darts'][1]]
    N, U = band(nbr, under, h1, h3, t)
    h2, d = nbr[h1], ((len(under), 1) if t else nbr[h3])
    assert N[d] == h1 and any(h2 in g and d in g for g in faces(N))
    N2, U2 = band(N, U, h2, d, 0)
    out = dict(dual_band_gives_K=isometric(spherogram.Link([list(x) for x in pd_code(N2, U2)]).exterior(), E_K))
    B = add_one_band(spherogram.Link([list(x) for x in pd]), ([h1, nbr[h3]], [], -t))
    out['add_one_band_components'] = len(B.link_components)
    out['add_one_band_alexander'] = alexander_seifert(B.seifert_matrix()) if len(B.link_components) == 1 else None
    out['add_one_band_isometric_to_J'] = isometric(B.exterior(), E_J) if hyperbolic and len(B.link_components) == 1 \
        else None
    return out


def slice_recheck(pd, hit, E_K):
    """Checks of a KnotInfo-slice hit: one component; Delta_J by Fox calculus and from a Seifert matrix equal to
    KnotInfo's; the exterior of J isometric to that of the KnotInfo diagram of X, with the orientation of the
    isometries; the Jones polynomial of J against that of the KnotInfo PD code of X (same convention) and its mirror
    image, and that of X against the polynomial KnotInfo lists; the dual band and spherogram's band code."""
    import spherogram
    P, X = [tuple(x) for x in hit['J_pd']], hit['J']
    L = spherogram.Link([list(x) for x in P])
    rec = dict(components=len(L.link_components), alexander_fox=alexander_fox(P),
               alexander_seifert=alexander_seifert(L.seifert_matrix()), alexander_knotinfo=G['alex'][X])
    E_J, E_X = L.exterior(), spherogram.Link(G['pd'][X]).exterior()
    hyperbolic = is_hyperbolic(E_J)
    rec['hyperbolic'] = hyperbolic
    rec['isometry_signs'] = isometry_signs(E_J, E_X) if hyperbolic else None
    VJ, VX = jones(P), jones([tuple(x) for x in G['pd'][X]])
    rec['jones_J'] = sorted(VJ.items())
    rec['jones'] = chirality(VJ, VX)
    rec['jones_X_equals_knotinfo_listed'] = VX == knotinfo_jones(G['jones'][X])
    rec.update(cocore_and_rebuild(pd, hit, hit['s'], E_K, E_J, hyperbolic))
    signs = rec['isometry_signs'] or []
    rec['chirality'] = ('J = mirror image of X' if rec['jones'] == 'mirror' else 'J = X' if rec['jones'] == 'same' else
                        'by the isometries: ' + ('J = X' if signs == [1] else 'J = mirror image of X' if signs == [-1]
                                                 else 'X amphichiral' if signs == [-1, 1] else 'undecided'))
    agree = (rec['jones'] == 'same' and 1 in signs) or (rec['jones'] == 'mirror' and -1 in signs) or \
        rec['jones'].startswith('equal')
    rec['passes'] = (rec['components'] == 1 and rec['alexander_fox'] == rec['alexander_seifert'] == rec['alexander_knotinfo']
                     and bool(signs) and agree and rec['jones_X_equals_knotinfo_listed']
                     and rec['dual_band_gives_K'] is True and rec['add_one_band_isometric_to_J'] is True)
    return rec


def certificate_recheck(pd, cert, E_K):
    """The bands of a Mobius band of K rebuilt by this script (J, and the unlink recognised again), the dual band and
    spherogram's band code."""
    import spherogram
    nbr, under, label = from_pd([tuple(x) for x in pd])
    f = faces(nbr)[cert['face']]
    P = pd_code(*band(nbr, under, f[cert['darts'][0]], f[cert['darts'][1]], cert['s']))
    rec = dict(J_rebuilt=[list(x) for x in P] == cert['J_pd'])
    E_J = spherogram.Link([list(x) for x in P]).exterior()
    hyperbolic = is_hyperbolic(E_J)
    if cert['kind'] == 'band to the unknot':
        rec['unknot'] = bool(unknot_check(P))
    else:
        s = cert.get('surface_of_J')
        D = [tuple(x) for x in (s['diagram_pd'] if s else P)]
        for b in (s['ribbon_bands'] if s else cert['ribbon_bands']):
            nb, un, lab = from_pd(D)
            g = faces(nb)[b['face']]
            Q = pd_code(*band(nb, un, g[b['darts'][0]], g[b['darts'][1]], b['r']))
            if b.get('after') is not None:
                rec.setdefault('after_rebuilt', []).append([list(x) for x in Q] == b['after'])
            D = Q
        rec['unlink'] = is_unlink(D, cert_components(cert))
    rec.update(cocore_and_rebuild(pd, cert, cert['s'], E_K, E_J, hyperbolic))
    ok_J = rec['add_one_band_isometric_to_J'] is True or (not hyperbolic and rec['add_one_band_components'] == 1)
    rec['passes'] = (rec['J_rebuilt'] and rec.get('unknot', rec.get('unlink')) is True
                     and all(rec.get('after_rebuilt', [True])) and rec['dual_band_gives_K'] is True and ok_J)
    return rec


def cert_components(cert):
    s = cert.get('surface_of_J')
    return (s or cert)['ribbon_bands'][-1]['components']


# ---------- diagrams of K

def new_diagrams(job, st):
    """New diagrams of K by spherogram's random Reidemeister moves from a base diagram: kind 'walk' (10 to 60 moves
    that do not decrease the crossing number, types II and III, on average 1 to 5 of type II), 'walk_level' (20 to 100
    moves of types I, II and III, then simplify('level')), 'walk_global' (20 to 200 such moves, then
    simplify('global')).  Kept if the signed Gauss code is new for the knot and there are c to c + X_MAX crossings;
    parent: the index of the base diagram."""
    import spherogram
    seen, bases, c = {tuple(map(tuple, k)) for k in job['seen']}, job['bases'], job['crossings']
    out, attempts = [], 0
    while len(out) < job['count'] and attempts < ATTEMPTS * job['count']:
        random.seed(job['seed'] + attempts)
        kind = job['kinds'][attempts % len(job['kinds'])]
        attempts += 1
        parent, base = bases[random.randrange(len(bases))]
        L = spherogram.Link([list(x) for x in base])
        if kind == 'walk':
            steps = random.randint(10, 60)
            L.backtrack(steps=steps, prob_type_1=0, prob_type_2=random.randint(1, 5) / steps)
        elif kind == 'walk_level':
            L.backtrack(steps=random.randint(20, 100), prob_type_1=0.05, prob_type_2=0.3); L.simplify('level')
        else:
            L.backtrack(steps=random.randint(20, 200), prob_type_1=0.1, prob_type_2=0.4); L.simplify('global')
        P = [tuple(x) for x in L.PD_code()]
        if len(L.link_components) != 1 or not c <= len(P) <= c + X_MAX:
            st['rejected_crossings'] += 1
            continue
        key = gauss_key(P)
        if key in seen:
            st['repeated_diagrams'] += 1
            continue
        seen.add(key)
        out.append(dict(source=kind, parent=parent, pd=[list(x) for x in P]))
    st['generation_attempts'] += attempts
    return out


# ---------- one round: a chunk of stored diagrams, or new diagrams

def run_round(job):
    """Scan the diagrams of one round of one knot (see the docstring); for a target the scan stops at the first Mobius
    band of K, for a control it goes on (any Mobius band is a contradiction)."""
    import spherogram
    t0, c0 = time.time(), time.process_time()
    name, n, pd0, control = job['knot'], job['det'], job['pd0'], job['group'] == 'control'
    st = Counter()
    if job['kind'] == 'stored':                 # checked against the KnotInfo exterior in mode mobius
        todo = [dict(index=i, source=src, pd=pd, crossings=len(pd), isometric_to_knotinfo=True)
                for i, src, pd in job['diagrams']]
    else:
        todo = new_diagrams(job, st)
        verify_diagrams(todo, pd0)
    E_K = spherogram.Link([list(x) for x in pd0]).exterior()
    scanned, certs, hits = [], [], []
    for pos, d in enumerate(todo):
        assert det_link(d['pd']) == n
        ds, dh = Counter(), []
        cert = search_diagram([tuple(x) for x in d['pd']], n, ds, dh)
        st.update(ds)
        d.update(counts=dict(ds), moebius_band=cert['kind'] if cert else None, knotinfo_slice=len(dh))
        scanned.append(d)
        where = dict(diagram=d['index']) if 'index' in d else dict(position=pos)     # position: in this round
        for h in dh:
            h.update(diagram_source=d['source'], diagram_pd=d['pd'], **where)
        hits += dh
        if cert:
            cert.update(diagram_source=d['source'], diagram_pd=d['pd'], verification='self-contained', **where)
            cert['recheck'] = certificate_recheck(d['pd'], cert, E_K)
            certs.append(cert)
            if not control:
                break
    # KnotInfo-slice J: a surface of J by our bands in related diagrams (at most SLICE_EXTRA distinct J per target)
    extra, done = [], set(job['slice_extra_done'])
    if not certs and not control:
        for h in hits:
            if h['J'] in done or len(done) >= SLICE_EXTRA:
                continue
            done.add(h['J']); extra.append(h['J'])
            s = surface_of_J([tuple(x) for x in h['J_pd']], h['J'], zlib.crc32((name + h['J']).encode()), st)
            if s:
                cert = {k: v for k, v in h.items() if k not in ('J', 'J_names', 'census')}
                cert.update(kind='band to a ribbon knot', surface_of_J=s, J_identified_as=h['J'],
                            verification='self-contained up to Reidemeister moves of J')
                cert['recheck'] = certificate_recheck(h['diagram_pd'], cert, E_K)
                certs.append(cert)
                break
    keep, per = [], Counter(job['slice_kept'])      # at most SLICE_PER_J hits per knot X are rechecked and stored
    for h in hits:
        per[h['J']] += 1
        if per[h['J']] <= SLICE_PER_J:
            h['recheck'] = slice_recheck(h['diagram_pd'], h, E_K)
            keep.append(h)
    return dict(knot=name, round=job['round'], kind=job['kind'], diagrams=scanned, certificates=certs,
                knotinfo_slice=keep, knotinfo_slice_count=len(hits),
                knotinfo_slice_J=dict(Counter(h['J'] for h in hits)), slice_extra=extra, counts=dict(st),
                seconds=round(time.time() - t0, 1), cpu=round(time.process_time() - c0, 1))


def task(job):
    """run_round in a worker; an exception becomes the record {error: traceback}."""
    try:
        return run_round(job)
    except Exception:
        import traceback
        return dict(knot=job['knot'], round=job['round'], kind=job['kind'], error=traceback.format_exc())


# ---------- the main process: rounds, checkpoints, progress

def params():
    return dict(S_max=S_MAX, R_max=R_MAX, S2_max=S2_MAX, unlink_tries=UNLINK_TRIES, chunk=CHUNK, new_rounds=NEW_ROUNDS, x_max=X_MAX,
                bases=BASES, attempts=ATTEMPTS, slice_extra=SLICE_EXTRA, j_variants=J_VARIANTS, x_diagrams=X_DIAGRAMS,
                x_found=X_FOUND, j_reach=J_REACH, slice_per_j=SLICE_PER_J, controls=CONTROLS,
                only=os.environ.get('MOS_ONLY'))


def knot_states(only=None):
    """Targets and controls with their stored diagrams (all diagrams of mode mobius, the KnotInfo diagram first)."""
    T = targets()
    old = data('klein_bottle_search_mobius.json')['knots']
    C = controls(old)
    names = [k for k in list(T) + C if not only or k in only]
    K = {}
    for k in names:
        o = old[k]
        assert (o['group'] == 'mobius' and not o['moebius_band_of_K']) if k in T or k in SETTLED else o['group'] == 'control'
        assert o['diagrams'][0]['source'] == 'knotinfo' and all(d['isometric_to_knotinfo'] for d in o['diagrams'])
        diags = [dict(index=i, source=d['source'], pd=d['pd'], crossings=len(d['pd']), stored=True)
                 for i, d in enumerate(o['diagrams'])]
        stored_rounds = 1 if k in C else -(-len(diags) // CHUNK)
        K[k] = dict(group='target' if k in T else 'control',
                    kind=T.get(k) or ('control (exact correction terms)' if k in SETTLED else f"control ({o['test']})"),
                    range=o['range'], h_witt=o['h_witt'], det=o['det'], crossings=knot_order(k)[0], pd0=diags[0]['pd'],
                    diagrams=diags, keys={gauss_key([tuple(x) for x in d['pd']]) for d in diags},
                    stored_rounds=stored_rounds, rounds=stored_rounds + (1 if k in C else NEW_ROUNDS), next_round=0,
                    certificates=[], knotinfo_slice=[], knotinfo_slice_count=0, knotinfo_slice_J=Counter(),
                    slice_extra_done=set(), errors=[], attempts=Counter(), cpu=0.0, done=False)
    return K


def bases(k):
    """Base diagrams of the random walks: the least crossing number first, then the next two, at most BASES."""
    c = min(d['crossings'] for d in k['diagrams'])
    ds = sorted((d for d in k['diagrams'] if d['crossings'] <= c + 2), key=lambda d: (d['crossings'], d['index']))
    return [(d['index'], d['pd']) for d in ds[:BASES]]


def make_job(name, k):
    r = k['next_round']
    job = dict(knot=name, round=r, det=k['det'], crossings=k['crossings'], pd0=k['pd0'], group=k['group'],
               slice_extra_done=sorted(k['slice_extra_done']), slice_kept=dict(Counter(h['J'] for h in k['knotinfo_slice'])))
    if r < k['stored_rounds']:
        return dict(job, kind='stored', diagrams=[(d['index'], d['source'], d['pd'])
                                                  for d in k['diagrams'] if d.get('stored')][r * CHUNK:(r + 1) * CHUNK])
    return dict(job, kind='new', count=CHUNK, kinds=['walk', 'walk_level', 'walk', 'walk_global'], bases=bases(k),
                seen=[[list(x) for x in key] for key in k['keys']],
                seed=zlib.crc32(f'{name} {r} {k["attempts"][r]}'.encode()))


def apply(K, res):
    """Take the result of a round into the state of its knot (live, or from the .partial file)."""
    k = K[res['knot']]
    k['attempts'][res['round']] += 1
    if 'error' in res:
        k['errors'].append(dict(round=res['round'], error=res['error'][-2000:]))
        if k['attempts'][res['round']] < 3:
            return
    else:
        k['cpu'] += res['cpu']
        index = {}
        for pos, d in enumerate(res['diagrams']):
            if 'index' in d:
                k['diagrams'][d['index']].update(d)
            else:
                d['index'] = index[pos] = len(k['diagrams'])
                k['diagrams'].append(d)
                k['keys'].add(gauss_key([tuple(x) for x in d['pd']]))
        for x in res['certificates'] + res['knotinfo_slice']:
            if 'position' in x:
                x['diagram'] = index[x.pop('position')]
        k['certificates'] += res['certificates']
        k['knotinfo_slice'] += res['knotinfo_slice']
        k['knotinfo_slice_count'] += res['knotinfo_slice_count']
        k['knotinfo_slice_J'].update(res['knotinfo_slice_J'])
        k['slice_extra_done'].update(res['slice_extra'])
    k['next_round'] = res['round'] + 1
    if (k['group'] == 'target' and k['certificates']) or k['next_round'] >= k['rounds']:
        k['done'] = True


def read_partial(K):
    """Records of the .partial file (the first line holds the parameters, which must be the current ones)."""
    n = 0
    with open(PARTIAL) as f:
        meta = json.loads(f.readline())
        assert meta['params'] == params(), ('parameters differ from those of the .partial file', meta['params'])
        for line in f:
            if line.strip():
                res = json.loads(line)
                if res['knot'] in K:
                    apply(K, res); n += 1
    return meta, n


def status_text(K, t_start, workers):
    T = {k: v for k, v in K.items() if v['group'] == 'target'}
    sc = [k for k, v in T.items() if v['certificates']]
    by = Counter(v['certificates'][0]['verification'] for v in T.values() if v['certificates'])
    ks = [k for k, v in T.items() if not v['certificates'] and v['knotinfo_slice_count']]
    done = sum(v['done'] for v in K.values())
    left = sum(v['rounds'] - v['next_round'] for v in K.values() if not v['done'])
    rounds = sum(v['next_round'] for v in K.values())
    cpu = sum(v['cpu'] for v in K.values())
    per = cpu / rounds if rounds else None
    eta = per * left / workers / 3600 if per else None
    ctl = {k: v for k, v in K.items() if v['group'] == 'control'}
    lines = [f'mobius_open_search  {time.strftime("%Y-%m-%d %H:%M")}  (started {time.strftime("%Y-%m-%d %H:%M", time.localtime(t_start))})',
             f'knots finished {done} / {len(K)} (targets {sum(v["done"] for v in T.values())} / {len(T)}, controls '
             f'{sum(v["done"] for v in ctl.values())} / {len(ctl)});  rounds done {rounds}, at most {left} to go',
             f'diagrams scanned {sum(len([d for d in v["diagrams"] if "counts" in d]) for v in K.values())}',
             f'CPU per round {per and round(per)} s;  remaining at most {eta and round(eta, 1)} h with {workers} workers',
             f'Mobius bands of targets (self-contained): {len(sc)} {dict(by)}: {sc}',
             f'targets with a KnotInfo-slice J only: {len(ks)}: ' +
             ', '.join(f'{k} ({",".join(sorted(K[k]["knotinfo_slice_J"]))})' for k in ks),
             f'controls with a Mobius band or a KnotInfo-slice J (contradictions): '
             f'{[k for k, v in ctl.items() if v["certificates"] or v["knotinfo_slice_count"]]}',
             f'errors: {[(k, e["round"]) for k, v in K.items() for e in v["errors"]]}']
    return '\n'.join(lines)


def summary(K):
    T = {k: v for k, v in K.items() if v['group'] == 'target'}
    rep = {k: v['certificates'][0] for k, v in T.items() if v['certificates']}
    ver = lambda x: [k for k, c in rep.items() if c['verification'] == x]
    ks = {k: dict(v['knotinfo_slice_J']) for k, v in T.items() if k not in rep and v['knotinfo_slice_count']}
    rest = [k for k in T if k not in rep and k not in ks]
    scanned = [d for v in K.values() for d in v['diagrams'] if 'counts' in d]
    hits = [h for v in K.values() for h in v['knotinfo_slice']]
    ctl = {k: v for k, v in K.items() if v['group'] == 'control'}
    nbands = lambda c: 1 + len((c.get('surface_of_J') or c).get('ribbon_bands', []))
    return dict(
        targets=len(T), moebius_band=len(rep), self_contained=ver('self-contained'),
        self_contained_up_to_reidemeister_moves_of_J=ver('self-contained up to Reidemeister moves of J'),
        knotinfo_slice_only=ks, remaining=rest,
        by_kind={g: dict(targets=sum(v['kind'] == g for v in T.values()),
                         moebius_band=sum(T[k]['kind'] == g for k in rep), knotinfo_slice_only=sum(T[k]['kind'] == g for k in ks))
                 for g in ('alternating', 'Montesinos', 'not alternating, not Montesinos')},
        representative=dict(kind=dict(Counter(c['kind'] for c in rep.values())),
                            bands=dict(Counter(nbands(c) for c in rep.values())),
                            s=dict(Counter(abs(c['s']) for c in rep.values())),
                            diagram_source=dict(Counter(c['diagram_source'] for c in rep.values())),
                            diagram_index={k: c['diagram'] for k, c in rep.items()}),
        certificates_failing_a_recheck=[[k, i] for k, v in K.items() for i, c in enumerate(v['certificates'])
                                        if not c['recheck']['passes']],
        knotinfo_slice_hits=dict(targets=sum(v['knotinfo_slice_count'] for v in T.values()), stored=len(hits),
                                 failing_a_recheck=[[h['J'], h['diagram']] for h in hits if not h['recheck']['passes']],
                                 chirality=dict(Counter(h['recheck']['chirality'] for h in hits)),
                                 J=dict(sum((v['knotinfo_slice_J'] for v in T.values()), Counter()).most_common())),
        controls=dict(knots=len(ctl), diagrams=sum(sum('counts' in d for d in v['diagrams']) for v in ctl.values()),
                      contradictions={k: dict(moebius_bands=len(v['certificates']), knotinfo_slice=dict(v['knotinfo_slice_J']))
                                      for k, v in ctl.items() if v['certificates'] or v['knotinfo_slice_count']}),
        diagrams=len(scanned), diagrams_by_source=dict(Counter(d['source'] for d in scanned)),
        diagrams_new=sum(not d.get('stored') for d in scanned),
        diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for d in scanned),
        counts=dict(sum((Counter(d['counts']) for d in scanned), Counter())),
        errors=[[k, e['round']] for k, v in K.items() for e in v['errors']])


def record(k):
    keep = ('index', 'source', 'parent', 'pd', 'crossings', 'isometric_to_knotinfo', 'counts', 'moebius_band', 'knotinfo_slice')
    rep = k['certificates'][0] if k['certificates'] else None
    return dict(group=k['group'], kind=k['kind'], range=k['range'], h_witt=k['h_witt'], det=k['det'], rounds=k['next_round'],
                diagrams=[{x: d[x] for x in keep if x in d} for d in k['diagrams'] if 'counts' in d],
                moebius_band=rep, verification=rep and rep['verification'], certificates=k['certificates'],
                knotinfo_slice=k['knotinfo_slice'], knotinfo_slice_count=k['knotinfo_slice_count'],
                knotinfo_slice_J=dict(k['knotinfo_slice_J']), slice_extra=sorted(k['slice_extra_done']),
                errors=k['errors'], cpu=round(k['cpu'], 1))


def run(workers, fresh):
    import snappy, spherogram, cypari, database_knotinfo
    from concurrent.futures.process import BrokenProcessPool
    only = os.environ.get('MOS_ONLY')
    K = knot_states(set(only.split(',')) if only else None)
    if fresh and os.path.exists(PARTIAL):
        os.replace(PARTIAL, PARTIAL + '.old')
    t_start = time.time()
    if os.path.exists(PARTIAL):
        meta, n = read_partial(K)
        t_start = meta['started']
        print(f'resumed from {PARTIAL}: {n} round records', flush=True)
    else:
        with open(PARTIAL, 'w') as f:
            f.write(json.dumps(dict(params=params(), started=t_start)) + '\n')
    init = knotinfo_data()
    order = {name: i for i, name in enumerate(K)}
    heap = [(k['next_round'], order[name], name) for name, k in K.items() if not k['done']]
    heapq.heapify(heap)
    running, crashes = {}, Counter()
    ex = ProcessPoolExecutor(workers, initializer=worker_init, initargs=(init,))

    def stop(signum, frame):
        raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, stop)
    try:
        while heap or running:
            while heap and len(running) < workers + 2:
                _, _, name = heapq.heappop(heap)
                running[ex.submit(task, make_job(name, K[name]))] = name
            finished, _ = wait(running, timeout=300, return_when=FIRST_COMPLETED)
            broken = False
            for fut in finished:
                name = running.pop(fut)
                try:
                    res = fut.result()
                except BrokenProcessPool:         # a worker died: run its round again in a new pool
                    broken = True
                    crashes[(name, K[name]['next_round'])] += 1
                    if crashes[(name, K[name]['next_round'])] < 5:
                        heapq.heappush(heap, (K[name]['next_round'], order[name], name))
                        continue
                    res = dict(knot=name, round=K[name]['next_round'], kind='?', error='worker died 5 times')
                with open(PARTIAL, 'a') as f:
                    f.write(json.dumps(res, separators=(',', ':'), default=str) + '\n')
                apply(K, res)
                k = K[name]
                print(f'{name} round {res["round"]} {res["kind"]}: ' + (
                    f'ERROR {res["error"][-300:]}' if 'error' in res else
                    f'diagrams {len(res["diagrams"])} cpu {res["cpu"]}s KnotInfo-slice J {res["knotinfo_slice_count"]} '
                    f'{res["knotinfo_slice_J"]} ' + ''.join(f'| MOBIUS BAND {c["kind"]} [{c["verification"]}] s {c["s"]} '
                                                            for c in res['certificates'])) +
                      (' -- done' if k['done'] else ''), flush=True)
                if not k['done']:
                    heapq.heappush(heap, (k['next_round'], order[name], name))
            if broken:
                for fut, name in running.items():
                    heapq.heappush(heap, (K[name]['next_round'], order[name], name))
                running = {}
                ex.shutdown(wait=False, cancel_futures=True)
                ex = ProcessPoolExecutor(workers, initializer=worker_init, initargs=(init,))
            with open(PROGRESS, 'w') as f:
                f.write(status_text(K, t_start, workers) + '\n')
    finally:
        for proc in list(getattr(ex, '_processes', {}).values()):
            proc.terminate()
        ex.shutdown(wait=False, cancel_futures=True)
    S = summary(K)
    meta = dict(params(), pari_stack=[PARI_STACK, PARI_SIZEMAX], started=t_start, seconds=round(time.time() - t_start),
                cpu_seconds=round(sum(k['cpu'] for k in K.values())),
                versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                              database_knotinfo=database_knotinfo.version(), cypari=cypari.__version__,
                              pari='.'.join(map(str, pari().version())), python=sys.version.split()[0]))
    with gzip.open(OUT, 'wt') as f:
        json.dump(dict(meta=meta, summary=S, knots={name: record(k) for name, k in K.items()}), f,
                  separators=(',', ':'), default=str)
    for f in (PARTIAL, PROGRESS):               # the output holds everything that the .partial file held
        os.remove(f)
    print('summary:', json.dumps({x: y for x, y in S.items() if x != 'representative'}, default=str))
    print('DONE: written', OUT, flush=True)


def finalize():
    """Write the output from the .partial file of a run that was stopped before all rounds were done: the knots without
    a certificate keep the rounds done so far (field rounds of each knot, rounds_done in meta.stopped)."""
    import snappy, spherogram, cypari, database_knotinfo
    only = os.environ.get('MOS_ONLY')
    K = knot_states(set(only.split(',')) if only else None)
    old, n = read_partial(K)
    S = summary(K)
    meta = dict(params(), pari_stack=[PARI_STACK, PARI_SIZEMAX], started=old['started'],
                seconds=round(os.path.getmtime(PARTIAL) - old['started']),
                cpu_seconds=round(sum(k['cpu'] for k in K.values())),
                versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                              database_knotinfo=database_knotinfo.version(), cypari=cypari.__version__,
                              pari='.'.join(map(str, pari().version())), python=sys.version.split()[0]),
                stopped=dict(records=n, knots_done=sum(k['done'] for k in K.values()), knots=len(K),
                             rounds_done=dict(sorted(Counter(k['next_round'] for k in K.values()).items()))))
    with gzip.open(OUT, 'wt') as f:
        json.dump(dict(meta=meta, summary=S, knots={name: record(k) for name, k in K.items()}), f,
                  separators=(',', ':'), default=str)
    print('summary:', json.dumps({x: y for x, y in S.items() if x != 'representative'}, default=str))
    print('stopped run:', json.dumps(meta['stopped']))
    print('written', OUT, flush=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if args[:1] == ['finalize']:
        finalize()
        sys.exit()
    if args[:1] == ['status']:
        only = os.environ.get('MOS_ONLY')
        K = knot_states(set(only.split(',')) if only else None)
        meta, n = read_partial(K)
        print(status_text(K, meta['started'], int(args[1]) if len(args) > 1 else 10))
        sys.exit()
    fresh = args[:1] == ['fresh']
    args = args[1:] if fresh else args
    run(int(args[0]) if args else 10, fresh)
