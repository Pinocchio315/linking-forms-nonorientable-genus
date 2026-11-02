#!/usr/bin/env python
"""[C] Proposition 3.2(a): correction terms of Sigma(K) from negative-definite plumbing trees, and the test of
[C] Proposition 2.6, for the open knots that are neither alternating nor Montesinos.

Targets.  Of the 375 knots with gamma_4 in {1, 2} that Proposition 3.1 and the first search of [C] Section 3.3 leave
open (data/klein_bottle_search_mobius.json: group 'mobius', no Moebius band of K), 316 are not alternating and have no
Montesinos notation in KnotInfo (data/open_gamma4_knotinfo_2026.10.5.json).  The follow-up paper on unknotting numbers
(arXiv:2610.25430; its repository negative-definite-plumbings is the argument of mode records) classified the double
branched covers of 62 of them.  This script classifies the other 254 by the same pipeline and applies the test of
Proposition 2.6 (rank1_dinv.test) wherever Sigma(K) or -Sigma(K) bounds a negative-definite plumbing tree with at most
two bad vertices.  All 375 have gamma_4 <= 2, so an obstructed knot has gamma_4 = 2.

Classes (those of the follow-up paper, Section 3.1):
  hyperbolic  Sigma(K) is verified hyperbolic;
  mixed       an incompressible torus cuts off a verified hyperbolic piece, so Sigma(K) has positive simplicial volume;
  nodef       a graph manifold, and the Neumann normal forms of Sigma(K) and -Sigma(K) both have p+ >= 1, so neither
              bounds a negative-definite plumbing tree (Neumann, Bull. Austral. Math. Soc. 40 (1989), Thm 1.2);
              nodef_alg when the two normal forms are read from an algebraic diagram;
  negdef      a negative-definite normal form (class (3)); the d-test applies when it has at most two bad vertices;
  undecided   otherwise; the record gives the reason.
Pipeline for the 254 (the stages of plumbing/classify47.py of the follow-up repository; the functions of that
repository that are used are copied here, each named in its docstring):
1. Sigma(K) is the 2-fold cyclic cover of the orbifold E(K)(2, 0) in SnapPy (never the (1, 0) filling of a cover of
   E(K)); its H_1 is compared with coker(V + V^T).  Regina's StandardTriangulation.recognise runs on randomised
   triangulations: 4 tries (reach_orb.classify), then 60 (recognise_more.one), then 150 for K and 150 for its mirror,
   every tenth simplified by simplifyExhaustive(1) (recognise_retry.one).  Regina names do not distinguish Sigma(K)
   from -Sigma(K) = Sigma(mirror of K).
2. (Sage) Covers without a name: verified hyperbolicity (verify_hyp.one).  A short geodesic is drilled and refilled,
   which gives another triangulation of Sigma(K) (H_1 checked), and SnapPy's verify_hyperbolicity certifies it.  When
   the 4 tries of verify_hyp.one certify nothing, 20 more are made, with the cover randomised until SnapPy finds a
   solution with shapes (drilling needs its dual curves).
3. Covers neither named nor certified: incompressible tori among the vertex normal surfaces in quadrilateral
   coordinates and the pieces cut along them (jsj.one; Jaco-Tollefson), and for covers without a geometric piece a
   search of the Pachner graph for a Regina name (pachner_recognise.one).  (Sage) The pieces with a geometric solution
   are verified hyperbolic (jsj_verify.one).
4. Normal forms (neumann/precheck.one).  For every collected Regina name of a star (SFS over S^2), a pair or a triple
   of Seifert pieces: the trees of both orientations (Hirzebruch-Jung legs and chains, node weights solved from
   |det| = D) whose name by the calibrated formula of the vendored neumann/formula.py is one of the Regina names of
   step 1 or of 16 more randomised recognitions of K and of its mirror.  APPLICABLE when one tree has p+ = 0, NOT
   APPLICABLE when both orientations give trees, all with p+ >= 1 (both trees are then listed, as in
   neumann/both_orientations.py), INCONCLUSIVE otherwise.
5. Algebraic diagrams (algebraic/; Lemma 2.7 of the follow-up paper).  For a cover left undecided by steps 1-4, an
   algebraic diagram is searched by random Reidemeister moves (alg_diagrams.search, seed 1), the plumbing tree of
   either orientation is read from it and reduced to its normal form (vendored algebraic_plumbing.py), and a second
   algebraic diagram, different from the first (seed 7, crosscheck_u.py), must give the same normal forms.  A knot
   whose KnotInfo diagram is algebraic gets the same computation on that diagram, as a check of steps 1-4.
6. Identification of every tree used for a class (neumann/identify_trees.one): the tree is triangulated as surgery on
   its plumbing link (vendored plumbing.py) and recognised by Regina, independently of formula.py; the name must be a
   name of Sigma(K), with 16 more recognitions of K and of its mirror when it is not among those collected.  Also the
   exact inertia (Descartes' rule on the characteristic polynomial), |det Q| = det K, coker Q = H_1(Sigma(K)), and the
   Casson-Walker identity L(Q) = +-R(K) (algebraic/cw_general.L_general against Mullins' R(K) = -D sigma/4 +
   D V'(-1)/(6 V(-1)) from the KnotInfo signature and Jones polynomial; graph_u1.R_of).
7. d-test.  Let Q be a negative-definite tree with |det Q| = det K and at most two bad vertices (-w(v) < deg(v)).  The
   plumbed 4-manifold X is simply connected, so lambda(+-Sigma(K)) = -Q^-1 (Gordon-Litherland, pp. 59-60), and
   d(dX, t) = max (c^2 + n)/4 over the characteristic covectors c of t (Ozsvath-Szabo, plumbed three-manifolds,
   Cor. 1.5, for at most two bad vertices).  The code of [C] applies with P = -Q: F = Form(P), test(F, (1, -1)),
   both orientations Z0 = +-dX.  This unsigned test is the criterion, since the Regina names fix Sigma(K) only up to
   orientation.  It is repeated by the independent implementation of montesinos_dinv_check.py (exact Fincke-Pohst
   minima; cyclic H_1 only), whose d-function must be the same multiset.
The 62 knots of the follow-up records keep the class, trees and identifications of those records (status() of the
follow-up's paper/make_tables.py) and get the same d-test.
Controls (mode controls).  (a) The d-test of every negative-definite normal form in the follow-up records; the 114
knots among them with gamma_4 = 1 must pass.  (b) The steps 1 (4 tries), 4, 5 (KnotInfo diagram only), 6 and 7 for
the non-alternating knots without Montesinos notation that have KnotInfo gamma_4 = 1 or a Moebius band of [C] and are
not in the follow-up classification.  No control may be obstructed.  (c) The knots of the records
rank1_dinv_targets/controls/mtargets/mcontrols.json whose KnotInfo diagram is algebraic: every negative-definite
algebraic normal form with at most two bad vertices must give the d-function of the filling used there (the
negative-definite Goeritz form; the star plumbing), up to sign, and the same verdict.
Obstructed knots (mode gap): the search of smooth_top_gap.py (a band to a knot with Delta = 1, hence a locally flat
Mobius band and gamma_4^top = 1 < gamma_4 = 2), with its rechecks.
Mode summary also compares the verdicts with those of an earlier, independent implementation of the test, recorded
in data/plumbing_dinv_routeP/: routeP_12.json for 12 of the targets and routeP_controls.json for the forms of
control (a) (fields pilot_routeP_12_reproduced and same_as_pilot).
Processes.  Every item runs in a fresh process, at most [workers] at a time, with a time limit per mode (TIMEOUT); the
record of an item that runs out of time or crashes says so.  Each mode writes its checkpoint
data/plumbing_dinv_<mode>.partial after every item and skips the items already there, so a mode can be resumed.
Usage: python scripts/plumbing_dinv.py records path/to/negative-definite-plumbings | recognise | jsj | normal | algebraic | identify | dtest | controls | gap [workers] | summary
       sage -python scripts/plumbing_dinv.py hyp | piece [workers]
  The modes in the order records, recognise, hyp, jsj, piece, normal, algebraic, identify, dtest, controls, gap,
  summary.  Sage needs SnapPy installed in it (README.md, Environment).
Output: data/plumbing_dinv.json (meta, summary, targets, controls, gap) and data/plumbing_dinv_gap.json (diagrams and
bands of the gap search), written by mode summary from the checkpoints.
"""
import itertools, json, math, os, re, sys, time, traceback, warnings
from collections import Counter
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
FU = os.path.join(ROOT, 'vendor', 'paper_followup', 'negative-definite-plumbings', 'plumbing')
sys.path[:0] = [HERE, FU, os.path.join(FU, 'neumann'), os.path.join(FU, 'algebraic')]
warnings.filterwarnings('ignore')

QUICK, MORE, RETRY, NAMES = 4, 60, 150, 16  # recognitions: reach_orb, recognise_more, recognise_retry, names_orb
HYP_TRIES, HYP_DRILLS, HYP_RANDOMISE = 4, 6, 8  # verify_hyp.one
HYP_MORE = 20                             # further tries when the four of verify_hyp.one certify nothing
JSJ_TRIANGULATIONS, PIECE_RANDOMISE = 12, 20  # jsj.smallest_sigma2, jsj_verify.one
PACHNER = ((2, 60000), (3, 150000))       # pachner_recognise.one: (height, triangulations visited)
NODE_RANGE = 40                           # precheck.R: node weights in [-R, R], the last one solved (|e| <= 3R)
ID_TRIES = (12, 48)                       # identify_trees.one
ALG_ITERS, ALG_SEEDS = 6000, (1, 7)       # alg_diagrams.ITERS; seeds of alg_diagrams.py and crosscheck_u.py
TIMEOUT = dict(recognise=1200, hyp=2400, jsj=3600, piece=2400, normal=3600, algebraic=3600, identify=2400, dtest=3600,
               controls=2400, gap=7200)
PARI_STACK, PARI_SIZEMAX = 2 ** 26, 2 ** 31
CLASSES = ('hyperbolic', 'mixed', 'nodef', 'nodef_alg', 'negdef', 'undecided')


# ---------- checkpoints and processes

def data(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


def partial(mode):
    return os.path.join(ROOT, 'data', f'plumbing_dinv_{mode}.partial')


def load(mode):
    f = partial(mode)
    return json.load(open(f)) if os.path.exists(f) else {}


def save(mode, obj):
    f = partial(mode)
    with open(f + '.tmp', 'w') as h:
        json.dump(obj, h, indent=0, default=str)
    os.replace(f + '.tmp', f)


def seed(name):
    """Seed Python's, SnapPy's and Regina's generators from the item name (spherogram also iterates over sets of objects
    hashed by address, so its random moves are not exactly reproducible)."""
    import random, zlib
    s = zlib.crc32(name.encode())
    random.seed(s)
    try:
        import snappy
        snappy.set_rand_seed(s)
    except Exception:
        pass
    try:
        import regina
        regina.RandomEngine.reseedWithDefault()
    except Exception:
        pass


def _child(fn, item, q):
    try:
        seed(item[0])
        out = fn(item)
    except BaseException:
        out = (item[0], dict(error=traceback.format_exc()[-3000:]))
    q.put(out)


def run(mode, fn, items, workers, show=()):
    """fn(item) -> (name, record) for the items (name = item[0]) not yet in the checkpoint of `mode`: one fresh process
    per item, at most `workers` at a time, the checkpoint rewritten after every item.  A process still running after
    TIMEOUT[mode] seconds is stopped; one that dies without a result is recorded as crashed."""
    import multiprocessing as mp, queue
    os.environ['PYTHONHASHSEED'] = '0'    # read by every new interpreter
    ctx = mp.get_context('spawn')
    out = load(mode)
    todo = [it for it in items if it[0] not in out]
    print(f'{mode}: {len(items)} items, {len(todo)} to do, {workers} workers', flush=True)
    q, running, t0 = ctx.Queue(), {}, time.time()

    def finish(name, rec):
        p, start = running.pop(name, (None, time.time()))
        rec.setdefault('seconds', round(time.time() - start, 1))
        out[name] = rec
        save(mode, out)
        keys = tuple(show) + ('error', 'timeout', 'crashed')
        print(f'[{time.time() - t0:6.0f}s] {mode} {name}: ' +
              json.dumps({k: rec[k] for k in keys if k in rec}, default=str)[:400], flush=True)
        return p
    while todo or running:
        while todo and len(running) < workers:
            it = todo.pop(0)
            p = ctx.Process(target=_child, args=(fn, it, q), daemon=True)
            p.start()
            running[it[0]] = (p, time.time())
        try:
            name, rec = q.get(timeout=1)
            p = finish(name, rec)
            if p is not None:
                p.join()
            continue
        except queue.Empty:
            pass
        for name, (p, start) in list(running.items()):
            if time.time() - start > TIMEOUT[mode]:
                p.kill(); p.join()
                finish(name, dict(timeout=TIMEOUT[mode]))
            elif p.exitcode not in (None, 0):
                finish(name, dict(crashed=p.exitcode))
    print(f'{mode}: done, {len(out)} records, {time.time() - t0:.0f}s', flush=True)
    return out


def pari_memory():
    from cypari import pari
    pari.allocatemem(PARI_STACK, PARI_SIZEMAX, silent=True)


# ---------- the double branched cover (step 1)

def mirror_pd(pd):
    import spherogram
    return [list(x) for x in spherogram.Link(pd).mirror().PD_code()]


def sigma2(E, randomisations=0):
    """Sigma(K): the 2-fold cyclic cover of the orbifold E(2, 0), E randomised first (reach_orb.classify)."""
    E2 = E.copy()
    for _ in range(randomisations):
        E2.randomize()
    E2.dehn_fill((2, 0))
    return E2.covers(2, cover_type='cyclic')[0]


def regina_name(C, exhaustive=False):
    """Regina's name of the closed manifold C, or None (recognise_retry.tries)."""
    import regina
    T = regina.Triangulation3(C.filled_triangulation()._to_string())
    T.intelligentSimplify()
    if exhaustive and T.size() <= 30:
        T.simplifyExhaustive(1)
    r = regina.StandardTriangulation.recognise(T)
    return r.manifold().name() if r and r.manifold() else None


def names_orb(pd, tries):
    """graph_u1.names_orb: Regina names of Sigma over `tries` triangulations (the first one not randomised)."""
    import snappy
    E = snappy.Link(pd).exterior()
    names = Counter()
    for t in range(tries):
        nm = regina_name(sigma2(E, 1 if t else 0))
        if nm:
            names[nm] += 1
    return dict(names)


def primary(invs):
    """Prime-power decomposition of a finite abelian group given by invariant factors."""
    out = []
    for x in invs:
        p = 2
        while x > 1:
            if x % p == 0:
                q = 1
                while x % p == 0:
                    x //= p; q *= p
                out.append(q)
            p += 1
    return sorted(out)


def recognise(item):
    """Step 1: reach_orb.classify (stop at the first name), recognise_more.one (60 tries; stop at try >= 8 once named),
    recognise_retry.one (150 tries for K, then 150 for its mirror if needed; stop after 3 names)."""
    n, pd, h1 = item
    import snappy
    E = snappy.Link(pd).exterior()
    names, homology, tries, geometric = Counter(), Counter(), Counter(), False
    for t in range(QUICK):
        C = sigma2(E, 1 if t else 0)
        homology[str(C.homology())] += 1
        tries['reach_orb'] += 1
        nm = regina_name(C)
        if nm:
            names[nm] += 1
            break
        st = C.solution_type()
        if ('positive' in st or 'negative' in st) and C.volume() > 0.5:
            geometric = True
    if not names:
        for t in range(MORE):
            C = sigma2(E, 1 if t else 0)
            tries['recognise_more'] += 1
            nm = regina_name(C)
            if nm:
                names[nm] += 1
                if t >= 8:
                    break
    if not names:
        for which, p in (('knot', pd), ('mirror', mirror_pd(pd))):
            E1 = snappy.Link(p).exterior()
            for t in range(RETRY):
                C = sigma2(E1, 1 + t % 3)
                tries['recognise_retry_' + which] += 1
                nm = regina_name(C, exhaustive=t % 10 == 9)
                if nm:
                    names[nm] += 1
                    if sum(names.values()) >= 3:
                        break
            if names:
                break
    C = sigma2(E)
    ok = primary(C.homology().elementary_divisors()) == primary(h1)
    return n, dict(names=dict(names), H1=dict(homology), H1_agrees_with_seifert=ok, geometric_solution=geometric,
                   tries=dict(tries), kinds=sorted({parse_name(x)[0] or 'other' for x in names}))


# ---------- verified hyperbolicity (step 2, Sage)

def hyp(item):
    """Step 2: verify_hyp.one (4 tries).  If these certify nothing, 20 more tries, in which the cover is also
    randomised until SnapPy finds a solution with shapes, since drilling needs its dual curves."""
    n, pd = item
    import sage.all  # noqa: F401  (interval arithmetic for SnapPy's verify module)
    import sage.version, snappy
    E = snappy.Link(pd).exterior()
    rec = dict(certified=False, versions=dict(sage=sage.version.version, snappy=snappy.__version__))
    for t in range(HYP_TRIES + HYP_MORE):
        C = sigma2(E, 1 if t else 0)
        if t >= HYP_TRIES:
            for s in range(HYP_RANDOMISE):
                if any(x in C.solution_type() for x in ('positively', 'negatively')):
                    break
                C.randomize()
        H = str(C.homology())
        for k in range(HYP_DRILLS):
            try:
                F = C.drill(k).filled_triangulation([0])
            except Exception:
                continue
            F.dehn_fill((1, 0))
            for s in range(HYP_RANDOMISE):
                if 'positively' in F.solution_type():
                    break
                F.randomize()
            if 'positively' not in F.solution_type() or str(F.homology()) != H:
                continue
            try:
                ok = F.verify_hyperbolicity()
                ok = bool(ok[0] if isinstance(ok, tuple) else ok)
            except Exception as e:
                rec['verify_error'] = repr(e)[:120]; ok = False
            if ok:
                rec.update(certified=True, tries=t + 1, phase=1 if t < HYP_TRIES else 2, drilled=k,
                           volume=float(F.volume()), H1=H)
                return n, rec
    return n, rec


# ---------- tori, pieces, Pachner graph (step 3)

def smallest_sigma2(pd, k=JSJ_TRIANGULATIONS):
    """jsj.smallest_sigma2."""
    import snappy, regina
    E = snappy.Link(pd).exterior()
    best = None
    for t in range(k):
        C = sigma2(E, 1 if t else 0)
        T = regina.Triangulation3(C.filled_triangulation()._to_string())
        T.intelligentSimplify()
        if best is None or T.size() < best.size():
            best = T
    return best


def describe_piece(P):
    """jsj.describe_piece."""
    import regina, snappy
    out = dict(tetrahedra=P.size(), boundary_components=P.countBoundaryComponents())
    try:
        out['solid_torus'] = bool(P.isSolidTorus())
    except Exception:
        out['solid_torus'] = None
    try:
        out['TxI'] = bool(P.isTxI()) if P.countBoundaryComponents() == 2 else False
    except Exception:
        out['TxI'] = None
    rec = regina.StandardTriangulation.recognise(P)
    out['regina'] = rec.manifold().name() if rec and rec.manifold() else None
    I = regina.Triangulation3(P)
    I.finiteToIdeal()
    I.intelligentSimplify()
    out['ideal_isosig'] = I.isoSig()
    try:
        M = snappy.Manifold(I.isoSig())
        st = None
        for t in range(10):
            st = M.solution_type()
            if 'positively' in st:
                break
            M.randomize()
        out['snappy_solution'] = st
        out['snappy_volume'] = float(M.volume())
        try:
            out['snappy_identify'] = [str(x) for x in M.identify()][:3]
        except Exception:
            out['snappy_identify'] = []
    except Exception as e:
        out['snappy_error'] = repr(e)[:80]
    return out


def explore(T, height, limit):
    """pachner_recognise.explore: the first Regina name among the triangulations within `height` of T."""
    import regina
    found, count = {}, [0]

    def action(sig, tri):
        count[0] += 1
        r = regina.StandardTriangulation.recognise(tri)
        if r and r.manifold():
            found[r.manifold().name()] = 1
            return True
        return count[0] >= limit
    T.retriangulate(height, 1, action)
    return found, count[0]


def jsj(item):
    """Step 3 (jsj.one, then pachner_recognise.one when no piece has a geometric solution)."""
    n, pd = item
    import regina
    T = smallest_sigma2(pd)
    rec = dict(tetrahedra=T.size(), H1=str(T.homology()))
    try:
        rec['irreducible'] = bool(T.isIrreducible())
    except Exception:
        rec['irreducible'] = None
    S = regina.NormalSurfaces(T, regina.NS_QUAD)
    tori = []
    for i in range(S.size()):
        s = S.surface(i)
        if not s.isCompact() or not s.isConnected() or s.eulerChar() != 0:
            continue
        if s.isOrientable() and s.isTwoSided():
            tori.append(('torus', i, s))
        elif not s.isOrientable() and not s.isTwoSided():
            tori.append(('klein-bottle double', i, s.doubleSurface()))
    rec['vertex_surfaces'], rec['candidate_tori'] = S.size(), len(tori)
    found = []
    for kind, i, s in tori:
        try:
            inc = bool(s.isIncompressible())
        except Exception:
            inc = None
        if not inc:
            continue
        cut = s.cutAlong()
        cut.intelligentSimplify()
        comps = cut.triangulateComponents() if hasattr(cut, 'triangulateComponents') else [cut]
        found.append(dict(kind=kind, surface=i, pieces=[describe_piece(regina.Triangulation3(c)) for c in comps]))
        if len(found) >= 3:
            break
    rec['incompressible_tori'] = found
    if not found:
        rec['verdict'] = 'atoroidal (no incompressible vertex torus)'
    elif any('positively' in (p.get('snappy_solution') or '') for f in found for p in f['pieces']):
        rec['verdict'] = 'toroidal, a piece with a geometric solution (likely hyperbolic piece: mixed)'
    else:
        rec['verdict'] = 'toroidal, no geometric piece found (graph manifold candidate)'
    if not rec['verdict'].startswith('toroidal, a piece'):
        T2 = smallest_sigma2(pd, 16)
        names, visited = {}, 0
        for height, limit in PACHNER:
            names, v = explore(T2, height, limit)
            visited += v
            if names:
                break
        rec['pachner'] = dict(names=sorted(names), visited=visited)
    return n, rec


def piece(item):
    """Step 3, verification of a geometric piece (jsj_verify.one, Sage)."""
    n, isosigs = item
    import sage.all  # noqa: F401
    import sage.version, snappy
    versions = dict(sage=sage.version.version, snappy=snappy.__version__)
    for sig in isosigs:
        try:
            M = snappy.Manifold(sig)
        except Exception:
            continue
        for t in range(PIECE_RANDOMISE):
            if 'positively' in M.solution_type():
                break
            M.randomize()
        if 'positively' not in M.solution_type():
            continue
        try:
            ok = M.verify_hyperbolicity()
            ok = bool(ok[0] if isinstance(ok, tuple) else ok)
        except Exception:
            ok = False
        if ok:
            return n, dict(certified_hyperbolic_piece=True, isosig=sig, volume=float(M.volume()), cusps=M.num_cusps(),
                           versions=versions)
    return n, dict(certified_hyperbolic_piece=False, versions=versions)


# ---------- normal forms from Regina names (step 4): general.py, invert.py, precheck.py, both_orientations.py

def cf(alpha, beta):
    """general.cf: Hirzebruch-Jung continued fraction of alpha/beta."""
    out = []
    while beta:
        a = -(-alpha // beta); out.append(a); alpha, beta = beta, a * beta - alpha
    return out


def parse_name(name):
    """general.parse_name: (kind, [(base, fibres)]) with kind star, pair, triple or None."""
    pieces = re.findall(r'SFS \[([^:]+): ([^\]]*)\]', name)
    out = []
    for base, fib in pieces:
        fibres = [tuple(int(x) for x in f.split(',')) for f in re.findall(r'\((-?\d+,-?\d+)\)', fib)]
        out.append((base.strip(), fibres))
    kind = 'star' if len(out) == 1 else ('pair' if len(out) == 2 else ('triple' if len(out) == 3 else None))
    if kind == 'star' and out[0][0] != 'S2':
        kind = None
    return kind, out


def plumbing_matrix(w, edges):
    Q = [[0] * len(w) for _ in w]
    for i, x in enumerate(w):
        Q[i][i] = x
    for a, b in edges:
        Q[a][b] = Q[b][a] = 1
    return Q


def leg_of(alpha, beta):
    """invert.leg_of."""
    return [-x for x in cf(alpha, beta % alpha)]


def build(nodes_w, legs, chains):
    """invert.build: the tree from the node weights (a row), the legs of each node and the chains between consecutive
    nodes."""
    w = list(nodes_w); edges = []; nodes = list(range(len(nodes_w)))

    def add(x):
        w.append(x); return len(w) - 1
    for v, ls in zip(nodes, legs):
        for leg in ls:
            prev = v
            for x in leg:
                y = add(x); edges.append((prev, y)); prev = y
    for i, ch in enumerate(chains):
        prev = nodes[i]
        for x in ch:
            y = add(x); edges.append((prev, y)); prev = y
        edges.append((prev, nodes[i + 1]))
    return w, edges, nodes


def det_of(w, edges):
    """invert.det_of (floating point; used only to propose node weights, every accepted tree is checked exactly)."""
    import numpy as np
    return round(np.linalg.det(np.array(plumbing_matrix(w, edges), dtype=float)))


TRIPLE_CONV = ('J', 'b', False, False, True, False)       # invert.TRIPLE_CONV


def name_of(kind, w, edges, nodes):
    """invert.name_of: the Regina name of the boundary of a tree by the calibrated formula (formula.py)."""
    import formula as F
    legs, chains = F.decompose(w, edges, nodes)
    if kind == 'star':
        return F.star_name(w[0], legs[0])
    if kind == 'pair':
        return F.pair_name(w[nodes[0]], legs[nodes[0]], w[nodes[1]], legs[nodes[1]], chains[0], ('J', 'b', 1, 1, False))
    import regina
    E, V, t0, t1, inv0, inv1 = TRIPLE_CONV
    M0 = F.chain_matrix(chains[0], E, V); M1 = F.chain_matrix(chains[1], E, V)
    d = M0[0][0] * M0[1][1] - M0[0][1] * M0[1][0]
    A = [[d * M0[1][1], -d * M0[0][1]], [-d * M0[1][0], d * M0[0][0]]]
    g = regina.GraphTriple(F.piece(w[nodes[0]], legs[nodes[0]], 1), F.piece(w[nodes[1]], legs[nodes[1]], 2),
                           F.piece(w[nodes[2]], legs[nodes[2]], 1), regina.Matrix2(*A[0], *A[1]),
                           regina.Matrix2(*M1[0], *M1[1]))
    return g.name()


def mats(name):
    """precheck.mats: the gluing matrices of a Regina name."""
    pattern = r'\[ (-?\d+),(-?\d+) \| (-?\d+),(-?\d+) \]'
    return [[[int(a), int(b)], [int(c), int(d)]] for a, b, c, d in re.findall(pattern, name)]


def chains_for(p):
    """precheck.chains_for: the chains between two nodes whose fibres meet p times."""
    if p == 1:
        return [[]]
    return [[-x for x in cf(p, q)] for q in range(1, p) if math.gcd(p, q) == 1]


def index(w, edges):
    """precheck.index (floating point; the inertia of every tree used is recomputed exactly by exact_inertia)."""
    import numpy as np
    ev = np.linalg.eigvalsh(np.array(plumbing_matrix(w, edges), dtype=float))
    return int((ev > 1e-9).sum()), int((ev < -1e-9).sum())


def realisations(name, names, D):
    """precheck.realisations: all trees, for every orientation of the pieces, whose formula name lies in `names`."""
    R = NODE_RANGE
    kind, pieces = parse_name(name)
    out = []
    if kind not in ('star', 'pair', 'triple'):
        return kind, out
    for orient in itertools.product((1, -1), repeat=len(pieces)):
        legs = [[leg_of(a, o * b) for a, b in f if a > 1] for o, (_, f) in zip(orient, pieces)]
        if kind == 'star':
            chain_sets, ranges = [[]], []
        elif kind == 'pair':
            chain_sets = [[c] for c in chains_for(abs(mats(name)[0][0][1]))]
            ranges = [range(-R, R + 1)]
        else:
            m, n = mats(name)
            pm, pn = abs(m[0][1]), abs(n[0][1])
            chain_sets = [list(c) for c in itertools.product(chains_for(pm), chains_for(pn))]
            chain_sets += [list(c) for c in itertools.product(chains_for(pn), chains_for(pm))]
            ranges = [range(-R, R + 1), range(-R, R + 1)]
        for chs in chain_sets:
            for part in itertools.product(*ranges):
                w0, e0, _ = build(list(part) + [-1], legs, chs)
                w1, e1, _ = build(list(part) + [-2], legs, chs)
                d0, d1 = det_of(w0, e0), det_of(w1, e1)
                a = d0 - d1
                for target in (D, -D):
                    if a == 0 or (target - d0) % a:
                        continue
                    e_last = (target - d0) // a - 1
                    if abs(e_last) > 3 * R:
                        continue
                    w, edges, nodes = build(list(part) + [e_last], legs, chs)
                    try:
                        nm = name_of(kind, w, edges, nodes)
                    except Exception:
                        continue
                    if nm in names:
                        out.append((orient, w, edges, nodes, index(w, edges)[0]))
    return kind, out


def precheck(D, names):
    """precheck.one: the best tree (least p+, then fewest vertices) of every orientation vector, and the verdict."""
    names = set(names)
    found, kinds = {}, set()
    for nm in sorted(names):
        kind, rs = realisations(nm, names, D)
        kinds.add(kind)
        for orient, w, edges, nodes, pplus in rs:
            found.setdefault(orient, []).append((pplus, w, edges, nodes))
    best = {o: min(v, key=lambda x: (x[0], len(x[1]))) for o, v in found.items()}
    neg = lambda o: tuple(-x for x in o)
    if any(b[0] == 0 for b in best.values()):
        verdict = 'APPLICABLE'
    elif any(neg(o) in best for o in best):
        verdict = 'NOT APPLICABLE'
    else:
        verdict = 'INCONCLUSIVE'
    key = lambda o: ''.join('+' if x > 0 else '-' for x in o)
    return dict(kinds=sorted(k for k in kinds if k), verdict=verdict, pplus={key(o): b[0] for o, b in best.items()},
                normal_forms={key(o): dict(weights=b[1], edges=[list(e) for e in b[2]], nodes=b[3])
                              for o, b in best.items()})


def canon(w, edges):
    """both_orientations.canon: canonical string of a weighted tree (minimum over all roots)."""
    adj = {i: [] for i in range(len(w))}
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    enc = lambda v, p: '(' + str(w[v]) + ''.join(sorted(enc(u, v) for u in adj[v] if u != p)) + ')'
    return min(enc(r, -1) for r in range(len(w)))


def in_normal_form(w, edges):
    """both_orientations.in_normal_form: every vertex of valency <= 2 has weight <= -2."""
    deg = Counter(v for e in edges for v in e)
    return all(w[i] <= -2 for i in range(len(w)) if deg[i] <= 2)


def all_trees(D, names):
    """both_orientations.all_trees: every accepted tree, up to isomorphism."""
    trees = {}
    for nm in sorted(names):
        for orient, w, edges, nodes, pplus in realisations(nm, set(names), D)[1]:
            trees.setdefault(canon(w, edges), dict(weights=w, edges=[list(e) for e in edges], pplus=pplus))
    return trees


def normal(item):
    """Step 4 for one knot: names of step 1 (or of the Pachner search) and 16 more of K and of its mirror."""
    n, pd, D, names0 = item
    names = Counter(names0)
    for p in (pd, mirror_pd(pd)):
        names.update(names_orb(p, NAMES))
    rec = precheck(D, sorted(names))
    rec['names'] = dict(names)
    trees = {}
    for v in rec['normal_forms'].values():
        trees.setdefault(canon(v['weights'], v['edges']), dict(weights=v['weights'], edges=v['edges']))
    rec['recomputed_all'] = False
    if rec['verdict'] == 'NOT APPLICABLE' and len(trees) < 2:
        trees = all_trees(D, sorted(names)); rec['recomputed_all'] = True
    for t in trees.values():
        t['pplus'] = index(t['weights'], t['edges'])[0]
        t['in_normal_form'] = in_normal_form(t['weights'], t['edges'])
    rec['trees'] = list(trees.values())
    return n, rec


# ---------- algebraic diagrams (step 5): alg_diagrams.py, crosscheck_u.py, u_plumbing.py

def alg_search(pd0, seed, avoid=None):
    """alg_diagrams.search: an algebraic diagram with the fewest crossings found by random Reidemeister moves."""
    import random, spherogram
    from basic_polyhedron import reduce_diagram
    random.seed(seed)
    K0 = spherogram.Link(pd0)
    found = []
    for it in range(ALG_ITERS):
        K = K0.copy()
        K.backtrack(random.randint(5, 40))
        for mode in ('level', 'basic', 'global'):
            L = K.copy()
            try:
                L.simplify(mode)
            except Exception:
                continue
            pd = [list(x) for x in L.PD_code()]
            if len(L.link_components) == 1 and reduce_diagram(pd) == 1 and not (avoid and avoid(pd)):
                found.append(dict(pd=pd, crossings=len(pd), iteration=it, mode=mode, seed=seed))
        if found and it >= 200:
            break
    return min(found, key=lambda d: (d['crossings'], d['iteration'])) if found else None


def face_sizes(pd):
    """crosscheck_u.face_sizes."""
    ends = {}
    for c, x in enumerate(pd):
        for i, e in enumerate(x):
            ends.setdefault(e, []).append((c, i))
    seen, sizes = set(), []
    for c0 in range(len(pd)):
        for i0 in range(4):
            if (c0, i0) in seen:
                continue
            c, i, k = c0, i0, 0
            while (c, i) not in seen:
                seen.add((c, i)); k += 1
                (c, i), = [p for p in ends[pd[c][(i + 1) % 4]] if p != (c, (i + 1) % 4)]
            sizes.append(k)
    assert len(sizes) == len(pd) + 2
    return sorted(sizes)


def diagram_form(pd):
    """crosscheck_u.diagram_form: canonical form of a knot diagram on the sphere."""
    m = {e: i for i, e in enumerate(sorted({e for x in pd for e in x}))}
    n2 = len(m)
    X = [[m[e] for e in x] for x in pd]
    assert n2 == 2 * len(pd) and all(x[2] == (x[0] + 1) % n2 for x in X)
    best = None
    for refl in (False, True):
        for rev in (False, True):
            for t in range(n2):
                T = []
                for x in X:
                    y = [x[0], x[3], x[2], x[1]] if refl else list(x)
                    if rev:
                        y = y[2:] + y[:2]
                    T.append(tuple((t - e) % n2 if rev else (e - t) % n2 for e in y))
                T = sorted(T)
                if best is None or T < best:
                    best = T
    return best


def differs(pd1, pd2):
    """crosscheck_u.differs."""
    if len(pd1) != len(pd2):
        return 'crossings'
    if face_sizes(pd1) != face_sizes(pd2):
        return 'face sizes'
    return 'diagram' if diagram_form(pd1) != diagram_form(pd2) else False


def exact_inertia(w, edges):
    """identify_trees.exact_inertia: (p+, p-, nullity) by Descartes' rule on the characteristic polynomial, and det."""
    import sympy
    n = len(w)
    Q = sympy.Matrix(plumbing_matrix(w, edges))
    x = sympy.Symbol('x')
    c = [int(t) for t in sympy.Poly(Q.charpoly(x).as_expr(), x).all_coeffs()]
    changes = lambda s: sum(1 for u, v in zip(s, s[1:]) if u * v < 0)
    nonzero = [t for t in c if t]
    zero = len(c) - 1 - max(i for i, t in enumerate(c) if t)
    neg = [t * (-1) ** (len(c) - 1 - i) for i, t in enumerate(c)]
    return dict(pplus=changes(nonzero), pminus=changes([t for t in neg if t]), nullity=zero, det=int(Q.det()))


def bad_vertices(w, edges):
    """u_plumbing.bad_vertices: -w(v) < deg(v)."""
    deg = Counter(v for e in edges for v in e)
    return [i for i in range(len(w)) if -w[i] < deg[i]]


def smith(w, edges):
    """Invariant factors > 1 of coker Q (exact)."""
    from bands import invariants
    return list(invariants(plumbing_matrix(w, edges)))


def parse_jones(s):
    """kinfo.parse_jones of the follow-up, on the KnotInfo string itself."""
    s = s.replace(' ', '').replace('(', '').replace(')', '')
    s = s.replace('^-', '^~').replace('/-', '/~')
    s = s.replace('-', '+-').replace('~', '-')
    out = {}
    for term in s.split('+'):
        if not term:
            continue
        m = re.fullmatch(r'(-?\d*)\*?(?:t(?:\^(-?\d+))?)?', term)
        if m:
            c, e = m.group(1), m.group(2)
            coef = int(c) if c not in ('', '-') else (-1 if c == '-' else 1)
            exp = int(e) if e is not None else (1 if 't' in term else 0)
            out[exp] = out.get(exp, 0) + coef
            continue
        m = re.fullmatch(r'(-?\d*)/t\^?(-?\d+)?', term)
        if m:
            c, e = m.group(1), m.group(2)
            coef = int(c) if c not in ('', '-') else (-1 if c == '-' else 1)
            exp = -(int(e) if e is not None else 1)
            out[exp] = out.get(exp, 0) + coef
            continue
        raise ValueError(f'cannot parse Jones term {term!r}')
    return {e: c for e, c in out.items() if c}


def R_of(jones, sigma, D):
    """graph_u1.R_of: Mullins' R(K) = -D sigma/4 + D V'(-1)/(6 V(-1)) (= -2 D lambda_CW(Sigma(K)))."""
    J = parse_jones(jones)
    V1 = sum(c * (-1) ** (e % 2) for e, c in J.items())
    dV1 = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in J.items())
    assert abs(V1) == D
    return Fraction(-D * sigma, 4) + Fraction(D * dV1, 6 * V1)


def L_general(w, edges):
    """cw_general.L_general: (|det Q|/12)(tr Q - 3 sigma(Q) + sum_v (2 - deg v)(Q^-1)_vv), the Casson-Walker invariant
    in the normalisation of R(K)."""
    import sympy
    Q = sympy.Matrix(plumbing_matrix(w, edges))
    n = Q.rows; H = abs(Q.det()); Qi = Q.inv()
    deg = Counter(v for e in edges for v in e)
    I = exact_inertia(w, edges)
    assert I['nullity'] == 0
    sig = I['pplus'] - I['pminus']
    val = H * (sum(w) - 3 * sig + sum((2 - deg[i]) * Qi[i, i] for i in range(n))) / 12
    num, den = sympy.fraction(sympy.nsimplify(val))
    return Fraction(int(num), int(den))


def tree_record(w, edges, D, h1, R):
    """The exact data of a tree: inertia, det, bad vertices, coker, Casson-Walker sign."""
    I = exact_inertia(w, edges)
    L = L_general(w, edges)
    return dict(weights=w, edges=[list(e) for e in edges], canon=canon(w, edges), vertices=len(w), **I,
                det_is_D=abs(I['det']) == D, bad=bad_vertices(w, edges), in_normal_form=in_normal_form(w, edges),
                coker=smith(w, edges), coker_is_H1=primary(smith(w, edges)) == primary(h1), L=str(L),
                cw_signs=[s for s in (1, -1) if s * L == R])


def algebraic(item):
    """Step 5 for one knot: the KnotInfo diagram if it is algebraic, and for an undecided cover two searched
    diagrams."""
    n, pd, D, h1, jones, sigma, search = item
    import spherogram
    from basic_polyhedron import reduce_diagram
    from algebraic_plumbing import plumbing_from_algebraic, normal_form
    R = R_of(jones, sigma, D)
    diagrams = []
    if reduce_diagram([list(x) for x in pd]) == 1:
        diagrams.append(dict(source='knotinfo', pd=pd))
    if search:
        first = alg_search(pd, ALG_SEEDS[0])
        if first:
            diagrams.append(dict(source='search', **first))
            second = alg_search(pd, ALG_SEEDS[1], avoid=lambda q: not differs(first['pd'], q))
            if second:
                diagrams.append(dict(source='second', differs_from_first=differs(first['pd'], second['pd']), **second))
    E0 = spherogram.Link(pd).exterior()
    for d in diagrams:
        d['forms'] = {}
        for s in (1, -1):
            w, e = plumbing_from_algebraic(d['pd'], s)
            wn, en = normal_form(w, e)
            d['forms'][str(s)] = dict(raw_vertices=len(w), raw_det=exact_inertia(w, e)['det'],
                                      **tree_record(wn, en, D, h1, R))
        if d['source'] != 'knotinfo':
            ok = None
            E1 = spherogram.Link(d['pd']).exterior()
            for t in range(8):
                try:
                    ok = bool(E0.is_isometric_to(E1)); break
                except Exception:
                    E1.randomize()
            d['isometric_to_knotinfo'] = ok
    canons = {tuple(sorted(f['canon'] for f in d['forms'].values())) for d in diagrams}
    rec = dict(R=str(R), diagrams=diagrams, same_normal_forms=len(canons) <= 1)
    if diagrams:
        f = diagrams[0]['forms']
        neg = [s for s in ('1', '-1') if f[s]['pplus'] == 0 and f[s]['pminus'] == f[s]['vertices']]
        rec['verdict'] = 'APPLICABLE' if neg else 'NOT APPLICABLE'
        rec['negdef_sign'] = neg
    else:
        rec['verdict'] = 'NO ALGEBRAIC DIAGRAM'
    return n, rec


# ---------- identification by triangulation (step 6)

def fast_regina_name(M, tries):
    """graph_u1.fast_regina_name: SnapPy randomisations and intelligentSimplify, the first Regina name."""
    import regina
    for t in range(tries):
        M2 = M.copy()
        if t:
            M2.randomize()
        T = regina.Triangulation3(M2.filled_triangulation()._to_string())
        T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec is not None and rec.manifold() is not None:
            return rec.manifold().name()
    return None


def regina_names(M, tries):
    """All Regina names of M over `tries` triangulations (the first not randomised).  Regina's names are not canonical
    (a small Seifert fibred space can also be named as a pair with a solid torus piece), but every name is a name of
    M."""
    import regina
    names = Counter()
    for t in range(tries):
        M2 = M.copy()
        if t:
            M2.randomize()
        T = regina.Triangulation3(M2.filled_triangulation()._to_string())
        T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec is not None and rec.manifold() is not None:
            names[rec.manifold().name()] += 1
    return dict(names)


def identify_tree(w, edges, pd, names, extra):
    """identify_trees.one for one tree: the first Regina name of the plumbed manifold (12, then 48 tries), compared with
    the names of Sigma(K), with 16 more recognitions of K and of its mirror when it is not among them.  If it still
    does not match, all names of 48 triangulations of the plumbed manifold are compared."""
    from plumbing import plumbed_manifold
    M = plumbed_manifold(w, [tuple(e) for e in edges])
    nm = None
    for tries in ID_TRIES:
        nm = fast_regina_name(M, tries)
        if nm:
            break
    matched = nm in names
    if not matched and extra is None:
        extra = set(names_orb(pd, NAMES)) | set(names_orb(mirror_pd(pd), NAMES))
    if nm and not matched:
        matched = nm in extra
    rec = dict(name=nm, matched=matched)
    if not matched:
        more = regina_names(M, ID_TRIES[1])
        rec.update(more_names=more, matched=bool(set(more) & (set(names) | extra)))
    return rec, extra


def identify(item):
    """Step 6 for the trees of one knot, with their exact data."""
    n, pd, D, h1, jones, sigma, names, trees = item
    names, extra = set(names), None
    R = R_of(jones, sigma, D)
    out = []
    for role, w, edges in trees:
        idn, extra = identify_tree(w, edges, pd, names, extra)
        out.append(dict(role=role, **idn, **tree_record(w, edges, D, h1, R)))
    return n, dict(R=str(R), extra_names_used=extra is not None, trees=out)


# ---------- the d-test (step 7)

def dtest_tree(D, w, edges):
    """The census d-test (rank1_dinv) on a negative-definite tree, and the independent implementation of
    montesinos_dinv_check.py."""
    from rank1_dinv import Form, test
    from owens_obstruction import positive_definite, det_int
    import montesinos_dinv_check as MDC
    Q = plumbing_matrix(w, edges)
    P = [[-v for v in row] for row in Q]
    bad = bad_vertices(w, edges)
    rec = dict(vertices=len(w), bad=len(bad), negative_definite=bool(positive_definite(P)), det=abs(det_int(P)))
    rec['det_is_D'] = rec['det'] == D
    if not (rec['negative_definite'] and rec['det_is_D']):
        rec['eligible'] = False
        return rec
    rec['eligible'] = len(bad) <= 2
    if not rec['eligible']:
        return rec
    F = Form(P)
    zero = tuple(0 for _ in F.mods)
    res = test(F, (1, -1))
    rec.update(h1=list(F.mods), d_spin=str(F.d[zero]), d=sorted(str(v) for v in F.d.values()), passing=res,
               obstructed=res is None)
    try:
        N, u, d = MDC.d_invariants(Q)
        ob, passing = MDC.obstructed(N, u, d)
        rec['independent'] = dict(obstructed=ob, d_spin=str(d[0]), passing=[list(p) for p in passing[:4]],
                                  same_d=sorted(d.values()) == sorted(F.d.values()),
                                  agrees=ob == rec['obstructed'] and d[0] == F.d[zero])
    except MDC.NotCyclic:
        rec['independent'] = 'H_1 not cyclic'
    return rec


def star_check(names, D):
    """For a Seifert fibred Sigma(K) over S^2: the route of montesinos_dinv_check.py, which shares no code with steps
    4-7 (the star plumbing of the Seifert invariants of a Regina name, its exact d-function, the rank-one test for
    cyclic H_1)."""
    import montesinos_dinv_check as MDC
    for nm in sorted(names):
        fib = MDC.parse(nm)
        if not fib:
            continue
        Q, e = MDC.plumbing(fib)
        try:
            N, u, d = MDC.d_invariants(Q)
        except MDC.NotCyclic:
            return dict(name=nm, result='H_1 not cyclic')
        ob, passing = MDC.obstructed(N, u, d)
        return dict(name=nm, det_is_D=N == D, rank=len(Q), obstructed=ob, d=sorted(str(v) for v in d.values()),
                    passing=[list(p) for p in passing[:4]])
    return None


def dtest(item):
    n, D, w, edges, names = item
    rec = dtest_tree(D, w, [tuple(e) for e in edges])
    sc = star_check(names, D)
    if sc is not None:
        if 'obstructed' in sc and rec.get('eligible'):
            dv = sorted(Fraction(x) for x in rec['d'])
            ds = sorted(Fraction(x) for x in sc['d'])
            sc['agrees'] = sc['obstructed'] == rec['obstructed'] and (dv == ds or dv == sorted(-x for x in ds))
        rec['star_check'] = sc
    return n, rec


# ---------- controls

def control(item):
    """Control (b): steps 1 (4 tries), 4, 5 (KnotInfo diagram), 6 and 7 for one knot with gamma_4 = 1."""
    n, pd, D, h1, jones, sigma = item
    import snappy
    from basic_polyhedron import reduce_diagram
    from algebraic_plumbing import plumbing_from_algebraic, normal_form
    rec, R = {}, R_of(jones, sigma, D)
    E = snappy.Link(pd).exterior()
    names = Counter()
    for t in range(QUICK):
        nm = regina_name(sigma2(E, 1 if t else 0))
        if nm:
            names[nm] += 1
            break
    trees = []
    if names:
        for p in (pd, mirror_pd(pd)):
            names.update(names_orb(p, NAMES))
        pre = precheck(D, sorted(names))
        rec.update(names=dict(names), verdict=pre['verdict'], pplus=pre['pplus'])
        if pre['verdict'] == 'APPLICABLE':
            o = [k for k, v in pre['pplus'].items() if v == 0][0]
            f = pre['normal_forms'][o]
            trees.append(('regina', f['weights'], [tuple(e) for e in f['edges']]))
    if reduce_diagram([list(x) for x in pd]) == 1:
        rec['algebraic'] = {}
        for s in (1, -1):
            wn, en = normal_form(*plumbing_from_algebraic(pd, s))
            I = exact_inertia(wn, en)
            rec['algebraic'][str(s)] = dict(canon=canon(wn, en), pplus=I['pplus'])
            if I['pplus'] == 0 and I['nullity'] == 0:
                trees.append(('algebraic', wn, [tuple(e) for e in en]))
    rec['trees'] = []
    for role, w, edges in trees:
        t = dict(role=role, **tree_record(w, edges, D, h1, R))
        if role == 'regina':
            idn, _ = identify_tree(w, edges, pd, set(names), None)
            t.update(idn)
        t['dtest'] = dtest_tree(D, w, edges)
        rec['trees'].append(t)
    canons = {t['canon'] for t in rec['trees']}
    rec['routes_agree'] = len(canons) <= 1
    rec['obstructed'] = any(t['dtest'].get('obstructed') for t in rec['trees'])
    return n, rec


def control_known(item):
    """Control (c): a knot of the rank1_dinv records whose KnotInfo diagram is algebraic, with the d-function of the
    filling used there (alternating: the negative-definite Goeritz form; Montesinos: the star plumbing of
    montesinos_u1).  Every negative-definite algebraic normal form with at most two bad vertices must give the same
    d-function up to sign (as a multiset; one of the two is +-Sigma(K)) and the recorded verdict."""
    n, pd, D, kind, notation, recorded = item
    from algebraic_plumbing import plumbing_from_algebraic, normal_form
    from rank1_dinv import Form
    from owens_obstruction import positive_definite
    if kind == 'alternating':
        from rank1_goeritz import goeritz_pair
        G = [G for G in goeritz_pair(pd) if positive_definite([[-v for v in row] for row in G])]
        assert len(G) == 1
        F0 = Form([[-v for v in row] for row in G[0]])
    else:
        sys.path.insert(0, os.path.join(ROOT, 'vendor', 'code', 'lower_bounds', 'montesinos'))
        from montesinos_u1 import plumbing_for
        Q, flipped = plumbing_for(notation, D)
        F0 = Form([[-int(v) for v in row] for row in Q])
    ref = sorted(F0.d.values())
    rec = dict(kind=kind, recorded_obstructed=recorded, forms={})
    for s in (1, -1):
        wn, en = normal_form(*plumbing_from_algebraic([list(x) for x in pd], s))
        I = exact_inertia(wn, en)
        f = dict(canon=canon(wn, en), vertices=len(wn), pplus=I['pplus'], bad=len(bad_vertices(wn, en)))
        if I['pplus'] == 0 and I['nullity'] == 0:
            t = dtest_tree(D, wn, [tuple(e) for e in en])
            f.update(eligible=t['eligible'], det_is_D=t['det_is_D'])
            if t['eligible']:
                d = sorted(Fraction(x) for x in t['d'])
                f.update(obstructed=t['obstructed'], same_d=d == ref, same_d_reversed=d == sorted(-x for x in ref),
                         independent=(t['independent'] if isinstance(t['independent'], str)
                                      else t['independent']['agrees']))
        rec['forms'][str(s)] = f
    used = [f for f in rec['forms'].values() if 'obstructed' in f]
    rec['dtests'] = len(used)
    rec['agrees'] = all(f['obstructed'] == recorded and (f['same_d'] or f['same_d_reversed']) for f in used)
    return n, rec


# ---------- the gap search

def gap(item):
    """smooth_top_gap.one for one obstructed knot, with its rechecks."""
    pari_memory()
    import smooth_top_gap as stg
    from database_knotinfo import link_list
    LL = link_list()[1:]
    stg.init_worker({r['name']: {k: r[k] for k in stg.KI_FIELDS} for r in LL},
                    {'K' + r['dt_name'].replace('_', ''): r['name'] for r in LL if r.get('dt_name')})
    return stg.one(tuple(item))


# ---------- modes

def followup_records(fu_root):
    """Classes of the follow-up records (status() of paper/make_tables.py) with the trees and identifications of the
    knots of class (2) and (3), and every negative-definite normal form in the records (control (a))."""
    import hashlib, subprocess
    T = os.path.join(fu_root, 'plumbing')
    used = []

    def J(*p):
        f = os.path.join(T, *p)
        used.append(f)
        return json.load(open(f))
    VH, JV, JS = J('verify_hyp.json'), J('jsj_verify.json'), J('jsj.json')
    p4, r13, fresh = J('neumann', 'run_pachner4.json'), J('neumann', 'run_retry13.json'), J('neumann', 'run_fresh.json')
    nf, nc3 = J('neumann', 'rederive_all.json'), J('neumann', 'normal_case3.json')
    ps0 = J('neumann', 'precheck_sigma0.json')
    star = {**J('star_u1.json'), **J('star_u1_orb.json')}
    gu, owt = J('graph_u1.json'), J('owens_u3_targets.json')
    pre = {**J('neumann', 'precheck_u12.json'), **J('neumann', 'precheck_u23.json')}
    UA, CU = J('algebraic', 'u_plumbing.json'), J('algebraic', 'crosscheck_u.json')
    BO, IT = J('neumann', 'both_orientations.json'), J('neumann', 'identify_trees.json')
    C47 = J('classify47_summary.json')['status']
    C47n, C47h, C47p = J('classify47_normal.json'), J('classify47_hyp.json'), J('classify47_piece.json')
    studied = {x['name'] for f in ('reach_u12_orb.json', 'reach_u23_orb.json') for x in J(f)}
    NF = (p4, r13, fresh, nf, nc3)

    def status(n):                        # paper/make_tables.status
        if n in VH and VH[n]['certified']:
            return 'hyperbolic'
        if n in JV and JV[n]['certified_hyperbolic_piece']:
            return 'mixed'
        for src in NF:
            if n in src:
                r = src[n]
                if r['applicability'] == 'NOT APPLICABLE':
                    return 'nodef'
                if r['applicability'] != 'APPLICABLE':
                    return 'other'
                if str(r.get('verdict', '')).startswith('OBSTRUCTED'):
                    return 'obstructed'
                return 'notest' if 'sigma' in str(r.get('verdict', '')) else 'tested'
        for src in (star, gu, owt):
            if n in src and src[n]:
                return 'nodef' if pre.get(n, {}).get('verdict') == 'NOT APPLICABLE' else 'other'
        if n in UA:
            return 'nodef_alg'
        if n in JS and JS[n]['verdict'].startswith('toroidal, no geometric'):
            return 'graph4'
        if n in ps0:
            return {'APPLICABLE': 'notest', 'NOT APPLICABLE': 'nodef'}.get(ps0[n]['verdict'], 'other')
        return 'other'

    def negdef_form(n):                   # paper/make_tables.negdef_form, and the file of the form
        for src in NF:
            if n in src:
                return src[n]['plumbing'], 'neumann/' + ['run_pachner4', 'run_retry13', 'run_fresh', 'rederive_all',
                                                        'normal_case3'][NF.index(src)] + '.json'
        r, f = (ps0[n], 'neumann/precheck_sigma0.json') if n in ps0 else (C47n[n], 'classify47_normal.json')
        return r['normal_forms'][[o for o, v in r['pplus'].items() if v == 0][0]], f

    def record(n):
        if n in studied:
            st = status(n)
        elif n in C47:
            st = {'negdef': 'negdef', 'nodef': 'nodef'}.get(C47[n], C47[n])
        else:
            return None
        cls = {'obstructed': 'negdef', 'tested': 'negdef', 'notest': 'negdef'}.get(st, st)
        rec = dict(followup_status=st, cls=cls if cls in CLASSES else 'undecided')
        if cls == 'hyperbolic':
            rec['evidence'] = VH[n] if n in VH else C47h[n]
        elif cls == 'mixed':
            rec['evidence'] = JV[n] if n in JV else C47p[n]
        elif cls == 'negdef':
            p, src = negdef_form(n)
            rec.update(tree=dict(weights=p['weights'], edges=[list(e) for e in p['edges']]), tree_source=src,
                       identification=IT.get(n))
        elif cls == 'nodef':
            rec.update(trees=BO[n]['trees'], identification=IT.get(n))
        elif cls == 'nodef_alg':
            rec.update(trees={s: UA[n][s] for s in ('1', '-1')}, second_diagram=CU[n].get('second_same_normal_forms'))
        return rec
    forms = {}                             # control (a): every negative-definite normal form of the records
    for src, f in zip(NF, ('run_pachner4', 'run_retry13', 'run_fresh', 'rederive_all', 'normal_case3')):
        for n, r in src.items():
            if n not in forms and r.get('applicability') == 'APPLICABLE':
                forms[n] = dict(weights=r['plumbing']['weights'], edges=r['plumbing']['edges'],
                                source=f'neumann/{f}.json')
    for n, r in ps0.items():
        if n not in forms and r.get('verdict') == 'APPLICABLE':
            o = [k for k, v in r['pplus'].items() if v == 0][0]
            forms[n] = dict(r['normal_forms'][o], source='neumann/precheck_sigma0.json')
    for n, s in C47.items():
        if s == 'negdef' and n not in forms:
            r = C47n[n]
            o = [k for k, v in r['pplus'].items() if v == 0][0]
            forms[n] = dict(r['normal_forms'][o], source='classify47_normal.json')
    try:
        git = lambda *a: subprocess.run(['git', '-C', fu_root, *a], capture_output=True, text=True).stdout.strip()
        commit, dirty = git('rev-parse', 'HEAD'), git('status', '--porcelain')
    except Exception:
        commit, dirty = None, None
    files = {os.path.relpath(f, fu_root): hashlib.sha256(open(f, 'rb').read()).hexdigest() for f in sorted(set(used))}
    classified = studied | set(C47)
    return record, forms, classified, dict(commit=commit, clean=dirty == '', files=files)


def mode_records(fu_root):
    """Mode records: the 316 targets with their KnotInfo data and the records of the follow-up repository at fu_root,
    the knots of control (b), and the forms of control (a)."""
    from database_knotinfo import link_list
    from bands import invariants
    from basic_polyhedron import reduce_diagram
    BY = {r['name']: r for r in link_list()[1:]}
    km = data('klein_bottle_search_mobius.json')['knots']
    left = sorted(k for k, v in km.items() if v['group'] == 'mobius' and not v['moebius_band_of_K'])
    assert len(left) == 375
    op, G4 = data('open_gamma4_knotinfo_2026.10.5.json'), data('gamma4_knotinfo_2026.10.5.json')
    nonalt_nonmont = lambda n: BY[n]['alternating'] == 'N' and not BY[n]['montesinos_notation'].startswith('K(')
    targets = [k for k in left if op[k]['alternating'] != 'Y' and not op[k]['montesinos_notation'].startswith('K(')]
    assert len(targets) == 316 and all(nonalt_nonmont(k) for k in targets)
    record, forms, classified, provenance = followup_records(fu_root)
    mob = {k for k, v in km.items() if v.get('moebius_band_of_K')}

    def info(n):
        r = BY[n]
        V = json.loads(r['seifert_matrix'])
        VV = [[V[i][j] + V[j][i] for j in range(len(V))] for i in range(len(V))]
        return dict(pd=json.loads(r['pd_notation']), D=int(r['determinant']), sigma=int(r['signature']),
                    jones=r['jones_polynomial'], h1=list(invariants(VV)),
                    u=(op[n] if n in op else r)['unknotting_number'].replace(' ', ''),   # KnotInfo 2026.10.5 if open
                    gamma4=G4.get(n), census_band=n in mob,
                    knotinfo_algebraic=reduce_diagram(json.loads(r['pd_notation'])) == 1)
    T = {}
    for n in targets:
        T[n] = info(n)
        fr = record(n)
        T[n]['followup'] = fr
    assert sum(r['followup'] is not None for r in T.values()) == 62
    controls = sorted(n for n in BY if re.fullmatch(r'\d+[an]?_\d+', n) and nonalt_nonmont(n)
                      and (G4.get(n) == '1' or n in mob) and n not in classified)
    C = {n: info(n) for n in controls}
    out = dict(targets=T, controls=C, followup_forms=forms, provenance=provenance)
    save('records', out)
    print('targets', len(T), '| in the follow-up records', sum(r['followup'] is not None for r in T.values()),
          Counter(r['followup']['cls'] for r in T.values() if r['followup']),
          '| new', sum(r['followup'] is None for r in T.values()), '| controls (b)', len(C),
          '| follow-up negative-definite normal forms', len(forms), '| follow-up repository', provenance['commit'],
          'clean' if provenance['clean'] else 'NOT CLEAN')


def new_targets():
    """The 254 targets that are not in the follow-up records."""
    return {n: r for n, r in load('records')['targets'].items() if r['followup'] is None}


def current_class(n, R, H, JJ, PV, N, A):
    """Class of a new target from the checkpoints, as classify47.status, completed by the algebraic diagrams."""
    if H.get(n, {}).get('certified'):
        return 'hyperbolic', 'verified hyperbolic'
    if PV.get(n, {}).get('certified_hyperbolic_piece'):
        return 'mixed', 'verified hyperbolic piece'
    v = N.get(n, {}).get('verdict')
    if v == 'APPLICABLE':
        return 'negdef', 'normal form from Regina names'
    if v == 'NOT APPLICABLE':
        return 'nodef', 'normal forms from Regina names'
    a = A.get(n, {})
    if a.get('verdict') in ('APPLICABLE', 'NOT APPLICABLE') and any(d['source'] != 'knotinfo' for d in a['diagrams']):
        return ('negdef' if a['verdict'] == 'APPLICABLE' else 'nodef_alg'), 'normal forms from algebraic diagrams'
    if a.get('verdict') in ('APPLICABLE', 'NOT APPLICABLE'):
        return ('negdef' if a['verdict'] == 'APPLICABLE' else 'nodef_alg'), 'normal forms from the KnotInfo diagram'
    reasons = []
    if n in N:
        reasons.append(f'normal form verdict {v}')
    if n in JJ:
        reasons.append(JJ[n]['verdict'])
    if n in A:
        reasons.append(a.get('verdict') or 'algebraic stage failed')
    if not R.get(n, {}).get('names') and not H.get(n, {}).get('certified'):
        reasons.append('no Regina name and not verified hyperbolic')
    return 'undecided', '; '.join(reasons) or 'not processed'


def main():
    mode = sys.argv[1]
    args = sys.argv[2:]
    workers = int(args[0]) if args and args[0].isdigit() else 10
    if mode == 'records':
        return mode_records(sys.argv[2])
    T = new_targets()
    if mode == 'recognise':
        items = [(n, r['pd'], r['h1']) for n, r in T.items()]
        run(mode, recognise, items, workers, show=('names', 'H1_agrees_with_seifert', 'tries'))
    elif mode == 'hyp':
        R = load('recognise')
        items = [(n, T[n]['pd']) for n in T if not R[n]['names']]
        run(mode, hyp, items, workers, show=('certified', 'volume'))
    elif mode == 'jsj':
        R, H = load('recognise'), load('hyp')
        items = [(n, T[n]['pd']) for n in T if not R[n]['names'] and not H.get(n, {}).get('certified')]
        run(mode, jsj, items, workers, show=('verdict', 'pachner'))
    elif mode == 'piece':
        JJ = load('jsj')
        items = [(n, [p['ideal_isosig'] for t in r['incompressible_tori'] for p in t['pieces']
                      if 'positively' in (p.get('snappy_solution') or '')]) for n, r in JJ.items() if 'verdict' in r]
        run(mode, piece, [x for x in items if x[1]], workers, show=('certified_hyperbolic_piece', 'volume'))
    elif mode == 'normal':
        R, JJ = load('recognise'), load('jsj')
        items = []
        for n, r in T.items():
            names = R[n]['names'] or {x: 1 for x in JJ.get(n, {}).get('pachner', {}).get('names', [])}
            if names:
                items.append((n, r['pd'], r['D'], names))
        run(mode, normal, items, workers, show=('verdict', 'pplus', 'kinds'))
    elif mode == 'algebraic':
        R, H, JJ, PV, N = load('recognise'), load('hyp'), load('jsj'), load('piece'), load('normal')
        items = []
        for n, r in T.items():
            cls, _ = current_class(n, R, H, JJ, PV, N, {})
            if cls == 'undecided' or r['knotinfo_algebraic']:
                items.append((n, r['pd'], r['D'], r['h1'], r['jones'], r['sigma'], cls == 'undecided'))
        run(mode, algebraic, items, workers, show=('verdict', 'negdef_sign', 'same_normal_forms'))
    elif mode == 'identify':
        rec_all = load('records')['targets']
        R, H, JJ, PV, N, A = (load(m) for m in ('recognise', 'hyp', 'jsj', 'piece', 'normal', 'algebraic'))
        items = []
        for n, r in rec_all.items():      # the negative-definite trees of the follow-up records, identified again
            if r['followup'] is not None and r['followup']['cls'] == 'negdef':
                t = r['followup']['tree']
                items.append((n, r['pd'], r['D'], r['h1'], r['jones'], r['sigma'], [],
                              [('followup', t['weights'], t['edges'])]))
        for n, r in T.items():
            cls, how = current_class(n, R, H, JJ, PV, N, A)
            names = sorted(set(R.get(n, {}).get('names', {})) | set(N.get(n, {}).get('names', {}))
                           | set(JJ.get(n, {}).get('pachner', {}).get('names', [])))
            trees = []
            if n in N and N[n]['verdict'] in ('APPLICABLE', 'NOT APPLICABLE'):
                for t in N[n]['trees']:
                    if N[n]['verdict'] == 'NOT APPLICABLE' or t['pplus'] == 0:
                        trees.append(('regina', t['weights'], t['edges']))
            for d in A.get(n, {}).get('diagrams', []):
                for s, f in d['forms'].items():
                    if cls in ('nodef_alg', 'negdef') and (cls == 'nodef_alg' or f['pplus'] == 0):
                        trees.append((f'algebraic {d["source"]} {s}', f['weights'], f['edges']))
            if trees:
                items.append((n, r['pd'], r['D'], r['h1'], r['jones'], r['sigma'], names, trees))
        run(mode, identify, items, workers, show=('trees',))
    elif mode == 'dtest':
        recs = load('records')['targets']
        R, H, JJ, PV, N, A, I = (load(m) for m in ('recognise', 'hyp', 'jsj', 'piece', 'normal', 'algebraic',
                                                    'identify'))
        items = []
        for n, r in recs.items():
            names = sorted(set(R.get(n, {}).get('names', {})) | set(N.get(n, {}).get('names', {}))
                           | {t['name'] for t in I.get(n, {}).get('trees', []) if t['matched'] and t['name']})
            if r['followup'] is not None:
                if r['followup']['cls'] == 'negdef':
                    t = r['followup']['tree']
                    items.append((n, r['D'], t['weights'], t['edges'], names))
                continue
            cls, how = current_class(n, R, H, JJ, PV, N, A)
            if cls != 'negdef':
                continue
            if N.get(n, {}).get('verdict') == 'APPLICABLE':
                t = next(t for t in N[n]['trees'] if t['pplus'] == 0)
            else:
                d = A[n]['diagrams'][0]
                t = d['forms'][A[n]['negdef_sign'][0]]
            items.append((n, r['D'], t['weights'], t['edges'], names))
        run(mode, dtest, items, workers, show=('bad', 'h1', 'd_spin', 'passing', 'obstructed'))
    elif mode == 'controls':
        from database_knotinfo import link_list
        from basic_polyhedron import reduce_diagram
        BY = {r['name']: r for r in link_list()[1:]}
        rec = load('records')
        forms = rec['followup_forms']
        items = [('followup:' + n, int(BY[n]['determinant']), f['weights'], f['edges'], [])
                 for n, f in sorted(forms.items())]
        run('controls', dtest, items, workers, show=('bad', 'passing', 'obstructed'))
        items = [(n, r['pd'], r['D'], r['h1'], r['jones'], r['sigma']) for n, r in sorted(rec['controls'].items())]
        run('controls', control, items, workers, show=('verdict', 'obstructed', 'routes_agree'))
        items = []
        for f, kind in (('rank1_dinv_targets.json', 'alternating'), ('rank1_dinv_controls.json', 'alternating'),
                        ('rank1_dinv_mtargets.json', 'montesinos'), ('rank1_dinv_mcontrols.json', 'montesinos')):
            for n, r in sorted(data(f).items()):
                pd = json.loads(BY[n]['pd_notation'])
                if r.get('definite') and reduce_diagram(pd) == 1:
                    items.append(('known:' + n, pd, int(BY[n]['determinant']), kind, BY[n]['montesinos_notation'],
                                  r['obstructed']))
        run('controls', control_known, items, workers, show=('kind', 'dtests', 'agrees'))
    elif mode == 'gap':
        D = load('dtest')
        R1, LG = data('rank1_targets.json'), data('linking_genus.json')
        recs = load('records')['targets']
        items = []
        for k in sorted(n for n, r in D.items() if r.get('obstructed')):
            r = R1[k]
            info = dict(det=r['det'], test='plumbing', h_witt=LG[k]['h_witt'], survivors=r['survivors'], s=r['s'],
                        jordan=r['jordan'], rank_one_full=any(x[0] == r['det'] for x in r['survivors']))
            items.append((k, info, recs[k]['pd'], None, []))
        run(mode, gap, items, workers, show=('certificate', 'J_found', 'contradictions'))
    elif mode == 'summary':
        summary()
    else:
        raise SystemExit(__doc__)


def sha256(path):
    import hashlib
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def tree_checks(t):
    """The exact checks of an identified tree (step 6)."""
    return dict(identified=t.get('matched') is True, det=t['det_is_D'], coker=t['coker_is_H1'],
                casson_walker=bool(t['cw_signs']), nullity=t['nullity'] == 0)


def summary():
    import numpy, snappy, spherogram, regina, sympy, cypari, database_knotinfo
    from cypari import pari
    rec = load('records')
    T, prov = rec['targets'], rec['provenance']
    modes = ('recognise', 'hyp', 'jsj', 'piece', 'normal', 'algebraic', 'identify', 'dtest', 'controls', 'gap')
    S = {m: load(m) for m in modes}
    R, H, JJ, PV, N, A, I, D, C, G = (S[m] for m in modes)
    G4 = data('gamma4_knotinfo_2026.10.5.json')
    km = data('klein_bottle_search_mobius.json')['knots']
    mob = {k for k, v in km.items() if v.get('moebius_band_of_K')}
    targets, problems = {}, []
    for n, r in sorted(T.items()):
        out = dict(D=r['D'], sigma=r['sigma'], u=r['u'], gamma4=r['gamma4'], h1=r['h1'],
                   knotinfo_algebraic=r['knotinfo_algebraic'])
        if r['followup'] is not None:
            f = r['followup']
            out.update(new=False, cls=f['cls'], how=f'follow-up records ({f["followup_status"]})', followup=f)
            if n in I:                    # the negative-definite tree of the records, identified again
                out['identify'] = I[n]
                if not all(all(tree_checks(x).values()) for x in I[n]['trees']):
                    problems.append((n, 'tree of the follow-up records not identified again'))
        else:
            cls, how = current_class(n, R, H, JJ, PV, N, A)
            out.update(new=True, cls=cls, how=how, **{m: S[m][n] for m in ('recognise', 'hyp', 'jsj', 'piece', 'normal',
                                                                           'algebraic', 'identify') if n in S[m]})
        if n in D:
            d = D[n]
            out['dtest'] = d
            out['verdict'] = 'obstructed' if d.get('obstructed') else 'passes' if d.get('eligible') else 'not eligible'
        targets[n] = out
    new = {n: t for n, t in targets.items() if t['new']}
    old = {n: t for n, t in targets.items() if not t['new']}
    # checks of the new classes
    chk = Counter()
    for n, t in new.items():
        chk['H1 of Sigma(K) = coker(V + V^T)'] += t['recognise']['H1_agrees_with_seifert'] is True
        if t['cls'] == 'negdef':
            trees = [x for x in t.get('identify', {}).get('trees', []) if x['pplus'] == 0]
            ok = trees and all(all(tree_checks(x).values()) and x['pminus'] == x['vertices'] for x in trees)
            chk['negdef: every negative-definite tree identified, |det Q| = D, coker Q = H_1, '
                'Casson-Walker'] += bool(ok)
            if not ok:
                problems.append((n, 'negdef tree check failed'))
        if t['cls'] == 'nodef':
            trees = [x for x in t['identify']['trees'] if x['role'] == 'regina']
            ok = (len({x['canon'] for x in trees}) == 2 and all(x['pplus'] >= 1 and x['in_normal_form'] for x in trees)
                  and all(all(tree_checks(x).values()) for x in trees))
            chk['nodef: two normal forms, both p+ >= 1, identified, |det Q| = D, coker Q = H_1, Casson-Walker'] += ok
            if not ok:
                problems.append((n, 'nodef tree check failed'))
        if t['cls'] == 'nodef_alg':
            a = t['algebraic']
            forms = [f for d in a['diagrams'] for f in d['forms'].values()]
            ok = (a['same_normal_forms'] and len(a['diagrams']) >= 2 and all(f['pplus'] >= 1 and f['det_is_D']
                  and f['coker_is_H1'] and f['cw_signs'] and f['in_normal_form'] for f in forms))
            chk['nodef_alg: two algebraic diagrams, same normal forms, both p+ >= 1, |det Q| = D, coker Q = H_1, '
                'Casson-Walker'] += ok
            if not ok:
                problems.append((n, 'nodef_alg check failed'))
    # the KnotInfo diagram, where it is algebraic, against the Regina route
    cross = {}
    for n, t in new.items():
        a = t.get('algebraic')
        if (not a or not a['diagrams'] or a['diagrams'][0]['source'] != 'knotinfo'
                or t['cls'] not in ('negdef', 'nodef')):
            continue
        forms = a['diagrams'][0]['forms']
        alg_negdef = {f['canon'] for f in forms.values() if f['pplus'] == 0}
        reg = {x['canon'] for x in t['identify']['trees'] if x['role'] == 'regina'}
        if t['cls'] == 'negdef':
            same = a['verdict'] == 'APPLICABLE' and reg <= alg_negdef
        else:
            same = a['verdict'] == 'NOT APPLICABLE' and reg == {f['canon'] for f in forms.values()}
        cross[n] = same
        if not same:
            problems.append((n, 'KnotInfo algebraic diagram disagrees with the Regina route'))
    # the d-test
    elig = {n: t for n, t in targets.items() if t.get('verdict') in ('obstructed', 'passes')}
    indep = {n: t['dtest']['independent'] for n, t in elig.items()}
    for n, x in indep.items():
        if isinstance(x, dict) and not (x['agrees'] and x['same_d']):
            problems.append((n, 'independent d-test disagrees'))
    def hf_red_twice(tree, d):
        """L(Q) - sum_t d(dX, t) = 2 sum_t rank HF_red(dX, t) (the Casson-Walker formula of the follow-up paper, Section
        3.4), so it must be a nonnegative even integer."""
        x = Fraction(tree['L']) - sum(Fraction(v) for v in d)
        return int(x) if x >= 0 and x.denominator == 1 and x.numerator % 2 == 0 else None
    cw_d = {}
    for n, t in elig.items():
        tree = next(x for x in t['identify']['trees'] if x['pplus'] == 0)
        cw_d[n] = hf_red_twice(tree, t['dtest']['d'])
        if cw_d[n] is None:
            problems.append((n, 'L(Q) - sum d is not a nonnegative even integer'))
    star = {n: t['dtest']['star_check'] for n, t in elig.items() if 'star_check' in t['dtest']}
    for n, x in star.items():
        if 'obstructed' in x and not (x.get('agrees') and x.get('det_is_D')):
            problems.append((n, 'star plumbing route disagrees'))
    # the verdicts of an earlier, independent implementation: routeP_12.json for 12 of the targets, routeP_controls.json
    # (read below) for the forms of control (a)
    ROUTEP = os.path.join(ROOT, 'data', 'plumbing_dinv_routeP')
    p12 = json.load(open(os.path.join(ROUTEP, 'routeP_12.json')))
    p12_same = all(targets[n].get('dtest', {}).get('obstructed') == p['obstructed'] for n, p in p12.items())
    if not p12_same:
        problems.append(('routeP_12', 'verdicts of routeP_12.json not reproduced'))
    # controls
    status = lambda n: '1 (KnotInfo)' if G4.get(n) == '1' else '1 (census band)' if n in mob else (G4.get(n) or '?')
    ca = {k.split(':', 1)[1]: v for k, v in C.items() if k.startswith('followup:')}
    cb = {k: v for k, v in C.items() if ':' not in k}
    cc = {k.split(':', 1)[1]: v for k, v in C.items() if k.startswith('known:')}
    pc = json.load(open(os.path.join(ROUTEP, 'routeP_controls.json')))
    table_a = Counter((status(n), str(v.get('obstructed', 'not eligible'))) for n, v in ca.items())
    controls_a = dict(forms=len(ca), table={f'gamma_4 {a}, obstructed {b}': c for (a, b), c in sorted(table_a.items())},
                      gamma4_one=sum(status(n).startswith('1') for n in ca),
                      gamma4_one_obstructed=sorted(n for n, v in ca.items()
                                                   if status(n).startswith('1') and v.get('obstructed')),
                      same_as_pilot=all(pc[n].get('obstructed') == ca[n].get('obstructed') for n in pc) if ca else None,
                      independent_disagrees=sorted(n for n, v in ca.items() if isinstance(v.get('independent'), dict)
                                                   and not (v['independent']['agrees'] and v['independent']['same_d'])))
    trees_b = [(n, x) for n, v in cb.items() for x in v.get('trees', [])]
    failed = lambda v: 'error' in v or 'timeout' in v or 'crashed' in v
    controls_b = dict(knots=len(cb), errors=sorted(n for n, v in cb.items() if failed(v)),
                      named=sum(bool(v.get('names')) for v in cb.values()),
                      verdicts=dict(Counter(v.get('verdict', 'no Regina name') for v in cb.values())),
                      knotinfo_algebraic=sum('algebraic' in v for v in cb.values()),
                      negdef_trees=dict(Counter(x['role'] for _, x in trees_b)),
                      dtests=sum(x['dtest'].get('eligible') is True for _, x in trees_b),
                      not_eligible_bad=sum(x['dtest'].get('eligible') is False for _, x in trees_b),
                      obstructed=sorted({n for n, x in trees_b if x['dtest'].get('obstructed')}),
                      regina_trees_identified=sum(x.get('matched') is True for _, x in trees_b
                                                  if x['role'] == 'regina'),
                      regina_trees=sum(x['role'] == 'regina' for _, x in trees_b),
                      routes_disagree=sorted(n for n, v in cb.items() if v.get('routes_agree') is False),
                      independent_disagrees=sorted({n for n, x in trees_b
                                                    if isinstance(x['dtest'].get('independent'), dict)
                                                    and not (x['dtest']['independent']['agrees']
                                                             and x['dtest']['independent']['same_d'])}),
                      casson_walker_d_fails=sorted({n for n, x in trees_b if x['dtest'].get('eligible')
                                                    and hf_red_twice(x, x['dtest']['d']) is None}),
                      tree_checks_fail=sorted({n for n, x in trees_b if not all(v for k, v in tree_checks(
                          dict(x, matched=x.get('matched', True))).items())}))
    used_c = [(n, sg, f) for n, v in cc.items() for sg, f in v.get('forms', {}).items() if 'obstructed' in f]
    controls_c = dict(knots=len(cc), by_kind=dict(Counter(v.get('kind') for v in cc.values())),
                      errors=sorted(n for n, v in cc.items() if failed(v)),
                      knots_with_dtest=sum(v.get('dtests', 0) > 0 for v in cc.values()), dtests=len(used_c),
                      not_eligible=sum(f.get('eligible') is False for v in cc.values()
                                       for f in v.get('forms', {}).values()),
                      disagree=sorted(n for n, v in cc.items() if v.get('agrees') is False),
                      verdicts=dict(Counter(f'{v["kind"]}, recorded obstructed {v["recorded_obstructed"]}'
                                            for v in cc.values() if v.get('dtests'))),
                      bad_vertices=dict(Counter(f'{cc[n]["kind"]}, {f["bad"]} bad' for n, _, f in used_c)),
                      alternating_orientation=dict(Counter(
                          f'sign {sg}: ' + ('symmetric d' if f['same_d'] and f['same_d_reversed']
                                            else 'same d' if f['same_d']
                                            else 'reversed d' if f['same_d_reversed'] else 'NEITHER')
                          for n, sg, f in used_c if cc[n]['kind'] == 'alternating')),
                      independent=dict(Counter(str(f['independent']) for _, _, f in used_c)))
    if controls_a['gamma4_one_obstructed'] or controls_b['obstructed']:
        problems.append(('controls', 'a control is obstructed'))
    if controls_c['disagree'] or controls_c['errors']:
        problems.append(('controls', 'a knot of the rank1_dinv records disagrees'))
    # gap search
    gap_summary, gap_full = {}, {}
    for k, g in sorted(G.items()):
        if 'diagrams' not in g:
            gap_summary[k] = g; continue
        gap_summary[k] = {x: g[x] for x in ('det', 'h_witt', 'rank_one_full', 'arcs', 'bands_tried', 'knot_results',
                                             'det1_knots', 'J_found', 'certificate', 'representative', 'contradictions',
                                             'seconds') if x in g}
        gap_summary[k].update(diagrams=len(g['diagrams']), hits=len(g['hits']),
                              hits_failing_a_recheck=sum(not h['passes'] for h in g['hits']),
                              diagrams_not_isometric=sum(d['isometric_to_knotinfo'] is not True for d in g['diagrams']),
                              fit_failures=len(g['fit_failures']))
        gap_full[k] = {x: g[x] for x in ('diagrams', 'hits', 'det1_not_hits', 'fit_failures', 'fox_det_checks')}
    count = lambda ts: dict(Counter(t['cls'] for t in ts.values()))
    obstructed = sorted(n for n, t in targets.items() if t.get('verdict') == 'obstructed')
    upper = {}                            # gamma_4 <= 2 for every obstructed knot
    for n in obstructed:
        r = km[n]
        if targets[n]['gamma4'] == '[1,2]':
            upper[n] = 'KnotInfo range [1,2]'
        elif r['range'] == '[1,3]' and (r.get('klein_bottle') or {}).get('certified') and r.get('self_contained'):
            upper[n] = 'punctured Klein bottle of the census (klein_bottle_search_mobius.json)'
        else:
            problems.append((n, 'no upper bound 2'))
    for n, t in targets.items():
        if t.get('verdict') == 'obstructed':
            t['upper_bound_2'] = upper.get(n)
    summ = dict(
        targets=len(targets), followup_records=len(old), new=len(new),
        classes=dict(all=count(targets), followup_records=count(old), new=count(new)),
        new_how=dict(Counter(t['how'] for t in new.values())),
        negdef=dict(knots=sorted(n for n, t in targets.items() if t['cls'] == 'negdef'),
                    bad_vertices=dict(Counter(str(t['dtest']['bad']) for t in targets.values() if 'dtest' in t)),
                    not_eligible=sorted(n for n, t in targets.items() if t.get('verdict') == 'not eligible')),
        obstructed=obstructed, obstructed_new=sorted(n for n in obstructed if targets[n]['new']),
        obstructed_by_range=dict(Counter(targets[n]['gamma4'] for n in obstructed)),
        obstructed_upper_bound=dict(Counter(upper.values())),
        passes={n: t['dtest']['passing'] for n, t in targets.items() if t.get('verdict') == 'passes'},
        undecided={n: t['how'] for n, t in targets.items() if t['cls'] == 'undecided'},
        checks=dict(chk, knots_new=len(new), knotinfo_algebraic_agrees=f'{sum(cross.values())} of {len(cross)}',
                    independent_dtest=dict(Counter('H_1 not cyclic' if isinstance(x, str) else
                                                   'agrees' if x['agrees'] and x['same_d'] else 'DISAGREES'
                                                   for x in indep.values())),
                    star_plumbing_route=dict(Counter('agrees' if x.get('agrees') and x.get('det_is_D') else
                                                     x.get('result', 'DISAGREES') for x in star.values())),
                    pair_or_triple_without_star_route=sorted(n for n in elig if n not in star),
                    casson_walker_d=dict(Counter('L(Q) - sum d in 2Z, >= 0' if v is not None else 'FAILS'
                                                 for v in cw_d.values())),
                    l_space=sorted(n for n, v in cw_d.items() if v == 0),
                    followup_negdef_trees_identified_again=sum(
                        all(all(tree_checks(x).values()) for x in t['identify']['trees'])
                        for t in old.values() if 'identify' in t),
                    pilot_routeP_12_reproduced=p12_same),
        controls=dict(followup_forms=controls_a, new_pipeline=controls_b, known_fillings=controls_c),
        gap=dict(knots=len(G), obstructed_not_searched=sorted(set(obstructed) - set(G)),
                 J_tau_epsilon={k: dict(Counter(f"tau {h['hfk']['tau']}, epsilon {h['hfk']['epsilon']}"
                                                for h in g['hits'] if 'hfk' in h))
                                for k, g in G.items() if g.get('certificate')},
                 certified_only_by_J_with_tau_zero=sorted(k for k, g in G.items() if g.get('certificate') and all(
                     h['hfk']['tau'] == 0 for h in g['hits'] if 'hfk' in h)),
                 errors=sorted(k for k, g in G.items() if 'diagrams' not in g),
                 certified={k: g['J_found'] for k, g in G.items() if g.get('certificate')},
                 not_certified=sorted(k for k, g in G.items() if not g.get('certificate')),
                 contradictions=sum(g.get('contradictions', 0) for g in G.values())),
        problems=problems)
    sage_versions = sorted({json.dumps(r['versions'], sort_keys=True) for M in (H, PV) for r in M.values()
                            if 'versions' in r})
    meta = dict(versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                              regina=regina.versionString(),
                              database_knotinfo=database_knotinfo.version(), cypari=cypari.__version__,
                              pari='.'.join(map(str, pari.version())), sympy=sympy.__version__, numpy=numpy.__version__,
                              python=sys.version.split()[0], sage=[json.loads(x) for x in sage_versions]),
                parameters=dict(quick=QUICK, more=MORE, retry=RETRY, names=NAMES,
                                hyp=[HYP_TRIES, HYP_MORE, HYP_DRILLS, HYP_RANDOMISE],
                                jsj_triangulations=JSJ_TRIANGULATIONS, piece_randomise=PIECE_RANDOMISE, pachner=PACHNER,
                                node_range=NODE_RANGE, identify_tries=ID_TRIES, algebraic=[ALG_ITERS, list(ALG_SEEDS)],
                                timeout=TIMEOUT),
                followup_repository=prov,
                vendored=dict((os.path.relpath(os.path.join(dp, f), ROOT), sha256(os.path.join(dp, f)))
                              for dp, _, fs in os.walk(FU) for f in sorted(fs) if f.endswith('.py')),
                inputs={f: sha256(os.path.join(ROOT, 'data', f)) for f in (
                    'klein_bottle_search_mobius.json', 'open_gamma4_knotinfo_2026.10.5.json',
                    'gamma4_knotinfo_2026.10.5.json', 'rank1_targets.json', 'linking_genus.json')},
                timeouts={m: sorted(n for n, r in S[m].items() if 'timeout' in r or 'crashed' in r or 'error' in r)
                          for m in S if any('timeout' in r or 'crashed' in r or 'error' in r for r in S[m].values())})
    json.dump(dict(meta=meta, summary=summ, targets=targets,
                   controls=dict(followup_forms=ca, new_pipeline=cb, known_fillings=cc), gap=gap_summary),
              open(os.path.join(ROOT, 'data', 'plumbing_dinv.json'), 'w'), indent=0, default=str)
    json.dump(gap_full, open(os.path.join(ROOT, 'data', 'plumbing_dinv_gap.json'), 'w'), default=str)
    print(json.dumps({k: v for k, v in summ.items() if k not in ('passes',)}, indent=1, default=str))


if __name__ == '__main__':
    main()
