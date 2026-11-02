#!/usr/bin/env python
"""[C] Proposition 3.2(b), with the method of [C] Section 2.2: exact V_j(C) from the knot Floer complex of C, and the
test of [C] Proposition 2.6 with exact correction terms, for the knots with unknotting number one that the bounds of
scripts/surgery_box.py leave open; mode gap: the bands of [C] Theorem B for the two knots that this test excludes.

Setting.  For each of the 50 knots of scripts/surgery_dinv.py and each of its routes, Sigma_K = S^3_{eD/q}(kappa) with
q = 2 (q = the slope of a drilling route), e = +-1, and a stored diagram of kappa; C = kappa (e = 1) or its mirror
(e = -1), and Y' = S^3_{D/q}(C) = e Sigma_K.  surgery_box.py tests every V between the bounds of Hom-Wu and Rasmussen;
it excludes 13 knots and the plumbing trees of plumbing_dinv.py 4 more.  The other 33 (paper_counts.py: ndleft & su)
are the targets here; the 17 excluded knots, the calibration knots of surgery_dinv.py and 11n_176 (part pilot of
data/surgery_dinv.json) are controls.  For every route that identified C (as in surgery_box.routes) we determine V_j(C)
from the UV = 0 complex.
1. The complex.  knot_floer_homology.pd_to_hfk(pd, complex=True) (Szabo's program) on the diagram of C gives a reduced
   model of CFK(C) over F_2[U,V]/(UV): generators x with (A, M), M = gr_U, and arrows, each U^a y (horizontal:
   A(y) = A(x) + a, M(y) - 2a = M(x) - 1) or V^b y (vertical: A(y) = A(x) - b, M(y) = M(x) - 1); the program records
   the coefficient only (decode()).  Checks: as many generators as the rank of HFK, no arrow of length 0, d^2 = 0
   over F_2[U,V]/(UV), the bigraded ranks of the HFK of C stored by surgery_dinv.py, Euler characteristic Delta_C.
2. V_s from H_*(A^_s) (every complex).  A_s^+ = C{max(i, j - s) >= 0} has one tower, with bottom in grading -2V_s
   (V_s >= 0), and A^_s = C{max(i, j - s) = 0} is the kernel of U, so H(A^_s) -> H(A_s^+) -U-> H(A_s^+) is exact and
   the bottom of the tower comes from H_{-2V_s}(A^_s).  A^_s is the Alexander grading s part of the UV = 0 complex
   (U^m y with A(y) - m = s in grading M - 2m, V^n y with A(y) + n = s in grading M): a term U^a V^b with a, b >= 1
   leaves the hook.  So -2V_s is one of the even gradings <= 0 of H_*(A^_s), and V_s is determined when there is
   only one (hat_homology(), field hat).  An L-space knot (H_*(A^_s) of rank one for all s, so S^3_N(C) is an L-space
   for N >= 2g - 1) has one such grading for each s, and V_s = t_s, the torsion coefficient (checked).
3. V_s from CFK over F_2[U,V] (second argument).  Homotopy-equivalent reduced complexes over F_2[U,V]/(UV) are
   isomorphic, and a bigraded matrix that is invertible mod (U, V) is invertible, so CFK(C) has a reduced model over
   F_2[U,V] that is the complex of step 1 mod UV, in the same basis.  Its differential has, in addition, only terms
   U^a V^b y with a, b >= 1 that both gradings allow: gr_U(y) - 2a = gr_U(x) - 1 and gr_V(y) - 2b = gr_V(x) - 1,
   gr_V = M - 2A (slots()).  We take every subset of these slots for which d^2 = 0 over F_2[U,V] and compute V_s from
   the tower of H_*(A_s^+) (tower_V(), the computation of surgery_dinv.V_model with grading bounds from the complex);
   V must be the same for all subsets, and equal to step 2.  Subsets with H_*(C|_{V=1}) = F_2[U] in grading 0, as for
   S^3, are counted too (information only).  Kinds of complexes:
   - L-space knot: a staircase has no slot (gr_U decreases and gr_V increases along it), so it is CFK itself.
   - thin (M = A - tau on all generators): no slot, since an arrow raising both gradings would need A(y) - A(x) >= 1
     and A(x) - A(y) >= 1; so step 1 gives CFK (Petkova's staircase plus boxes is not used).
   - no diagonal slot: the same conclusion for a complex that is not thin.
   - diagonal slots: the gradings allow some terms U^a V^b; the subsets decide whether they can change V.
   If V is determined by neither step, the test of step 6 is run on every V allowed by both (and by V_j >= V_{j+1} >=
   V_j - 1); a verdict common to all of them is still a verdict.
4. The chirality of C.  surgery_dinv.py fixes it by an orientation-preserving isometry of exteriors (isometry) or by
   Regina's names, which give Sigma_K = +-S^3_{eD/q}(kappa) only up to orientation (unsigned test only).  Otherwise:
   - khovanov: if the reduced Khovanov homology of K over F_2 (KnotInfo) has rank D, Sigma_K is an L-space
     (Ozsvath-Szabo, branched double covers), so C is an L-space knot with D/q >= 2g - 1; exactly one chirality is.
   - covers: exterior_to_link (careful_perturbation) guarantees that the diagram L has E(L) = +-E(kappa) with the
     meridian, so kappa = L or its mirror and Sigma_K = S^3_{eD/q}(L) or -S^3_{-eD/q}(L).  The connected covers of
     degree n (conjugacy classes of index-n subgroups of pi_1, SnapPy) with their H_1 are invariants of the unoriented
     manifold; at the first n <= COVER_DEGREES where those of S^3_{eD/q}(L) and S^3_{-eD/q}(L) differ, those of
     Sigma_K (from the KnotInfo PD code) agree with exactly one of them, which fixes kappa.  For the route that
     decides each knot, GAP (LowIndexSubgroupsFpGroup, AbelianInvariants, under Sage) recomputes the covers of that
     degree.
   - otherwise both chiralities are tested, and a verdict needs both to agree.
   The other chirality has the dual complex (gradings negated, arrows reversed); for a route without a fixed chirality
   it is also computed by the program on the mirror diagram, and steps 2-3 must give the same data.
5. Correction terms: surgery_dinv.surgery_form (Ni-Wu, Prop. 1.6 and Rem. 2.10, with label i_s + qk for t_s + k mu,
   2 i_s = q - 1 mod p, lambda(mu, mu) = -q/p), Z = -e Y' as in surgery_dinv.py.  Checks: Casson-Walker against
   Mullins' R(K), e (sum_i d(L(D, q), i) - 2q sum_k t_k) = R; lambda_Z = (V + V^T)^-1 (KnotInfo Seifert matrix); for an
   L-space Y' also sum_t d(Y', t) = e R; V_j <= g - j (Rasmussen, g_4 <= g), V_j >= 1 for j < max(tau, nu) (Hom-Wu),
   V_j >= V_{j+1} >= V_j - 1; when nu >= 0 (nu from the program), rank HF^(Y') = D + 2 max(0, (2 nu - 1) q - D)
   + q sum_s (rank H_*(A^_s) - 1) (Ozsvath-Szabo, rational surgeries, Prop. 9.5, for p/q > 0) is at most the rank of
   the reduced Khovanov homology of K; the d-functions of the routes of a knot are isometric (up to sign when a route
   has only Regina's names).  The field epsilon records whether V = 0 when epsilon(C) = 0; it is not asserted, as no
   such theorem is used.
6. The test of Proposition 2.6 (rank1_dinv.test(Z, (1, -1)), repeated by montesinos_dinv_check.obstructed); the
   signed variant of rank1_dinv (s = sigma + 4 Arf = 2 or 6) is recorded when the orientation is known.  A knot is
   excluded (gamma_4 = 2, as gamma_4 <= 2 by its KnotInfo range [1, 2]) when the test fails on its exact d-function.
   For an excluded knot, surgery_box.box is also run with V_min from the chirality fixed here (box_with_chirality).
Controls.  calibration (surgery_dinv.CALIBRATION, alternating, Sigma_K an L-space): the d-function must be isometric
to the Goeritz d-function of rank1_dinv, with the same verdict, and no knot with KnotInfo value 1 is excluded.  others
(the 17 excluded knots): each is excluded again, V lies in the box of surgery_box.py, and for the six with a star
plumbing (surgery_dinv.py, field plumbing) the d-function is isometric to that of the plumbing up to sign.  pilot:
11n_176 must pass.  tensor: connected sums A # B of two torus knots, whose CFK is the tensor product of two staircases
(the program's complexes of A and B): steps 2-3 on the program's complex of a diagram of A # B must give the V of the
tensor product, and the H_*(A^_s) of the two complexes must agree.
Mode gap: for the excluded targets, the band search of smooth_top_gap.py (surgery_dinv.gap_records): one band from K to
a knot J with Delta_J = 1 gives a locally flat Mobius band (Freedman), so gamma_4^top = 1 < gamma_4 = 2; every hit
is rechecked as in smooth_top_gap.recheck, and J must not be known to be smoothly slice.
Usage: python scripts/surgery_cfk.py [workers] | gap [workers]
Output: data/surgery_cfk.json (checkpoint data/surgery_cfk.partial); mode gap: data/surgery_cfk_gap.json (checkpoint
data/surgery_cfk_gap.partial)
"""
import hashlib, json, os, shutil, subprocess, sys, time, warnings
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from rank1_dinv import Form, test
import surgery_dinv as sd

CFK_SECONDS = 3600                        # time limit of one run of Szabo's program
MAX_SLOTS = 16                            # the subsets of step 3 are enumerated up to this many slots
COVER_DEGREES = 7                         # covers of degree 2..COVER_DEGREES in step 4
TENSOR = [('3_1', 1, '3_1', 1), ('5_1', 1, '3_1', -1), ('8_19', 1, '3_1', 1), ('8_19', 1, '3_1', -1),
          ('8_19', 1, '5_1', -1), ('5_1', 1, '8_19', -1), ('10_124', 1, '8_19', -1), ('7_1', 1, '8_19', -1),
          ('9_1', 1, '10_124', -1)]       # (A, sign, B, sign): A # B or A # mirror(B), KnotInfo diagrams


# ---------- step 1: the complex

def pd_key(pd):
    """The key of a diagram in the jobs of Szabo's program and in the field complexes of the output."""
    return hashlib.sha1(json.dumps([[int(v) for v in c] for c in pd]).encode()).hexdigest()[:16]


def mirror_pd(pd):
    """The PD code of the mirror image."""
    import spherogram
    return [list(c) for c in spherogram.Link([list(c) for c in pd]).mirror().PD_code()]


def _cfk_child(pd, q):
    try:
        import knot_floer_homology as kfh
        h = kfh.pd_to_hfk(pd, complex=True)
        q.put(('ok', dict(generators=[[int(x), int(A), int(M)] for x, (A, M) in sorted(h['generators'].items())],
                          differentials=[[int(x), int(y), int(c)] for (x, y), c in sorted(h['differentials'].items())],
                          ranks=sorted([int(A), int(M), int(r)] for (A, M), r in h['ranks'].items()),
                          **{k: h[k] for k in ('tau', 'nu', 'epsilon', 'seifert_genus', 'fibered', 'L_space_knot',
                                               'total_rank', 'modulus')})))
    except BaseException as ex:
        q.put(('error', repr(ex)))


def run_jobs(jobs, workers, partial, limit=CFK_SECONDS, log=print):
    """Szabo's program on every diagram of jobs ({key: pd}), each in its own process with a time limit; results are
    checkpointed to `partial` (JSON) as they come."""
    import multiprocessing as mp
    import queue as Q
    ctx = mp.get_context('spawn')
    done = json.load(open(partial)) if os.path.exists(partial) else {}
    todo = sorted((k for k in jobs if k not in done), key=lambda k: len(jobs[k]))
    running, t_start = {}, time.time()
    while todo or running:
        while todo and len(running) < workers:
            k = todo.pop(0)
            q = ctx.Queue()
            p = ctx.Process(target=_cfk_child, args=(jobs[k], q))
            p.start()
            running[k] = (p, q, time.time())
        time.sleep(0.2)
        for k in list(running):
            p, q, t0 = running[k]
            try:
                res = q.get(timeout=0.05)
            except Q.Empty:
                if p.is_alive() and time.time() - t0 <= limit:
                    continue
                try:
                    res = q.get(timeout=2)
                except Q.Empty:
                    res = ('timeout', limit) if p.is_alive() else ('error', 'exit code %s' % p.exitcode)
                p.kill()
            p.join(5)
            del running[k]
            status, val = res
            done[k] = dict(status=status, seconds=round(time.time() - t0, 1), **(val if status == 'ok' else dict(error=val)))
            json.dump(done, open(partial + '.tmp', 'w'))
            os.replace(partial + '.tmp', partial)
            log('%s %s %.0fs (%d left, %d running, %.0fs)' % (k, status, time.time() - t0, len(todo), len(running),
                                                             time.time() - t_start))
    return done


def decode(cx):
    """(gens, arrows) of a record of Szabo's program: gens {x: (A, M)}, arrows (x, y, a, b) for the term U^a V^b y of
    d x, with a, b read off the gradings (step 1); an arrow of length 0 or of neither form stops the script."""
    gens = {x: (A, M) for x, A, M in cx['generators']}
    arrows = []
    for x, y, c in cx['differentials']:
        assert c % 2, 'coefficient %d' % c
        (Ax, Mx), (Ay, My) = gens[x], gens[y]
        if My == Mx - 1 and Ax > Ay:
            arrows.append((x, y, 0, Ax - Ay))
        elif Ay > Ax and My - 2 * (Ay - Ax) == Mx - 1:
            arrows.append((x, y, Ay - Ax, 0))
        else:
            raise AssertionError('arrow %d -> %d is neither horizontal nor vertical' % (x, y))
    return gens, arrows


def dual(gens, arrows):
    """The complex of the mirror image: gradings negated, arrows reversed with the same powers."""
    return {x: (-A, -M) for x, (A, M) in gens.items()}, [(y, x, a, b) for x, y, a, b in arrows]


def tensor(c1, c2):
    """The tensor product of two complexes (CFK of a connected sum)."""
    (g1, a1), (g2, a2) = c1, c2
    idx = {(x, y): i for i, (x, y) in enumerate((x, y) for x in g1 for y in g2)}
    gens = {i: (g1[x][0] + g2[y][0], g1[x][1] + g2[y][1]) for (x, y), i in idx.items()}
    arrows = [(idx[x, z], idx[y, z], a, b) for x, y, a, b in a1 for z in g2]
    arrows += [(idx[z, x], idx[z, y], a, b) for x, y, a, b in a2 for z in g1]
    return gens, arrows


def alexander(gens):
    """Delta_C as {exponent: coefficient}, the Euler characteristic of HFK."""
    out = Counter()
    for A, M in gens.values():
        out[A] += 1 if M % 2 == 0 else -1
    return {e: c for e, c in out.items() if c}


def rank_f2(rows):
    """Rank over F_2 of vectors given as sets of basis elements."""
    index, piv, r = {}, {}, 0
    for row in rows:
        v = 0
        for e in row:
            v ^= 1 << index.setdefault(e, len(index))
        while v:
            h = v.bit_length() - 1
            if h not in piv:
                piv[h] = v; r += 1
                break
            v ^= piv[h]
    return r


def d_squared(gens, arrows, mod_uv=False):
    """Pairs (x, z) on which d^2 has an odd coefficient, over F_2[U,V] (or F_2[U,V]/(UV)); the monomial of a pair is
    fixed by the gradings, so the coefficient counts the paths x -> y -> z."""
    nxt = defaultdict(list)
    for x, y, a, b in arrows:
        nxt[x].append((y, a, b))
    bad = []
    for x in gens:
        cnt = Counter()
        for y, a, b in nxt[x]:
            for z, a2, b2 in nxt[y]:
                if not (mod_uv and a + a2 and b + b2):
                    cnt[z] += 1
        bad += [(x, z) for z, c in cnt.items() if c % 2]
    return bad


# ---------- steps 2 and 3: V_s

def hat_homology(gens, arrows, s):
    """{grading: dimension} of H_*(A^_s) (step 2); terms carrying both U and V vanish."""
    nxt = defaultdict(list)
    for x, y, a, b in arrows:
        nxt[x].append((y, a, b))
    basis = {}
    for y, (A, M) in gens.items():
        if A >= s:
            basis[(y, A - s, 0)] = M - 2 * (A - s)
        if A <= s:
            basis[(y, 0, s - A)] = M
    def bd(e):
        y, m, n = e
        out = set()
        for z, a, b in nxt[y]:
            if not (m + a and n + b):
                assert (z, m + a, n + b) in basis
                out ^= {(z, m + a, n + b)}
        return out
    by = defaultdict(list)
    for e, gr in basis.items():
        by[gr].append(e)
    rk = {gr: rank_f2([bd(e) for e in es]) for gr, es in by.items()}
    return {gr: h for gr, es in by.items() for h in [len(es) - rk[gr] - rk.get(gr + 1, 0)] if h}


def tower_V(gens, arrows, s):
    """V_s of the complex over F_2[U,V] with these arrows (of any kind): -2V_s is the bottom grading of the tower of
    H_*(A_s^+), where [x, i] (j = i + A(x)) has grading M(x) + 2i.  The tower class is taken in an even grading
    > max M, where A_s^+ agrees with CFK^infinity, and multiplied by U until it is a boundary."""
    nxt = defaultdict(list)
    for x, y, a, b in arrows:
        nxt[x].append((y, a))
    inside = lambda x, i: max(i, i + gens[x][0] - s) >= 0
    elems = lambda gr: [(x, (gr - M) // 2) for x, (A, M) in gens.items()
                        if (gr - M) % 2 == 0 and inside(x, (gr - M) // 2)]
    def bd(e):
        out = set()
        for y, a in nxt[e[0]]:
            if inside(y, e[1] - a):
                out ^= {(y, e[1] - a)}
        return out
    def boundary(z, gr):
        B = [bd(e) for e in elems(gr + 1)]
        return rank_f2(B + [z]) == rank_f2(B)
    gr = max(M for A, M in gens.values()) + 1
    gr += gr % 2
    z, piv = None, []                     # a cycle of grading gr that is not a boundary (H_gr(A_s^+) = F_2)
    for e in elems(gr):
        src, im = {e}, bd(e)
        for p, ps, pi in piv:
            if p in im:
                im ^= pi; src ^= ps
        if im:
            piv.append((next(iter(im)), src, im))
        elif not boundary(src, gr):
            z = src
            break
    assert z is not None, 'no tower class'
    g = max(abs(A) for A, M in gens.values())
    last = None
    while not boundary(z, gr):
        last = gr
        z = {(x, i - 1) for x, i in z if inside(x, i - 1)}
        gr -= 2
        assert gr >= -2 * g - 4, 'the tower does not end'
    return -last // 2


def slots(gens):
    """The diagonal slots (x, y, a, b), a, b >= 1: the terms U^a V^b y of d x that both gradings allow (step 3)."""
    out = []
    for x, (Ax, Mx) in gens.items():
        for y, (Ay, My) in gens.items():
            if (My - Mx + 1) % 2 == 0:
                a = (My - Mx + 1) // 2
                b = Ax - Ay + a
                if a >= 1 and b >= 1:
                    out.append((x, y, a, b))
    return out


def s3_homology(gens, arrows):
    """True if H_*(C|_{V=1}) is F_2[U] with top grading 0 (that of S^3); below min M - 2 the chain groups are
    U-periodic, so the gradings down to min M - 4 suffice."""
    nxt = defaultdict(list)
    for x, y, a, b in arrows:
        nxt[x].append((y, a))
    Ms = [M for A, M in gens.values()]
    elems = lambda k: [(x, (M - k) // 2) for x, (A, M) in gens.items() if (M - k) % 2 == 0 and M >= k]
    bd = lambda e: {(y, e[1] + a) for y, a in nxt[e[0]]}
    for k in range(max(Ms) + 1, min(Ms) - 5, -1):
        h = len(elems(k)) - rank_f2([bd(e) for e in elems(k)]) - rank_f2([bd(e) for e in elems(k + 1)])
        if h != (1 if k <= 0 and k % 2 == 0 else 0):
            return False
    return True


def lifts(gens, arrows, g, limit=MAX_SLOTS):
    """Step 3: the subsets of the slots with d^2 = 0 over F_2[U,V], with the V of each and the S^3 condition."""
    sl = slots(gens)
    if len(sl) > limit:
        return dict(slots=[list(x) for x in sl], enumerated=False)
    Vs, s3 = Counter(), Counter()
    for m in range(2 ** len(sl)):
        arr = arrows + [sl[i] for i in range(len(sl)) if m >> i & 1]
        if d_squared(gens, arr):
            continue
        V = tuple(tower_V(gens, arr, s) for s in range(g + 1))
        Vs[V] += 1
        s3[V] += s3_homology(gens, arr)
    return dict(slots=[list(x) for x in sl], enumerated=True, subsets=2 ** len(sl), d2_zero=sum(Vs.values()),
                V=[dict(V=list(v), subsets=c, with_S3_homology=s3[v]) for v, c in sorted(Vs.items())])


def analyse(gens, arrows, tau, g):
    """Steps 2 and 3 for one complex: its kind, V (None if neither step determines it) and both arguments."""
    hat = {s: hat_homology(gens, arrows, s) for s in range(-g - 1, g + 2)}
    lspace = all(sum(H.values()) == 1 for H in hat.values())
    thin = all(M == A - tau for A, M in gens.values())
    lf = lifts(gens, arrows, g)
    assert not lf['enumerated'] or lf['d2_zero'], 'no subset of the slots gives d^2 = 0'
    kind = 'L-space knot' if lspace else 'thin' if thin else 'no diagonal slot' if not lf['slots'] else 'diagonal slots'
    assert kind not in ('L-space knot', 'thin') or not lf['slots'], 'a slot in a staircase or a thin complex'
    even = {s: sorted(k for k in hat[s] if k % 2 == 0 and k <= 0) for s in range(g + 1)}
    V_hat = [-even[s][0] // 2 if len(even[s]) == 1 else None for s in range(g + 1)]
    V_lifts = lf['V'][0]['V'] if lf['enumerated'] and len(lf['V']) == 1 else None
    out = dict(kind=kind, hat=[dict(s=s, rank=sum(hat[s].values()), even_nonpositive=even[s]) for s in range(g + 1)],
               hat_rank_excess=sum(sum(H.values()) - 1 for H in hat.values()), V_hat=V_hat, lifts=lf, V_lifts=V_lifts)
    if None not in V_hat:
        assert V_lifts in (None, V_hat), 'steps 2 and 3 disagree'
        out.update(V=V_hat, argument='hat and lifts' if V_lifts else 'hat')
    else:
        out.update(V=V_lifts, argument='lifts' if V_lifts else None)
    if lspace:
        out['torsion'] = sd.torsion(alexander(gens))
        assert out['V'] == out['torsion'], 'V of an L-space knot is not t_s'
    return out


def candidates(an, g):
    """Every V allowed by steps 2 and 3 (and V_j >= V_{j+1} >= V_j - 1, V_g = 0) when V is not determined."""
    from itertools import product
    cand = [[v] if v is not None else [-k // 2 for k in h['even_nonpositive']] for v, h in zip(an['V_hat'], an['hat'])]
    allowed = {tuple(x['V']) for x in an['lifts']['V']} if an['lifts']['enumerated'] else None
    return [list(V) for V in product(*cand) if V[-1] == 0 and all(V[j] >= V[j + 1] >= V[j] - 1 for j in range(g))
            and (allowed is None or V in allowed)]


# ---------- step 4: chirality

def khovanov_rank(BY, name):
    """rank of the reduced Khovanov homology of K over F_2 (KnotInfo), as in surgery_dinv.khovanov_check."""
    import sympy
    poly = BY[name].get('khovanov_reduced_mod2_polynomial') or ''
    if not poly:
        return None
    t, q = sympy.symbols('t q')
    return int(sympy.sympify(poly.replace('^', '**'), locals={'t': t, 'q': q}).subs({t: 1, q: 1}))


def abelian_invariants(divisors):
    """The invariants (prime powers, and 0 for Z) of the abelian group with these elementary divisors."""
    import sympy
    out = []
    for d in divisors:
        out += [0] if d == 0 else [p ** k for p, k in sympy.factorint(d).items()]
    return sorted(out)


def cover_profile(M, n):
    """The connected covers of degree n (conjugacy classes of index-n subgroups of pi_1), by their H_1."""
    return sorted(abelian_invariants(C.homology().elementary_divisors()) for C in M.covers(n))


def presentation(M):
    """The generators and relators of SnapPy's presentation of pi_1(M), for gap_check."""
    G = M.fundamental_group()
    return dict(generators=G.generators(), relators=G.relators())


def chirality_by_covers(SigK, pd_L, e, D, q, degrees=COVER_DEGREES):
    """Step 4 (covers): 'construction' (kappa = L), 'mirror' (kappa = mirror of L) or None."""
    import snappy, spherogram
    Y = {}
    for sgn, nm in ((1, 'plus'), (-1, 'minus')):
        E = spherogram.Link([list(c) for c in pd_L]).exterior()
        E.dehn_fill((sgn * e * D, q))
        Y[nm] = snappy.Manifold(E.filled_triangulation())
    rec = dict(slopes=dict(plus=[e * D, q], minus=[-e * D, q]), counts={})
    for n in range(2, degrees + 1):
        P = dict(sigma=cover_profile(SigK, n), plus=cover_profile(Y['plus'], n), minus=cover_profile(Y['minus'], n))
        rec['counts'][str(n)] = {k: len(v) for k, v in P.items()}
        if P['plus'] != P['minus']:
            dec = 'construction' if P['sigma'] == P['plus'] else 'mirror' if P['sigma'] == P['minus'] else None
            assert dec, 'the covers of Sigma_K agree with neither surgery'
            rec.update(degree=n, decided=dec, profiles=P,
                       presentations={k: presentation(M) for k, M in (('sigma', SigK), ('plus', Y['plus']),
                                                                       ('minus', Y['minus']))})
            return rec
    rec['decided'] = None
    return rec


GAP_CHILD = r'''
import json, sys
from sage.all import libgap
out = {}
for job in json.load(sys.stdin):
    gens = job['generators']
    F = libgap.FreeGroup(len(gens))
    fg = F.GeneratorsOfGroup()
    def word(w):
        x = F.One()
        for ch in w:
            i = gens.index(ch.lower())
            x = x * (fg[i] if ch.islower() else fg[i] ** -1)
        return x
    H = F / [word(r) for r in job['relators']]
    L = [K for K in libgap.LowIndexSubgroupsFpGroup(H, job['degree']) if int(libgap.Index(H, K)) == job['degree']]
    out[job['id']] = sorted(sorted(int(a) for a in K.AbelianInvariants()) for K in L)
print(json.dumps(out))
'''


def gap_check(decisions):
    """The covers of the deciding degree recomputed by GAP under Sage: {key: agrees}, or None without Sage."""
    sage = shutil.which('sage')
    if not sage or not decisions:
        return None
    jobs = [dict(id='%s|%s' % (k, nm), degree=d['degree'], **d['presentations'][nm])
            for k, d in decisions.items() for nm in ('sigma', 'plus', 'minus')]
    r = subprocess.run([sage, '-python', '-c', GAP_CHILD], input=json.dumps(jobs), capture_output=True, text=True,
                       timeout=7200)
    res = json.loads(r.stdout.strip().splitlines()[-1])
    return {k: all(res['%s|%s' % (k, nm)] == d['profiles'][nm] for nm in ('sigma', 'plus', 'minus'))
            for k, d in decisions.items()}


# ---------- steps 5 and 6

def dtest(D, q, e, V, s, oriented):
    """The test on Z = -e Y', Y' = S^3_{D/q}(C); the signed variant only when the orientation is known."""
    import montesinos_dinv_check as MC
    FY = sd.surgery_form(D, q, V)
    FZ = sd.negate(FY) if e == 1 else FY
    passing = test(FZ, (1, -1))
    ob2, _ = MC.obstructed(FZ.D, FZ.adj[0][0], {k: FZ.d[(k,)] for k in range(FZ.D)})
    assert ob2 == (passing is None), 'the two implementations of the test differ'
    rec = dict(obstructed=passing is None, passing=passing, obstructed_independent_code=ob2)
    if oriented:
        allowed = {2: (1,), 6: (-1,)}.get(s, (1, -1))
        signed = passing if passing and passing['sign'] in allowed else test(FZ, allowed)
        rec.update(passing_signed=signed, obstructed_signed=signed is None)
    return FY, FZ, rec


def checks(D, q, e, c, info, FY):
    """Step 5 checks for one chirality c of one route."""
    V, g, an = c['analysis']['V'], c['genus'], c['analysis']
    t = sd.torsion(alexander(c['gens']))
    dl = sum(sd.d_lens(D, q, i) for i in range(D))
    out = dict(casson_walker=e * (dl - 2 * q * (t[0] + 2 * sum(t[1:]))) == info['R'],
               linking_form=sd.lens_form_check(info['seifert'], D, e * q)['plus'],
               monotone=V[-1] == 0 and all(V[j] >= V[j + 1] >= V[j] - 1 for j in range(len(V) - 1)),
               rasmussen=all(v <= max(g - j, 0) for j, v in enumerate(V)),
               hom_wu=all(V[j] >= 1 for j in range(max(c['tau'], c['nu'] if c['nu'] is not None else 0, 0))),
               epsilon=c['epsilon'] != 0 or not any(V))   # V = 0 when epsilon(C) = 0: recorded, not a theorem used here
    if an['kind'] == 'L-space knot' and D >= q * (2 * g - 1):
        out['L_space_sum_d'] = sum(FY.d.values()) == e * info['R']
    return out


def run_test(D, q, e, s, c, info, oriented):
    """Steps 5 and 6 for one chirality: (record, Form of Z or None)."""
    an = c['analysis']
    if an['V'] is None:
        vs = candidates(an, c['genus'])
        ob = {test(sd.negate(sd.surgery_form(D, q, V)) if e == 1 else sd.surgery_form(D, q, V), (1, -1)) is None for V in vs}
        return dict(V=None, candidates=vs, obstructed=ob.pop() if len(ob) == 1 else None), None
    FY, FZ, rec = dtest(D, q, e, an['V'], s, oriented)
    rec['checks'] = checks(D, q, e, c, info, FY)
    assert all(v for k, v in rec['checks'].items() if k != 'epsilon'), rec['checks']   # epsilon is only recorded
    rec.update(V=an['V'], d_Z={str(k[0]): str(v) for k, v in FZ.d.items()}, lambda_Z_mu_mu='%d/%d' % (FZ.adj[0][0], D))
    return rec, FZ


# ---------- one route, one knot

ORIENTED = ('isometry', 'khovanov', 'covers')


def routes(rec):
    """The routes of surgery_dinv.py that identified C, without duplicates (surgery_box.routes), with labels."""
    out = [('diagram %d crossing %d' % (r['diagram_index'], r['crossing']) if 'diagram_index' in r else
            'crossing %d' % r['crossing'], r)
           for r in rec.get('routes', []) if r.get('kappa_hat') and not r.get('duplicate_of')]
    dr = [r for r in rec.get('drill_routes', []) if r.get('kappa_hat') and not r.get('duplicate_of')]
    return out + [('drilling %d' % k, r) for k, r in enumerate(dr)]


def route_item(lab, r):
    """Slope, sign, chirality status and the diagrams of C (and of its mirror) for one route."""
    q = r['slope_p_q'][1] if r.get('kind') == 'drilling' else 2
    ori = r['kappa']['diagram'].get('orientation')
    pd0 = r['kappa']['hfk_diagram']['pd']             # kappa, or its mirror for 'mirror (Regina names)'
    flip = (r['e'] == -1) != (ori == 'mirror (Regina names)')
    status = 'isometry' if ori in ('diagram', 'mirror', 'amphichiral') else 'regina names' if ori else None
    it = dict(label=lab, q=q, e=r['e'], orientation=ori, status=status, pd_C=mirror_pd(pd0) if flip else pd0,
              expected_ranks=r['kappa_hat']['hfk']['ranks'], kappa_identify=r['kappa'].get('identify'))
    if status is None:                                # the diagram of exterior_to_link: kappa = L or its mirror
        it.update(pd_L=pd0, pd_mirror=pd0 if flip else mirror_pd(pd0))
    return it


def complex_record(cx):
    """gens, arrows, tau, nu, epsilon, genus of a record of Szabo's program, with the checks of step 1."""
    gens, arrows = decode(cx)
    al = alexander(gens)
    ok = dict(reduced=len(gens) == cx['total_rank'], d2_mod_uv=not d_squared(gens, arrows, True),
              ranks=sorted([A, M, r] for (A, M), r in Counter(gens.values()).items()) == sorted(cx['ranks']),
              alexander=sum(al.values()) == 1 and all(al.get(-k) == c for k, c in al.items()))
    assert all(ok.values()), ok
    g = cx['seifert_genus']
    return dict(gens=gens, arrows=arrows, tau=cx['tau'], nu=cx['nu'], epsilon=cx['epsilon'], genus=g, checks=ok,
                analysis=analyse(gens, arrows, cx['tau'], g))


def same_hat(c1, c2):
    """The two complexes have the same H_*(A^_s) for |s| <= g + 1."""
    g = c1['genus']
    return all(hat_homology(c1['gens'], c1['arrows'], s) == hat_homology(c2['gens'], c2['arrows'], s)
               for s in range(-g - 1, g + 2))


def one(name, rec_sd, info, cxs, covers):
    """Steps 2-6 for every route of one knot; covers caches the decisions of step 4 by covers."""
    D, s = rec_sd['det'], rec_sd['s']
    out = dict(name=name, det=D, s=s, R=str(info['R']), khovanov_rank=info['khovanov'],
               sigma_L_space=info['khovanov'] == D, routes=[])
    forms = []
    for lab, r in routes(rec_sd):
        it = route_item(lab, r)
        q, e = it['q'], it['e']
        rr = dict(route=lab, slope=[e * D, q], e=e, orientation=it['orientation'], kappa_identify=it['kappa_identify'],
                  key=pd_key(it['pd_C']), pd=it['pd_C'])
        out['routes'].append(rr)
        cx = cxs.get(rr['key'])
        if not cx or cx['status'] != 'ok':
            rr['flag'] = 'complex not computed: %s' % (cx or {}).get('status', 'no job')
            continue
        C = complex_record(cx)
        assert sorted(cx['ranks']) == sorted(map(list, it['expected_ranks'])), 'HFK differs from surgery_dinv.py'
        rr['complex_checks'] = C['checks']
        cand = {'C': C}
        status, which = it['status'], 'C'
        if status is None:                            # step 4: the other chirality, by duality and by the program
            gm, am = dual(C['gens'], C['arrows'])
            m = dict(gens=gm, arrows=am, tau=-C['tau'], nu=None, epsilon=-C['epsilon'], genus=C['genus'],
                     analysis=analyse(gm, am, -C['tau'], C['genus']))
            cm = cxs.get(pd_key(it['pd_mirror']))
            if cm and cm['status'] == 'ok':
                P = complex_record(cm)
                rr['mirror_by_program'] = dict(key=pd_key(it['pd_mirror']), same_hat=same_hat(m, P),
                                               same_V=P['analysis']['V'] == m['analysis']['V'], nu=P['nu'])
                assert rr['mirror_by_program']['same_hat'] and rr['mirror_by_program']['same_V'], 'dual and program differ'
                m['nu'] = P['nu']
            cand['mirror'] = m
            if out['sigma_L_space']:
                ok = {k: c['analysis']['kind'] == 'L-space knot' and D >= q * (2 * c['genus'] - 1) for k, c in cand.items()}
                assert sum(ok.values()) == 1, 'Sigma_K is an L-space; L-space candidates %s' % ok
                status, which = 'khovanov', 'C' if ok['C'] else 'mirror'
            else:
                ck = '%s|%s|%d|%d' % (name, pd_key(it['pd_L']), e, q)
                if ck not in covers:
                    covers[ck] = chirality_by_covers(info['sigma_cover'](), it['pd_L'], e, D, q)
                cv = covers[ck]
                rr['covers'] = dict(key=ck, **{k: v for k, v in cv.items() if k != 'presentations'})
                which = {'construction': 'C', 'mirror': 'mirror', None: None}[cv['decided']]
                status = 'covers' if which else None
        rr.update(chirality=status, C_is='kappa_hat' if which == 'C' else 'mirror of kappa_hat' if which else None)
        if out['sigma_L_space'] and which:
            c = cand[which]
            assert c['analysis']['kind'] == 'L-space knot' and D >= q * (2 * c['genus'] - 1), 'Sigma_K is an L-space'
        rr['candidates'] = {}
        for k, c in cand.items():
            t, FZ = run_test(D, q, e, s, c, info, status in ORIENTED and k == which)
            hf = (D + 2 * max(0, (2 * c['nu'] - 1) * q - D) + q * c['analysis']['hat_rank_excess']
                  if c['nu'] is not None and c['nu'] >= 0 else None)          # rank HF^(Y'), Prop. 9.5 (nu >= 0)
            rr['candidates'][k] = dict(kind=c['analysis']['kind'], tau=c['tau'], nu=c['nu'], epsilon=c['epsilon'],
                                       genus=c['genus'], generators=len(c['gens']), arrows=len(c['arrows']),
                                       hf_rank=hf, analysis=c['analysis'], test=t)
            if k == which and hf is not None and info['khovanov'] is not None:
                assert hf <= info['khovanov'], 'rank HF^(Y\') exceeds the Khovanov rank of K'
            if k == which and FZ is not None:
                forms.append((lab, status, FZ))
        if which:
            rr['verdict'] = rr['candidates'][which]['test']['obstructed']
        else:                                         # chirality not fixed: both chiralities must agree
            v = {c['test']['obstructed'] for c in rr['candidates'].values()}
            rr['verdict'] = v.pop() if len(v) == 1 and None not in v else None
    verds = {r['verdict'] for r in out['routes'] if r.get('verdict') is not None}
    assert len(verds) <= 1, 'routes disagree'
    out['verdict'] = 'excluded' if verds == {True} else 'passes' if verds == {False} else 'undetermined'
    if forms:                                         # the d-functions of the routes are isometric
        lab0, st0, F0 = next((f for f in forms if f[1] in ORIENTED), forms[0])
        out['routes_agree'] = {lab: bool(sd.same_form(F0, F) or (not (st in ORIENTED and st0 in ORIENTED)
                                                             and sd.same_form(F0, sd.negate(F))))
                               for lab, st, F in forms}
        assert all(out['routes_agree'].values()), 'the d-functions of two routes differ'
        out['route_used'], out['form'] = lab0, F0
    return out


# ---------- controls and the run

def self_test():
    """tower_V against surgery_dinv.V_model on the model complexes of surgery_dinv.model_checks, and the hat argument
    on them."""
    out = {}
    tests = {'T(2,3)': {1: 1, 0: -1, -1: 1}, 'T(3,4)': {3: 1, 2: -1, 0: 1, -2: -1, -3: 1},
             'T(3,5)': {4: 1, 3: -1, 1: 1, 0: -1, -1: 1, -3: -1, -4: 1}}
    for nm, al in tests.items():
        g = max(al)
        gl, ar = sd.staircase(al)
        for tag, (G, A) in (('staircase', (gl, ar)), ('dual', sd.dual_complex(gl, ar))):
            gens = dict(enumerate(G))
            V = [tower_V(gens, A, j) for j in range(g + 1)]
            an = analyse(gens, A, g if tag == 'staircase' else -g, g)
            out['%s %s' % (nm, tag)] = V == [sd.V_model(G, A, j, g) for j in range(g + 1)] == an['V']
    return out


def tensor_controls(cxs, BY):
    """Connected sums A # B of torus knots: steps 2-3 on the program's complex of a diagram of A # B against the
    V of the tensor product of the program's complexes of A and B (staircases, so CFK itself)."""
    out = {}
    for a, sa, b, sb in TENSOR:
        keys, _ = tensor_keys(BY, a, sa, b, sb)
        nm = '%s%s # %s%s' % ('' if sa == 1 else '-', a, '' if sb == 1 else '-', b)
        if any(cxs.get(k, {}).get('status') != 'ok' for k in keys):
            out[nm] = dict(flag='complex not computed')
            continue
        A, B, S = (complex_record(cxs[k]) for k in keys)
        assert not A['analysis']['lifts']['slots'] and not B['analysis']['lifts']['slots']   # (dual) staircases
        Tg, Ta = tensor((A['gens'], A['arrows']), (B['gens'], B['arrows']))
        g = A['genus'] + B['genus']
        assert S['genus'] == g
        V = [tower_V(Tg, Ta, j) for j in range(g + 1)]
        T = dict(gens=Tg, arrows=Ta, genus=g)
        out[nm] = dict(
            tau=S['tau'], genus=g, kind=S['analysis']['kind'], slots=len(S['analysis']['lifts']['slots']),
            V_tensor=V, V_hat=S['analysis']['V_hat'], V_lifts=S['analysis']['V_lifts'], V=S['analysis']['V'],
            argument=S['analysis']['argument'], same_hat=same_hat(S, T), agrees=S['analysis']['V'] == V)
    return out


def tensor_keys(BY, a, sa, b, sb):
    import spherogram
    A = spherogram.Link(json.loads(BY[a]['pd_notation']))
    B = spherogram.Link(json.loads(BY[b]['pd_notation']))
    A, B = (A if sa == 1 else A.mirror()), (B if sb == 1 else B.mirror())
    pds = [[list(c) for c in K.PD_code()] for K in (A, B, A.connected_sum(B))]
    return [pd_key(p) for p in pds], pds


def knot_inputs(BY, name, rec_sd):
    """KnotInfo data of K: Mullins' R(K), the Seifert matrix, the Khovanov rank, and Sigma_K (computed when needed)."""
    r = BY[name]
    D = rec_sd['det']
    assert int(r['determinant']) == D
    R = sd.mullins_R(D, int(r['signature']), r['jones_polynomial'])
    assert str(R) == rec_sd['R'], 'R(K) differs from surgery_dinv.py'
    cache = {}
    def sigma():
        if not cache:
            cache['S'] = sd.sigma_cover(json.loads(r['pd_notation']))
            assert str(cache['S'].homology()) == 'Z/%d' % D
        return cache['S']
    return dict(R=R, seifert=json.loads(r['seifert_matrix']), khovanov=khovanov_rank(BY, name), sigma_cover=sigma)


def jobs_of(parts, BY):
    """The diagrams for Szabo's program: C for every route, its mirror where step 4 needs it, and the tensor controls."""
    jobs = {}
    for part, recs in parts.items():
        for rec in recs.values():
            for lab, r in routes(rec):
                it = route_item(lab, r)
                jobs[pd_key(it['pd_C'])] = it['pd_C']
                if it['status'] is None and part != 'calibration':
                    jobs[pd_key(it['pd_mirror'])] = it['pd_mirror']
    for a, sa, b, sb in TENSOR:
        keys, pds = tensor_keys(BY, a, sa, b, sb)
        jobs.update(zip(keys, pds))
    return jobs


def chosen(r):
    """The record of the chirality of C fixed for route r (None if it is not fixed)."""
    k = {'kappa_hat': 'C', 'mirror of kappa_hat': 'mirror'}.get(r.get('C_is'))
    return r['candidates'][k] if k else None


def row(rec):
    """The route that decides a knot, in short."""
    r = next((r for r in rec['routes'] if r['route'] == rec.get('route_used')), None) or \
        next((r for r in rec['routes'] if r.get('verdict') is not None), None)
    if r is None:
        return dict(verdict=rec['verdict'], reason=sorted({x.get('flag', 'chirality not fixed') for x in rec['routes']}))
    c = chosen(r) or r['candidates']['C']
    return dict(verdict=rec['verdict'], route=r['route'], slope='%d/%d' % tuple(r['slope']), chirality=r['chirality'],
                kind=c['kind'], argument=c['analysis']['argument'], slots=len(c['analysis']['lifts']['slots']),
                V=c['analysis']['V'], passing=c['test'].get('passing'), routes=len(rec['routes']),
                routes_with_verdict=sum(x.get('verdict') is not None for x in rec['routes']))


def gap(BY, workers):
    """Mode gap: smooth_top_gap.one (by surgery_dinv.gap_records) for the targets that the test excludes."""
    import smooth_top_gap as stg
    C = json.load(open(os.path.join(ROOT, 'data', 'surgery_cfk.json')))
    names = C['summary']['excluded']
    partial = os.path.join(ROOT, 'data', 'surgery_cfk_gap.partial')
    recs = {k: v for k, v in sd.gap_records(names, BY, workers, partial).items() if k in names}
    per = {k: sd.gap_summary(v) for k, v in sorted(recs.items())}
    summary = dict(knots=names, certified=sorted(k for k, v in per.items() if v.get('certificate')),
                   not_certified=sorted(k for k, v in per.items() if not v.get('certificate')),
                   contradictions=sum(v.get('contradictions') or 0 for v in per.values()), knots_summary=per)
    parameters = dict(t_max=stg.T_MAX, extra=stg.EXTRA, minimal_target=stg.MINIMAL_TARGET,
                      nonminimal_tries=stg.NONMINIMAL_TRIES)
    out = dict(meta=dict(versions=sd.versions(), parameters=parameters), summary=summary, records=dict(sorted(recs.items())))
    json.dump(out, open(os.path.join(ROOT, 'data', 'surgery_cfk_gap.json'), 'w'), default=str)
    print(json.dumps(out['summary'], indent=1, default=str))
    assert not summary['contradictions'] and not any(v.get('error') for v in per.values())


if __name__ == '__main__':
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}
    if sys.argv[1:2] == ['gap']:
        gap(BY, int(sys.argv[2]) if len(sys.argv) > 2 else 2)
        sys.exit()
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    load = lambda f: json.load(open(os.path.join(ROOT, 'data', f)))
    SD, SB, PL = load('surgery_dinv.json'), load('surgery_box.json'), load('plumbing_dinv.json')
    G4 = load('gamma4_knotinfo_2026.10.5.json')
    targets = sorted(set(SD['targets']) - set(PL['summary']['obstructed']) - set(SB['summary']['box_obstructed']))
    others = sorted(set(SD['targets']) - set(targets))
    assert len(targets) == 33 and len(others) == 17 and all(G4[n] == '[1,2]' for n in targets)
    parts = dict(targets={n: SD['targets'][n] for n in targets}, others={n: SD['targets'][n] for n in others},
                 pilot=SD['pilot'], calibration=SD['calibration'])
    partial = os.path.join(ROOT, 'data', 'surgery_cfk.partial')
    jobs = jobs_of(parts, BY)
    cxs = run_jobs(jobs, workers, partial)
    covers, res = {}, {}
    for part, recs in parts.items():
        res[part] = {}
        for n, rec in sorted(recs.items()):
            res[part][n] = one(n, rec, knot_inputs(BY, n, rec), cxs, covers) if routes(rec) else \
                dict(name=n, det=rec['det'], routes=[], verdict='no route')
            print(part, n, res[part][n]['verdict'], flush=True)
    decided = {k: v for k, v in covers.items() if v.get('decided')}
    used = {r['covers']['key'] for part in res.values() for v in part.values() for r in v['routes']
            if r['route'] == v.get('route_used') and r.get('chirality') == 'covers'}
    gap_ok = gap_check({k: decided[k] for k in sorted(used)})     # GAP for the route that decides each knot
    assert gap_ok is None or all(gap_ok.values()), 'GAP and SnapPy covers differ'
    # controls
    cal = {}
    for n, v in res['calibration'].items():
        if 'form' not in v:
            cal[n] = dict(verdict=v['verdict'])
            continue
        G = sd.goeritz_P(json.loads(BY[n]['pd_notation']))
        c = sd.calibrate(G, v['form'])
        cal[n] = dict(verdict=v['verdict'], goeritz_agrees=c['agrees'], units=c['units'],
                      goeritz_excluded=test(Form(G), (1, -1)) is None, knotinfo=G4.get(n))
        cal[n]['same_verdict'] = cal[n]['goeritz_excluded'] == (v['verdict'] == 'excluded')
    oth = {}
    for n, v in res['others'].items():
        assert len(v['routes']) == len(SB['targets'][n]['routes'])
        inbox = []
        for r, b in zip(v['routes'], SB['targets'][n]['routes']):    # the box of surgery_box.py, route by route
            c = chosen(r)
            if c and c['analysis']['V'] is not None:
                V = c['analysis']['V'] + [0] * (len(b['V_max']) - len(c['analysis']['V']))
                inbox.append(all(lo <= x <= hi for x, lo, hi in zip(V, b['V_min'], b['V_max'])))
        o = dict(verdict=v['verdict'], V_in_box=all(inbox) and bool(inbox))
        fib = (SD['targets'][n].get('plumbing') or {}).get('fibres')
        if fib and 'form' in v:
            import montesinos_dinv_check as MC
            Q, _ = MC.plumbing(fib)
            N, u, dd = MC.d_invariants(Q)
            FP = sd.Cyclic(N, u, [dd[k] for k in range(N)])
            o['plumbing_agrees'] = bool(sd.same_form(v['form'], FP) or sd.same_form(v['form'], sd.negate(FP)))
        oth[n] = o
    import surgery_box as SBX
    for n, v in res['targets'].items():               # the box of surgery_box.py, with V_min from the fixed chirality
        if v['verdict'] == 'excluded':
            v['box_with_chirality'] = {r['route']: SBX.box(v['det'], r['e'], chosen(r)['tau'], b['g_alg'])['obstructed']
                                       for r, b in zip(v['routes'], SB['targets'][n]['routes'])
                                       if r.get('chirality') in ORIENTED}
    ten = tensor_controls(cxs, BY)
    meta = dict(versions=sd.versions(), lens_labelling_check=sd.lens_check(), self_test=self_test(),
                parameters=dict(cfk_seconds=CFK_SECONDS, max_slots=MAX_SLOTS, cover_degrees=COVER_DEGREES, tensor=TENSOR),
                jobs=len(jobs), jobs_not_ok={k: v['status'] for k, v in cxs.items() if k in jobs and v['status'] != 'ok'})
    T = res['targets']
    problems = dict(
        meta=[k for k, v in list(meta['lens_labelling_check'].items()) + list(meta['self_test'].items()) if not v],
        calibration=[n for n, c in cal.items() if c['verdict'] != 'no route'
                     and not (c.get('goeritz_agrees') and c.get('same_verdict'))],
        calibration_gamma4_one_excluded=[n for n, c in cal.items() if c['verdict'] == 'excluded' and G4.get(n) == '1'],
        others=[n for n, o in oth.items()
                if o['verdict'] != 'excluded' or not o['V_in_box'] or o.get('plumbing_agrees') is False],
        pilot=[n for n, v in res['pilot'].items() if v['verdict'] != 'passes'],
        tensor=[k for k, t in ten.items() if not (t.get('agrees') and t.get('same_hat'))],
        gap=None if gap_ok is None else [k for k, x in gap_ok.items() if not x])
    summary = dict(
        targets=len(T), excluded=sorted(n for n, v in T.items() if v['verdict'] == 'excluded'),
        excluded_by_box_with_chirality=sorted(n for n, v in T.items() if v.get('box_with_chirality')
                                              and all(v['box_with_chirality'].values())),
        passes=sorted(n for n, v in T.items() if v['verdict'] == 'passes'),
        undetermined={n: row(v) for n, v in T.items() if v['verdict'] == 'undetermined'},
        table={n: row(v) for n, v in sorted(T.items())},
        kinds=dict(Counter(row(v).get('kind') for v in T.values())),
        chirality=dict(Counter(row(v).get('chirality') for v in T.values())),
        covers=dict(decided=len(decided), gap_rechecked=None if gap_ok is None else len(gap_ok),
                    routes=sum(1 for v in T.values() for r in v['routes'] if r.get('chirality') == 'covers')),
        controls=dict(calibration=dict(knots=len(cal), with_route=sum(c['verdict'] != 'no route' for c in cal.values()),
                                       goeritz_agrees=sum(bool(c.get('goeritz_agrees')) for c in cal.values()),
                                       excluded=sorted(n for n, c in cal.items() if c['verdict'] == 'excluded')),
                      others=dict(knots=len(oth), excluded=sum(o['verdict'] == 'excluded' for o in oth.values()),
                                  plumbing_agrees=sorted(n for n, o in oth.items() if o.get('plumbing_agrees'))),
                      pilot={n: row(v) for n, v in res['pilot'].items()},
                      tensor={k: dict(V=t.get('V'), V_tensor=t.get('V_tensor'), kind=t.get('kind'), slots=t.get('slots'),
                                      agrees=t.get('agrees'), same_hat=t.get('same_hat')) for k, t in ten.items()}),
        problems=problems)
    for part in res.values():
        for v in part.values():
            v.pop('form', None)
    out = dict(meta=meta, summary=summary, targets=res['targets'], others=res['others'], pilot=res['pilot'],
               calibration=res['calibration'], controls=dict(calibration=cal, others=oth, tensor=ten),
               covers=dict(covers),
               complexes={k: dict(generators=v['generators'], differentials=[[x, y] for x, y, c in v['differentials']])
                          for k, v in cxs.items() if k in jobs and v['status'] == 'ok'})
    json.dump(out, open(os.path.join(ROOT, 'data', 'surgery_cfk.json'), 'w'), indent=0, default=str)
    keys = ('excluded', 'excluded_by_box_with_chirality', 'passes', 'undetermined', 'kinds', 'chirality', 'covers',
            'problems')
    print(json.dumps({k: summary[k] for k in keys}, indent=1, default=str))
    assert not any(v for v in problems.values() if v), problems
