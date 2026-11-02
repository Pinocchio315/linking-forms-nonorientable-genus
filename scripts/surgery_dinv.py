#!/usr/bin/env python
"""[C] Proposition 3.2(b) and [C] Section 3.8 (11n_176): correction terms of Sigma(K) from surgery descriptions for
knots with unknotting number one (the Montesinos trick), and the test of [C] Proposition 2.6.

Setting.  If switching one crossing c of a diagram of K gives the unknot U, then Sigma(K) = S^3_{e D/2}(kappa) for a
sign e and a knot kappa in S^3 = Sigma(U), the lift of the crossing arc; D = det K (Montesinos trick).  The script
finds c, kappa and e, computes d(Sigma(K), t) by the formula of Ni and Wu, and applies the d-test of Proposition 2.6
(rank1_dinv.test).  Mode targets treats the 50 knots of Proposition 3.2(b), the open knots with KnotInfo u = 1 that are
neither alternating nor listed as Montesinos; mode pilot treats 11n_176 (Section 3.8), and mode calibration the 32
alternating knots of CALIBRATION.  KnotInfo gives no Montesinos notation for 13 crossings, and Regina recognises some of
these Sigma(K) as Seifert fibred (see 7).

Pipeline for one knot (route()); every step is checked and the record of the knot keeps the checks.
1. Unknotting crossing.  Switch each crossing of the KnotInfo PD code (switch(): [a,b,c,d] becomes [b,c,d,a] or
   [d,a,b,c], whichever starts with an incoming strand).  Certificate: spherogram simplify('global'), which uses only
   Reidemeister moves, reduces the switched diagram to no crossings.
2. Crossing circle C around the two edges that run into c (with_circle()), so lk(U, C) = +-2.  Checks: filling C by
   its meridian gives E(K); in the switched diagram U and C reduce to no crossings; filling C by (1, f) gives E(K) for
   exactly one f in {1, -1} (is_isometric_to; True is rigorous).
3. Double branched cover over U.  N = the 2-fold cover of E(U u C)(2, 0) with two cusps over C; N is the exterior of the
   two lifts C1, C2 in S^3, parallel copies of kappa.  SnapPy's peripheral curves on the cusps of a cover need not be
   lifts of those of C, so slopes s1, s2 with N(s1, s2) = S^3 are searched, certified by Regina's 3-sphere
   recognition.  Lemma: if kappa is nontrivial, these pairs are exactly (m1 + k h, m2 - k h), k in Z (m_i the
   meridians, h the fibre of (pair of pants) x S^1 = N(kappa) - C1 - C2), because S^3 has no incompressible torus and
   E(kappa) has incompressible boundary, so the pants side filled along s1, s2 is a solid torus with meridian
   mu_kappa.  The twist along the annulus between C1 and C2 is an orientation-preserving self-homeomorphism of N
   realising the shift by k.  Hence X = N(s2) is E(kappa) with meridian s1, oriented.  (kappa is nontrivial:
   otherwise Sigma(K) is a lens space and K is 2-bridge, Hodgson-Rubinstein; also checked.)
4. The sign e, and a diagram of kappa.  X gets the basis (s1, homological longitude); X(1, 0) = S^3 (Regina).  The
   sign e is decided by every source that decides it, and they must agree (decide_sign()): the Casson-Walker invariant
   (step 6), the linking form, drilling Sigma(K) along a short geodesic (an orientation-preserving isometry to X
   sending the drilled meridian to +-(e D, 2)), and closed isometries Sigma(K) = X(e D, 2).
   A diagram of kappa: exterior_to_link (Dunfield-Obeidin-Rudd; with careful_perturbation the exterior and the
   meridian are guaranteed), then an isometry from its exterior to X with cusp map of determinant +1 sending the
   meridian to +-meridian (else the mirror).  When X is not hyperbolic: Regina names (regina_chirality()); if the
   name of S^3_{e D/2}(diagram) is a name of Sigma(K), then Sigma(K) = +-S^3_{e D/2}(diagram), which decides the
   unsigned test (no signed verdict is then recorded); otherwise the plumbing (7), or else the verdict must be the
   same for both chiralities.
5. V_j.  Y' = S^3_{D/2}(kappa_hat) = e Sigma(K), kappa_hat = kappa (e = 1) or its mirror (e = -1).  HFK of kappa by
   knot_floer_homology (Szabo's program; girth <= 20, otherwise a Reidemeister variant; in a separate process with a
   time limit); HFK of the mirror by negating both gradings.  V_from(): HFK of staircase shape with tau = g > 0 (an
   L-space knot): V_j = t_j (torsion coefficients); the dual staircase (tau = -g, the mirror of an L-space knot):
   V_j = 0; thin HFK (Maslov = Alexander - tau): V_j = ceil((tau - j)/2) for 0 <= j <= tau, V_j = 0 if tau <= 0
   (Petkova: staircase of T(2, 2 tau + 1) plus acyclic boxes).  Each case is checked against V_model(), an exact F_2
   computation of H_*(A_s^+) of the model complex.  Otherwise V_j is not determined: the knot Floer homology is
   recorded, and bounds_test() runs the d-test for every V allowed by V_j >= V_{j+1} >= V_j - 1,
   V_j <= ceil((g - j)/2) and V_j >= 1 for j < nu (nu^+ >= nu, Hom-Wu); the knot stays flagged.  This sharper
   upper bound, from Rasmussen's later work and not checked here against the original, enters only the field
   obstructed_by_bounds, which [C] does not use: the counts of [C] come from surgery_box.py, with the bound
   V_j <= g - j of Rasmussen's thesis (Cor. 7.4).
   Delta of kappa: Fox calculus on the diagram (PARI) against the Euler characteristic of HFK.
6. d-invariants.  Ni-Wu: d(S^3_{p/q}(K), i) = d(L(p, q), i) - 2 max(V_{floor(i/q)}, V_{floor((p+q-1-i)/q)}), p, q > 0.
   Labels: t_s + k mu has label i_s + q k with 2 i_s = q - 1 mod p, and lambda(mu, mu) = -q/p (the Gordon-Litherland
   sign of [C]), checked on lens spaces against the sharp linear plumbing (lens_check).
   Casson-Walker: -2 p lambda_CW(S^3_{p/q}(K)) = sum_i d(L(p, q), i) - 2 q sum_{k in Z} t_k(K), and Mullins gives
   -2 D lambda_CW(Sigma(K)) = R(K) = -D sigma/4 + D V'(-1)/(6 V(-1)) (KnotInfo signature and Jones polynomial), so
   e (sum_i d(L(D, 2), i) - 4 sum t_k) = R(K).  This uses only Delta(kappa), so it does not check the chirality of
   kappa, and for an L-space knot with D/2 >= 2g - 1 the equality sum_t d(Y', t) = e R(K) follows from it.
7. Other routes and checks.  If no unknotting crossing of the KnotInfo diagram gives V_j: other diagrams of K
   (other_diagrams(), checked isometric to the KnotInfo exterior) and their unknotting crossings, each new kappa once.
   Drilling routes (drill_route()): Sigma(K) drilled along a short geodesic, an S^3 filling among the slopes of length
   <= 6 (6-theorem), so Sigma(K) = S^3_{p/q}(kappa'); same steps 4-6 with that p/q.  Plumbing (plumbing_check(), the
   code of montesinos_dinv_check.py): if Regina recognises Sigma(K) as a Seifert fibred space over S^2, the
   negative-definite star plumbing gives the d-invariants of +-Sigma(K) independently of kappa; its d-function must be
   isometric to +-that of the routes, and it decides the chirality of a diagram of a non-hyperbolic X.  The
   d-functions of all routes must be isometric (same_form()).
8. d-test.  Z = -Sigma(K) = -e Y' is the Z of [C], whose linking form is (V + V^T)^-1 (KnotInfo Seifert matrix;
   checked by Jordan data).  Unsigned test: test(F, (1, -1)), the criterion; it is repeated by the independent
   implementation for cyclic H_1 in montesinos_dinv_check.obstructed.  Signed test as in rank1_dinv.one
   (s = sigma + 4 Arf mod 8 = 2 allows only Z0 = Z, s = 6 only Z0 = -Z).  For a knot that passes, the first passing
   (sign, l, |M|) is recorded.
9. Obstructed knots (mode gap): the search of smooth_top_gap.py, a band to a knot with Delta = 1, hence a locally flat
   Mobius band and gamma_4^top = 1 (smooth_top_gap.one with its rechecks).
Verdicts: obstructed (the unsigned test fails for the computed d-function); passes (the first passing (sign, l, |M|) is
recorded); obstructed_by_bounds (V_j not determined, but every V allowed by the bounds of step 5 is obstructed: a
result conditional on those inequalities); passes_for_all_V (every V allowed by the chirality-free bounds passes);
flagged (with the reason).
Mode summary also applies bounds_verdict() to the flagged knots and tries torus_sum_route(): if HFK of kappa_hat is
that of a connected sum A # B of two torus knots and Regina gives S^3_{D/2}(A # B) a name of Sigma(K), then
Sigma(K) = +-S^3_{D/2}(A # B) and V_j comes exactly from the tensor product of the two complexes.  It uses the rank of
the reduced Khovanov homology of K (khovanov_check(): rank D makes Sigma(K) an L-space), recomputes, for every
obstructed knot whose kappa is a census knot, HFK from a diagram of the census triangulation (census_recheck()), and
the bounds with the independent implementation (bounds_recheck()), and checks that the stored d-functions give the
stored verdicts (stored_check()).
Calibration (mode calibration, [C] Section 3.5): alternating knots with u = 1, not Montesinos, with hyperbolic
Sigma(K); the d-function must agree with the Goeritz d-function of rank1_dinv (an isomorphism of groups carrying lambda
and d, with Z = -e Y'), and no knot with KnotInfo value 1 may be obstructed.
Usage: python scripts/surgery_dinv.py calibration | targets | pilot | gap [workers] | summary
  The modes calibration, targets and pilot first; gap searches the obstructed knots of the checkpoint of targets;
  summary last.  All modes need Regina.
Output: data/surgery_dinv.json (meta, summary, calibration, pilot (the record of 11n_176), targets, gap) and
data/surgery_dinv_gap.json (the diagrams and bands of the gap search), written by mode summary from the checkpoints
data/surgery_dinv_<mode>.partial of the other modes.
"""
import json, os, random, re, sys, time, traceback, warnings
from fractions import Fraction as Fr
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from rank1_dinv import Form, test
from rank1_linking import invariants, legendre, primes_of
from rank1_goeritz import goeritz_pair
from smooth_top_gap import orientation, alexander_fox, alexander_hfk, pari

PARI_STACK, PARI_SIZEMAX = 2 ** 26, 2 ** 31
SIMPLIFY_TRIES = 20                       # simplify('global') attempts for an unknotting certificate
SLOPE_BOX = 4                             # |a|, |b| <= SLOPE_BOX for the slopes s1, s2 on the lifts of C
ISO_TRIES = 30                            # randomised retries of a closed is_isometric_to
DRILL_SEGMENTS = 8                        # max_segments of the dual curves of Sigma(K) that are drilled
KNOT_SECONDS = 3600                       # soft time budget per knot (no new route after a third or a half of it);
                                          # a knot still running after twice this is stopped and recorded as timed out
HFK_SECONDS = 900                         # time limit of one knot Floer homology computation
REGINA_SECONDS = 60                       # time limit of one 3-sphere recognition of a candidate filling
MAX_GIRTH, GIRTH_TRIES = 20, 30           # Szabo's program needs girth <= 20; Reidemeister variants tried


# ---------- step 1: an unknotting crossing

def switch(pd, c):
    """The PD code with crossing c switched; the new entry starts with its incoming strand."""
    entering, _ = orientation([tuple(x) for x in pd])
    a, b, x, d = pd[c]
    new = [list(X) for X in pd]
    new[c] = [b, x, d, a] if (c, 1) in entering else [d, a, b, x]
    return new


def simplifies_to_unknot(pd, tries=SIMPLIFY_TRIES):
    """Number of simplify('global') runs until the diagram has no crossings, or None."""
    import spherogram
    for t in range(tries):
        L = spherogram.Link([list(x) for x in pd])
        L.simplify('global')
        if len(L.crossings) == 0:
            return t + 1
    return None


def unknotting_crossings(pd):
    """The crossings of pd whose switch gives a diagram that simplify('global') reduces to no crossings."""
    out = []
    for c in range(len(pd)):
        S = switch(pd, c)
        t = simplifies_to_unknot(S)
        if t:
            out.append(dict(crossing=c, entry=list(pd[c]), switched_entry=S[c], simplify_runs=t))
    return out


# ---------- step 2: the crossing circle

def with_circle(pd, c):
    """K u C, C around the two edges that run into crossing c (slot 0 and the incoming over-slot)."""
    entering, _ = orientation([tuple(x) for x in pd])
    f = 0 if (c, 1) in entering else 3
    g = (f + 1) % 4
    a, b = pd[c][f], pd[c][g]
    assert list(pd[c]).count(a) == 1 and list(pd[c]).count(b) == 1 and a != b
    big = max(max(X) for X in pd) + 1
    e0a, e0b, e0c, e1a, e1b, e1c, c1, c2, c3, c4 = range(big, big + 10)
    new = []
    for i, X in enumerate(pd):
        X = list(X)
        if i == c:
            X[f], X[g] = e0a, e1a
        else:
            X = [e0c if v == a else e1c if v == b else v for v in X]
        new.append(X)
    new += [[c1, e0c, c2, e0b],           # C under the edge at slot f
            [e0b, c4, e0a, c1],           # C over it
            [e1b, c3, e1a, c4],           # C over the edge at slot f + 1
            [c2, e1c, c3, e1b]]           # C under it
    return new, (f, g)


def isometric(A, B, tries=8):
    """True (rigorous), False (SnapPy's answer, not rigorous) or None (no answer)."""
    for _ in range(tries):
        try:
            return bool(A.is_isometric_to(B))
        except Exception:
            A.randomize(); B.randomize()
    return None


def circle_checks(pd, c):
    """The checks of step 2 for the crossing circle at crossing c: the record, U u C, and the indices of C and U."""
    import spherogram
    KC_pd, slots = with_circle(pd, c)
    KC = spherogram.Link(KC_pd)
    comps = KC.link_components
    ci = [i for i, comp in enumerate(comps) if len(comp) == 4]
    assert len(comps) == 2 and len(ci) == 1, 'crossing circle: wrong components'
    ci = ci[0]; ki = 1 - ci
    EK = spherogram.Link([list(x) for x in pd]).exterior()
    E = KC.exterior(); E.dehn_fill((1, 0), ci)
    rec = dict(slots=list(slots), KC_pd=KC_pd, C_index=ci, meridian_filling_is_K=isometric(E, EK))
    UC_pd = switch(KC_pd, c)
    UC = spherogram.Link(UC_pd)
    assert [len(x) for x in UC.link_components] == [len(x) for x in comps]
    rec['UC_pd'] = UC_pd
    rec['linking_number'] = UC.linking_matrix()[0][1]
    rec['U_unknot_runs'] = simplifies_to_unknot([list(x) for x in UC.sublink([ki]).PD_code()]) if len(UC.sublink([ki]).crossings) else 0
    rec['C_unknot_runs'] = simplifies_to_unknot([list(x) for x in UC.sublink([ci]).PD_code()]) if len(UC.sublink([ci]).crossings) else 0
    EU = UC.exterior()
    fill = {}
    for f in (1, -1):
        G = EU.copy(); G.dehn_fill((1, f), ci)
        fill[f] = isometric(G, EK.copy())
    rec['C_filling_gives_K'] = {str(f): v for f, v in fill.items()}
    rec['ok'] = (rec['meridian_filling_is_K'] is True and rec['U_unknot_runs'] is not None and rec['C_unknot_runs'] is not None
                 and abs(rec['linking_number']) == 2 and sorted(v is True for v in fill.values()) == [False, True])
    return rec, UC, ci, 1 - ci


# ---------- step 3: the double branched cover over U, and an S^3 filling of the two lifts of C

def ab(word, gens):
    """The exponent sums of a word of SnapPy's fundamental group (lower case: a generator, upper case: its inverse)."""
    v = dict.fromkeys(gens, 0)
    for ch in word:
        if ch.islower():
            v[ch] += 1
        else:
            v[ch.lower()] -= 1
    return [v[g] for g in gens]


def regina_sphere(M):
    """Regina's 3-sphere recognition of the Dehn filling M (all cusps filled): exact."""
    import regina
    T = regina.SnapPeaTriangulation(M._to_string()).filledAll()
    T.simplify()
    return bool(T.isSphere())


def _sphere_child(text, q):
    try:
        import regina
        T = regina.SnapPeaTriangulation(text).filledAll()
        T.simplify()
        q.put(bool(T.isSphere()))
    except BaseException as ex:
        q.put(repr(ex))


def regina_sphere_timed(M, limit=REGINA_SECONDS):
    """regina_sphere in a separate process: True, False, or None when the time limit is reached."""
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    q = ctx.Queue()
    p = ctx.Process(target=_sphere_child, args=(M._to_string(), q))
    p.start()
    t0 = time.time()
    while time.time() - t0 < limit:
        try:
            r = q.get(timeout=0.5)
            p.join(5)
            return r if isinstance(r, bool) else None
        except Exception:
            if not p.is_alive():
                return None
    p.kill(); p.join()
    return None


def s3_candidate(M):
    """For a Dehn filling M with trivial H_1: True if Regina recognises S^3.  Fillings that SnapPy finds hyperbolic
    (positively oriented solution, volume > 0.5) are skipped as certainly not S^3 (this only discards candidates;
    every S^3 used is certified by Regina), since 3-sphere recognition can be very slow on other homology spheres;
    Regina runs in a separate process with a time limit, and a candidate it does not certify in time is skipped."""
    if M.fundamental_group().num_generators() == 0:
        return regina_sphere_timed(M, 4 * REGINA_SECONDS) is True
    H = improved(M, 10)
    if hyperbolic(H) and H.volume() > 0.5:
        return False
    return regina_sphere_timed(M) is True


def slope_list(box):
    """The slopes (a, b), b >= 0, with |a|, b <= box, ordered by |a| + |b|, (1, 0) first."""
    out = [(a, b) for a in range(-box, box + 1) for b in range(0, box + 1) if gcd(a, b) == 1 and (b > 0 or a > 0)]
    return sorted(out, key=lambda s: (abs(s[0]) + abs(s[1]), s != (1, 0), s))


def sphere_pair(N, j1, j2, box=SLOPE_BOX):
    """Slopes (s1, s2) on cusps j1, j2 with N(s1, s2) = S^3; trivial H_1 first, then pi_1 and Regina."""
    S = slope_list(box)
    pairs = sorted(((s, t) for s in S for t in S), key=lambda st: (sum(map(abs, st[0])) + sum(map(abs, st[1])), st))
    tried = 0
    for s1, s2 in pairs:
        M = N.copy(); M.dehn_fill(s1, j1); M.dehn_fill(s2, j2)
        if M.homology().order() != 1:
            continue
        tried += 1
        G = M.fundamental_group()
        trivial = G.num_generators() == 0
        if s3_candidate(M):
            return dict(s1=list(s1), s2=list(s2), pi1_trivial_presentation=trivial, regina_sphere=True,
                        homology_spheres_tried=tried)
    return None


def double_cover(UC, ui):
    """The double cover N of E(U u C)(2, 0) branched over U with two cusps over C, and those two cusps."""
    E = UC.exterior(); E.dehn_fill((2, 0), ui)
    covs = E.covers(2)
    good = []
    for N in covs:
        fl = N.cusp_info('filling')
        openc = [i for i, f in enumerate(fl) if f == (0.0, 0.0)]
        closed = [i for i, f in enumerate(fl) if f != (0.0, 0.0)]
        if len(openc) == 2 and len(closed) == 1 and fl[closed[0]] == (1.0, 0.0):
            good.append((N, openc))
    assert len(good) == 1, 'double branched cover: %d candidates' % len(good)
    return good[0]


def h1_free_rank_one(M):
    """For M with H_1 = Z: the map of the peripheral curves of cusp 0 to H_1 = Z (exact, from pi_1)."""
    import sympy
    G = M.fundamental_group()
    gens = G.generators()
    R = sympy.Matrix([ab(r, gens) for r in G.relators()]) if G.relators() else sympy.zeros(0, len(gens))
    ker = R.nullspace() if R.rows else [sympy.eye(len(gens))[:, i] for i in range(len(gens))]
    assert len(ker) == 1, 'H_1 has rank %d' % len(ker)
    v = ker[0]
    den = sympy.ilcm(*[x.q for x in v])
    v = [int(x * den) for x in v]
    g = 0
    for x in v:
        g = gcd(g, x)
    v = [x // g for x in v]
    m, l = G.peripheral_curves()[0]
    return sum(a * b for a, b in zip(v, ab(m, gens))), sum(a * b for a, b in zip(v, ab(l, gens)))


def kappa_exterior(N, j2, pair):
    """X = N(s2 on j2), with peripheral basis (s1, homological longitude) on the remaining cusp."""
    import snappy
    s1, s2 = pair['s1'], pair['s2']
    M = N.copy(); M.dehn_fill(s2, j2)
    X = snappy.Manifold(M.filled_triangulation())
    assert X.num_cusps() == 1
    a, b = s1
    # complete s1 to a basis with determinant +1
    c, d = next((c, d) for c in range(-abs(b) - 1, abs(b) + 2) for d in range(-abs(a) - 1, abs(a) + 2) if a * d - b * c == 1)
    X.set_peripheral_curves([(a, b), (c, d)], 0)
    im_m, im_l = h1_free_rank_one(X)
    assert abs(im_m) == 1, 'the meridian does not generate H_1'
    x = -im_l * im_m
    X.set_peripheral_curves([(1, 0), (x, 1)], 0)
    im_m2, im_l2 = h1_free_rank_one(X)
    Y = X.copy(); Y.dehn_fill((1, 0))
    rec = dict(H1=str(X.homology()), meridian_image=im_m2, longitude_image=im_l2, filling_10_sphere=regina_sphere(Y),
               solution_type=X.solution_type(), volume=round(float(X.volume()), 8))
    try:
        rec['identify'] = [str(m) for m in X.identify()]
    except Exception:
        rec['identify'] = []
    G = X.fundamental_group()
    rec['pi1_generators'] = G.num_generators()
    rec['ok'] = (rec['H1'] == 'Z' and im_l2 == 0 and abs(im_m2) == 1 and rec['filling_10_sphere'] and G.num_generators() >= 2)
    return X, rec


# ---------- step 4: the sign e, and a diagram of kappa

def sigma_cover(pd):
    """Sigma(K): the 2-fold cyclic cover of the orbifold E(K)(2, 0), as a closed manifold."""
    import snappy, spherogram
    E = spherogram.Link([list(x) for x in pd]).exterior()
    E.dehn_fill((2, 0))
    S = E.covers(2, cover_type='cyclic')
    assert len(S) == 1
    return snappy.Manifold(S[0].filled_triangulation())


def closed_isometric(A, B, tries=ISO_TRIES):
    import snappy
    A = snappy.Manifold(A); B = snappy.Manifold(B)
    for t in range(tries):
        try:
            return bool(A.is_isometric_to(B)), t + 1
        except Exception:
            A.randomize(); B.randomize()
    return None, tries


def hyperbolic(M):
    return M.solution_type() == 'all tetrahedra positively oriented'


def improved(M, tries=30):
    """A copy of M retriangulated until SnapPy finds a solution with positively oriented tetrahedra (if it can)."""
    import snappy
    M = snappy.Manifold(M) if not isinstance(M, snappy.Manifold) else M.copy()
    for _ in range(tries):
        if hyperbolic(M):
            break
        M.randomize()
    return M


def slope_checks(SigK, X, D):
    """Closed isometries Sigma(K) = X(e D, 2) (unoriented), and isometries from Sigma(K) drilled along a short geodesic
    to X (cusp maps: the first column is the image of the meridian of the drilled curve)."""
    import snappy
    S = improved(SigK)
    out = dict(sigma_solution_type=S.solution_type(), sigma_volume=round(float(S.volume()), 8))
    for e in (1, -1):
        Y = X.copy(); Y.dehn_fill((e * D, 2))
        F = improved(snappy.Manifold(Y.filled_triangulation()))
        r = dict(volume=round(float(F.volume()), 8), solution_type=F.solution_type(), H1=str(Y.homology()))
        if hyperbolic(S) and hyperbolic(F):
            r['isometric'], r['tries'] = closed_isometric(F, S)
        out['filling_%+d' % e] = r
    drill = []
    if hyperbolic(S):
        try:
            curves = S.dual_curves(max_segments=DRILL_SEGMENTS)
        except Exception:
            curves = []
        Xh = improved(X)
        for i, g in enumerate(curves):
            try:
                Dr = improved(S.drill(i, max_segments=DRILL_SEGMENTS), 10)
                iso = Dr.is_isometric_to(Xh, return_isometries=True) if hyperbolic(Dr) else []
            except Exception:
                continue
            if not iso:
                continue
            maps = []
            for h in iso:
                cm = h.cusp_maps()[0]
                m = [[int(cm[0, 0]), int(cm[0, 1])], [int(cm[1, 0]), int(cm[1, 1])]]
                maps.append(dict(cusp_map=m, det=m[0][0] * m[1][1] - m[0][1] * m[1][0], image_of_meridian=[m[0][0], m[1][0]]))
            drill.append(dict(curve=i, curve_data=str(g), isometries=maps))
            break
    out['drilling'] = drill
    return out


def decide_sign(cw, lf, slope, D):
    """The sign e with Sigma(K) = S^3_{e D/2}(kappa), from every source that decides it; they must agree.
    Casson-Walker (exact; decides unless both signs fit), the linking form (exact; decides unless <2/D> = <-2/D>),
    drilling (an orientation-preserving isometry, rigorous), closed isometry (True for one sign, False for the other
    with a different volume; by the Montesinos trick the true sign is among those that are isometric)."""
    src = {}
    if cw['1'] != cw['-1']:
        src['casson_walker'] = 1 if cw['1'] else -1
    if lf['plus'] != lf['minus']:
        src['linking_form'] = 1 if lf['plus'] else -1
    dr = set()
    for d in slope['drilling']:
        for m in d['isometries']:
            if m['det'] == 1 and m['image_of_meridian'] in ([D, 2], [-D, -2]):
                dr.add(1)
            elif m['det'] == 1 and m['image_of_meridian'] in ([-D, 2], [D, -2]):
                dr.add(-1)
    if len(dr) == 1:
        src['drilling'] = dr.pop()
    elif len(dr) == 2:
        src['drilling'] = 0
    f1, f2 = slope['filling_+1'], slope['filling_-1']
    if abs(f1['volume'] - f2['volume']) > 1e-4:
        if f1.get('isometric') is True and f2.get('isometric') is not True:
            src['closed_isometry'] = 1
        elif f2.get('isometric') is True and f1.get('isometric') is not True:
            src['closed_isometry'] = -1
    vals = set(src.values())
    if not cw['1'] and not cw['-1']:
        return None, src, 'Casson-Walker check fails for both signs'
    if len(vals) != 1 or 0 in vals:
        return None, src, 'sign sources: %s' % src
    return vals.pop(), src, None


def kappa_diagram(X):
    """A diagram of kappa: exterior_to_link (Dunfield-Obeidin-Rudd) with careful_perturbation, which guarantees the
    exterior and the meridian; its own final check (check_answer, up to 1000 high-precision retriangulations) is
    replaced by an isometry from the exterior of the diagram to X with cusp map of determinant +1 sending the meridian
    to +-meridian (the mirror diagram if only orientation-reversing ones exist)."""
    import snappy
    snappy.set_rand_seed(0); random.seed(0)
    L = X.exterior_to_link(check_answer=False)
    for _ in range(3):
        L.simplify('global')
    rec = dict(crossings=len(L.crossings), orientation=None)
    det = lambda c: c[0, 0] * c[1, 1] - c[0, 1] * c[1, 0]
    Xh = improved(X)
    if hyperbolic(Xh):
        found = {}
        for name, K in (('diagram', L), ('mirror', L.mirror())):
            E = improved(K.exterior())
            iso = []
            for _ in range(10):
                try:
                    iso = E.is_isometric_to(Xh, return_isometries=True)
                    break
                except Exception:
                    E.randomize()
            cms = [h.cusp_maps()[0] for h in iso] if iso else []
            found[name] = [c for c in cms if det(c) == 1 and abs(c[0, 0]) == 1 and c[1, 0] == 0]
            if found[name]:
                rec['cusp_maps_%s' % name] = [[[int(c[0, 0]), int(c[0, 1])], [int(c[1, 0]), int(c[1, 1])]] for c in cms]
        if found['diagram'] and not found['mirror']:
            rec['orientation'] = 'diagram'
        elif found['mirror'] and not found['diagram']:
            rec['orientation'] = 'mirror'; L = L.mirror()
        elif found['diagram'] and found['mirror']:
            rec['orientation'] = 'amphichiral'
    rec['pd'] = [list(x) for x in L.PD_code()]
    return L, rec


# ---------- step 5: knot Floer homology and V_j

def torsion(alex):
    """t_j for j >= 0 of a symmetric Delta = {exponent: coefficient}, Delta(1) = 1."""
    g = max(alex)
    return [sum(k * alex.get(j + k, 0) for k in range(1, g + 1)) for j in range(g + 1)]


def symmetric(coeffs):
    """Coefficient list from the lowest degree (normalize) -> symmetric {exponent: coefficient}."""
    n = len(coeffs) - 1
    assert n % 2 == 0 and coeffs == coeffs[::-1] and sum(coeffs) == 1
    return {e - n // 2: c for e, c in enumerate(coeffs) if c}


def _hfk_child(pd, q):
    try:
        import knot_floer_homology as kfh
        h = kfh.pd_to_hfk(pd)
        q.put(('ok', {k: (v if k != 'ranks' else [[a, m, r] for (a, m), r in v.items()]) for k, v in h.items()}))
    except BaseException as ex:
        q.put(('error', repr(ex)))


def hfk_timed(pd, limit):
    """knot_floer_homology.pd_to_hfk in a separate process with a time limit: (status, result)."""
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    q = ctx.Queue()
    p = ctx.Process(target=_hfk_child, args=(pd, q))
    p.start()
    t0 = time.time()
    while time.time() - t0 < limit:
        try:
            res = q.get(timeout=1)
            p.join(5)
            return res
        except Exception:
            if not p.is_alive():
                return ('error', 'exit code %s' % p.exitcode)
    p.kill(); p.join()
    return ('timeout', limit)


def girth(L):
    import knot_floer_homology as kfh
    return kfh.pd_to_morse(L.PD_code())['girth']


def hfk_record(K, limit=HFK_SECONDS):
    """HFK of the knot of diagram K by Szabo's program, which needs a Morse presentation of girth <= MAX_GIRTH.  If K
    has none, K is changed by random Reidemeister moves (backtrack, then simplify('global')), which keep the knot
    and its chirality.  Returns (record or None, the diagram used, notes)."""
    best = (girth(K), K)
    for t in range(GIRTH_TRIES):
        if best[0] <= MAX_GIRTH:
            break
        M = K.copy()
        M.backtrack(steps=20 + 10 * t, prob_type_1=0, prob_type_2=0.5)
        M.simplify('global')
        g = girth(M)
        if g < best[0]:
            best = (g, M)
    g, K = best
    notes = dict(girth=g, crossings=len(K.crossings))
    if g > MAX_GIRTH:
        return None, K, dict(notes, failure='girth %d > %d' % (g, MAX_GIRTH))
    status, h = hfk_timed([list(x) for x in K.PD_code()], limit)
    if status != 'ok':
        return None, K, dict(notes, failure='%s %s' % (status, h))
    ranks = {(a, m): r for a, m, r in h['ranks']}
    tau = h['tau']
    rec = dict(tau=tau, nu=h['nu'], epsilon=h['epsilon'], seifert_genus=h['seifert_genus'], fibered=h['fibered'],
               L_space_knot=h['L_space_knot'], total_rank=h['total_rank'], thin=all(m == a - tau for (a, m) in ranks),
               ranks=sorted([a, m, r] for (a, m), r in ranks.items()), alexander_hfk=alexander_hfk(ranks))
    return rec, K, notes


def staircase(alex):
    """The staircase complex determined by Delta (as for an L-space knot): generators (A, M) and arrows
    (from, to, horizontal length, vertical length); x_0 has A = g, M = 0; x_odd -> x_odd-1 horizontal,
    x_odd -> x_odd+1 vertical."""
    ns = sorted((e for e, c in alex.items() if c), reverse=True)
    if not all(alex[n] == (-1) ** k for k, n in enumerate(ns)):
        return None, None
    M = [0]
    for k in range(1, len(ns)):
        M.append(M[k - 1] + 1 - 2 * (ns[k - 1] - ns[k]) if k % 2 else M[k - 1] - 1)
    arrows = []
    for k in range(1, len(ns), 2):
        arrows += [(k, k - 1, ns[k - 1] - ns[k], 0), (k, k + 1, 0, ns[k] - ns[k + 1])]
    return [(n, m) for n, m in zip(ns, M)], arrows


def dual_complex(gens, arrows):
    """The model complex of the mirror image: gradings negated, arrows reversed."""
    return [(-a, -m) for a, m in gens], [(y, x, a, b) for x, y, a, b in arrows]


def V_model(gens, arrows, s, g):
    """V_s of a model CFK^infinity over F_2: the tower of H_*(A_s^+), A_s^+ = C{max(i, j - s) >= 0}, has bottom
    grading -2 V_s; [x, i, j] has j - i = A(x) and grading M(x) + 2 i (exact linear algebra over F_2)."""
    lo, hi = -2 * g - 6, 2 * g + 8
    def elems(gr):
        return [(x, (gr - M) // 2) for x, (A, M) in enumerate(gens)
                if (gr - M) % 2 == 0 and max((gr - M) // 2, (gr - M) // 2 + A - s) >= 0]
    def boundary(e):
        x, i = e
        out = set()
        for (u, y, a, b) in arrows:
            if u == x and max(i - a, i + gens[x][0] - b - s) >= 0:
                out ^= {(y, i - a)}
        return out
    def in_span(vecs, t):
        red, t = [], set(t)
        for r in vecs:
            r = set(r)
            for p, rr in red:
                if p in r:
                    r ^= rr
            if r:
                red.append((min(r), r))
        for p, rr in red:
            if p in t:
                t ^= rr
        return not t
    def cycles(gr):
        red, basis = [], []
        for e in elems(gr):
            src, im = {e}, set(boundary(e))
            for p, rs, ri in red:
                if p in im:
                    im ^= ri; src ^= rs
            if im:
                red.append((min(im), src, im))
            else:
                basis.append(src)
        return basis
    gr = hi - hi % 2
    z = next(c for c in cycles(gr) if not in_span([boundary(e) for e in elems(gr + 1)], c))
    last = None
    while gr >= lo and not in_span([boundary(e) for e in elems(gr + 1)], z):
        last = gr
        z = {(x, i - 1) for (x, i) in z if max(i - 1, i - 1 + gens[x][0] - s) >= 0}
        gr -= 2
    return -last // 2


def mirror_hfk(h):
    """HFK of the mirror: both gradings negated (HFK_d(-K, a) = HFK_{-d}(K, -a)), tau and epsilon negated."""
    return dict(tau=-h['tau'], epsilon=-h['epsilon'], seifert_genus=h['seifert_genus'], fibered=h['fibered'],
                total_rank=h['total_rank'], thin=h['thin'], ranks=sorted([-a, -m, r] for a, m, r in h['ranks']),
                alexander_hfk=h['alexander_hfk'], note='negated gradings of the HFK of kappa')


def V_from(hk, alex):
    """(V_j for j >= 0, source, check) or (None, reason, None).
    - HFK of staircase shape (rank one at each exponent of Delta, with the Maslov gradings of staircase(Delta)) and
      tau = g > 0: CFK^infinity is that staircase (an L-space knot), V_j = t_j; check: V_model of the staircase.
    - The same after negating both gradings, tau = -g: the mirror of an L-space knot, CFK^infinity is the dual
      staircase, V_j = 0; check: V_model of the dual staircase.
    - Thin (Maslov = Alexander - tau on all of HFK): CFK^infinity is the staircase of T(2, 2|tau| + 1) (its dual if
      tau < 0) plus acyclic boxes (Petkova), V_j = ceil((tau - j)/2) for 0 <= j <= tau if tau >= 0, V_j = 0 if
      tau <= 0; check: V_model of that staircase."""
    tau, g = hk['tau'], hk['seifert_genus']
    ranks = sorted(tuple(r) for r in hk['ranks'])
    if g == 0:
        return [0], 'unknot', True
    gens, arrows = staircase(alex)
    if gens is not None:
        if ranks == sorted((a, m, 1) for a, m in gens) and tau == g:
            V = torsion(alex)
            return V, 'L-space knot', [V_model(gens, arrows, s, g) for s in range(g + 1)] == V
        dg, da = dual_complex(gens, arrows)
        if ranks == sorted((a, m, 1) for a, m in dg) and tau == -g:
            return [0], 'mirror of an L-space knot', [V_model(dg, da, s, g) for s in range(g + 1)] == [0] * (g + 1)
    if hk['thin']:
        t = abs(tau)
        if t == 0:
            return [0], 'thin', True
        sg, sa = staircase({k: (-1) ** (t - k) for k in range(-t, t + 1)})
        if tau < 0:
            sg, sa = dual_complex(sg, sa)
        model = [V_model(sg, sa, s, t) for s in range(t + 1)]
        V = [-((j - tau) // 2) for j in range(tau + 1)] + [0] if tau > 0 else [0]
        return V, 'thin', model == (V[:t + 1] if tau > 0 else [0] * (t + 1))
    return None, 'neither of staircase shape nor thin', None


# ---------- step 6: correction terms

def d_lens(p, q, i):
    """d(L(p,q), i) for L(p,q) = p/q surgery on the unknot (Ozsvath-Szabo recursion, Ni-Wu labels)."""
    if p == 1:
        return Fr(0)
    return Fr(-1, 4) + Fr((2 * i + 1 - p - q) ** 2, 4 * p * q) - d_lens(q, p % q, i % q)


def d_surgery(p, q, V, i):
    """d(S^3_{p/q}(K), i) for p, q > 0 (Ni-Wu), V_j = 0 beyond the list."""
    v = lambda j: V[j] if 0 <= j < len(V) else 0
    return d_lens(p, q, i) - 2 * max(v(i // q), v((p + q - 1 - i) // q))


class Cyclic(Form):
    """H = Z/p generated by mu, D lambda(k mu, l mu) = a k l mod p, d(t_s + k mu) given."""

    def __init__(self, p, a, d):
        self.D, self.mods, self.n = p, (p,), 1
        self.lift = {(k,): [k] for k in range(p)}
        self.adj = [[a % p]]
        self.d = {(k,): d[k] for k in range(p)}


def surgery_form(p, q, V):
    """Y = S^3_{p/q}(K), p odd: lambda(mu, mu) = -q/p; t_s + k mu has the Ni-Wu label i_s + q k, 2 i_s = q - 1 mod p."""
    i_s = (q - 1) * pow(2, -1, p) % p
    return Cyclic(p, -q, [d_surgery(p, q, V, (i_s + q * k) % p) for k in range(p)])


def negate(F):
    """The Form of -Y: lambda and d negated."""
    return Cyclic(F.D, -F.adj[0][0], [-F.d[(k,)] for k in range(F.D)])


def hj(p, q):
    """p/q = a1 - 1/(a2 - 1/...), a_i >= 2."""
    out, x = [], Fr(p, q)
    while True:
        a = -((-x.numerator) // x.denominator)
        out.append(a)
        if x == a:
            return out
        x = 1 / (a - x)


def lens_check(pairs=((5, 2), (7, 2), (9, 2), (11, 2), (63, 2), (7, 3), (13, 5), (25, 7), (63, 4), (175, 2), (223, 2), (27, 2))):
    """The Ni-Wu labelling with mu = e_1 against the sharp linear plumbing of -p/q (bounded by -L(p, q))."""
    out = {}
    for p, q in pairs:
        a = hj(p, q)
        n = len(a)
        Q = [[(-a[i] if i == j else 1 if abs(i - j) == 1 else 0) for j in range(n)] for i in range(n)]
        Fz = Form([[-v for v in row] for row in Q])          # boundary of X_Q = S^3_{-p/q}(U) = -L(p, q); mu = e_1
        Fl = surgery_form(p, q, [])
        e1 = Fz.cm.coords([1] + [0] * (n - 1))
        order = [Fz.mul(k, e1) for k in range(p)]
        assert len(set(order)) == p
        ok = all(Fz.d[order[k]] == -Fl.d[(k,)] for k in range(p)) and Fz.lam(e1, e1) % p == q % p
        out['%d/%d' % (p, q)] = ok
    return out


def jones_at_minus_one(s):
    import sympy
    t = sympy.Symbol('t')
    V = sympy.sympify(s.replace('^', '**'), locals={'t': t})
    return sympy.Rational(V.subs(t, -1)), sympy.Rational(sympy.diff(V, t).subs(t, -1))


def mullins_R(D, sigma, jones):
    """Mullins' R(K) = -D sigma/4 + D V'(-1)/(6 V(-1)) = -2 D lambda_CW(Sigma(K)), from the KnotInfo signature and Jones
    polynomial."""
    v, dv = jones_at_minus_one(jones)
    assert abs(v) == D
    return Fr(-D * sigma, 4) + Fr(int(D * dv.p), int(6 * v.p * dv.q))


def lens_form_check(V_K, D, a):
    """<a/D> and <-a/D> against lambda_K = (V + V^T)^-1 on coker(V + V^T), by Jordan data (exact)."""
    n = len(V_K)
    A = [[V_K[i][j] + V_K[j][i] for j in range(n)] for i in range(n)]
    inv = invariants(A, D)
    def iso(a):
        for p, J in inv.items():
            assert len(J) == 1, 'H_1 is not cyclic'
            v, u = J[0]
            w = a * (D // p ** v) % p
            if legendre(w * u, p) != 1:
                return False
        return True
    return dict(plus=iso(a), minus=iso(-a))


# ---------- the pipeline for one knot

MAX_ROUTES = 3                            # unknotting crossings carried to the end (their d-functions must agree)
MAX_DRILL_CURVES, MAX_DRILL_FORMS = 15, 2  # short geodesics of Sigma(K) tried, and drilled descriptions carried to the end
DIAGRAMS_MINIMAL, DIAGRAMS_NONMINIMAL = 30, 60   # other diagrams of K searched when no route gives V_j
BOUNDS_MAX = 20000                        # largest number of admissible V vectors tested


def regina_names(M, tries=60, enough=5):
    """Names of the closed manifold M (a SnapPy filling) by Regina's recognition of standard triangulations, over
    randomised triangulations (as montesinos_dinv_check.sfs_name).  Two triangulations with the same name are
    homeomorphic, orientation included when the names are not invariant under reflection."""
    import regina
    from collections import Counter
    names = Counter()
    for t in range(tries):
        F = M.filled_triangulation()
        if t:
            F.randomize()
        T = regina.Triangulation3(F._to_string())
        T.intelligentSimplify(); T.intelligentSimplify()
        st = regina.StandardTriangulation.recognise(T)
        if st is not None and st.manifold() is not None:
            names[st.manifold().name()] += 1
        if sum(names.values()) >= enough:
            break
    return dict(names)


def regina_chirality(L, p, q, sigma_names):
    """Which of the diagram L and its mirror gives S^3_{p/q} with a Regina name of Sigma(K).  Equal names describe the
    same manifold, so a match proves Sigma(K) = +-S^3_{p/q}(that diagram) (the names do not see the orientation here);
    that is all the unsigned d-test needs.  A name that does not match is not a proof of non-homeomorphism."""
    out = {}
    for nm, K in (('diagram', L), ('mirror', L.mirror())):
        E = K.exterior(); E.dehn_fill((p, q))
        names = regina_names(E)
        out[nm] = dict(names=names, matches=bool(set(names) & set(sigma_names)))
    return out


def kappa_key(xr):
    """A key for kappa when X is hyperbolic (SnapPy identification, else the volume); None otherwise."""
    if not xr['solution_type'].startswith(('all tetrahedra positively', 'contains negatively')) or xr['volume'] < 0.5:
        return None
    return tuple(sorted(xr['identify'])) or ('volume', round(xr['volume'], 5))


def route(pd, c, info, SigK, R, seen=None, sigma_names=None):
    """Steps 2-6 for the unknotting crossing c: the record, and the Form of Z = -Sigma(K) (or None).  With `seen`
    (keys of kappa already treated), a kappa met before is not treated again."""
    D = info['det']
    rec = dict(crossing=c)
    stage = 'crossing circle'
    try:
        circ, UC, _, ui = circle_checks(pd, c)
        rec['circle'] = circ
        assert circ['ok'], 'crossing circle checks failed'
        stage = 'double branched cover'
        N, (j1, j2) = double_cover(UC, ui)
        rec['cover'] = dict(H1=str(N.homology()), fillings=[list(map(float, f)) for f in N.cusp_info('filling')],
                            open_cusps=[j1, j2])
        pair = sphere_pair(N, j1, j2)
        assert pair, 'no S^3 filling of the two lifts found'
        rec['cover'].update(pair)
        stage = 'kappa exterior'
        X, xr = kappa_exterior(N, j2, pair)
        rec['kappa'] = xr
        assert xr['ok'], 'kappa exterior checks failed'
        key = kappa_key(xr)
        if seen is not None and key is not None:
            if key in seen:
                rec['duplicate_of'] = seen[key]
                return rec, None
            seen[key] = 'crossing %d' % c
        stage = 'slope'
        rec['slope'] = slope_checks(SigK, X, D)
        stage = 'kappa diagram'
        L, kd = kappa_diagram(X)
        rec['kappa']['diagram'] = kd
        stage = 'knot Floer homology'
        alex_fox = alexander_fox([tuple(x) for x in kd['pd']]) if kd['pd'] else [1]
        rec['kappa']['alexander_fox'] = alex_fox
        alex = symmetric(alex_fox)
        hk, L, notes = hfk_record(L, info.get('hfk_seconds', HFK_SECONDS))
        rec['kappa']['hfk'] = hk
        rec['kappa']['hfk_diagram'] = dict(notes, pd=[list(x) for x in L.PD_code()])
        if hk is None:
            rec['flag'] = 'knot Floer homology of kappa not computed: ' + notes['failure']
            return rec, None
        assert hk['alexander_hfk'] == alex_fox, 'Delta(kappa): Fox calculus and HFK differ'
        tk = torsion(alex)
        sum_t = tk[0] + 2 * sum(tk[1:])
        stage = 'sign'
        dl = sum(d_lens(D, 2, i) for i in range(D))
        cw = {str(e): e * (dl - 4 * sum_t) == R for e in (1, -1)}
        rec['casson_walker'] = dict(R=str(R), sum_d_lens=str(dl), sum_torsion=sum_t, holds=cw)
        lf = lens_form_check(info['seifert'], D, 2)
        rec['linking_form'] = dict(plus=lf['plus'], minus=lf['minus'])  # <2/D>, <-2/D> against (V + V^T)^-1
        e, src, why = decide_sign(cw, rec['linking_form'], rec['slope'], D)
        rec['sign_sources'] = src
        assert e is not None, why
        assert rec['linking_form']['plus' if e == 1 else 'minus'], 'lambda_Z is not (V + V^T)^-1'
        rec['e'] = e
        rec['surgery'] = 'Sigma(K) = S^3_{%d/2}(kappa)' % (e * D)
        if kd['orientation'] is None and sigma_names:     # X is not hyperbolic: Regina names decide the chirality
            stage = 'chirality by Regina'
            ch = regina_chirality(L, e * D, 2, sigma_names)
            kd['regina_chirality'] = ch
            if ch['diagram']['matches'] and not ch['mirror']['matches']:
                kd['orientation'] = 'diagram (Regina names)'
            elif ch['mirror']['matches'] and not ch['diagram']['matches']:
                kd['orientation'] = 'mirror (Regina names)'
                L = L.mirror(); hk = mirror_hfk(hk)
                kd['pd'] = [list(x) for x in L.PD_code()]
        stage = 'V_j'
        hh = hk if e == 1 else mirror_hfk(hk)
        rec['kappa_hat'] = dict(which='kappa' if e == 1 else 'mirror of kappa', hfk=hh)
        V, src, check = V_from(hh, alex)
        rec['V'], rec['V_source'], rec['V_model_check'] = V, src, check
        if V is not None:
            assert check, 'V_j against the model complex'
        if V is None:
            rec['flag'] = 'V_j not determined: ' + src
            rec['bounds'] = bounds_test(D, e, hh, kd['orientation'] is not None)
            return rec, None
        FZ_alt = None
        if kd['orientation'] is None:     # X is not hyperbolic: also the other chirality of the diagram
            rec['kappa']['diagram']['note'] = 'chirality from the construction of exterior_to_link only'
            V2, src2, check2 = V_from(mirror_hfk(hh), alex)
            rec['other_chirality'] = dict(V=V2, V_source=src2, V_model_check=check2)
            if V2 is not None:
                F2 = surgery_form(D, 2, V2)
                FZ_alt = negate(F2) if e == 1 else F2
        stage = 'd-invariants'
        FY = surgery_form(D, 2, V)       # Y' = S^3_{D/2}(kappa_hat) = e Sigma(K)
        sd = sum(FY.d.values())
        rec['sum_d'] = str(sd)
        g = hh['seifert_genus']
        rec['L_space'] = src in ('L-space knot', 'unknot') and D >= 2 * (2 * g - 1)   # D/2 >= 2g - 1
        rec['sum_d_equals_eR'] = sd == e * R
        if rec['L_space']:
            assert rec['sum_d_equals_eR'], 'Y\' is an L-space but sum d != e R(K)'
        rec['d_spin'] = str(FY.d[(0,)])
        FZ = negate(FY) if e == 1 else FY  # Z = -e Y'
        rec['d_Z'] = {str(k[0]): str(v) for k, v in FZ.d.items()}
        rec['lambda_Z_mu_mu'] = '%d/%d' % (FZ.adj[0][0], D)
        if kd['orientation'] is None:
            rec['unverified_chirality'] = True
            return rec, (FZ, FZ_alt)
        return rec, FZ
    except AssertionError as ex:
        rec.update(error=str(ex), stage=stage)
    except Exception as ex:
        rec.update(error=repr(ex), stage=stage, traceback=traceback.format_exc()[-2000:])
    return rec, None


def drill_route(S, i, info, R, box=6):
    """Route through the knot kappa' whose exterior is Sigma(K) drilled along the i-th short geodesic."""
    import snappy
    D = info['det']
    rec = dict(kind='drilling', curve=i)
    stage = 'drilling'
    try:
        Dr = S.drill(i, max_segments=DRILL_SEGMENTS)
        X = snappy.Manifold(Dr.filled_triangulation())
        assert X.num_cusps() == 1
        Z1 = X.copy(); Z1.dehn_fill((1, 0))
        assert str(Z1.homology()) == 'Z/%d' % D
        if str(X.homology()) != 'Z':
            rec['flag'] = 'the drilled manifold has H_1 = %s' % X.homology()
            return rec, None
        # an S^3 filling of a hyperbolic X is exceptional, so its slope has length <= 6 (6-theorem)
        try:
            cand = [tuple(t) for t in improved(X).short_slopes(length=6.01)[0]]
        except Exception:
            cand = []
        s1 = None
        for t in cand + [t for t in slope_list(box) if t not in cand]:
            M = X.copy(); M.dehn_fill(t)
            if M.homology().order() == 1 and s3_candidate(M):
                s1 = t
                break
        if s1 is None:
            rec['flag'] = 'no S^3 filling among the slopes of length <= 6 and in the box'
            return rec, None
        rec['s3_slope'] = list(s1)
        a, b = s1
        c, d = next((c, d) for c in range(-abs(b) - 1, abs(b) + 2) for d in range(-abs(a) - 1, abs(a) + 2) if a * d - b * c == 1)
        # (1,0) of the drilled cusp in the basis (s1, t1): solve (1,0) = u s1 + v t1
        u, v = d, -b                                        # inverse of [[a, b], [c, d]] with determinant 1: rows -> (d,-b; -c,a)
        assert (u * a + v * c, u * b + v * d) == (1, 0)
        X.set_peripheral_curves([(a, b), (c, d)], 0)
        im_m, im_l = h1_free_rank_one(X)
        assert abs(im_m) == 1
        x = -im_l * im_m
        X.set_peripheral_curves([(1, 0), (x, 1)], 0)       # longitude (x, 1) in the basis (s1, t1)
        # the drilled meridian u s1 + v t1 = (u - v x) s1 + v (x s1 + t1)
        p, q = u - v * x, v
        if q < 0:
            p, q = -p, -q
        assert abs(p) == D and q > 0 and gcd(p, q) == 1
        Y = X.copy(); Y.dehn_fill((p, q))
        assert str(Y.homology()) == 'Z/%d' % D
        rec['slope_p_q'] = [p, q]
        Xs = X.copy(); Xs.dehn_fill((1, 0))
        rec['kappa'] = dict(H1=str(X.homology()), filling_10_sphere=regina_sphere(Xs), solution_type=X.solution_type(),
                            volume=round(float(X.volume()), 8), longitude_image=h1_free_rank_one(X)[1])
        try:
            rec['kappa']['identify'] = [str(m) for m in improved(X).identify()]
        except Exception:
            rec['kappa']['identify'] = []
        isom, _ = closed_isometric(Y.filled_triangulation(), S)
        rec['filling_isometric_to_sigma'] = isom
        stage = 'kappa diagram'
        L, kd = kappa_diagram(X)
        rec['kappa']['diagram'] = kd
        stage = 'knot Floer homology'
        alex_fox = alexander_fox([tuple(x) for x in kd['pd']]) if kd['pd'] else [1]
        rec['kappa']['alexander_fox'] = alex_fox
        alex = symmetric(alex_fox)
        hk, L, notes = hfk_record(L, info.get('hfk_seconds', HFK_SECONDS))
        rec['kappa']['hfk'] = hk
        rec['kappa']['hfk_diagram'] = dict(notes, pd=[list(x) for x in L.PD_code()])
        if hk is None:
            rec['flag'] = 'knot Floer homology of kappa not computed: ' + notes['failure']
            return rec, None
        assert hk['alexander_hfk'] == alex_fox, 'Delta(kappa): Fox calculus and HFK differ'
        tk = torsion(alex)
        sum_t = tk[0] + 2 * sum(tk[1:])
        e = 1 if p > 0 else -1                              # Sigma(K) = S^3_{p/q}(kappa'); p < 0: Sigma(K) = -S^3_{|p|/q}(mirror)
        P = abs(p)
        dl = sum(d_lens(P, q, i) for i in range(P))
        rec['casson_walker'] = dict(R=str(R), sum_d_lens=str(dl), sum_torsion=sum_t, holds=e * (dl - 2 * q * sum_t) == R)
        assert rec['casson_walker']['holds'], 'Casson-Walker check of the drilled description fails'
        lf = lens_form_check(info['seifert'], D, e * q)     # lambda_Z(mu, mu) = e q / D
        rec['linking_form'] = lf
        assert lf['plus'], 'lambda_Z is not (V + V^T)^-1'
        rec['e'] = e
        rec['surgery'] = 'Sigma(K) = S^3_{%d/%d}(kappa\')' % (p, q)
        if kd['orientation'] is None:
            rec['flag'] = 'the chirality of the diagram of kappa\' is not checked by an isometry'
            return rec, None
        hh = hk if e == 1 else mirror_hfk(hk)
        rec['kappa_hat'] = dict(which='kappa\'' if e == 1 else 'mirror of kappa\'', hfk=hh)
        V, src, check = V_from(hh, alex)
        rec['V'], rec['V_source'], rec['V_model_check'] = V, src, check
        if V is None:
            rec['flag'] = 'V_j not determined: ' + src
            return rec, None
        assert check, 'V_j against the model complex'
        FY = surgery_form(P, q, V)
        rec['sum_d'] = str(sum(FY.d.values()))
        rec['d_spin'] = str(FY.d[(0,)])
        FZ = negate(FY) if e == 1 else FY
        rec['d_Z'] = {str(k[0]): str(v) for k, v in FZ.d.items()}
        rec['lambda_Z_mu_mu'] = '%d/%d' % (FZ.adj[0][0], D)
        return rec, FZ
    except AssertionError as ex:
        rec.update(error=str(ex), stage=stage)
    except Exception as ex:
        rec.update(error=repr(ex), stage=stage, traceback=traceback.format_exc()[-2000:])
    return rec, None


def admissible_count(g, lower):
    """The number of vectors admissible_V(g, lower) would list (dynamic programming over j)."""
    cap = lambda j: -((j - g) // 2)
    ways = {v: 1 for v in range(0, cap(0) + 1) if (0 >= lower or v >= 1)}
    for j in range(1, g + 1):
        new = {}
        for v, w in ways.items():
            for u in (v, v - 1):
                if 0 <= u <= cap(j) and (j >= lower or u >= 1):
                    new[u] = new.get(u, 0) + w
        ways = new
    return ways.get(0, 0)


def admissible_V(g, lower):
    """Every V = (V_0, ..., V_g) with V_j >= V_{j+1} >= V_j - 1 (Ni-Wu, Lemma 2.4), 0 <= V_j <= ceil((g - j)/2)
    (the sharper bound V_j <= ceil((g_4 - j)/2) with g_4 <= g, the Seifert genus from HFK; see the module docstring),
    V_g = 0, and V_j >= 1 for
    j < lower (lower = tau or nu: tau <= nu <= nu^+ = min{j : V_j = 0}, Hom-Wu)."""
    out = []
    def grow(V):
        j = len(V)
        if j == g + 1:
            if V[-1] == 0:
                out.append(list(V))
            return
        for v in ((V[-1], V[-1] - 1) if V else range(0, -(-g // 2) + 1)):
            if 0 <= v <= -((j - g) // 2) and (j >= lower or v >= 1):
                grow(V + [v])
    grow([])
    return out


def bounds_test(D, e, hh, chirality_checked):
    """The d-test for every V allowed by the bounds when V_j is not determined.  Variants: 'none' (monotonicity and
    the sharper upper bound only; independent of the chirality of kappa), 'tau' and 'nu' (also V_j >= 1 for j < tau or nu,
    which needs the chirality of kappa_hat)."""
    g = hh['seifert_genus']
    variants = [('none', 0)]
    if chirality_checked:
        variants.append(('tau', max(hh['tau'], 0)))
        if hh.get('nu') is not None:
            variants.append(('nu', max(hh['nu'], 0)))
    out = {}
    for label, lower in variants:
        n = admissible_count(g, lower)
        if n > BOUNDS_MAX:
            out[label] = dict(lower=lower, admissible=n, skipped=True)
            continue
        obs, passing = 0, None
        for V in admissible_V(g, lower):
            FY = surgery_form(D, 2, V)
            FZ = negate(FY) if e == 1 else FY
            t = test(FZ, (1, -1))
            if t is None:
                obs += 1
            elif passing is None:
                passing = dict(V=V, passing=t)
        out[label] = dict(lower=lower, admissible=n, obstructed=obs, all_obstructed=n > 0 and obs == n, first_passing=passing)
    return out


def flag_text(rec):
    """The reasons of the routes that did not give a d-function (routes met before are left out)."""
    rs = [r for r in rec.get('routes', []) + rec.get('drill_routes', []) if not r.get('duplicate_of')]
    return '; '.join(sorted({str(r.get('flag') or r.get('error')) for r in rs if r.get('flag') or r.get('error')})) or 'no route'


def bounds_verdict(rec):
    """For a knot whose V_j is not determined (verdict flagged): 'obstructed_by_bounds' if every V allowed by the
    bounds is obstructed (weakest variant first: none, tau, nu), 'passes_for_all_V' if every V allowed by the
    chirality-free bounds passes.  For routes whose chirality is unchecked, the tau-bounds under the chirality of the
    construction are recorded (information only)."""
    routes = rec.get('routes', [])
    for label in ('none', 'tau', 'nu'):
        hit = [r for r in routes if r.get('bounds', {}).get(label, {}).get('all_obstructed')]
        if hit:
            rec.update(verdict='obstructed_by_bounds', bounds_variant=label,
                       bounds_route='crossing %d' % hit[0]['crossing'] + (' of diagram %d' % hit[0]['diagram_index']
                                                                          if 'diagram_index' in hit[0] else ''))
            return rec
    free = [r for r in routes if r.get('bounds', {}).get('none', {}).get('admissible') and r['bounds']['none'].get('obstructed') == 0]
    if free:
        rec.update(verdict='passes_for_all_V', bounds_variant='none',
                   bounds_route='crossing %d' % free[0]['crossing'] + (' of diagram %d' % free[0]['diagram_index']
                                                                       if 'diagram_index' in free[0] else ''))
        return rec
    info = {}
    for r in routes:
        if r.get('bounds') and 'tau' not in r['bounds'] and r.get('kappa_hat'):
            key = 'crossing %d' % r['crossing'] + (' of diagram %d' % r['diagram_index'] if 'diagram_index' in r else '')
            b = bounds_test(rec['det'], r['e'], r['kappa_hat']['hfk'], True).get('tau')
            info[key] = b
    if info:
        rec['bounds_construction_chirality'] = info
    return rec


def same_form(F1, F2):
    """Units u with d_2(u k) = d_1(k) and lambda_2(u mu, u mu) = lambda_1(mu, mu) (cyclic Forms)."""
    D = F1.D
    return [u for u in range(1, D) if gcd(u, D) == 1 and (u * u * F2.adj[0][0] - F1.adj[0][0]) % D == 0
            and all(F2.d[((u * k) % D,)] == F1.d[(k,)] for k in range(D))]


def other_diagrams(pd0, seed, minimal=DIAGRAMS_MINIMAL, nonminimal=DIAGRAMS_NONMINIMAL):
    """Diagrams of K other than the KnotInfo one (spherogram many_diagrams, and random Reidemeister II/III moves from
    them, as in smooth_top_gap.diagrams_of), each checked isometric to the KnotInfo exterior (answer True)."""
    import spherogram
    from smooth_top_gap import gauss_key
    E0 = spherogram.Link([list(x) for x in pd0]).exterior()
    keys, out = {gauss_key([tuple(x) for x in pd0])}, []
    random.seed(seed)
    found = spherogram.Link([list(x) for x in pd0]).many_diagrams(target=minimal, tries=150, method='backtrack')
    c = min(len(D.crossings) for D in found)
    for i in range(nonminimal):
        random.seed(seed + 1 + i)
        D = found[i % len(found)].copy()
        D.backtrack(steps=1 + i % 4, prob_type_1=0, prob_type_2=0.5)
        if c < len(D.crossings) <= c + 4:
            found.append(D)
    for D in found:
        P = [list(x) for x in D.PD_code()]
        k = gauss_key([tuple(x) for x in P])
        if k in keys:
            continue
        keys.add(k)
        if isometric(spherogram.Link(P).exterior(), E0) is True:
            out.append(P)
    return out


def plumbing_check(name, names=None):
    """Independent d-invariants when Regina recognises Sigma(K) as a Seifert fibred space over S^2 (the code of
    montesinos_dinv_check.py: negative-definite star plumbing, exact closest-vector enumeration, its own rank-one
    test).  Returns the record and the Form of the plumbed boundary Y_P (= +-Sigma(K)), or (record, None)."""
    import snappy, regina
    import montesinos_dinv_check as MC
    if names is None:
        snappy.set_rand_seed(0); regina.RandomEngine.reseedWithDefault()
        names = dict(MC.sfs_name(name))
    rec = dict(regina=dict(names))
    fib = MC.parse(max(names, key=names.get)) if names else None
    if not fib:
        return rec, None
    Q, _ = MC.plumbing(fib)
    try:
        N, u, d = MC.d_invariants(Q)
    except MC.NotCyclic:
        rec['note'] = 'H_1 not cyclic'
        return rec, None
    ob, passing = MC.obstructed(N, u, d)
    rec.update(fibres=fib, plumbing_rank=len(Q), det=N, obstructed=ob, passing=[list(x) for x in passing[:3]])
    return rec, Cyclic(N, u, [d[k] for k in range(N)])


def verdicts(FZ, s):
    """The first passing configuration of the unsigned test and of the signed test (s = 2: Z0 = Z only, s = 6: Z0 = -Z
    only), or None."""
    passing = test(FZ, (1, -1))
    allowed = {2: (1,), 6: (-1,)}.get(s, (1, -1))
    signed = passing if passing and passing['sign'] in allowed else test(FZ, allowed)
    return passing, signed


def one(item):
    """item = (name, pd, info); info has det, sigma, jones, seifert, s, and for calibration the Goeritz matrix."""
    name, pd, info = item
    t0 = time.time()
    import snappy
    pari().allocatemem(PARI_STACK, PARI_SIZEMAX, silent=True)
    snappy.set_rand_seed(zlib_seed(name)); random.seed(zlib_seed(name))
    D = info['det']
    rec = dict(name=name, det=D, sigma=info['sigma'], s=info['s'], arf=info['arf'])
    try:
        uc = unknotting_crossings(pd)
        rec['unknotting_crossings'] = uc
        SigK = sigma_cover(pd)
        assert str(SigK.homology()) == 'Z/%d' % D, 'H_1(Sigma(K)) = %s' % SigK.homology()
        R = mullins_R(D, info['sigma'], info['jones'])
        rec['R'] = str(R)
        import montesinos_dinv_check as MC  # Regina names, the plumbing and the second implementation of the test
        snappy.set_rand_seed(0)
        sigma_names = dict(MC.sfs_name(name))           # Regina names of Sigma(K) (from the KnotInfo PD code)
        rec['sigma_regina_names'] = sigma_names
        routes, forms, seen = [], [], {}
        for u in uc:                      # the Montesinos trick at each unknotting crossing
            if len(forms) == MAX_ROUTES or time.time() - t0 > KNOT_SECONDS / 3:
                break
            r, FZ = route(pd, u['crossing'], info, SigK, R, seen, sigma_names)
            routes.append(r)
            if FZ is not None:
                forms.append(('crossing %d' % u['crossing'], FZ))
        if not forms:                     # other diagrams of K: other unknotting crossings, possibly other kappa
            searched = []
            for di, P in enumerate(other_diagrams(pd, zlib_seed(name))):
                if forms or time.time() - t0 > KNOT_SECONDS / 3:
                    break
                ucs = unknotting_crossings(P)
                searched.append(dict(diagram=di, crossings=len(P), unknotting_crossings=[u['crossing'] for u in ucs]))
                for u in ucs:
                    r, FZ = route(P, u['crossing'], info, SigK, R, seen, sigma_names)
                    if 'duplicate_of' in r:
                        continue
                    r['diagram'], r['diagram_index'] = P, di
                    routes.append(r)
                    if FZ is not None:
                        forms.append(('diagram %d crossing %d' % (di, u['crossing']), FZ))
                        break
            rec['diagram_search'] = dict(diagrams=len(searched), searched=searched, kappas_seen=len(seen))
        # routes whose diagram of kappa has an unchecked chirality: (FZ, FZ of the other chirality)
        unverified = {k: F for k, F in forms if isinstance(F, tuple)}
        forms = [(k, F) for k, F in forms if not isinstance(F, tuple)] + [(k, F) for k, F in unverified.items()]
        rec['routes'] = routes
        drills = []                       # surgery descriptions from short geodesics of Sigma(K)
        S = improved(SigK, 100)
        if S.solution_type() in ('all tetrahedra positively oriented', 'contains negatively oriented tetrahedra') \
                and time.time() - t0 < KNOT_SECONDS / 2:
            try:
                curves = S.dual_curves(max_segments=DRILL_SEGMENTS)
            except Exception:
                curves = []
            nd = 0
            for i in range(min(len(curves), MAX_DRILL_CURVES)):
                if nd == MAX_DRILL_FORMS or time.time() - t0 > KNOT_SECONDS / 2:
                    break
                r, FZ = drill_route(S, i, info, R)
                drills.append(r)
                if FZ is not None:
                    forms.append(('geodesic %d' % i, FZ)); nd += 1
        rec['drill_routes'] = drills
        pr, FP = plumbing_check(name, sigma_names)
        rec['plumbing'] = pr
        if not forms and FP is None:
            rec.update(verdict='flagged', flag=flag_text(rec))
            bounds_verdict(rec)
            return name, rec
        checked = [(k, F) for k, F in forms if not isinstance(F, tuple)]
        if FP is not None:                # the plumbing decides the chirality of an unchecked diagram
            for k, (F1, F2) in list(unverified.items()):
                ok1 = bool(same_form(F1, FP) or same_form(F1, negate(FP)))
                ok2 = F2 is not None and bool(same_form(F2, FP) or same_form(F2, negate(FP)))
                pr.setdefault('chirality', {})[k] = dict(construction=ok1, other=ok2)
                assert ok1 or ok2, 'the plumbing d-function differs from both chiralities of route %s' % k
                if ok1 and not ok2:
                    checked.append((k, F1))
                elif ok2 and not ok1:
                    checked.append((k, F2))
        rec['unverified_chirality_routes'] = sorted(unverified)
        if not checked and unverified and FP is None:
            # no check of the chirality: the verdict must be the same for both chiralities of every such route
            out = {}
            for k, (F1, F2) in unverified.items():
                v1 = verdicts(F1, info['s'])[0] is None
                v2 = verdicts(F2, info['s'])[0] is None if F2 is not None else None
                out[k] = dict(construction_obstructed=v1, other_obstructed=v2)
            rec['chirality_robustness'] = out
            vals = {x for o in out.values() for x in (o['construction_obstructed'], o['other_obstructed'])}
            if len(vals) != 1 or None in vals:
                rec.update(verdict='flagged', flag='chirality of kappa unchecked and the verdict depends on it')
                return name, rec
            rec.update(verdict='obstructed' if vals.pop() else 'passes', route_used=sorted(unverified)[0],
                       note='verdict robust to the unchecked chirality of kappa', passing_signed=None, obstructed_signed=None)
            k0 = sorted(unverified)[0]
            rec['passing'] = verdicts(unverified[k0][0], info['s'])[0]
            rec['obstructed'] = rec['passing'] is None
            F0 = unverified[k0][0]
            rec['obstructed_independent_code'] = MC.obstructed(F0.D, F0.adj[0][0], {k: F0.d[(k,)] for k in range(F0.D)})[0]
            assert rec['obstructed_independent_code'] == rec['obstructed'], 'the two implementations of the d-test differ'
            if 'goeritz' in info:         # calibration: the Goeritz d-function decides the chirality
                rec['calibration'] = dict(construction=calibrate(info['goeritz'], unverified[k0][0]),
                                          other=calibrate(info['goeritz'], unverified[k0][1]) if unverified[k0][1] is not None else None)
                rec['calibration']['agrees'] = rec['calibration']['construction']['agrees']
            return name, rec
        forms = checked
        if forms:
            src, FZ = forms[0]
            rec['route_used'] = src
            rec['routes_agree'] = {k: bool(same_form(FZ, F)) for k, F in forms[1:]}
            assert all(rec['routes_agree'].values()), 'the d-functions of two routes differ'
            if FP is not None:
                pr['agrees_with_route'] = dict(plus=bool(same_form(FZ, FP)), minus=bool(same_form(FZ, negate(FP))))
                assert pr['agrees_with_route']['plus'] or pr['agrees_with_route']['minus'], 'plumbing d differs'
            def label(r):
                if r.get('kind') == 'drilling':
                    return 'geodesic %d' % r['curve']
                return ('diagram %d crossing %d' % (r['diagram_index'], r['crossing'])) if 'diagram_index' in r else 'crossing %d' % r['crossing']
            r0 = next(r for r in routes + drills if label(r) == src)
            rec.update(e=r0['e'], surgery=r0['surgery'], V=r0['V'], V_source=r0['V_source'], d_Z=r0['d_Z'],
                       lambda_Z_mu_mu=r0['lambda_Z_mu_mu'])
        else:                             # only the plumbing: Y_P = +-Sigma(K); the unsigned test needs no orientation
            src, FZ = 'plumbing', FP
            rec['route_used'] = src
            rec['d_plumbing'] = {str(k[0]): str(v) for k, v in FP.d.items()}
        rec['passing'], rec['passing_signed'] = verdicts(FZ, info['s'])
        # Regina names identify Sigma(K) with +-S^3_{eD/2}(diagram) without the orientation: the unsigned test is decided,
        # the signed one is not; the plumbing alone does not fix the orientation either
        no_orientation = src == 'plumbing' or 'Regina' in str((r0 if forms else {}).get('kappa', {}).get('diagram', {}).get('orientation'))
        if no_orientation:
            rec['passing_signed'] = None
        rec['obstructed'] = rec['passing'] is None
        rec['obstructed_signed'] = None if no_orientation else rec['passing_signed'] is None
        ob2, _ = MC.obstructed(FZ.D, FZ.adj[0][0], {k: FZ.d[(k,)] for k in range(FZ.D)})
        rec['obstructed_independent_code'] = ob2
        assert ob2 == rec['obstructed'], 'the two implementations of the d-test differ'
        if FP is not None:
            assert pr['obstructed'] == rec['obstructed'], 'plumbing verdict differs'
        rec['verdict'] = 'obstructed' if rec['obstructed'] else 'passes'
        if 'goeritz' in info:
            rec['calibration'] = calibrate(info['goeritz'], FZ)
    except AssertionError as ex:
        rec.update(verdict='error', error=str(ex))
    except Exception as ex:
        rec.update(verdict='error', error=repr(ex), traceback=traceback.format_exc()[-2000:])
    finally:
        rec['seconds'] = round(time.time() - t0, 1)
    return name, rec


def zlib_seed(name):
    """The seed of the random generators for one knot."""
    import zlib
    return zlib.crc32(name.encode()) % 2 ** 31


def calibrate(G, FZ):
    """The Goeritz Form G (rank1_dinv, Z = boundary of X_G) against the surgery Form FZ: units u with
    d_G(k g) = d_Z(u k) for all k and lambda_G(g, g) = u^2 lambda_Z(mu, mu), for a generator g of H_G."""
    F = Form(G)
    D = F.D
    assert D == FZ.D
    zero = tuple(0 for _ in F.mods)
    g = next(x for x in F.lift if all(F.mul(D // p, x) != zero for p in primes_of(D)))
    pts = [zero]
    for k in range(1, D):
        pts.append(F.add(pts[-1], g))
    units = [u for u in range(1, D) if gcd(u, D) == 1
             and F.lam(g, g) % D == (u * u * FZ.adj[0][0]) % D
             and all(F.d[pts[k]] == FZ.d[((u * k) % D,)] for k in range(D))]
    return dict(units=units, agrees=bool(units), d_spin_goeritz=str(F.d[zero]),
                sum_d_goeritz=str(sum(F.d.values())))


# ---------- running knots in separate processes, with a time limit and checkpoints

def _child(item, q):
    try:
        q.put(one(item))
    except BaseException as ex:            # never leave the parent waiting
        q.put((item[0], dict(name=item[0], verdict='error', stage='process', error=repr(ex))))


def run_items(items, workers, partial, limit=2 * KNOT_SECONDS, log=print):
    """one() for every item, each in its own process; results are checkpointed to `partial` (JSON) as they come."""
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    out = json.load(open(partial)) if os.path.exists(partial) else {}
    todo = [it for it in items if it[0] not in out]
    if out:
        log('resuming: %d done, %d to do' % (len(out), len(todo)))
    running = {}
    while todo or running:
        while todo and len(running) < workers:
            it = todo.pop(0)
            q = ctx.Queue()
            p = ctx.Process(target=_child, args=(it, q))
            p.start()
            running[it[0]] = (p, q, time.time())
        time.sleep(0.5)
        for name in list(running):
            p, q, t0 = running[name]
            res = None
            try:
                res = q.get_nowait()
            except Exception:
                pass
            if res is not None:
                p.join(5)
                out[res[0]] = res[1]
            elif not p.is_alive():
                out[name] = dict(name=name, verdict='error', stage='process', error='exit code %s' % p.exitcode)
            elif time.time() - t0 > limit:
                p.kill(); p.join()
                out[name] = dict(name=name, verdict='flagged', flag='time limit %d s' % limit, seconds=limit)
            else:
                continue
            del running[name]
            r = out[name]
            log('%s: %s %s | V %s | e %s | passing %s | routes %s | %ss' % (
                name, r.get('verdict'), r.get('flag') or r.get('error') or '', r.get('V_source'), r.get('e'), r.get('passing'),
                len(r.get('routes', [])), r.get('seconds')), flush=True)
            json.dump(out, open(partial + '.tmp', 'w'))
            os.replace(partial + '.tmp', partial)
    return out


# ---------- inputs

def knot_info(BY, name, s=None):
    """KnotInfo data of K: det, the signature (KnotInfo, and of V + V^T), Arf from det mod 8, s = sigma + 4 Arf mod 8
    (s if given, as in data/rank1_targets.json), the Jones polynomial and the Seifert matrix V."""
    from rank1_linking import signature
    r = BY[name]
    V = json.loads(r['seifert_matrix'])
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    D = int(r['determinant'])
    sig = signature(A)
    arf = 0 if D % 8 in (1, 7) else 1
    info = dict(det=D, sigma=int(r['signature']), sigma_seifert=sig, arf=arf, s_seifert=(sig + 4 * arf) % 8,
                jones=r['jones_polynomial'], seifert=V, s=s if s is not None else (sig + 4 * arf) % 8)
    return info


def convention_check(BY):
    """KnotInfo 3_1 (sigma = -2) is the right-handed trefoil for HFK (tau = 1), and SnapPy's (5, 1) filling of its
    exterior is a lens space: SnapPy's (p, q) filling is p/q surgery."""
    import spherogram
    T = spherogram.Link(json.loads(BY['3_1']['pd_notation']))
    out = dict(sigma=int(BY['3_1']['signature']), tau=T.knot_floer_homology()['tau'])
    for s in ((5, 1), (-5, 1)):
        E = T.exterior(); E.dehn_fill(s)
        out['filling %s lens space' % (s,)] = E.fundamental_group().num_generators() == 1
    out['ok'] = out['sigma'] == -2 and out['tau'] == 1 and out['filling (5, 1) lens space'] and not out['filling (-5, 1) lens space']
    return out


def versions():
    import snappy, spherogram, cypari, database_knotinfo, regina, sympy
    import knot_floer_homology as kfh
    return dict(snappy=snappy.__version__, spherogram=spherogram.__version__, database_knotinfo=database_knotinfo.version(),
                cypari=cypari.__version__, pari='.'.join(map(str, pari().version())), regina=regina.versionString(),
                sympy=sympy.__version__, knot_floer_homology=getattr(kfh, '__version__', '1.2.2'), python=sys.version.split()[0])


def targets():
    """The 50 knots of [C] Proposition 3.2(b): those of the 316 knots of Proposition 3.2 (the targets of
    data/plumbing_dinv.json, which are neither alternating nor Montesinos) whose KnotInfo unknotting number is 1."""
    P = json.load(open(os.path.join(ROOT, 'data', 'plumbing_dinv.json')))['targets']
    T = [k for k, v in P.items() if str(v['u']) == '1']
    assert len(T) == 50
    return T


# Calibration ([C] Section 3.5): alternating knots with KnotInfo u = 1, not Montesinos, whose Sigma(K) is hyperbolic.
# For each crossing number 9-13, the first knots in KnotInfo order with distinct determinants and a SnapPy solution
# with positively oriented tetrahedra (within 10 randomisations), and further knots whose determinant has a square
# factor.  The outcome of the SnapPy test varies between runs, so the list is fixed here.
CALIBRATION = ['9_34', '10_102', '10_104', '10_114', '10_118', '10_119', '11a_73', '11a_101', '11a_114', '11a_136',
               '11a_168', '11a_189', '11a_209', '11a_283', '11a_300', '11a_326', '12a_209', '12a_256', '12a_261',
               '12a_287', '12a_301', '12a_309', '12a_310', '12a_340', '12a_358', '12a_388', '12a_451', '12a_799',
               '13a_713', '13a_773', '13a_885', '13a_1085']


def calibration_knots(BY):
    for n in CALIBRATION:
        r = BY[n]
        assert r['alternating'] == 'Y' and r['unknotting_number'] == '1' and not r['montesinos_notation'].startswith('K(')
    return list(CALIBRATION)


def goeritz_P(pd):
    """P = -G for the negative-definite Goeritz matrix G of an alternating diagram (rank1_goeritz.goeritz_pair)."""
    from owens_obstruction import positive_definite
    Gs = goeritz_pair(pd)
    neg = [G for G in Gs if positive_definite([[-v for v in row] for row in G])]
    assert len(neg) == 1
    return [[-v for v in row] for row in neg[0]]


# ---------- step 9: locally flat Mobius bands for the obstructed knots (smooth_top_gap.py)

def _gap_init(ki, dt):
    import smooth_top_gap as stg
    stg.init_worker(ki, dt)
    pari().allocatemem(PARI_STACK, PARI_SIZEMAX, silent=True)


def _gap_one(item):
    import smooth_top_gap as stg
    try:
        return stg.one(item)
    except BaseException as ex:            # a PARI error must not kill the pool
        return item[0], dict(error=repr(ex), traceback=traceback.format_exc()[-2000:])


def gap_records(names, BY, workers, partial, log=print):
    """smooth_top_gap.one (face bands to a knot with Delta = 1) for the given knots; returns the full records."""
    import smooth_top_gap as stg
    from multiprocessing import get_context
    LL = list(BY.values())
    ki = {r['name']: {k: r[k] for k in stg.KI_FIELDS} for r in LL}
    dt = {'K' + r['dt_name'].replace('_', ''): r['name'] for r in LL if r.get('dt_name')}
    R1 = json.load(open(os.path.join(ROOT, 'data', 'rank1_targets.json')))
    LG = json.load(open(os.path.join(ROOT, 'data', 'linking_genus.json')))
    out = json.load(open(partial)) if os.path.exists(partial) else {}
    items = []
    for k in names:
        if k in out:
            continue
        r = R1[k]
        info = dict(det=r['det'], test='surgery', h_witt=LG[k]['h_witt'], survivors=r['survivors'], s=r['s'],
                    jordan=r['jordan'], rank_one_full=any(x[0] == r['det'] for x in r['survivors']))
        items.append((k, info, json.loads(BY[k]['pd_notation']), None, []))
    with get_context('spawn').Pool(workers, initializer=_gap_init, initargs=(ki, dt)) as p:
        for name, rec in p.imap_unordered(_gap_one, items):
            out[name] = rec
            log('%s: certificate %s | J %s | diagrams %s | det1 %s | hits %s | contradictions %s | %ss %s' % (
                name, rec.get('certificate'), rec.get('J_found'), len(rec.get('diagrams', [])), rec.get('det1_knots'),
                len(rec.get('hits', [])), rec.get('contradictions'), rec.get('seconds'), rec.get('error', '')), flush=True)
            json.dump(out, open(partial + '.tmp', 'w'), default=str)
            os.replace(partial + '.tmp', partial)
    return out


def gap_summary(rec):
    """The part of a smooth_top_gap record kept in data/surgery_dinv.json (the diagrams go to the gap file)."""
    keep = ('det', 'h_witt', 'rank_one_full', 'arcs', 'bands_tried', 'bands_split', 'bands_outside_range', 'knot_results',
            'det1_knots', 'J_found', 'certificate', 'representative', 'contradictions', 'seconds', 'fox_det_checks', 'error')
    out = {k: rec.get(k) for k in keep if k in rec}
    out['diagrams'] = len(rec.get('diagrams', []))
    out['diagrams_not_isometric'] = sum(d.get('isometric_to_knotinfo') is not True for d in rec.get('diagrams', []))
    out['hits'] = len(rec.get('hits', []))
    out['hits_failing_a_recheck'] = sum(not h.get('passes') for h in rec.get('hits', []))
    out['fit_failures'] = len(rec.get('fit_failures', []))
    return out


def census_recheck(rec):
    """For a decided knot whose kappa (the route used) is a SnapPy census knot: a second diagram from the census
    triangulation, independent of the cover; its HFK must be that of our diagram when the isometry from our exterior to
    the census manifold preserves orientation, and its mirror when it reverses it."""
    import snappy, spherogram
    r = next((r for r in rec.get('routes', []) if (r.get('V') is not None or r.get('bounds')) and r.get('kind') != 'drilling'
              and r.get('kappa', {}).get('hfk')), None)
    if r is None:
        return None
    names = [n.split('(')[0] for n in r['kappa'].get('identify', [])]
    cen = next((n for n in names if re.fullmatch(r'K\d+_\d+', n)), None)
    if cen is None:
        return None
    M = snappy.Manifold(cen)
    L0 = spherogram.Link(r['kappa']['diagram']['pd'])
    iso = improved(L0.exterior()).is_isometric_to(improved(M), return_isometries=True)
    dets = sorted({int(h.cusp_maps()[0].det()) for h in iso})
    snappy.set_rand_seed(1)
    L1 = M.exterior_to_link(check_answer=False)
    h1, _, notes = hfk_record(L1)
    if h1 is None:
        return dict(census=cen, failure=notes.get('failure'))
    h0 = r['kappa']['hfk']                      # HFK of the diagram before a mirror chosen by Regina names
    if str(r['kappa']['diagram'].get('orientation')).startswith('mirror (Regina'):
        h0 = mirror_hfk(h0)
    same = sorted(map(tuple, h1['ranks'])) == sorted(map(tuple, h0['ranks']))
    mirror = sorted(map(tuple, h1['ranks'])) == sorted((-a, -m, x) for a, m, x in h0['ranks'])
    ok = (dets == [1] and same) or (dets == [-1] and mirror) or (dets == [-1, 1] and same and mirror)
    return dict(census=cen, cusp_map_dets=dets, hfk_equal=same, hfk_mirror=mirror, consistent=ok)


def bounds_recheck(rec):
    """For a knot obstructed by the bounds: every admissible V again, with the independent implementation of the test
    (montesinos_dinv_check.obstructed); all must be obstructed."""
    import montesinos_dinv_check as MC
    lab = rec['bounds_variant']
    out = []
    for r in rec['routes']:
        b = r.get('bounds', {}).get(lab)
        if not b or not b.get('all_obstructed'):
            continue
        hh, D, e = r['kappa_hat']['hfk'], rec['det'], r['e']
        Vs = admissible_V(hh['seifert_genus'], b['lower'])
        ok = 0
        for V in Vs:
            FY = surgery_form(D, 2, V)
            FZ = negate(FY) if e == 1 else FY
            ob, _ = MC.obstructed(FZ.D, FZ.adj[0][0], {k: FZ.d[(k,)] for k in range(FZ.D)})
            ok += bool(ob)
        out.append(dict(route=r['crossing'], diagram=r.get('diagram_index'), admissible=len(Vs), obstructed_independent=ok))
    return out


# Torus knots of KnotInfo (positive: tau = g) for the connected-sum route
TORUS = {'3_1': 'T(2,3)', '5_1': 'T(2,5)', '7_1': 'T(2,7)', '8_19': 'T(3,4)', '9_1': 'T(2,9)', '10_124': 'T(3,5)'}
_TORUS_DATA = {}


def tensor(c1, c2):
    """The tensor product of two model complexes (CFK^infinity of a connected sum, Ozsvath-Szabo)."""
    g1, a1 = c1; g2, a2 = c2
    n2 = len(g2)
    gens = [(A1 + A2, M1 + M2) for (A1, M1) in g1 for (A2, M2) in g2]
    idx = lambda i, j: i * n2 + j
    arrows = [(idx(x, j), idx(y, j), a, b) for (x, y, a, b) in a1 for j in range(n2)]
    arrows += [(idx(i, x), idx(i, y), a, b) for (x, y, a, b) in a2 for i in range(len(g1))]
    return gens, arrows


def torus_sum_route(rec, BY):
    """For a knot whose V_j is not determined: if HFK of kappa_hat equals HFK of a connected sum A # B of two torus
    knots (either chirality), and Regina gives Y = S^3_{D/2}(A # B) a name of Sigma(K), then Sigma(K) = +-Y (equal
    names describe the same manifold; the names here do not see the orientation), and V_j(A # B) is computed exactly
    from the tensor product of the two (dual) staircases.  The unsigned d-test does not depend on the orientation, so it
    is decided; the signed test needs the orientation and is recorded only when exactly one of +-Y has the name."""
    import spherogram
    from itertools import combinations_with_replacement
    sig = rec.get('sigma_regina_names') or {}
    if not sig:
        return None
    for r in rec.get('routes', []):
        hh = r.get('kappa_hat', {}).get('hfk')
        if not hh or r.get('V') is not None:
            continue
        ranks = sorted(tuple(x) for x in hh['ranks'])
        if not _TORUS_DATA:
            for nm in TORUS:
                K = spherogram.Link(json.loads(BY[nm]['pd_notation']))
                h0 = hfk_record(K)[0]
                assert h0['tau'] == h0['seifert_genus'] > 0               # positive torus knot
                st = staircase(symmetric(alexander_hfk({(x, m): c for x, m, c in h0['ranks']})))
                assert sorted(tuple(x) for x in h0['ranks']) == sorted((x, m, 1) for x, m in st[0])
                _TORUS_DATA[nm] = (K, st, h0['seifert_genus'])
        for a, b in combinations_with_replacement(sorted(TORUS), 2):
            for sa in (1, -1):
                for sb in (1, -1):
                    pieces, cplx, genus, tau = [], [], 0, 0
                    for nm, sg in ((a, sa), (b, sb)):
                        K, st, g0 = _TORUS_DATA[nm]
                        pieces.append(K if sg == 1 else K.mirror())
                        cplx.append(st if sg == 1 else dual_complex(*st))
                        genus += g0; tau += sg * g0
                    if genus != hh['seifert_genus'] or tau != hh['tau']:
                        continue
                    L = pieces[0].connected_sum(pieces[1])
                    hl = hfk_record(L)[0]
                    if hl is None or sorted(tuple(x) for x in hl['ranks']) != ranks:
                        continue
                    D, e = rec['det'], r['e']
                    E = L.exterior(); E.dehn_fill((D, 2))          # Y = S^3_{D/2}(A # B)
                    names = regina_names(E)
                    Em = L.mirror().exterior(); Em.dehn_fill((-D, 2))  # -Y = S^3_{-D/2}(mirror of A # B)
                    names_m = regina_names(Em)
                    plus = bool(set(names) & set(sig))             # Sigma(K) = Y would match e = 1
                    minus = bool(set(names_m) & set(sig))          # Sigma(K) = -Y would match e = -1
                    out = dict(route=r['crossing'], diagram=r.get('diagram_index'),
                               kappa_hat='%s%s # %s%s' % ('' if sa == 1 else '-', TORUS[a], '' if sb == 1 else '-', TORUS[b]),
                               hfk_equal=True, regina_names=names, regina_names_reversed=names_m, name_of_Y_is_a_name_of_sigma=plus,
                               name_of_minus_Y_is_a_name_of_sigma=minus,
                               orientation=('Sigma(K) = Y' if plus and not minus else 'Sigma(K) = -Y' if minus and not plus else None))
                    if not (plus or minus):
                        return out
                    g = genus
                    gens, arrows = tensor(*cplx)
                    V = [V_model(gens, arrows, j, g) for j in range(g + 1)]
                    FY = surgery_form(D, 2, V)
                    eo = 1 if out['orientation'] == 'Sigma(K) = Y' else -1 if out['orientation'] == 'Sigma(K) = -Y' else e
                    FZ = negate(FY) if eo == 1 else FY             # Z = -Sigma(K); the unsigned test does not use eo
                    import montesinos_dinv_check as MC
                    passing = test(FZ, (1, -1))
                    ob2, _ = MC.obstructed(FZ.D, FZ.adj[0][0], {k: FZ.d[(k,)] for k in range(FZ.D)})
                    assert ob2 == (passing is None)
                    allowed = {2: (1,), 6: (-1,)}.get(rec['s'], (1, -1))
                    signed = None
                    if out['orientation'] is not None:
                        signed = passing if passing and passing['sign'] in allowed else test(FZ, allowed)
                    out.update(V=V, d_Z={str(k[0]): str(x) for k, x in FZ.d.items()}, lambda_Z_mu_mu='%d/%d' % (FZ.adj[0][0], D),
                               passing=passing, passing_signed=signed, obstructed=passing is None, obstructed_independent_code=ob2,
                               within_bounds=V in admissible_V(g, max(hh['tau'], 0)))
                    return out
    return None


def stored_check(rec):
    """The stored record is self-consistent: d_Z is the Ni-Wu d-function of the stored V (with Z = -e Y'), and the
    unsigned test on the stored d_Z and lambda_Z gives the stored verdict."""
    if rec['verdict'] not in ('obstructed', 'passes') or not rec.get('d_Z') or rec.get('lambda_Z_mu_mu') is None:
        return None
    D = rec['det']
    a = int(rec['lambda_Z_mu_mu'].split('/')[0])
    FZ = Cyclic(D, a, [Fr(rec['d_Z'][str(k)]) for k in range(D)])
    ok_test = (test(FZ, (1, -1)) is None) == (rec['verdict'] == 'obstructed')
    ok_d = None
    if rec.get('V') is not None and rec.get('e') in (1, -1) and str(rec.get('route_used', '')).startswith(('crossing', 'diagram')):
        FY = surgery_form(D, 2, rec['V'])
        F2 = negate(FY) if rec['e'] == 1 else FY
        ok_d = all(F2.d[(k,)] == FZ.d[(k,)] for k in range(D)) and F2.adj[0][0] == a % D
    return dict(test=ok_test, d_from_V=ok_d)


def khovanov_check(rec, BY):
    """rank Kh~(K; F_2) (KnotInfo, reduced, mod 2) = det K implies that Sigma(K) is an L-space: the spectral sequence of
    Ozsvath-Szabo from reduced Khovanov homology to HF-hat(Sigma(K)) gives rank HF-hat <= rank Kh~, and rank HF-hat >= det.
    Then kappa_hat (Y' = S^3_{D/2}(kappa_hat) = +-Sigma(K), D/2 > 0) must be a positive L-space knot: every route that
    gave V_j must say so, and a route with an unchecked chirality is decided by it."""
    import sympy
    poly = BY[rec['name']].get('khovanov_reduced_mod2_polynomial') or ''
    if not poly:
        return None
    t, q = sympy.symbols('t q')
    rk = int(sympy.sympify(poly.replace('^', '**'), locals={'t': t, 'q': q}).subs({t: 1, q: 1}))
    out = dict(rank=rk, det=rec['det'], sigma_is_L_space=rk == rec['det'])
    if not out['sigma_is_L_space']:
        return out
    routes = [r for r in rec.get('routes', []) if r.get('V') is not None]
    out['routes_with_L_space_Yprime'] = sorted({bool(r.get('L_space')) for r in routes})
    assert all(r.get('L_space') for r in routes), 'Sigma(K) is an L-space, but a route gives a kappa_hat that is not'
    if rec.get('chirality_robustness'):     # the construction chirality gives an L-space Y', the other one does not
        r0 = next(r for r in routes if r.get('unverified_chirality'))
        oc = r0.get('other_chirality', {})
        assert r0['V_source'] == 'L-space knot' and oc.get('V_source') != 'L-space knot'
        D = rec['det']
        FZ = Cyclic(D, int(r0['lambda_Z_mu_mu'].split('/')[0]), [Fr(r0['d_Z'][str(k)]) for k in range(D)])
        passing = test(FZ, (1, -1))
        import montesinos_dinv_check as MC
        ob2, _ = MC.obstructed(FZ.D, FZ.adj[0][0], {k: FZ.d[(k,)] for k in range(FZ.D)})
        assert ob2 == (passing is None)
        allowed = {2: (1,), 6: (-1,)}.get(rec['s'], (1, -1))
        signed = passing if passing and passing['sign'] in allowed else test(FZ, allowed)
        out['chirality'] = 'construction (the other chirality gives a non-L-space Y\')'
        rec.update(verdict_before_khovanov=rec['verdict'], verdict='obstructed' if passing is None else 'passes',
                   route_used='crossing %d' % r0['crossing'], e=r0['e'], surgery=r0['surgery'], V=r0['V'], V_source=r0['V_source'],
                   d_Z=r0['d_Z'], lambda_Z_mu_mu=r0['lambda_Z_mu_mu'], passing=passing, passing_signed=signed,
                   obstructed=passing is None, obstructed_signed=signed is None, obstructed_independent_code=ob2)
        rec.pop('flag', None)
    return out


def model_checks():
    """V_j formulas against the model complexes: staircases (torus knots and the kappa of 11n_176), their duals, and
    the thin staircases of T(2, 2 tau + 1)."""
    tests = {'T(2,3)': {1: 1, 0: -1, -1: 1}, 'T(3,4)': {3: 1, 2: -1, 0: 1, -2: -1, -3: 1},
             'T(3,5)': {4: 1, 3: -1, 1: 1, 0: -1, -1: 1, -3: -1, -4: 1},
             'kappa(11n_176)': {e - 9: c for e, c in enumerate([1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 1, -1, 1, 0, -1, 1, 0, -1, 1]) if c}}
    out = {}
    for nm, al in tests.items():
        g = max(al)
        gens, arrows = staircase(al)
        dg, da = dual_complex(gens, arrows)
        out[nm] = dict(staircase=[V_model(gens, arrows, j, g) for j in range(g + 1)] == torsion(al),
                       dual=[V_model(dg, da, j, g) for j in range(g + 1)] == [0] * (g + 1))
    for tau in range(1, 7):
        gens, arrows = staircase({k: (-1) ** (tau - k) for k in range(-tau, tau + 1)})
        out['thin tau=%d' % tau] = [V_model(gens, arrows, j, tau) for j in range(tau + 1)] == [-((j - tau) // 2) for j in range(tau + 1)]
    return out


def summary(cal, tgt, pilot, gap):
    """The field summary of data/surgery_dinv.json."""
    from collections import Counter
    G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))
    ver = Counter(v['verdict'] for v in tgt.values())
    obs = sorted(k for k, v in tgt.items() if v['verdict'] == 'obstructed')
    out = dict(
        calibration=dict(knots=len(cal), verdicts=dict(Counter(v['verdict'] for v in cal.values())),
                         agree_with_goeritz=sorted(k for k, v in cal.items() if v.get('calibration', {}).get('agrees') is True),
                         disagree=sorted(k for k, v in cal.items() if v.get('calibration') and not v['calibration']['agrees']),
                         not_completed={k: v.get('flag') or v.get('error') for k, v in cal.items() if v['verdict'] not in ('passes', 'obstructed')},
                         obstructed_with_knotinfo_value_1=sorted(k for k, v in cal.items() if v['verdict'] == 'obstructed' and G4.get(k) == '1')),
        pilot={k: dict(verdict=v['verdict'], passing=v.get('passing'), surgery=v.get('surgery'), V=v.get('V')) for k, v in pilot.items()},
        targets=dict(knots=len(tgt), verdicts=dict(ver),
                     obstructed=obs, gamma4_two=[k for k in obs if G4.get(k) == '[1,2]'],
                     obstructed_by_bounds={k: v.get('bounds_variant') for k, v in sorted(tgt.items()) if v['verdict'] == 'obstructed_by_bounds'},
                     passes_for_all_V=sorted(k for k, v in tgt.items() if v['verdict'] == 'passes_for_all_V'),
                     bounds_construction_chirality={k: {r: (b or {}).get('all_obstructed') for r, b in v['bounds_construction_chirality'].items()}
                                                    for k, v in sorted(tgt.items()) if v.get('bounds_construction_chirality')},
                     obstructed_signed_only=sorted(k for k, v in tgt.items() if v['verdict'] == 'passes' and v.get('obstructed_signed')),
                     passes={k: dict(unsigned=v.get('passing'), signed=v.get('passing_signed')) for k, v in sorted(tgt.items()) if v['verdict'] == 'passes'},
                     flagged={k: v.get('flag') for k, v in sorted(tgt.items()) if v['verdict'] == 'flagged'},
                     errors={k: v.get('error') for k, v in sorted(tgt.items()) if v['verdict'] == 'error'},
                     route_used=dict(Counter((v.get('route_used') or '').split(' ')[0] for v in tgt.values() if v['verdict'] in ('passes', 'obstructed'))),
                     V_source=dict(Counter(v.get('V_source') for v in tgt.values() if v['verdict'] in ('passes', 'obstructed'))),
                     plumbing_checked=sorted(k for k, v in tgt.items() if v.get('plumbing', {}).get('fibres')),
                     sigma_L_space_by_khovanov=sorted(k for k, v in tgt.items() if (v.get('khovanov') or {}).get('sigma_is_L_space')),
                     chirality_by_khovanov=sorted(k for k, v in tgt.items() if (v.get('khovanov') or {}).get('chirality'))),
        gap=dict(knots=sorted(gap), certified=sorted(k for k, v in gap.items() if v.get('certificate')),
                 not_certified=sorted(k for k, v in gap.items() if not v.get('certificate')),
                 contradictions=sum(v.get('contradictions') or 0 for v in gap.values())))
    return out


if __name__ == '__main__':
    from database_knotinfo import link_list
    mode = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    BY = {r['name']: r for r in link_list()[1:]}
    R1 = json.load(open(os.path.join(ROOT, 'data', 'rank1_targets.json')))
    P = lambda m: os.path.join(ROOT, 'data', 'surgery_dinv_%s.partial' % m)
    if mode in ('calibration', 'targets', 'pilot'):
        names = ['11n_176'] if mode == 'pilot' else targets() if mode == 'targets' else calibration_knots(BY)
        items = []
        for n in names:
            info = knot_info(BY, n, R1[n]['s'] if n in R1 else None)
            if mode == 'calibration':
                info['goeritz'] = goeritz_P(json.loads(BY[n]['pd_notation']))
                info['hfk_seconds'] = 3 * HFK_SECONDS      # some calibration knots have large kappa
            items.append((n, json.loads(BY[n]['pd_notation']), info))
        res = run_items(items, workers, P(mode))
        print(json.dumps({k: (v.get('verdict'), v.get('passing')) for k, v in sorted(res.items())}))
    elif mode == 'gap':
        T = json.load(open(P('targets')))
        names = sorted(k for k, v in T.items() if v.get('verdict') in ('obstructed', 'obstructed_by_bounds'))
        gap_records(names, BY, workers, P('gap'))
    elif mode == 'summary':
        cal, tgt, pil = (json.load(open(P(m))) for m in ('calibration', 'targets', 'pilot'))
        gp = json.load(open(P('gap'))) if os.path.exists(P('gap')) else {}
        gp = {k: v for k, v in gp.items() if tgt.get(k, {}).get('verdict') in ('obstructed', 'obstructed_by_bounds')}
        for k, v in tgt.items():                 # the linking-form test excludes none of the targets
            assert R1[k]['obstructed'] is False, k
        meta = dict(versions=versions(), convention_check=convention_check(BY), lens_labelling_check=lens_check(),
                    V_model_check=model_checks(),
                    parameters=dict(simplify_tries=SIMPLIFY_TRIES, slope_box=SLOPE_BOX, iso_tries=ISO_TRIES,
                                    drill_segments=DRILL_SEGMENTS, knot_seconds=KNOT_SECONDS, hfk_seconds=HFK_SECONDS,
                                    max_girth=MAX_GIRTH, girth_tries=GIRTH_TRIES, max_routes=MAX_ROUTES,
                                    max_drill_curves=MAX_DRILL_CURVES, max_drill_forms=MAX_DRILL_FORMS,
                                    diagrams=[DIAGRAMS_MINIMAL, DIAGRAMS_NONMINIMAL]))
        assert meta['convention_check']['ok'] and all(meta['lens_labelling_check'].values()) and all(
            all(v.values()) if isinstance(v, dict) else v for v in meta['V_model_check'].values())
        for part in (cal, tgt, pil):             # Khovanov rank = det: Sigma(K) is an L-space
            for k, v in part.items():
                v['khovanov'] = khovanov_check(v, BY)
        for k, v in tgt.items():                 # flagged knots: the reasons, and the verdicts of bounds_verdict()
            if v['verdict'] == 'flagged':
                if not v.get('chirality_robustness'):
                    v['flag'] = flag_text(v)
                bounds_verdict(v)
        for k, v in tgt.items():                 # no signed verdict when Regina names gave Sigma(K) = +-Y, as in one()
            r0 = next((r for r in v.get('routes', []) if v.get('route_used') in ('crossing %d' % r.get('crossing', -1),
                       'diagram %s crossing %s' % (r.get('diagram_index'), r.get('crossing')))), None)
            if r0 and 'Regina' in str(r0.get('kappa', {}).get('diagram', {}).get('orientation')):
                v['passing_signed'] = None; v['obstructed_signed'] = None
        for k, v in tgt.items():                 # V_j exactly when kappa_hat is a connected sum of two torus knots
            if v['verdict'] in ('flagged', 'obstructed_by_bounds', 'passes_for_all_V'):
                ts = torus_sum_route(v, BY)
                if ts:
                    v['torus_sum_route'] = ts
                    if ts.get('V') is not None:
                        assert ts['within_bounds'], k
                        v.update(verdict_before_torus_sum=v['verdict'], verdict='obstructed' if ts['obstructed'] else 'passes',
                                 route_used='connected sum %s (Regina)' % ts['kappa_hat'], V=ts['V'], V_source='tensor product',
                                 d_Z=ts['d_Z'], lambda_Z_mu_mu=ts['lambda_Z_mu_mu'], passing=ts['passing'],
                                 passing_signed=ts['passing_signed'], obstructed=ts['obstructed'],
                                 obstructed_signed=(ts['passing_signed'] is None) if ts['orientation'] else None,
                                 obstructed_independent_code=ts['obstructed_independent_code'])
        for k, v in tgt.items():                 # independent HFK of kappa from the census triangulation
            if v['verdict'] in ('obstructed', 'obstructed_by_bounds'):
                v['census_recheck'] = census_recheck(v)
                assert v['census_recheck'] is None or v['census_recheck'].get('consistent', True), k
            if v['verdict'] == 'obstructed_by_bounds':   # the bounds again, with the independent implementation
                v['bounds_recheck'] = bounds_recheck(v)
                assert v['bounds_recheck'] and all(x['obstructed_independent'] == x['admissible'] for x in v['bounds_recheck']), k
        for part in (cal, tgt, pil):             # the stored d-functions reproduce the verdicts
            for k, v in part.items():
                v['stored_check'] = stored_check(v)
                assert v['stored_check'] is None or (v['stored_check']['test'] and v['stored_check']['d_from_V'] is not False), k
        out = dict(meta=meta, summary=summary(cal, tgt, pil, {k: gap_summary(v) for k, v in gp.items()}),
                   calibration=dict(sorted(cal.items())), pilot=pil, targets=dict(sorted(tgt.items())),
                   gap={k: gap_summary(v) for k, v in sorted(gp.items())})
        json.dump(out, open(os.path.join(ROOT, 'data', 'surgery_dinv.json'), 'w'), indent=0, default=str)
        if gp:
            json.dump({k: dict(diagrams=[dict(source=d['source'], pd=d['pd'], isometric_to_knotinfo=d.get('isometric_to_knotinfo'))
                                         for d in v.get('diagrams', [])],
                               hits=[{x: h.get(x) for x in ('diagram', 'face', 'darts', 'edges', 't', 'J_label', 'J_pd', 'passes',
                                                            'contradiction')} for h in v.get('hits', [])])
                       for k, v in sorted(gp.items())},
                      open(os.path.join(ROOT, 'data', 'surgery_dinv_gap.json'), 'w'), default=str)
        print(json.dumps(out['summary'], indent=1, default=str))
