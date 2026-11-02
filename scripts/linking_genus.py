#!/usr/bin/env python
"""h_lambda and h^oplus of the 12,965 prime knots with at most 13 crossings ([C] Proposition 3.1(a), and the
congruences of [T] Lemma 2.5 in [C] Section 3.5).

- h_witt is h_lambda ([T] Definition 4.1): the least rank h of an odd or empty integral symmetric form Q of odd
  determinant with [lambda_Q] = [lambda_K] in the Witt group.  It is a lower bound for gamma_4^top ([T] Theorem 4.2):
  if K bounds a locally flat F with b_1(F) = h, then W = Sigma_2(B^4, F) has [lambda_K] = [lambda_{Q_W}] ([T] Lemma
  3.1), rank Q_W = h, and Q_W is odd when F is nonorientable ([T] Lemma 3.3).  It is computed by [T] Theorem 4.4
  (h_witt below), with the 2-adic condition (R_h) evaluated through the oddities of [T] Lemma 2.4.
- h_summand is h^oplus ([T] Definition 4.8): the least rank h of an odd Q of odd determinant with
  lambda_K = lambda_Q + metabolic, orthogonally.  It would be a lower bound if the splitting lemma (GL Lemma E.1) held,
  but that lemma fails in general ([T] Example 3.5), and h^oplus is not a lower bound for gamma_4^top ([T] Theorem 5.1).
  Always h_witt <= h_summand.  By [T] Proposition 4.9 they differ only at an exceptional prime (p = 3 mod 4, all
  odd-level Jordan blocks of rank one in one unit class, 3 mod 4 of them; this needs v_p(det K) >= 9), and then
  h_summand is 3 or 4 (h_summand_closed below).  No knot up to 14 crossings has one, and the two agree on all 12,965.
- sigma8_from_lambda is sigma mod 8 computed from lambda by the first congruence of [T] Lemma 2.5 (minus the sum of
  the p-excesses), and sigma_check compares it with sigma = sign(V + V^T).

The local test of h_summand for a given h (passes below).  Let D = det K, and for each odd p | D take the Jordan data
of lambda_p (level v, unit u).  A form Q with det Q = tau l (tau = +-1, D / l a square) exists iff (1) and (2) hold.
(1) At every p | D, Q_p is (a sub-form of each level-v block of lambda_p: rank s_v, discriminant e_v) + (a unimodular
    part of rank h - sum s_v), with
    - complement Witt-trivial at the odd levels, and
    - if sum s_v = h, Legendre(prod e_v, p) = Legendre(tau D_p', p), where D_p' = D / p^{v_p(D)}.
    A sub-form of rank s_v < r_v may have either discriminant; s_v = r_v forces e_v = d_v.
(2) There is an odd unimodular rank-h Z_2-lattice with det = tau D (mod 8) and oddity t, and a signature S with
    (-1)^{(h-S)/2} = tau and S = t + sigma (mod 8).
    This is the oddity formula: the p-excess of Q equals that of lambda, since metabolic forms have excess 0, and
    sigma(K) = -sum_p excess(lambda_p) (mod 8), because V + V^T is even.  Conway-Sloane: local symbols plus the oddity
    formula give a global genus.
Usage: python scripts/linking_genus.py [workers]
Output: data/linking_genus.json
"""
import json, os, sys
from itertools import product
from multiprocessing import Pool
from rank1_linking import invariants, legendre, signature
from rank1_goeritz import det
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')


def excess(J, p):
    """p-excess (mod 8) of Jordan data [(level, unit)]: sum over constituents of (p^v - 1) + 4[v odd, unit nonsquare]."""
    return sum((p ** v - 1) + (4 if v % 2 and legendre(u, p) == -1 else 0) for v, u in J) % 8


def blocks(J, p):
    L = {}
    for v, u in J:
        r, d = L.get(v, (0, 1))
        L[v] = (r + 1, d * u % p)
    return L


def local_ok(J, p, h, Dp, tau):
    n = next(x for x in range(2, p) if legendre(x, p) == -1)
    L = blocks(J, p)
    options = []
    for v, (r, d) in L.items():
        opts = [(0, 1)] + [(s, e) for s in range(1, r) for e in (1, n)] + [(r, d)]
        options.append([(v, r, d, s, e) for s, e in opts])
    for choice in product(*options):
        s_tot = sum(c[3] for c in choice)
        if s_tot > h:
            continue
        rank, disc = 0, 1
        for v, r, d, s, e in choice:
            if v % 2:
                rank += r - s; disc = disc * d * e % p
        if rank % 2 or legendre((-1) ** (rank // 2) * disc, p) != 1:
            continue
        if s_tot < h:
            return True
        prod_e = 1
        for c in choice:
            prod_e = prod_e * c[4] % p
        if legendre(prod_e, p) == legendre(tau * Dp, p):
            return True
    return False


ODDITIES = {}
def oddities(h, d8):
    """Oddities of odd unimodular rank-h Z_2-lattices with det = d8 (mod 8)."""
    key = (h, d8)
    if key not in ODDITIES:
        out = set()
        for us in product((1, 3, 5, 7), repeat=h):
            pr = 1
            for u in us:
                pr = pr * u % 8
            if pr == d8 % 8:
                out.add(sum(us) % 8)
        ODDITIES[key] = out
    return ODDITIES[key]


def passes(inv, D, sigma8, h):
    if h == 0:
        return all(local_ok(J, p, 0, 1, 1) for p, J in inv.items())
    for tau in (1, -1):
        # the 2-adic condition (R_h) of [T] Theorem 4.4; for h = 1 it follows from the local conditions ([T] Cor. 4.5)
        sig_ok = any((-1) ** ((h - S) // 2) == tau and (S - sigma8) % 8 in oddities(h, tau * D % 8)
                     for S in range(-h, h + 1, 2))
        if not sig_ok:
            continue
        if all(local_ok(J, p, h, D // p ** sum(v for v, _ in J), tau) for p, J in inv.items()):
            return True
    return False


def h_summand(inv, D, sigma8, hmax=6):
    """h^oplus by a direct search over orthogonal summands (the local test passes() for h = 0, 1, ...)."""
    for h in range(hmax + 1):
        if passes(inv, D, sigma8, h):
            return h
    return None


def anisotropic(J, p):
    """Rank a (0, 1, 2) and discriminant of the anisotropic representative of the Witt class of lambda_p."""
    odd = [u for v, u in J if v % 2]
    r, d = len(odd), 1
    for u in odd:
        d = d * u % p
    if r % 2:
        return 1, (-1) ** ((r - 1) // 2) * d % p
    if legendre((-1) ** (r // 2) * d, p) == 1:
        return 0, 1
    return 2, (-1) ** (r // 2 - 1) * d % p


def exceptional_primes(inv):
    """Primes p = 3 mod 4 whose odd-level Jordan blocks all have rank one and one unit class, and are 3 mod 4 in number
    ([T] Definition 3.11).  Only at such primes do the summand-level and the Witt-level tests differ."""
    out = []
    for p, J in inv.items():
        odd = {}
        for v, u in J:
            if v % 2:
                odd.setdefault(v, []).append(u)
        if (p % 4 == 3 and odd and all(len(us) == 1 for us in odd.values()) and len(odd) % 4 == 3
                and len({legendre(us[0], p) for us in odd.values()}) == 1):
            out.append(p)
    return out


def h_summand_closed(inv, D, sigma8):
    """h^oplus by [T] Proposition 4.9: h_witt without exceptional primes; otherwise 3 if the symbols (u_p D'_p / p) of
    the exceptional primes agree and 4 if not."""
    T = exceptional_primes(inv)
    if not T:
        return h_witt(inv, D, sigma8)
    symbols = {legendre(next(u for v, u in inv[p] if v % 2) * (D // p ** sum(v for v, _ in inv[p])), p) for p in T}
    return 3 if len(symbols) == 1 else 4


def h_witt(inv, D, sigma8):
    """h_lambda by [T] Theorem 4.4: the least h such that some odd rank-h Q of odd determinant has
    [lambda_Q] = [lambda].  At p the form Q_p needs a_p odd-level constituents; if a_p = h, it has no free unit, and
    disc A_p = tau D_p' (mod squares)."""
    A = {p: anisotropic(J, p) for p, J in inv.items()}
    for h in range(0, 4):
        if any(a > h for a, _ in A.values()):
            continue
        for tau in (1, -1):
            if h and not any((-1) ** ((h - S) // 2) == tau and (S - sigma8) % 8 in oddities(h, tau * D % 8)
                             for S in range(-h, h + 1, 2)):
                continue
            if all(a < h or legendre(d, p) == legendre(tau * (D // p ** sum(v for v, _ in inv[p])), p)
                   for p, (a, d) in A.items() if a):
                return h
    return None


def one(item):
    name, V = item
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    D = abs(det(A))
    inv = invariants(A, D) if D > 1 else {}
    sig = signature(A)
    sigma8 = (-sum(excess(J, p) for p, J in inv.items())) % 8
    return name, dict(det=D, sigma=sig, sigma8_from_lambda=sigma8, sigma_check=(sig - sigma8) % 8 == 0,
                      h_summand=h_summand(inv, D, sigma8), h_witt=h_witt(inv, D, sigma8))


if __name__ == '__main__':
    from collections import Counter
    from database_knotinfo import link_list
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    BY = {r['name']: r for r in link_list()[1:]}
    G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))
    with Pool(workers) as p:
        out = dict(p.imap_unordered(one, [(k, json.loads(BY[k]['seifert_matrix'])) for k in G4], chunksize=20))
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', 'linking_genus.json'), 'w'), indent=0)
    print('knots', len(out), '| sigma = -sum excess (mod 8):', sum(r['sigma_check'] for r in out.values()),
          '| h_lambda (= h_witt) distribution:', sorted(Counter(r['h_witt'] for r in out.values()).items()),
          '| h_summand equals h_witt:', sum(r['h_summand'] == r['h_witt'] for r in out.values()))
    R1 = {}
    for f in ('rank1_targets.json', 'rank1_controls.json'):
        R1.update(json.load(open(os.path.join(ROOT, 'data', f))))
    R2 = {}
    for f in ('rank2_targets.json', 'rank2_controls.json'):
        R2.update(json.load(open(os.path.join(ROOT, 'data', f))))
    print('rank 1 vs h_witt >= 2:', Counter((R1[k]['obstructed'], out[k]['h_witt'] >= 2) for k in R1))
    print('rank 2 (odd) vs h_witt >= 3:', Counter((R2[k]['obstructed_odd'], out[k]['h_witt'] >= 3) for k in R2))
    print('h_witt >= 3 among all knots:', sorted(k for k, r in out.items() if r['h_witt'] >= 3))
    print('h_witt vs KnotInfo value for h_witt >= 3:', Counter((r['h_witt'], G4[k]) for k, r in out.items() if r['h_witt'] >= 3))
