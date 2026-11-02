#!/usr/bin/env python
"""Jabuka's rational Witt span ws(K) [Jabuka, Def. 1.3] from the rational Witt class of V + V^T alone, compared with
h_lambda and with [T] Theorem 4.6 ([C] Section 3.6 and the column ws of [C] Table 3).

ws(K) is the least rank of a form B over Q with d_p[B] = d_p W(K) for every prime p, where d_p: W(Q) -> W(F_p) is the
second residue homomorphism [Jabuka, (8) and Sec. 2.3.3]: d_p<p^v u> = <u mod p> for v odd and 0 for v even (for p = 2,
W(F_2) = Z/2 and d_2[B] is the parity of v_2(det B)).  The script uses no linking form, Jordan decomposition or lattice:
- V + V^T is diagonalised over Q in exact fractions, and d_p W(K) is read off from the entries of odd p-adic valuation;
  its anisotropic representative over F_p has rank a_p and discriminant delta_p ([T] (1));
- rank 0 and rank 1 are decided directly;
- ranks 2 and 3 are decided by the existence of a rational form with prescribed localisations: for a determinant d in
  Q^* and forms B_v over Q_v of rank n and determinant d at all places v, with Hasse invariants
  e_v = prod_{i<j} (a_i, a_j)_v equal to 1 for almost all v, there is a rational form B with B (x) Q_v = B_v for all v
  if and only if prod_v e_v = 1 (Conway-Sloane, Ch. 15: the product formula is the sum formula (15), Thm. 5 is the
  existence, Sec. 5.3 relates the Hasse invariant to the p-excess, Sec. 6.1 classifies p-adic binary forms by
  determinant and p-excess).  The Hilbert symbols are the standard ones; at p = 2 this is [Jabuka, eq. (7)], and
  tests() checks [Jabuka, Ex. 2.14 and 2.15] and the product formula.
At a prime q not dividing 2D a form with unit determinant and d_q = 0 has e_q = 1: either all entries are units, or two
of them have odd valuation and span a hyperbolic plane.  So only the places infinity, 2 and the primes dividing D = det K
matter.  A prime at which d has odd valuation forces a residue of odd rank, so d runs over the products of a sign and
primes of D.
The result is compared with the closed form of [T] Theorem 4.6(c), (d) (computed from a_p, delta_p, sigma and D), with
h_lambda of [T] Theorem 4.4 and Corollary 4.5 (recomputed from the same residues), and with data/linking_genus.json and
data/hthree_certificates.json.  For the knots with ws = 2 < h_lambda = 3 an even binary form of determinant -l with the
same residues is found, as in the proof of [T] Theorem 4.6(c).
For the 196 knots with 14 crossings and h_lambda = 3 (data/rank1_ls14.json), d_p W(K) is read off from a Goeritz matrix
of the SnapPy diagram instead: the Gordon-Litherland form of any spanning surface has the same residues as V + V^T,
since the moves relating spanning surfaces change its rational Witt class by [<+-1>] or not at all [Jabuka, proof of
Prop. 3.3].  There sigma mod 8 is taken from the residues by [T] Lemma 2.5 (and compared with data/rank2_ls14.json);
the Goeritz matrix may present the mirror image, which changes neither ws nor the conditions of [T] Theorem 4.6.
Usage: python scripts/witt_span.py [workers]
Output: data/witt_span.json (both parts; the 14-crossing part needs SnapPy)
"""
import json, os, sys
from fractions import Fraction
from itertools import product
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')


def vp(x, p):
    """p-adic valuation of a nonzero integer or Fraction."""
    x = Fraction(x)
    n, d, v = abs(x.numerator), x.denominator, 0
    while n % p == 0:
        n //= p; v += 1
    while d % p == 0:
        d //= p; v -= 1
    return v


def unit_part(x, p):
    """x / p^{v_p(x)} as a Fraction."""
    return Fraction(x) / Fraction(p) ** vp(x, p)


def unit_mod(x, p):
    """The residue mod p of the unit part of x (p odd)."""
    u = unit_part(x, p)
    return u.numerator * pow(u.denominator, -1, p) % p


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def primes_of(n):
    n, out, q = abs(n), [], 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def diagonalise(A):
    """Exact symmetric elimination over Q: a list of nonzero Fractions d_i with A congruent to diag(d_i) over Q."""
    M = [[Fraction(x) for x in row] for row in A]
    out = []
    while M:
        n = len(M)
        k = next((i for i in range(n) if M[i][i] != 0), None)
        if k is None:                       # zero diagonal: e_i -> e_i + e_j makes the (i, i) entry 2 M_ij != 0
            i, j = next((i, j) for i in range(n) for j in range(n) if M[i][j] != 0)
            for c in range(n):
                M[i][c] += M[j][c]
            for r in range(n):
                M[r][i] += M[r][j]
            k = i
        piv = M[k][k]
        out.append(piv)
        M = [[M[r][c] - M[r][k] * M[k][c] / piv for c in range(n) if c != k] for r in range(n) if r != k]
    return out


def aniso(units, p):
    """(a, delta) for the F_p-space <u_1, ..., u_r>: rank and discriminant class (Legendre symbol) of the anisotropic
    space in its Witt class ([T] (1)); delta = 0 when a = 0."""
    r, d = len(units), 1
    for u in units:
        d = d * u % p
    if r % 2:
        return 1, legendre((-1) ** ((r - 1) // 2) * d, p)
    if legendre((-1) ** (r // 2) * d, p) == 1:
        return 0, 0
    return 2, legendre((-1) ** (r // 2 - 1) * d, p)


def residue(entries, p):
    """d_p of the diagonal form diag(entries) over Q_p, as (a, delta)."""
    return aniso([unit_mod(x, p) for x in entries if vp(x, p) % 2], p)


def hilbert(a, b, v):
    """Hilbert symbol (a, b)_v for nonzero rationals a, b; v is a prime or 'inf'.  For a = p^al u, b = p^be w with units
    u, w: (-1)^(al be (p-1)/2) (u/p)^be (w/p)^al for odd p, and (-1)^(e(u) e(w) + al o(w) + be o(u)) for p = 2 with
    e(u) = (u-1)/2 and o(u) = (u^2-1)/8 mod 2 [Jabuka, eq. (7)]."""
    a, b = Fraction(a), Fraction(b)
    if v == 'inf':
        return -1 if a < 0 and b < 0 else 1
    al, be = vp(a, v), vp(b, v)
    if v == 2:
        u, w = unit_part(a, 2), unit_part(b, 2)
        u = u.numerator * u.denominator % 8                 # the class mod 8 of a 2-adic unit u (1/u = u mod 8)
        w = w.numerator * w.denominator % 8
        eps = lambda x: (x - 1) // 2 % 2
        om = lambda x: (x * x - 1) // 8 % 2
        return (-1) ** (eps(u) * eps(w) + al * om(w) + be * om(u))
    u, w = unit_mod(a, v), unit_mod(b, v)
    return (-1) ** (al * be * ((v - 1) // 2) % 2) * legendre(u, v) ** (be % 2) * legendre(w, v) ** (al % 2)


def hasse(entries, v):
    e = 1
    for i in range(len(entries)):
        for j in range(i + 1, len(entries)):
            e *= hilbert(entries[i], entries[j], v)
    return e


def square_classes(v):
    if v == 'inf':
        return [1, -1]
    if v == 2:
        return [1, 3, 5, 7, 2, 6, 10, 14]
    n = next(x for x in range(2, v) if legendre(x, v) == -1)
    return [1, n, v, v * n]


def det_classes(P):
    """Squarefree odd integers built from a sign and the primes P."""
    out = []
    for bits in product((0, 1), repeat=len(P)):
        m = 1
        for q, b in zip(P, bits):
            m *= q ** b
        out += [m, -m]
    return out


def rank_at_most(target, P, n):
    """Is there a rational form B of rank n with d_p[B] = target[p] for every prime p?  (target[p] = (a_p, delta_p) for
    p in P, the odd primes of D; d_p = 0 elsewhere, including p = 2.)"""
    if n == 0:
        return all(t[0] == 0 for t in target.values())
    for d in det_classes(P):
        if n == 1:
            if all(residue([d], p) == target[p] for p in P):
                return True
            continue
        sets = []
        for v in ['inf', 2] + list(P):
            E = set()
            for xs in product(square_classes(v), repeat=n - 1):
                rest = Fraction(d)
                for x in xs:
                    rest *= x
                B = list(xs) + [rest]                       # determinant d times a square
                if v not in ('inf', 2) and residue(B, v) != target[v]:
                    continue
                E.add(hasse(B, v))
            if not E:
                break
            sets.append(E)
        else:
            if any(len(E) == 2 for E in sets):
                return True
            prod_e = 1
            for E in sets:
                prod_e *= next(iter(E))
            if prod_e == 1:
                return True
    return False


def theorem(target, P, sigma, D):
    """The values given by [T] Theorems 4.4 and 4.6 and Corollary 4.5(c), from a_p, delta_p, sigma and D."""
    Dp = {p: D // p ** vp(D, p) for p in P}
    eps = {p: delta * legendre(Dp[p], p) for p, (a, delta) in target.items() if a}
    a2 = [p for p in P if target[p][0] == 2]
    L1 = lambda tau: all(target[p][0] <= 1 for p in P) and all(eps[p] == legendre(tau, p) for p in P if target[p][0] == 1)
    yas = (sigma % 8, D % 8) in ((0, 5), (4, 1))           # sigma + 4 Arf = 4 (mod 8), [T] (5)
    c1 = any(p % 4 == 1 and eps[p] == -1 for p in a2)
    c2 = len({eps[p] for p in a2 if p % 4 == 3}) == 2
    neg3 = any(p % 4 == 3 and eps[p] == -1 for p in a2)
    c3 = neg3 and yas
    c3_41 = neg3 and (sigma % 8, D % 8) == (4, 1)
    crit = [name for name, c in (('(i)', c1), ('(ii)', c2), ('(iii)', c3)) if c]
    if all(t[0] == 0 for t in target.values()):
        h = 0
    elif any(L1(tau) and sigma % 8 == tau * (1 - D) % 8 for tau in (1, -1)):
        h = 1
    elif crit:
        h = 3
    else:
        h = 2
    if all(t[0] == 0 for t in target.values()):
        ws = 0
    elif any(L1(tau) for tau in (1, -1)):
        ws = 1
    elif c1 or c2 or c3_41:
        ws = 3
    else:
        ws = 2
    return dict(h_theorem=h, ws_theorem=ws, criteria=crit, L1_signs=[tau for tau in (1, -1) if L1(tau)],
                eps={str(p): e for p, e in eps.items() if target[p][0] == 2})


def even_binary(target, P, D):
    """An even binary form [[2a, b], [b, 2c]] of determinant -l, l = prod_p p^{a_p}, with the residues target, found by
    a search over a ([T] Theorem 4.6(c) proves that one exists in the case (iii) with (sigma, D) = (0, 5) mod 8)."""
    l = 1
    for p in P:
        l *= p ** target[p][0]
    for m in range(1, 400):
        for a in (m, -m):
            for b in range(-abs(2 * a) + 1, abs(2 * a), 2):
                if (b * b - l) % (4 * a) == 0:
                    c = (b * b - l) // (4 * a)
                    Q = [[2 * a, b], [b, 2 * c]]
                    diag = diagonalise(Q)
                    if all(residue(diag, p) == target[p] for p in P):
                        return Q
    return None


def one(item):
    name, V = item
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    diag = diagonalise(A)
    det = Fraction(1)
    for x in diag:
        det *= x
    assert det.denominator == 1 and det.numerator % 2
    D = abs(det.numerator)
    sigma = sum(1 if x > 0 else -1 for x in diag)
    P = primes_of(D)
    target = {p: residue(diag, p) for p in P}
    ws = next(r for r in range(4) if rank_at_most(target, P, r))
    rec = dict(det=D, sigma=sigma, a={str(p): target[p][0] for p in P}, delta={str(p): target[p][1] for p in P}, ws=ws,
               **theorem(target, P, sigma, D))
    return name, rec


def one_ls14(name):
    import snappy
    from rank1_goeritz import goeritz_pair
    G1, _ = goeritz_pair([list(c) for c in snappy.Link(name).PD_code()])
    diag = diagonalise(G1)
    det = Fraction(1)
    for x in diag:
        det *= x
    assert det.denominator == 1 and det.numerator % 2
    D = abs(det.numerator)
    P = primes_of(D)
    target = {p: residue(diag, p) for p in P}
    sigma8 = -sum(a * (p - 1) + (4 if delta == -1 else 0) for p, (a, delta) in target.items()) % 8   # [T] Lemma 2.5
    ws = next(r for r in range(4) if rank_at_most(target, P, r))
    rec = dict(det=D, sigma8=sigma8, a={str(p): target[p][0] for p in P}, delta={str(p): target[p][1] for p in P}, ws=ws,
               **theorem(target, P, sigma8, D))
    return name, rec


def tests():
    """Hilbert symbols: the examples of [Jabuka, Ex. 2.14 and 2.15] and the product formula on a grid."""
    assert hilbert(13, 24, 2) == -1 and hilbert(-1, -1, 2) == -1
    for a in range(-40, 41):
        for b in range(-40, 41):
            if a and b:
                e = hilbert(a, b, 'inf') * hilbert(a, b, 2)
                for q in primes_of(2 * a * b * 3):
                    if q > 2:
                        e *= hilbert(a, b, q)
                assert e == 1, (a, b)


if __name__ == '__main__':
    from collections import Counter
    from database_knotinfo import link_list
    tests()
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    BY = {r['name']: r for r in link_list()[1:]}
    G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))
    with Pool(workers) as pool:
        out = dict(pool.imap_unordered(one, [(k, json.loads(BY[k]['seifert_matrix'])) for k in G4], chunksize=20))
    lg = json.load(open(os.path.join(ROOT, 'data', 'linking_genus.json')))
    cert = json.load(open(os.path.join(ROOT, 'data', 'hthree_certificates.json')))['knots']
    for k, r in out.items():
        r['h_lambda'] = lg[k]['h_witt']
        assert r['sigma'] == lg[k]['sigma'] and r['det'] == lg[k]['det'], k
    bad_ws = sorted(k for k, r in out.items() if r['ws'] != r['ws_theorem'])
    bad_h = sorted(k for k, r in out.items() if r['h_lambda'] != r['h_theorem'])
    gap = sorted(k for k, r in out.items() if r['ws'] != r['h_lambda'])
    for k in gap:
        r = out[k]
        P = [int(p) for p in r['a']]
        r['even_binary_form'] = even_binary({p: (r['a'][str(p)], r['delta'][str(p)]) for p in P}, P, r['det'])
    summary = dict(
        knots=len(out),
        ws_distribution=dict(sorted(Counter(r['ws'] for r in out.values()).items())),
        h_lambda_distribution=dict(sorted(Counter(r['h_lambda'] for r in out.values()).items())),
        ws_differs_from_theorem=bad_ws, h_lambda_differs_from_theorem=bad_h,
        ws_less_than_h_lambda=gap,
        pairs=dict(sorted(Counter(f"{r['ws']},{r['h_lambda']}" for r in out.values()).items())),
        gap_all_iii_05=all(out[k]['criteria'] == ['(iii)'] and (out[k]['sigma'] % 8, out[k]['det'] % 8) == (0, 5)
                           for k in gap),
        gap_even_forms_found=all(out[k]['even_binary_form'] is not None for k in gap),
        ws_three=sorted(k for k, r in out.items() if r['ws'] == 3),
        hthree_criteria_agree=all(out[k]['criteria'] == cert[k]['criteria'] and
                                  out[k]['eps'] == {str(p): e for p, e in cert[k]['eps']} for k in cert),
        jabuka_examples={k: out[k]['ws'] for k in ('9_40', '7_5', '9_49', '11n_133', '12a_561', '12a_664', '12a_780',
                                                   '12a_907', '12n_276', '3_1', '5_1', '7_1')},
    )
    print(json.dumps(summary, indent=1))
    ls = json.load(open(os.path.join(ROOT, 'data', 'rank1_ls14.json')))
    r2ls = json.load(open(os.path.join(ROOT, 'data', 'rank2_ls14.json')))
    names = sorted(k for k, v in ls.items() if v['h_witt'] == 3)
    with Pool(min(workers, 4)) as pool:
        out14 = dict(pool.imap_unordered(one_ls14, names, chunksize=4))
    arf = lambda D: 0 if D % 8 in (1, 7) else 1
    for k, r in out14.items():
        r['h_lambda'] = ls[k]['h_witt']
        assert r['det'] == ls[k]['det'] == r2ls[k]['det'], k
        assert r2ls[k]['s'] in ((r['sigma8'] + 4 * arf(r['det'])) % 8, (-r['sigma8'] + 4 * arf(r['det'])) % 8), k
    summary14 = dict(
        knots=len(out14),
        ws_distribution=dict(sorted(Counter(r['ws'] for r in out14.values()).items())),
        ws_differs_from_theorem=sorted(k for k, r in out14.items() if r['ws'] != r['ws_theorem']),
        h_lambda_differs_from_theorem=sorted(k for k, r in out14.items() if r['h_lambda'] != r['h_theorem']),
        ws_less_than_h_lambda=sorted(k for k, r in out14.items() if r['ws'] < r['h_lambda']),
        criteria=dict(sorted(Counter(' '.join(r['criteria']) + (' (4,1)' if r['criteria'] == ['(iii)'] and
                                                                (r['sigma8'], r['det'] % 8) == (4, 1) else '')
                                     for r in out14.values()).items())),
    )
    print(json.dumps(summary14, indent=1))
    json.dump(dict(summary=summary, knots=dict(sorted(out.items())), summary_ls14=summary14,
                   knots_ls14=dict(sorted(out14.items()))), open(os.path.join(ROOT, 'data', 'witt_span.json'), 'w'), indent=0)
