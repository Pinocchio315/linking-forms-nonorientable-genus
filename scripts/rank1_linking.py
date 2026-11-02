#!/usr/bin/env python
"""Rank-one linking-form test: [C] Proposition 3.1(b), and its checks in [C] Section 3.5.

lambda = lambda_K is the linking form lambda(x, y) = x^T (V + V^T)^{-1} y mod 1 on G = coker(V + V^T), for the Seifert
matrix V of the KnotInfo package; it is minus the linking form of Sigma_K ([T] (6)), and the test is symmetric under
lambda -> -lambda.  The script finds all (n, eps), with eps = +-1 and det K / n a square, for which
lambda = <eps/n> + mu with mu metabolic (orthogonal sum); K is obstructed when there is none.

This orthogonal-sum condition is the one that the splitting lemma (Gilmer, Lemma 1; GL Lemma E.1) would give, and that
lemma fails in general ([T] Example 3.5).  What [T] proves is the Witt-level condition [lambda] = [<eps/l>], with
Q_W = (q), l = |q| and det K / l a square ([T] Lemma 3.1 and Theorem 4.2), which is condition (b) below alone; it
obstructs K iff h_lambda(K) >= 2 ([C] Corollary 2.3(b)).  Condition (a) can exclude more only at an exceptional prime
([T] Definition 3.11 and Proposition 4.9), and survivors() refuses inputs with an exceptional prime.  An exceptional
prime p needs p^9 | det K, and every knot with at most 14 crossings has det K <= 1,145 < 3^9, so every obstruction
found here is a Witt-level obstruction (computed directly by linking_genus.py for 13 crossings and by rank1_ls14.py
for 14).

Method.  |G| = det K is odd, so lambda splits over the odd primes p | det K, and a form on a p-group is determined by
its Jordan data: diagonalising V + V^T over Z_(p) gives entries p^v u, and the p-part of lambda is the orthogonal sum of
<u^{-1}/p^v> over v >= 1 ([T] (3)).  Then
  (a) <eps/n> is an orthogonal summand iff for every p | n, with a = v_p(n), the level-a block is nonempty and, if it
      has rank one, its unit is in the square class of eps * n / p^a mod p (forms of rank >= 2 over F_p are universal);
  (b) the complement is metabolic iff the Witt class of lambda - <eps/n> vanishes: for every p, the F_p-form of all
      units at odd levels, together with -eps * n / p^a when a is odd, has even rank r and (-1)^{r/2} * disc a square.

Sign refinement (field obstructed_signed).  Let s = sign(V + V^T) + 4 Arf(K) mod 8.  If W is the double branched cover
of B^4 along a Mobius band and Q_W = (q), then eps = Sign W = sign q.  This follows from the Witt-level statement
[lambda_Y] = [-lambda_{Q_W}] ([T] Lemma 3.1), compared with the simply connected cover along a pushed-in Seifert
surface, whose intersection form is +-(V + V^T) and whose boundary has linking form minus its inverse
(Gordon-Litherland, pp. 59-60; the same sign as in GL Thm. 6).  GL Cor. 8 with Brown invariant +-1 then gives eps = +1
if s = 2 and eps = -1 if s = 6; the global sign convention cancels.  Check: T(2,3) with s = 2 has lambda = <1/3>,
T(2,7) with s = 2 has lambda = <1/7>.

Modes: validate, every knot of determinant at most 200, by the test and by brute(), an enumeration of G ([C] Section
3.5); controls, the knots with KnotInfo value 1, none of which may be obstructed ([C] Section 3.5); targets, the knots
whose KnotInfo range starts at 1 ([C] Proposition 3.1(b)).
Usage: python scripts/rank1_linking.py validate | controls | targets [workers]
Output: data/rank1_<mode>.json
"""
import json, os, sys
from fractions import Fraction
from math import gcd, isqrt
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')


def vp(x, p):
    v, num, den = 0, x.numerator, x.denominator
    while num % p == 0:
        num //= p; v += 1
    while den % p == 0:
        den //= p; v -= 1
    return v


def unit_mod_p(x, p):
    num, den = x.numerator, x.denominator
    while num % p == 0:
        num //= p
    while den % p == 0:
        den //= p
    return num * pow(den, -1, p) % p


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def primes_of(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def jordan(A, p):
    """Diagonalise the symmetric integer matrix A over Z_(p), p odd: list of (valuation, unit mod p)."""
    n = len(A)
    M = [[Fraction(x) for x in row] for row in A]
    active, out = list(range(n)), []
    while active:
        best = None
        for i in active:
            for j in active:
                if M[i][j] != 0:
                    key = (vp(M[i][j], p), i != j)
                    if best is None or key < best[0]:
                        best = (key, i, j)
        (v, offdiag), i, j = best
        if offdiag:                       # every diagonal entry has larger valuation: replace e_i by e_i + e_j
            for c in range(n):
                M[i][c] += M[j][c]
            for r in range(n):
                M[r][i] += M[r][j]
        piv = M[i][i]
        assert vp(piv, p) == v
        for r in active:
            if r != i and M[r][i] != 0:
                f = M[r][i] / piv         # p-integral, since v_p(M[r][i]) >= v
                for c in range(n):
                    M[r][c] -= f * M[i][c]
                for s in range(n):
                    M[s][r] -= f * M[s][i]
        out.append((v, unit_mod_p(piv, p)))
        active.remove(i)
    return out


def signature(A):
    """Exact signature of a nonsingular symmetric integer matrix: congruence diagonalisation over Q."""
    n = len(A)
    M = [[Fraction(x) for x in row] for row in A]
    active, sig = list(range(n)), 0
    while active:
        i = next((i for i in active if M[i][i] != 0), None)
        if i is None:                     # zero diagonal: replace e_i by e_i + e_j
            i, j = next((i, j) for i in active for j in active if M[i][j] != 0)
            for c in range(n):
                M[i][c] += M[j][c]
            for r in range(n):
                M[r][i] += M[r][j]
        piv = M[i][i]
        for r in active:
            if r != i and M[r][i] != 0:
                f = M[r][i] / piv
                for c in range(n):
                    M[r][c] -= f * M[i][c]
                for s in range(n):
                    M[s][r] -= f * M[s][i]
        sig += 1 if piv > 0 else -1
        active.remove(i)
    return sig


def invariants(A, D):
    inv = {}
    for p in primes_of(D):
        J = [(v, u) for v, u in jordan(A, p) if v > 0]
        assert sum(v for v, _ in J) == max(k for k in range(64) if D % p ** k == 0)
        inv[p] = J
    return inv


def witt_zero(units, p):
    r, d = len(units), 1
    for u in units:
        d = d * u % p
    return r % 2 == 0 and legendre((-1) ** (r // 2) * d, p) == 1


def survivors(inv, D, allow_exceptional=False):
    """All (n, eps) with lambda = <eps/n> + metabolic.  With an exceptional prime this summand-level test can exclude
    more than the proven Witt-level test ([T] Proposition 4.9), so it is refused unless allow_exceptional."""
    if not allow_exceptional:
        from linking_genus import exceptional_primes
        bad = exceptional_primes(inv)
        if bad:
            raise ValueError(f'exceptional primes {bad}: the orthogonal-sum test is not a proven obstruction here; '
                             'use linking_genus.h_witt')
    out = []
    choices = [[]]
    for p, J in inv.items():
        e = sum(v for v, _ in J)
        choices = [c + [(p, a)] for c in choices for a in range(e % 2, e + 1, 2)]
    for c in choices:
        n = 1
        for p, a in c:
            n *= p ** a
        for eps in (1, -1):
            ok = True
            for p, a in c:
                J = inv[p]
                cof = eps * (n // p ** a)
                if a:
                    block = [u for v, u in J if v == a]
                    if not block or (len(block) == 1 and legendre(block[0] * cof, p) != 1):
                        ok = False; break
                odd = [u for v, u in J if v % 2] + ([-cof % p] if a % 2 else [])
                if not witt_zero(odd, p):
                    ok = False; break
            if ok:
                out.append((n, eps))
    return out


def brute(A):
    """The same condition by enumeration of G (small groups only).  Metabolisers are searched among subgroups with two
    generators, which is complete only when their order is below 27 (the least odd order needing three generators);
    larger cases raise an error instead of returning an incomplete answer."""
    sys.path.insert(0, os.path.join(ROOT, 'vendor', 'paper_followup', 'negative-definite-plumbings', 'rank_two'))
    from bf2 import LinkingForm
    L = LinkingForm(A)
    D, el, Bx = L.D, L.elems, L.Bx
    key = {tuple(b): k for k, b in enumerate(Bx)}
    lam = lambda x, by: sum(a * b for a, b in zip(x, by)) % D          # D * lambda(x, y) mod D
    out = set()
    for c, bc in zip(el, Bx):
        g = D
        for t in bc:
            g = gcd(g, t)
        n = D // g                                                     # order of c
        s = lam(c, bc)
        for eps in (1, -1):
            if s != (eps * (D // n)) % D or isqrt(D // n) ** 2 != D // n:
                continue
            comp = [k for k, b in enumerate(Bx) if lam(c, b) == 0]
            m = len(comp)
            if m == 1:
                out.add((n, eps)); continue
            iso = [k for k in comp if lam(el[k], Bx[k]) == 0]
            half = isqrt(m)
            if half >= 27:
                raise ValueError(f'metaboliser of order {half}: two generators need not suffice')
            found = False
            for k1 in iso:                     # subgroups generated by at most two isotropic elements
                for k2 in iso:
                    if lam(el[k1], Bx[k2]) != 0:
                        continue
                    span = {tuple([0] * len(bc))}
                    gens = [tuple(Bx[k1]), tuple(Bx[k2])]
                    frontier = list(span)
                    while frontier:
                        nxt = []
                        for x in frontier:
                            for gv in gens:
                                y = tuple((a + b) % D for a, b in zip(x, gv))
                                if y not in span:
                                    span.add(y); nxt.append(y)
                        frontier = nxt
                    if len(span) == half:
                        found = True; break
                if found:
                    break
            if found:
                out.add((n, eps))
    return sorted(out)


def one(item):
    name, V = item
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    import sympy
    D = abs(int(sympy.Matrix(A).det())) if n else 1
    inv = invariants(A, D) if D > 1 else {}
    surv = survivors(inv, D) if D > 1 else [(1, 1), (1, -1)]
    sig, arf = signature(A), 0 if D % 8 in (1, 7) else 1
    s = (sig + 4 * arf) % 8
    eps = {0: (1, -1), 2: (1,), 6: (-1,)}.get(s, ())  # s = 4: no Mobius band (Yasuhara)
    return name, dict(det=D, jordan={str(p): J for p, J in inv.items()}, survivors=surv, obstructed=not surv,
                      sigma=sig, arf=arf, s=s, obstructed_signed=not [e for _, e in surv if e in eps])


def one_validate(item):
    name, V = item
    rec = one(item)[1]
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    b = brute(A) if rec['det'] > 1 else [(1, -1), (1, 1)]
    return name, dict(det=rec['det'], invariant=sorted(map(tuple, rec['survivors'])), brute=b,
                      agree=sorted(map(tuple, rec['survivors'])) == b)


if __name__ == '__main__':
    from collections import Counter
    from database_knotinfo import link_list
    mode = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    BY = {r['name']: r for r in link_list()[1:]}
    G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))
    sm = lambda n: json.loads(BY[n]['seifert_matrix'])
    if mode == 'validate':           # every knot of determinant at most 200, both methods
        names = [n for n in G4 if int(BY[n]['determinant']) <= 200]
        f = one_validate
    elif mode == 'controls':         # knots with gamma_4 = 1: the test must not obstruct any of them
        names = [n for n, v in G4.items() if v == '1']
        f = one
    else:                            # open ranges with lower end 1
        names = [n for n, v in G4.items() if v.startswith('[1,')]
        f = one
    with Pool(workers) as p:
        out = dict(p.imap_unordered(f, [(n, sm(n)) for n in names], chunksize=10))
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', f'rank1_{mode}.json'), 'w'), indent=0)
    if mode == 'validate':
        print('knots', len(out), '| agree', sum(r['agree'] for r in out.values()),
              '| disagreements', [n for n, r in out.items() if not r['agree']][:10])
    else:
        print(mode, 'knots', len(out), '| obstructed', sum(r['obstructed'] for r in out.values()),
              '| obstructed with sign', sum(r['obstructed_signed'] for r in out.values()))
        print('Arf mismatches with KnotInfo:', [n for n, r in out.items() if r['arf'] != int(BY[n]['arf_invariant'])],
              '| sign(V + V^T) = KnotInfo signature:', sum(r['sigma'] == int(BY[n]['signature']) for n, r in out.items()),
              '| s:', Counter(r['s'] for r in out.values()))
        if mode == 'targets':
            print('obstructed by range:', Counter(G4[n] for n, r in out.items() if r['obstructed']))
            print('obstructed with sign by range:', Counter(G4[n] for n, r in out.items() if r['obstructed_signed']))
