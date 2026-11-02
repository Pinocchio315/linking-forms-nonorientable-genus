#!/usr/bin/env python
"""Rank-two test: can K bound a punctured Klein bottle F in B^4?  The check of h_lambda >= 3 by binary forms ([C]
Proposition 3.1(a) and Section 3.5, Linking forms), and the binary forms of [C] Section 3.6.

W = Sigma_2(B^4, F) has b_2 = 2 ([T] Lemma 3.3) and an integral binary intersection form Q = Q_W with |det Q| = l and
det K / l a square, and [lambda_Y] = [-lambda_Q] for Y = Sigma_K ([T] Lemma 3.1).  As lambda_K, presented by V + V^T,
is isometric to -lambda_Y, this reads [lambda_K] = [lambda_Q].
- The linking-form test uses the orthogonal condition lambda_K = lambda_Q + metabolic (summand_ok), the conclusion of
  the splitting lemma (GL Lemma E.1), which fails in general ([T] Example 3.5).  It can exclude more than the
  Witt-level condition only at an exceptional prime ([T] Definition 3.11 and Proposition 4.9), which needs
  v_p(det K) >= 9 and does not occur up to 14 crossings; linking_signs() refuses such inputs.  On the knots with at most
  13 crossings the result with odd Q agrees with h_lambda >= 3 (linking_genus.py), so the obstructions found here are
  Witt-level obstructions.
- Brown invariants of punctured Klein bottles lie in {0, +-2}.  GL Cor. 8 (s = sigma + 4 Arf = Sign W - beta mod 8) then
  allows Sign W in {2, 0, -2} for s = 0, {2, 0} for s = 2, {2, -2} for s = 4 and {0, -2} for s = 6 (ALLOWED).

A sign S is excluded when one of two things holds.
- Linking form: no integral binary form Q of signature S with |det Q| = l (det K / l a square) has
  lambda_K = lambda_Q + metabolic.  For odd p this is checked level by level on Jordan data:
  - the level-v block of lambda_Q is an orthogonal summand of that of lambda_K (smaller rank, or equal rank and
    discriminant);
  - the complement is metabolic (Witt class zero at the odd levels).
  Binary forms are enumerated with a representative of every class:
  - definite: reduced forms, a <= sqrt(4l/3);
  - indefinite, not isotropic: in a basis whose first vector takes the least nonzero |value| |a| and with |2b| <= |a|
    (a shear), |c| >= |a|, so ac < 0 (else l = b^2 - ac < 0) and l = b^2 + |ac| >= a^2: |a| <= sqrt(l);
  - and the isotropic forms when l is a square.
- Donaldson (alternating knots, S = +-2): the closed definite manifold forces the definite Goeritz lattice to embed in
  Z^{rank+2}; S = 2 needs G_pos and S = -2 needs -G_neg (Jabuka-Kelly Thm 2.11, jk_lattice.py conventions).
K bounds no punctured Klein bottle (gamma_4 >= 3) when every allowed sign is excluded (field obstructed).  The field
obstructed_linking_only uses the linking form alone, and obstructed_odd and obstructed_odd_signed use it with odd Q
only, of any signature or of an allowed one.
Modes: controls, the knots whose KnotInfo upper bound is at most 2 (value 1, value 2 or range [1,2]); targets, the
knots whose range starts at 2, and those whose range starts at 1 but is not [1,2] that rank1_linking.py or
rank1_dinv.py obstructs; checks, three knots with gamma_4 = 3: 8_18 (Jabuka-Kelly Prop. 3.1) and 4_1 # 5_1,
4_1 # -5_1 (GL Thm. A).
Usage: python scripts/rank2.py controls | targets | checks [workers]
Output: data/rank2_<mode>.json
"""
import json, os, sys
from math import isqrt
from multiprocessing import Pool
from rank1_linking import jordan, invariants, legendre, signature
from rank1_goeritz import goeritz_pair, det
from jk_lattice import embeds
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, os.path.join(ROOT, 'vendor', 'code', 'lower_bounds', 'owens'))
from owens_obstruction import positive_definite
ALLOWED = {0: (2, 0, -2), 2: (2, 0), 4: (2, -2), 6: (0, -2)}


def binary_forms(l):
    """(Q, signature) with det Q = l (definite) or -l (indefinite): a representative of every class."""
    out = []
    for a in range(1, isqrt(4 * l // 3) + 2):
        for b in range(-(a // 2), a // 2 + 1):
            if (l + b * b) % a == 0 and (l + b * b) // a >= a:
                c = (l + b * b) // a
                out += [([[a, b], [b, c]], 2), ([[-a, -b], [-b, -c]], -2)]
    m = isqrt(l) + 1                                       # |a| <= sqrt(l), see the docstring
    for a in [x for x in range(-m, m + 1) if x]:
        for b in range(-(abs(a) // 2), abs(a) // 2 + 1):
            if (b * b - l) % a == 0:
                out.append(([[a, b], [b, (b * b - l) // a]], 0))
    r = isqrt(l)
    if r * r == l:
        out += [([[0, r], [r, c]], 0) for c in range(2 * r)]
    return out


def levels(J, p):
    L = {}
    for v, u in J:
        r, d = L.get(v, (0, 1))
        L[v] = (r + 1, d * u % p)
    return L


def summand_ok(JA, JQ, p):
    """lambda_A = lambda_Q + metabolic on the p-parts, from Jordan data (valuation >= 1, unit mod p)."""
    LA, LQ = levels(JA, p), levels(JQ, p)
    rank, disc = 0, 1
    for v in set(LA) | set(LQ):
        ra, da = LA.get(v, (0, 1)); rq, dq = LQ.get(v, (0, 1))
        if rq > ra or (rq == ra and legendre(da * dq, p) != 1):
            return False
        if v % 2:
            rank += ra - rq; disc = disc * da * dq % p
    return rank % 2 == 0 and legendre((-1) ** (rank // 2) * disc, p) == 1


def linking_signs(inv, D):
    """Pairs (S, parity) for which some binary Q of signature S passes the linking-form test; parity is 'odd' or
    'even' (Q even iff both diagonal entries are even)."""
    from linking_genus import exceptional_primes
    if exceptional_primes(inv):
        raise ValueError('exceptional prime: the orthogonal-sum test is not a proven obstruction '
                         '([T] Definition 3.11 and Proposition 4.9)')
    ok = set()
    for n in range(1, isqrt(D) + 1):
        if D % (n * n):
            continue
        l = D // (n * n)
        for Q, S in binary_forms(l):
            key = (S, 'odd' if Q[0][0] % 2 or Q[1][1] % 2 else 'even')
            if key in ok:
                continue
            JQ = {p: [(v, u) for v, u in jordan(Q, p) if v > 0] for p in inv}
            if all(summand_ok(inv[p], JQ[p], p) for p in inv):
                ok.add(key)
        if len(ok) == 6:
            break
    return ok


def one(item):
    name, V, pd = item
    n = len(V)
    A = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    D = abs(det(A))
    inv = invariants(A, D) if D > 1 else {}
    s = (signature(A) + 4 * (0 if D % 8 in (1, 7) else 1)) % 8
    pairs = linking_signs(inv, D)
    link = {S for S, _ in pairs}
    odd = {S for S, par in pairs if par == 'odd'}
    don = {}
    if pd is not None:
        Gs = goeritz_pair(pd)
        neg = [G for G in Gs if positive_definite([[-v for v in row] for row in G])]
        pos = [G for G in Gs if positive_definite(G)]
        if len(neg) == 1 and len(pos) == 1:
            for S, L in ((2, pos[0]), (-2, [[-v for v in row] for row in neg[0]])):
                if S in link and S in ALLOWED[s]:
                    don[S] = embeds(L, len(L) + 2)
    excluded = {S: S not in link or don.get(S) is False for S in (2, 0, -2)}
    return name, dict(det=D, s=s, linking_signs=sorted(link), odd_signs=sorted(odd), donaldson=don,
                      obstructed=all(excluded[S] for S in ALLOWED[s]),
                      obstructed_linking_only=all(S not in link for S in ALLOWED[s]),
                      obstructed_odd=not odd, obstructed_odd_signed=all(S not in odd for S in ALLOWED[s]))


def block_sum(*Vs):
    out, k = [], sum(len(V) for V in Vs)
    off = 0
    for V in Vs:
        for row in V:
            out.append([0] * off + list(row) + [0] * (k - off - len(V)))
        off += len(V)
    return out


if __name__ == '__main__':
    from collections import Counter
    from database_knotinfo import link_list
    mode = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    BY = {r['name']: r for r in link_list()[1:]}
    G4 = json.load(open(os.path.join(ROOT, 'data', 'gamma4_knotinfo_2026.10.5.json')))
    sm = lambda k: json.loads(BY[k]['seifert_matrix'])
    pdn = lambda k: json.loads(BY[k]['pd_notation']) if BY[k]['alternating'] == 'Y' else None
    mirror = lambda V: [[-x for x in row] for row in zip(*V)]
    if mode == 'checks':                  # gamma_4 = 3: 8_18 (JK Prop. 3.1), 4_1 # 5_1 and 4_1 # -5_1 (GL Thm. A)
        items = [('8_18', sm('8_18'), pdn('8_18')), ('4_1#5_1', block_sum(sm('4_1'), sm('5_1')), None),
                 ('4_1#-5_1', block_sum(sm('4_1'), mirror(sm('5_1'))), None)]
    else:
        lb2 = set()
        for f in ('rank1_targets.json', 'rank1_dinv_targets.json', 'rank1_dinv_mtargets.json'):
            lb2 |= {k for k, r in json.load(open(os.path.join(ROOT, 'data', f))).items() if r.get('obstructed')}
        if mode == 'controls':            # upper bound at most 2: a punctured Klein bottle is known
            names = [k for k, v in G4.items() if v in ('1', '2', '[1,2]')]
        else:                             # gamma_4 >= 2 by KnotInfo or by the rank-one tests, upper bound at least 3
            names = [k for k, v in G4.items() if v.startswith('[2,') or (v.startswith('[1,') and v != '[1,2]' and k in lb2)]
        items = [(k, sm(k), pdn(k)) for k in names]
    with Pool(workers) as p:
        out = dict(p.imap_unordered(one, items, chunksize=8))
    json.dump(dict(sorted(out.items())), open(os.path.join(ROOT, 'data', f'rank2_{mode}.json'), 'w'), indent=0)
    print(mode, 'knots', len(out), '| obstructed', sum(r['obstructed'] for r in out.values()),
          '| by linking form alone', sum(r['obstructed_linking_only'] for r in out.values()),
          '| no odd Q at all (topological)', sum(r['obstructed_odd'] for r in out.values()),
          '| no odd Q with allowed sign', sum(r['obstructed_odd_signed'] for r in out.values()))
    if mode == 'checks':
        for k, r in out.items():
            print('  ', k, r)
    else:
        print('s:', Counter(r['s'] for r in out.values()), '| linking signs:', Counter(tuple(r['linking_signs']) for r in out.values()).most_common(8))
        if mode == 'targets':
            ob = sorted(k for k, r in out.items() if r['obstructed'])
            print('obstructed by range:', Counter(G4[k] for k in ob), '|', ob[:40])
