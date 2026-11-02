#!/usr/bin/env python
"""[C] Theorem B and [C] Section 3.8: an independent check of the test of [C] Proposition 2.6 for the three Montesinos
knots of Theorem B (12n_25, 12n_198, 12n_550) and, with the argument literature, for 11n_57.

No code of rank1_dinv.py or of vendor/ is used:
- Sigma_2(K) from the KnotInfo PD code by the orbifold route (fill the exterior (2,0), take the 2-fold cyclic cover);
- Regina recognises it as a Seifert fibred space over S^2 (randomised triangulations), so the KnotInfo Montesinos
  notation is not used;
- the negative-definite star plumbing of that space (normal form; it has at most one bad vertex) is built here;
- d(Y, t) = max (K^2 + n)/4 over the characteristic vectors of t (Ozsvath-Szabo, plumbed three-manifolds, Cor. 1.5),
  by an exact closest-vector enumeration (Fincke-Pohst in rational arithmetic, with an exact LDL^T factorisation);
- the rank-one test of Proposition 2.6 for Z = +-Y, for cyclic H_1 (decided by the adjugate, not by the standard basis).
Controls: the non-alternating Montesinos knots with KnotInfo value 1 (data/rank1_dinv_mcontrols.json) whose double
branched cover Regina names as a Seifert fibred space over S^2 with cyclic H_1; none may be obstructed.  The d-invariant
of the spin structure is compared with that of rank1_dinv.py where that has one (data/rank1_dinv_mtargets.json,
data/rank1_dinv_mcontrols.json).
With the argument literature the script checks instead the knots of LITERATURE, for which a published source gives
gamma_4 = 1 although rank1_dinv.py excludes a smooth Mobius band (Section 3.8).
The randomised triangulations of SnapPy and Regina are seeded, so the output is reproducible.  sfs_name, parse,
plumbing, d_invariants and obstructed are also the independent implementation of the test in surgery_dinv.py,
surgery_cfk.py, surgery_drill.py and plumbing_dinv.py.
Usage: python scripts/montesinos_dinv_check.py [literature]
Output: data/montesinos_dinv_check.json; with literature, data/montesinos_dinv_check_literature.json
"""
import json, math, os, re, sys, warnings
from fractions import Fraction as Fr
from collections import Counter
warnings.filterwarnings('ignore')
import snappy, regina
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
LITERATURE = ['11n_57']        # Fairchild, Prop. 3.2, and the table of Lee-Sabloff give gamma_4 = 1


def sfs_name(k, tries=40):
    """Regina's names of Sigma_2(K), with their counts, over at most `tries` triangulations of the cover (the first as
    SnapPy gives it, the others randomised), until five are named."""
    E = snappy.Link(json.loads(BY[k]['pd_notation'])).exterior(); E.dehn_fill((2, 0))
    C = E.covers(2, cover_type='cyclic')[0]
    names = Counter()
    for t_ in range(tries):
        F = C.filled_triangulation()
        if t_: F.randomize()
        t = regina.Triangulation3(F._to_string()); t.intelligentSimplify()
        s = regina.StandardTriangulation.recognise(t)
        if s is not None and s.manifold() is not None:
            names[s.manifold().name()] += 1
        if sum(names.values()) >= 5:
            break
    return names


def parse(name):
    """The Seifert invariants [(a, b), ...] of a Regina name 'SFS [S2: (a,b) ...]', or None for any other name."""
    m = re.fullmatch(r'SFS \[S2: (.*)\]', name)
    return [tuple(map(int, f.strip('()').split(','))) for f in m.group(1).split(') (')] if m else None


def hj(p, q):
    """The Hirzebruch-Jung continued fraction p/q = a_1 - 1/(a_2 - ...), a_i >= 2, for p > q > 0."""
    out = []
    while q:
        a = -(-p // q); out.append(a); p, q = q, a * q - p
    return out


def plumbing(fibres):
    """The star-shaped plumbing of the Seifert fibred space over S^2 with invariants (a_i, b_i), in the orientation with
    Euler number e = sum b_i/a_i <= 0 (all b_i negated if e > 0): the central weight sum floor(b_i/a_i), and for each
    fibre with r_i = b_i - floor(b_i/a_i) a_i > 0 a leg with the weights -c of the expansion a_i/r_i = [c_1, c_2, ...].
    Returns the matrix Q and e."""
    e = sum(Fr(b, a) for a, b in fibres)
    if e > 0:
        fibres = [(a, -b) for a, b in fibres]; e = -e
    e0, legs = 0, []
    for a, b in fibres:
        n = b // a; w = b - n * a; e0 += n
        if w:
            legs.append(hj(a, w))
    verts = [e0]; edges = []
    for leg in legs:
        prev = 0
        for c in leg:
            verts.append(-c); edges.append((prev, len(verts) - 1)); prev = len(verts) - 1
    n = len(verts)
    Q = [[0] * n for _ in range(n)]
    for i, w in enumerate(verts):
        Q[i][i] = w
    for i, j in edges:
        Q[i][j] = Q[j][i] = 1
    return Q, e


def inv(Q):
    """Q^-1 over Q (Gauss-Jordan, exact)."""
    n = len(Q); A = [[Fr(Q[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        piv = A[c][c]; A[c] = [x / piv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def det(Q):
    """det Q (exact)."""
    n = len(Q); A = [[Fr(x) for x in row] for row in Q]; d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return 0
        if p != c: A[c], A[p] = A[p], A[c]; d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return int(d)


def quad(Q, y):
    """y^T Q y."""
    n = len(Q); return sum(y[i] * Q[i][j] * y[j] for i in range(n) for j in range(n))


class NotCyclic(Exception):
    """coker Q is not cyclic."""


def generator(Q, Qi):
    """A generator of coker Q, or NotCyclic.  coker Q is cyclic iff the entries of adj(Q), the (n-1)-minors, are coprime
    (they generate the ideal s_1 ... s_{n-1} of the Smith form).  Then for each prime p | det Q some column j_p of adj(Q)
    is nonzero mod p, and x = e_{j_p} mod p for all p (Chinese remainders) has adj(Q) x nonzero mod every such p, so
    Q^{-1} x = adj(Q) x / det Q has order |det Q|."""
    n = len(Q); D = det(Q); N = abs(int(D))
    adj = [[int(Qi[i][j] * D) for j in range(n)] for i in range(n)]
    if math.gcd(*[a for row in adj for a in row]) != 1:
        raise NotCyclic
    primes, m = [], N
    for q in range(2, N + 1):
        if q * q > m:
            break
        if m % q == 0:
            primes.append(q)
            while m % q == 0:
                m //= q
    if m > 1:
        primes.append(m)
    x, mod = [0] * n, 1
    for p in primes:
        j = next(j for j in range(n) if any(adj[i][j] % p for i in range(n)))
        for i in range(n):
            x[i] += mod * ((int(i == j) - x[i]) * pow(mod, -1, p) % p)
        mod *= p
    return x


def ldl(P):
    """P = sum_i d_i (x_i + sum_{j>i} mu[i][j] x_j)^2 exactly, for P symmetric positive definite."""
    n = len(P); A = [[Fr(x) for x in row] for row in P]
    d, mu = [Fr(0)] * n, [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        d[i] = A[i][i]
        assert d[i] > 0
        for j in range(i + 1, n):
            mu[i][j] = A[i][j] / d[i]
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j][k] -= d[i] * mu[i][j] * mu[i][k]
    return d, mu


def cv_min(P, y0):
    """min over integer z of (y0 + 2z)^T P (y0 + 2z), exactly.  With c = -y0/2 this is 4 q(z - c) for q(x) = x^T P x, and
    q(z - c) = sum_i d_i (z_i - center_i)^2 with center_i = c_i - sum_{j>i} mu_ij (z_j - c_j).  Fincke-Pohst enumeration
    from the last coordinate down, in rational arithmetic: every z with q(z - c) <= best/4 is visited."""
    n = len(P); d, mu = ldl(P)
    c = [-Fr(v) / 2 for v in y0]
    value = lambda z: quad(P, [y0[i] + 2 * z[i] for i in range(n)])
    best = [value([round(x) for x in c])]
    z = [0] * n
    def rec(i, partial):
        center = c[i] - sum(mu[i][j] * (z[j] - c[j]) for j in range(i + 1, n))
        rem = best[0] / 4 - partial
        if rem < 0:
            return
        r2 = rem / d[i]                                        # (z_i - center)^2 <= r2
        s = math.isqrt(r2.numerator // r2.denominator) + 1     # s > sqrt(r2)
        for zi in range(math.floor(center - s), math.ceil(center + s) + 1):
            t = zi - center
            val = partial + d[i] * t * t
            if val > best[0] / 4:
                continue
            z[i] = zi
            if i == 0:
                v = value(z)
                if v < best[0]:
                    best[0] = v
            else:
                rec(i - 1, val)
    rec(n - 1, Fr(0))
    return best[0]


def d_invariants(Q):
    """For the boundary Y of the negative-definite plumbing Q with cyclic coker Q (otherwise NotCyclic): N = |det Q|,
    u with lambda_Y(x0, x0) = u/N for the generator x0 of generator(), and d with d[x] = d(Y, t_s + x x0), the maximum
    of (K^2 + n)/4, K^2 = K^T Q^-1 K, over the characteristic covectors K with [K] = 2 x x0 in coker Q."""
    n = len(Q); Qi = inv(Q); N = abs(det(Q))
    P = [[-x for x in row] for row in Q]                       # P positive definite
    order = lambda v: math.lcm(*[Fr(sum(Qi[i][j] * v[j] for j in range(n))).denominator for i in range(n)])
    x0 = generator(Q, Qi)                                     # raises NotCyclic
    assert order(x0) == N
    in_im = lambda v: all(sum(Qi[i][j] * v[j] for j in range(n)).denominator == 1 for i in range(n))
    Kb = [Q[i][i] % 2 for i in range(n)]
    c = next(c for c in range(N) if in_im([Kb[i] - c * x0[i] for i in range(n)]))
    u = Fr(-quad(Qi, x0)) % 1 * N                               # lambda_Y(x0, x0) = u / N, lambda = -Q^-1
    d = {}
    for j in range(N):
        K = [Kb[i] + 2 * j * x0[i] for i in range(n)]
        y0 = [sum(Qi[i][t] * K[t] for t in range(n)) for i in range(n)]
        best = cv_min(P, y0)
        label = (c * (N + 1) // 2 + j) % N                     # x(t) = [K]/2 in Z/N, coordinate of x0 = 1
        d[label] = (-best + n) / 4
    return N, int(u), d


def obstructed(N, u, d):
    """The test of [C] Proposition 2.6 for H_1(Y) = Z/N, lambda_Y(x0, x0) = u/N and d[x] = d(Y, t_s + x x0): for
    Z0 = eps Y, each l = N/m^2 (M = <N/m>, of order m, is isotropic; P = M^perp = <m>) and each g = m s with
    lambda_Z0(g, g) = 1/l, it checks 4 d(Z0, t_s + j g + y) >= 1 - (2j - l)^2/l for 0 <= j < l and y in M.  Returns
    (no (eps, l, s) passes, the list of those that pass)."""
    out = []
    for eps in (1, -1):                                        # Z0 = Y or -Y
        for m in range(1, math.isqrt(N) + 1):
            if N % (m * m): continue
            l = N // (m * m)
            for s in range(1, l + 1):
                if math.gcd(s, l) != 1 or (eps * u * s * s - 1) % l: continue
                g = (m * s) % N
                ok = all(4 * eps * d[(j * g + (N // m) * t) % N] >= 1 - Fr((2 * j - l) ** 2, l)
                         for j in range(l) for t in range(m))
                if ok:
                    out.append((eps, l, s))
    return not out, out


if __name__ == '__main__':
    snappy.set_rand_seed(0); regina.RandomEngine.reseedWithDefault()   # the randomised triangulations are reproducible
    ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    literature = sys.argv[1:] == ['literature']
    targets = LITERATURE if literature else ['12n_25', '12n_198', '12n_550']
    controls = [] if literature else sorted(json.load(open(os.path.join(ROOT, 'data', 'rank1_dinv_mcontrols.json'))))
    D = {}
    for f in ('rank1_dinv_mtargets.json', 'rank1_dinv_mcontrols.json'):
        D.update(json.load(open(os.path.join(ROOT, 'data', f))))
    out, skipped = {}, {}
    for k in targets + controls:
        names = sfs_name(k)
        fib = parse(names.most_common(1)[0][0]) if names else None
        if not fib:
            skipped[k] = dict(names); continue
        Q, e = plumbing(fib)
        try:
            N, u, d = d_invariants(Q)
        except NotCyclic:
            skipped[k] = 'H_1 not cyclic'; continue
        assert N == int(BY[k]['determinant']), k
        ob, passing = obstructed(N, u, d)
        out[k] = dict(role='target' if k in targets else 'control', regina=names.most_common(1)[0][0], plumbing_rank=len(Q),
                      det=N, u=u, d={str(x): str(v) for x, v in sorted(d.items())}, d_spin=str(d[0]),
                      d_spin_rank1_dinv=D.get(k, {}).get('d_spin'), obstructed=ob,
                      obstructed_rank1_dinv=D.get(k, {}).get('obstructed'), passing=[list(p) for p in passing[:4]])
        print(k, out[k]['role'], out[k]['regina'], 'obstructed', ob, 'd_spin', d[0], flush=True)
    summary = dict(targets_obstructed=all(out[k]['obstructed'] for k in targets),
                   controls_checked=sum(v['role'] == 'control' for v in out.values()),
                   controls_obstructed=[k for k, v in out.items() if v['role'] == 'control' and v['obstructed']],
                   d_spin_disagreements=[k for k, v in out.items() if v['d_spin_rank1_dinv'] is not None
                                         and Fr(v['d_spin']) != Fr(v['d_spin_rank1_dinv'])],
                   skipped=len(skipped))
    json.dump(dict(summary=summary, knots=out, skipped=skipped),
              open(os.path.join(ROOT, 'data', 'montesinos_dinv_check_literature.json' if literature
                                else 'montesinos_dinv_check.json'), 'w'), indent=1)
    print(summary)
