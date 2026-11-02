#!/usr/bin/env python
"""Independent check of the certificates of scripts/donaldson_surfaces.py ([C] Proposition 3.2(c), Checks).

No code of the search (donaldson_surfaces.py, rank1_goeritz.py, jk_lattice.py, owens_obstruction.py) is used.  For each
certificate (a diagram of K, a colour and the lattice L = -eps G_V, eps = +1 for s = 2 and -1 for s = 6):
1. s: sigma(K) as the signature of V + V^T of the KnotInfo Seifert matrix (exact, LDL over Q) and Arf(K) from det K
   mod 8 (Levine-Murasugi), s = sigma + 4 Arf mod 8.
2. The diagram is a diagram of K: SnapPy finds an isometry from its exterior to that of the KnotInfo PD code that
   extends to the knots and whose cusp map has determinant +1 (so the diagram is not that of the mirror image).
3. Goeritz matrices: the faces of the PD code are traced again, coloured, and the matrix of Gordon and Litherland is
   formed for each colour (incidence number eta = +1 at a crossing X[a,b,c,d] for the colour of the corners between a
   and b and between c and d, -1 for the other colour; the convention is checked on 3_1, whose Seifert surface has
   signature -2); -eps times the matrix of the colour must be isometric to L (PARI qfisom).
4. L is positive definite (LDL over Q) and |det L| = det K.
5. The linking form presented by L is that presented by -eps (V + V^T) (Jordan data of rank1_linking.invariants at
   every odd prime; this compares the Goeritz matrix with the Seifert matrix, so it also checks the sign); the record
   says whether the opposite sign would also pass.
6. L does not embed in Z^{rank L + 1}: a search over all integer vectors of the required norms, rows placed in an order
   with many earlier neighbours, and partial matrices kept only in the canonical form under permutations and sign
   changes of columns (first nonzero entry of each column positive, columns in non-increasing lexicographic order).
Also recorded: for every knot searched, whether its KnotInfo diagram has a checkerboard surface with a definite form
(summary field knotinfo_diagram_definite).
Usage: python scripts/donaldson_check.py [workers]
Output: data/donaldson_check.json
"""
import json, os, sys, time, warnings
from fractions import Fraction
from math import isqrt, prod
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from rank1_linking import invariants, legendre

BY = None


def init():
    global BY
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}


# ---------- exact linear algebra

def ldl_signs(A):
    """(positive, negative, zero) counts of the diagonal of an LDL^T decomposition over Q with symmetric pivoting."""
    M = [[Fraction(x) for x in row] for row in A]
    n, pos, neg = len(M), 0, 0
    rows = list(range(n))
    while rows:
        p = next((i for i in rows if M[i][i] != 0), None)
        if p is None:
            i = rows[0]; j = next((j for j in rows if M[i][j] != 0), None)
            if j is None:
                rows.remove(i); continue
            for k in range(n):                             # replace row/column i by i + j
                M[i][k] += M[j][k]
            for k in range(n):
                M[k][i] += M[k][j]
            continue
        d = M[p][p]
        pos += d > 0; neg += d < 0
        rows.remove(p)
        for i in rows:
            f = M[i][p] / d
            for j in rows:
                M[i][j] -= f * M[p][j]
    return pos, neg, n - pos - neg


def positive_definite(A):
    return ldl_signs(A) == (len(A), 0, 0)


def definite(A):
    """Whether A is nonempty and positive or negative definite."""
    pos, neg, zero = ldl_signs(A)
    return bool(A) and zero == 0 and (pos == 0 or neg == 0)


def determinant(A):
    M = [[Fraction(x) for x in row] for row in A]
    n, d = len(M), Fraction(1)
    for k in range(n):
        p = next((r for r in range(k, n) if M[r][k] != 0), None)
        if p is None:
            return 0
        if p != k:
            M[k], M[p] = M[p], M[k]; d = -d
        d *= M[k][k]
        for r in range(k + 1, n):
            f = M[r][k] / M[k][k]
            for c in range(k, n):
                M[r][c] -= f * M[k][c]
    return int(d)


# ---------- Goeritz matrices of Gordon and Litherland

def goeritz(pd):
    """The matrices of the two colours of a PD code (slots counterclockwise from the incoming under-strand), each on the
    regions of that colour with the first region removed, as [colour of the corners (a,b) and (c,d) of the first
    crossing, the other colour]."""
    ends = {}
    for c, X in enumerate(pd):
        for i, e in enumerate(X):
            ends.setdefault(e, []).append((c, i))
    other = {}
    for e, v in ends.items():
        assert len(v) == 2
        other[v[0]], other[v[1]] = v[1], v[0]
    region, nreg = {}, 0
    for c in range(len(pd)):                              # region of the corner between slots i and i + 1 of c
        for i in range(4):
            if (c, i) in region:
                continue
            x = (c, i)
            while x not in region:
                region[x] = nreg
                d, j = other[(x[0], (x[1] + 1) % 4)]       # leave along slot i + 1; continue at the far end
                x = (d, j)
            nreg += 1
    assert nreg == len(pd) + 2
    col = {region[(0, 0)]: 0}                             # opposite corners: same colour; adjacent corners: different
    changed = True
    while changed:
        changed = False
        for c in range(len(pd)):
            for i in range(4):
                a, b, o = region[(c, i)], region[(c, (i + 1) % 4)], region[(c, (i + 2) % 4)]
                for x, y, same in ((a, o, True), (a, b, False)):
                    if x in col and y not in col:
                        col[y] = col[x] if same else 1 - col[x]; changed = True
                    elif y in col and x not in col:
                        col[x] = col[y] if same else 1 - col[y]; changed = True
                    elif x in col and y in col:
                        assert (col[x] == col[y]) == same, 'no checkerboard colouring'
    assert len(col) == nreg
    out = []
    for colour in (0, 1):
        regs = sorted(r for r in range(nreg) if col[r] == colour)
        ix = {r: k for k, r in enumerate(regs)}
        g = [[0] * len(regs) for _ in regs]
        for c in range(len(pd)):
            first = 0 if col[region[(c, 0)]] == colour else 1     # the colour sits at the corners first, first + 2
            eta = 1 if first == 0 else -1
            a, b = region[(c, first)], region[(c, first + 2)]
            if a != b:
                g[ix[a]][ix[b]] -= eta; g[ix[b]][ix[a]] -= eta
        for i in range(len(regs)):
            g[i][i] = -sum(g[i][j] for j in range(len(regs)) if j != i)
        out.append([row[1:] for row in g[1:]])
    return out


def isometric(A, B):
    """Whether the positive definite forms A and B are isometric (PARI qfisom); False if B is not positive definite."""
    from cypari import pari
    if len(A) != len(B) or not positive_definite(B):
        return False
    if not A:
        return True
    res = pari.qfisom(pari.matrix(len(A), len(A), [x for r in A for x in r]),
                      pari.matrix(len(B), len(B), [x for r in B for x in r]))
    return str(pari.type(res)) == 't_MAT'


# ---------- the embedding search

def vectors(norm, N):
    out = []

    def rec(i, rem, v):
        if i == N:
            if rem == 0:
                out.append(tuple(v))
            return
        r = isqrt(rem)
        for x in range(-r, r + 1):
            v.append(x); rec(i + 1, rem - x * x, v); v.pop()
    rec(0, norm, [])
    return out


def embeds(L):
    n = len(L); N = n + 1
    order = [max(range(n), key=lambda i: L[i][i])]
    while len(order) < n:
        rest = [i for i in range(n) if i not in order]
        order.append(max(rest, key=lambda i: (sum(L[i][j] != 0 for j in order), L[i][i])))
    L = [[L[a][b] for b in order] for a in order]
    V = {a: vectors(a, N) for a in {L[i][i] for i in range(n)}}
    rows = []

    def canonical(v):
        for c in range(N):
            if v[c] < 0 and all(r[c] == 0 for r in rows):
                return False
            if c + 1 < N and v[c] < v[c + 1] and all(r[c] == r[c + 1] for r in rows):
                return False
        return True

    def place(i):
        if i == n:
            return True
        for v in V[L[i][i]]:
            if all(sum(a * b for a, b in zip(v, rows[j])) == L[i][j] for j in range(i)) and canonical(v):
                rows.append(v)
                if place(i + 1):
                    return True
                rows.pop()
        return False
    return place(0)


# ---------- one certificate

def linking_class(A, D):
    """For each odd p | D and level v of the Jordan data of A: (rank, Legendre symbol of the discriminant)."""
    out = {}
    for p, J in invariants(A, D).items():
        lv = {}
        for v, u in J:
            lv.setdefault(v, []).append(u)
        out[p] = {v: (len(us), legendre(prod(us) % p, p)) for v, us in lv.items()}
    return out


def chirality(pd, pd0, tries=8):
    import snappy, spherogram            # spherogram's exterior() needs snappy imported
    A = spherogram.Link([list(x) for x in pd]).exterior()
    B = spherogram.Link([list(x) for x in pd0]).exterior()
    for _ in range(tries):
        try:
            return sorted({int(m.cusp_maps()[0].det()) for m in A.is_isometric_to(B, return_isometries=True)
                           if m.extends_to_link()})
        except Exception:
            A.randomize(); B.randomize()
    return None


def one(item):
    name, cert = item
    t0 = time.time()
    r = BY[name]
    D = int(r['determinant'])
    V = json.loads(r['seifert_matrix'])
    A = [[V[i][j] + V[j][i] for j in range(len(V))] for i in range(len(V))]
    pos, neg, _ = ldl_signs(A)
    sigma = pos - neg
    arf = 0 if D % 8 in (1, 7) else 1
    s = (sigma + 4 * arf) % 8
    eps = {2: 1, 6: -1}.get(s)
    pd, L = [tuple(x) for x in cert['pd']], cert['lattice']
    out = dict(knot=name, det=D, sigma=sigma, sigma_knotinfo=int(r['signature']), s=s, rank=len(L))
    out['s_ok'] = eps is not None and sigma == int(r['signature'])
    out['chirality'] = chirality(pd, json.loads(r['pd_notation']))
    out['diagram_ok'] = out['chirality'] == [1]
    G = goeritz(pd)
    out['goeritz_isometric'] = [isometric(L, [[-eps * x for x in row] for row in M]) for M in G] if eps else None
    out['goeritz_ok'] = bool(out['goeritz_isometric']) and sum(out['goeritz_isometric']) >= 1
    out['definite'] = positive_definite(L)
    out['det_ok'] = abs(determinant(L)) == D
    if eps:
        target = linking_class([[-eps * x for x in row] for row in A], D)
        opposite = linking_class([[eps * x for x in row] for row in A], D)
        own = linking_class(L, D)
        out['linking_form_ok'] = own == target
        out['linking_form_sign_decides'] = target != opposite
    out['embeds'] = embeds(L)
    out['passes'] = bool(out['s_ok'] and out['diagram_ok'] and out['goeritz_ok'] and out['definite'] and out['det_ok']
                         and out.get('linking_form_ok') and out['embeds'] is False)
    out['seconds'] = round(time.time() - t0, 1)
    return out


def convention_check():
    """On the KnotInfo PD code of 3_1 (sigma = -2): the matrix of rank 2 belongs to the Seifert surface (genus one) and
    must be negative definite, the one of rank 1 is that of the Mobius band; and the embedding search on small lattices."""
    G = goeritz([tuple(x) for x in json.loads(BY['3_1']['pd_notation'])])
    two = [M for M in G if len(M) == 2]
    tests = [([[2, 1], [1, 2]], True), ([[3]], False), ([[5]], True), ([[7]], False), ([[2, 1, 0], [1, 2, 1], [0, 1, 2]], True)]
    return dict(matrices=G, sigma_3_1=int(BY['3_1']['signature']),
                ok=int(BY['3_1']['signature']) == -2 and len(two) == 1 and ldl_signs(two[0]) == (0, 2, 0)
                and all(embeds(A) == e for A, e in tests))


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    init()
    conv = convention_check()
    assert conv['ok'], conv
    T = json.load(open(os.path.join(ROOT, 'data', 'donaldson_surfaces.json')))['targets']['knots']
    items = [(k, r['certificate']) for k, r in sorted(T.items()) if r['certificate']]
    minimal = sorted(k for k in T
                     if any(definite(M) for M in goeritz([tuple(x) for x in json.loads(BY[k]['pd_notation'])])))
    t0 = time.time()
    with Pool(workers, initializer=init) as P:
        res = P.map(one, items)
    recs = {r['knot']: r for r in res}
    summary = dict(certificates=len(recs), passes=sorted(k for k, r in recs.items() if r['passes']),
                   failing=sorted(k for k, r in recs.items() if not r['passes']),
                   linking_form_sign_decides=sum(bool(r.get('linking_form_sign_decides')) for r in recs.values()),
                   searched=len(T), knotinfo_diagram_definite=minimal)
    import snappy, cypari
    meta = dict(convention_check=conv, seconds=round(time.time() - t0),
                versions=dict(snappy=snappy.__version__, cypari=cypari.__version__))
    json.dump(dict(meta=meta, summary=summary, knots=recs), open(os.path.join(ROOT, 'data', 'donaldson_check.json'), 'w'),
              indent=1)
    print(json.dumps({k: (v if not isinstance(v, list) or len(v) < 8 else len(v)) for k, v in summary.items()}))
