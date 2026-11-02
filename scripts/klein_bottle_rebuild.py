#!/usr/bin/env python
"""Rebuild of the surfaces of the first search of [C] Section 3.3 (Proposition 3.3 and Section 3.8), stored by
scripts/klein_bottle_search.py, from the KnotInfo PD codes.  It shares no code with the search.

Certificates.  Each is a chain of face bands that ends in the unknot or the two-component unlink:
  K -> J -> unknot                  moebius_band_of_J, kind 'band to the unknot'          b_1 = 2
  K -> J -> J' -> unlink            moebius_band_of_J, kind 'band to a ribbon knot'       b_1 = 2
  K -> J -> J' -> unknot            klein_bottle_of_J (J with KnotInfo gamma_4 2)         b_1 = 3
  K -> unknot, K -> J' -> unlink    moebius_band_of_K (group upper of mode ranges,        b_1 = 1
                                    and group mobius of mode mobius)
Every band but the ribbon band at the end joins a knot to a knot, so the first band makes the surface (the bands at
successive levels of a collar of S^3, and a disc in B^4 for each component of the final unlink) nonorientable, with
one boundary component, chi = (discs) - (bands) and b_1 = 1 - chi.  A chain is self-contained when its surface of J
lies in the diagram of J made by the first band; the chain then starts from the KnotInfo PD code of K
(database_knotinfo pd_notation), which must be the stored diagram of K, or, for a further diagram of the search (stage
B), from the stored diagram, whose exterior must be that of the KnotInfo diagram (SnapPy is_isometric_to; by
Gordon-Luecke the diagram is then K or its mirror image, and gamma_4 is mirror-invariant).  For such a diagram the
script also records whether Reidemeister moves of spherogram (simplification, or random moves followed by
simplification) turn it into the KnotInfo diagram or its reflection, or both into a common diagram (reaches, meets),
which links it to the KnotInfo PD code without SnapPy.  The other certificates (surface of J in the KnotInfo diagram of
J, or in another stored diagram of J) are rebuilt in two parts and reported apart; for them the script records whether
Reidemeister moves turn the diagram of J made by the band into the diagram of the surface of J (or its reflection).  A
certificate is diagrammatic when all of it is linked to the KnotInfo PD code of K by band moves and Reidemeister moves
found here (no SnapPy, no KnotInfo value of gamma_4).
Rebuild.  Nothing of scripts/bands.py or scripts/klein_bottle_search.py is used.  A band is given by the
number of its face, the two positions i, j on the face and the number t of half twists (t > 0: the strand from SW to
NE is over, the band drawn from W to E; [C] Figure 1).  In klein_bottle_search a face is the list of the ends
(crossing, slot) it leaves from, face on the left, from its least end, and the faces are numbered by their least ends.
This script traverses the same faces in reverse (after arriving at slot k of a crossing it leaves by slot k + 1, slots
counterclockwise), so it recovers the numbering; the edges of the face and the two edges of the band must agree with
those stored.  It finds the ends A, B of edge i and C, D of edge j (A, B and C, D in the order of the face).  The band
move joins A to D and B to C (t = 0), or runs them through |t| new crossings.  It is done twice: by spherogram's band
code (links/bands/core.add_one_band of Dunfield-Gong, crossing strands A and D, num_twist = -t; it checks planarity),
and by this script (planarity checked by V - E + F = 2 for every connected piece).  After every band the two diagrams
must be isomorphic (PD codes equal up to renumbering crossings and edges, each crossing turned by 0 or 2 slots: an
isomorphism of oriented planar 4-valent maps keeping over and under, found by this script), isomorphic to the diagram
stored in the certificate when there is one, and have the stated number of components.  The final diagram (both
versions) must simplify to no crossings by spherogram's Reidemeister simplification alone (simplify('global'), at
most TRIES attempts), leaving one unknotted component or two unlinked ones; knot Floer homology is not used.
Usage: python scripts/klein_bottle_rebuild.py [workers] [search_output.json ...]
The default search outputs are data/klein_bottle_search_ranges.json and data/klein_bottle_search.json; the output of
mode mobius is rebuilt alone: python scripts/klein_bottle_rebuild.py 10 data/klein_bottle_search_mobius.json.
Environment, for tests: REBUILD_OUT=path replaces the output file.
Output: data/klein_bottle_rebuild.json, for data/klein_bottle_search_mobius.json data/klein_bottle_rebuild_mobius.json
"""
import json, os, sys, time, warnings
from collections import Counter, defaultdict
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
warnings.filterwarnings('ignore')
TRIES = 20                                # attempts of simplify('global') on a final diagram
REACH_TRIES = 20000                       # attempts to turn one stored diagram into another (reaches)


# ---------- planar diagrams: PD codes (KnotTheory convention: slots counterclockwise from the incoming under-strand)

def ends(P):
    """Edge label -> its two ends (crossing, slot)."""
    o = defaultdict(list)
    for c, X in enumerate(P):
        for k, e in enumerate(X):
            o[e].append((c, k))
    assert all(len(v) == 2 for v in o.values()), 'not a PD code'
    return o


def joins(P):
    """End (crossing, slot) -> the other end of its edge."""
    nb = {}
    for p, q in ends(P).values():
        nb[p] = q; nb[q] = p
    return nb


def faces(nb, ncross):
    """Faces as lists of steps (departure, arrival): after arriving at (c, k), leave by (c, k + 1)."""
    seen, out = set(), []
    for c in range(ncross):
        for k in range(4):
            if (c, k) in seen:
                continue
            steps, p = [], (c, k)
            while p not in seen:
                seen.add(p)
                dep = (p[0], (p[1] + 1) % 4)
                p = nb[dep]
                steps.append((dep, p))
            assert p == (c, k)
            out.append(steps)
    return out


def pieces(nb, ncross):
    parent = list(range(ncross))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (c, _), (d, _) in nb.items():
        parent[find(c)] = find(d)
    return len({find(c) for c in range(ncross)})


def planar(nb, ncross):
    return ncross - len(nb) // 2 + len(faces(nb, ncross)) == 2 * pieces(nb, ncross)


def face_by_index(P, index):
    """The edges of face number index in the numbering of klein_bottle_search: there a face is the list of the ends it
    leaves from (here: arrives at), in its direction (here: reversed), from its least end (crossing, slot); faces are
    numbered by their least end."""
    nb = joins(P)
    out = []
    for steps in faces(nb, len(P)):
        darts = [a for _, a in steps][::-1]
        r = darts.index(min(darts))
        darts = darts[r:] + darts[:r]
        out.append((darts[0], [P[c][k] for c, k in darts]))
    return sorted(out)[index][1]


def locate(P, face_edges, darts):
    """Ends A, B of the edge face_edges[i] and C, D of face_edges[j], in the order of the face of klein_bottle_search
    (this traversal runs the other way: it arrives at A and C and leaves from B and D)."""
    nb = joins(P)
    assert len(set(face_edges)) == len(face_edges)
    found = []
    for steps in faces(nb, len(P)):
        labels = [P[d[0]][d[1]] for d, _ in steps]
        if sorted(labels) == sorted(face_edges):
            found.append((labels, steps))
    assert len(found) == 1, 'face not found, or not unique'
    labels, steps = found[0]
    labels, steps = labels[::-1], steps[::-1]
    r = labels.index(face_edges[0])
    assert labels[r:] + labels[:r] == list(face_edges), 'the face has another cyclic order'
    at = {lab: st for lab, st in zip(labels, steps)}
    (b, a), (d, c) = at[face_edges[darts[0]]], at[face_edges[darts[1]]]
    return a, b, c, d


def band_move(P, A, B, C, D, t):
    """This script's band move (see the docstring): PD code of the result and its number of components."""
    n = len(P)
    nb = joins(P)
    for p in (A, B, C, D):
        del nb[p]

    def join(p, q):
        nb[p] = q; nb[q] = p
    under = {c: 0 for c in range(n)}      # the under-strand occupies slots under[c] and under[c] + 2
    if t == 0:
        join(A, D); join(B, C)
    else:                                 # new crossings, slots 0 NE, 1 NW, 2 SW, 3 SE (counterclockwise)
        xs = list(range(n, n + abs(t)))
        for x in xs:
            under[x] = 1 if t > 0 else 0  # t > 0: NW-SE under, so the strand from SW to NE is over
        join(A, (xs[0], 1)); join(B, (xs[0], 2))
        for x, y in zip(xs, xs[1:]):
            join((x, 0), (y, 1)); join((x, 3), (y, 2))
        join((xs[-1], 0), D); join((xs[-1], 3), C)
    m = n + abs(t)
    assert len(nb) == 4 * m and planar(nb, m), 'not planar'
    lab, arriving, seen, k, comps = {}, set(), set(), 0, 0
    for c in range(m):
        for s in range(4):
            if (c, s) in seen:
                continue
            comps += 1
            p = (c, s)
            while p not in seen:
                q = nb[p]
                seen.add(p); seen.add(q)
                k += 1; lab[p] = lab[q] = k; arriving.add(q)
                p = (q[0], (q[1] + 2) % 4)
    out = []
    for c in range(m):
        u = under[c]
        s0 = u if (c, u) in arriving else u + 2
        assert ((c, s0) in arriving) and ((c, (s0 + 2) % 4) not in arriving)
        out.append(tuple(lab[(c, (s0 + j) % 4)] for j in range(4)))
    return out, comps


def isomorphism(P, Q):
    """Ends of P -> ends of Q for an isomorphism of the diagrams (see the docstring), or None.  A connected piece is
    determined by the image of one crossing and its turn; the pieces are matched by backtracking."""
    if len(P) != len(Q):
        return None
    nP, nQ = joins(P), joins(Q)

    def grow(c, c2, r, m, rot, used):
        m[c], rot[c] = c2, r; used.add(c2)
        stack = [c]
        while stack:
            c = stack.pop()
            for k in range(4):
                d, l = nP[(c, k)]
                d2, l2 = nQ[(m[c], (k + rot[c]) % 4)]
                r = (l2 - l) % 4
                if d in m:
                    if (m[d], rot[d]) != (d2, r):
                        return False
                elif r not in (0, 2) or d2 in used:
                    return False
                else:
                    m[d], rot[d] = d2, r; used.add(d2); stack.append(d)
        return True

    def solve(m, rot, used):
        free = [c for c in range(len(P)) if c not in m]
        if not free:
            return m, rot
        for c2 in range(len(Q)):
            for r in (0, 2):
                if c2 not in used:
                    m2, rot2, used2 = dict(m), dict(rot), set(used)
                    if grow(free[0], c2, r, m2, rot2, used2):
                        found = solve(m2, rot2, used2)
                        if found:
                            return found
        return None
    found = solve({}, {}, set())
    if found is None:
        return None
    m, rot = found
    return {(c, k): (m[c], (k + rot[c]) % 4) for c in range(len(P)) for k in range(4)}


# ---------- chains

def simplifies(L, ncomp):
    """Attempts of simplify('global') until no crossing is left with ncomp unlinked unknotted components (or None)."""
    for i in range(1, TRIES + 1):
        S = L.copy(); S.simplify('global')
        if not S.crossings and S.unlinked_unknot_components == ncomp:
            return i
    return None


def run_chain(start, steps):
    """Apply the bands of steps to the diagram start, twice (spherogram and this script); see the docstring."""
    import snappy, spherogram
    from spherogram.links.bands.core import add_one_band
    rebuilt = [tuple(x) for x in start]
    sph = [tuple(x) for x in start]
    nbm = joins(rebuilt)
    rec = dict(start_crossings=len(start), start_planar=planar(nbm, len(rebuilt)), steps=[])
    ok = rec['start_planar'] and spherogram.Link([list(x) for x in start]).is_planar()
    stored = [tuple(x) for x in start]
    for st in steps:
        fe = face_by_index(stored, st['face'])               # the face from its number; stored edges must agree
        consistent = (st.get('face_edges') in (None, fe)) and [fe[i] for i in st['darts']] == st['edges']
        a, b, c, d = locate(stored, fe, st['darts'])
        f_rebuilt, f_sph = isomorphism(stored, rebuilt), isomorphism(stored, sph)
        if f_rebuilt is None or f_sph is None:
            rec['steps'].append(dict(located=False)); return dict(rec, ok=False)
        rebuilt, comps = band_move(rebuilt, *(f_rebuilt[p] for p in (a, b, c, d)), st['t'])
        L = add_one_band(spherogram.Link([list(x) for x in sph]), ([f_sph[a], f_sph[d]], [], -st['t']))
        sph = [tuple(x) for x in L.PD_code()]
        s = dict(t=st['t'], edges=st['edges'], face_consistent=consistent, crossings=len(rebuilt), components_rebuilt=comps,
                 components_spherogram=len(L.link_components), spherogram_equals_rebuilt=isomorphism(rebuilt, sph) is not None)
        s['components_ok'] = s['components_rebuilt'] == s['components_spherogram'] == st['components']
        if st.get('after') is not None:
            s['stored_equals_rebuilt'] = isomorphism([tuple(x) for x in st['after']], rebuilt) is not None
            stored = [tuple(x) for x in st['after']]
        if st.get('crossings') is not None:
            s['stored_crossings_ok'] = st['crossings'] == len(rebuilt)
        ok = ok and consistent and s['components_ok'] and s['spherogram_equals_rebuilt'] \
            and s.get('stored_equals_rebuilt', True) and s.get('stored_crossings_ok', True)
        rec['steps'].append(s)
        last = L
    if steps[-1].get('final'):
        n = steps[-1]['components']
        rec['simplify_tries_spherogram'] = simplifies(last, n)
        rec['simplify_tries_rebuilt'] = simplifies(spherogram.Link([list(x) for x in rebuilt]), n)
        ok = ok and rec['simplify_tries_spherogram'] is not None and rec['simplify_tries_rebuilt'] is not None
        rec['discs'] = n
        rec['bands'] = len(steps)
    return dict(rec, ok=ok)


def step(r, twist, **kw):
    """A band of a certificate: face number, positions on the face, its edges (and the face's edges if stored)."""
    return dict(face=r['face'], darts=r['darts'], edges=r['edges'], face_edges=r.get('face_edges'), t=r[twist], **kw)


def surface_steps(s):
    """Steps of the surface of J (or of K, for moebius_band_of_K)."""
    if 'unknotting_band_of_J1' in s:                         # klein_bottle_of_J
        u = s['unknotting_band_of_J1']
        return [step(s, 't', after=s['J1_pd'], components=1),
                step(u, 's', components=1, crossings=u['crossings'], final=True)]
    if 'ribbon_band_of_J1' in s:                             # band to a ribbon knot
        r = s['ribbon_band_of_J1']
        return [step(s, 's', after=s['J1_pd'], components=1),
                step(r, 'r', components=2, crossings=r['crossings'], final=True)]
    return [step(s, 's', components=1, crossings=s['crossings'], final=True)]


def job(item):
    """One certificate: (file, knot, index or 'K', record, KnotInfo PD of K, KnotInfo PD of J or None, representative)."""
    import snappy, spherogram
    src, name, idx, c, pd_K, pd_J, is_rep = item
    t0 = time.time()
    out = dict(file=src, knot=name, certificate=idx, diagram_source=c['diagram_source'])
    if c['diagram_source'] == 'knotinfo':                    # start from KnotInfo's PD code itself
        out['K_is_knotinfo_pd'] = [list(x) for x in c['diagram_pd']] == pd_K
        start = pd_K
    else:
        start = c['diagram_pd']
        E = spherogram.Link([list(x) for x in start]).exterior()
        E0 = spherogram.Link(pd_K).exterior()
        out['K_isometric_to_knotinfo'] = isometric(E, E0)
        out['K_simplifies_to'] = reaches(start, pd_K) if is_rep else None     # the claimed certificate only
    if idx == 'K':                                           # moebius_band_of_K
        chains = [('K', start, surface_steps(c))]
        out['self_contained'] = True
    else:
        out['J'] = c['J']
        first = step(c, 't', after=c['J_pd'], components=1)
        s = c.get('moebius_band_of_J') or c.get('klein_bottle_of_J')
        out['surface'] = None if not s else 'klein_bottle_of_J' if 'unknotting_band_of_J1' in s else s['kind']
        out['self_contained'] = bool(s) and s['diagram'] == 'J made by the band'
        if out['self_contained']:
            chains = [('K', start, [first] + surface_steps(s))]
        else:
            chains = [('K', start, [first])]
            if s:
                D = pd_J if s['diagram'] == 'KnotInfo diagram of J' else s['diagram_pd']
                chains.append(('J', D, surface_steps(s)))
                if is_rep:                                    # the claimed certificate only (slow when it fails)
                    out['J_made_reaches_surface_diagram'] = reaches(c['J_pd'], D)
    out['chains'] = []
    for part, st, steps in chains:
        try:
            out['chains'].append(dict(part=part, **run_chain(st, steps)))
        except (AssertionError, IndexError, KeyError, ValueError) as e:     # a malformed certificate fails
            out['chains'].append(dict(part=part, ok=False, error=f'{type(e).__name__}: {e}'))
    head = out['chains'][0]
    out['ok'] = all(ch['ok'] for ch in out['chains']) and bool(out.get('K_is_knotinfo_pd') or
                                                                 out.get('K_isometric_to_knotinfo'))
    if out['self_contained'] and head['ok']:
        out['b1'] = 1 - (head['discs'] - head['bands'])
    # without SnapPy: the diagram of K is KnotInfo's or turns into it, and the chain is self-contained or the diagram of
    # J made by the band turns into the diagram of the surface of J (a reflection: mirror the surface of J)
    out['diagrammatic'] = (all(ch['ok'] for ch in out['chains'])
                           and bool(out.get('K_is_knotinfo_pd') or out.get('K_simplifies_to'))
                           and (out['self_contained'] or bool(out.get('J_made_reaches_surface_diagram'))))
    if out['diagrammatic'] and not out['self_contained']:
        tail = out['chains'][1]
        out['b1'] = 1 - (tail['discs'] - tail['bands'] - 1)
    out['seconds'] = round(time.time() - t0, 2)
    return out


def reaches(start, P0):
    """Whether Reidemeister moves of spherogram turn the diagram start into the diagram P0 (the KnotInfo diagram of K,
    or a stored diagram of J) or its reflection (a diagram of the mirror image), up to isomorphism: then start is a
    diagram of the same knot or of its mirror image without SnapPy.  Attempts (variant): simplify ('level', 'global'),
    then random Reidemeister moves followed by simplify('global'); at most REACH_TRIES attempts.  Returns ['diagram' or
    'reflection', attempt]; otherwise the result of meets (a diagram that both turn into), or None."""
    import spherogram
    targets = {'diagram': [tuple(x) for x in P0], 'reflection': reflection(P0)}
    for i in range(REACH_TRIES):
        Q = variant(start, i)
        for label, T in targets.items():
            if len(Q) == len(T) and isomorphism(Q, T) is not None:
                return [label, i + 1]
    return meets(start, P0)


def reflection(P):
    """PD code of the reflected diagram (counterclockwise order reversed): a diagram of the mirror image."""
    return [(a, d, c, b) for a, b, c, d in P]


def variant(P, i):
    """Attempt i of reaches and meets: simplify (attempt 0: 'level'), or random Reidemeister moves followed by
    simplify('global'): in turn 100 moves as in spherogram's many_diagrams, and 1 to 12 moves that do not decrease the
    crossing number."""
    import spherogram
    L = spherogram.Link([list(x) for x in P])
    if i >= 2 and i % 2:
        L.backtrack(100)
    elif i >= 2:
        L.backtrack(steps=1 + i % 12, prob_type_1=0.2, prob_type_2=0.4)
    L.simplify('level' if i == 0 else 'global')
    return [tuple(x) for x in L.PD_code()]


def meets(A, B):
    """Fallback of reaches: whether Reidemeister moves of spherogram turn A and B (or the reflection of B) into a common
    diagram, up to isomorphism (variants of both, at most REACH_TRIES attempts each); ['common diagram' or 'common
    diagram with the reflection', attempt] or None."""
    def shape(P):                                            # bucket: crossings and face sizes
        return len(P), tuple(sorted(len(f) for f in faces(joins(P), len(P))))
    seen = {'A': defaultdict(list), 'B': defaultdict(list), 'R': defaultdict(list)}
    other = {'A': (('B', 'common diagram'), ('R', 'common diagram with the reflection')), 'B': (('A', 'common diagram'),),
             'R': (('A', 'common diagram with the reflection'),)}
    for i in range(REACH_TRIES):
        VA, VB = variant(A, i), variant(B, i)
        for side, V in (('A', VA), ('B', VB), ('R', reflection(VB))):
            if not V:
                continue
            key = shape(V)
            if any(isomorphism(V, W) is not None for W in seen[side][key]):
                continue
            seen[side][key].append(V)
            for o, label in other[side]:
                if any(isomorphism(V, W) is not None for W in seen[o][key]):
                    return [label, i + 1]
    return None


def isometric(A, B, tries=8):
    for _ in range(tries):
        try:
            return bool(A.is_isometric_to(B))
        except Exception:
            A.randomize(); B.randomize()
    return None


def same_band(c, r):
    return all(c[k] == r[k] for k in ('diagram', 'face', 'edges', 't', 'J'))


if __name__ == '__main__':
    import snappy, spherogram, database_knotinfo
    from database_knotinfo import link_list
    args = sys.argv[1:]
    workers = int(args[0]) if args and args[0].isdigit() else 10
    files = [a for a in args if not a.isdigit()] or [os.path.join(ROOT, 'data', f) for f in
                                                      ('klein_bottle_search_ranges.json', 'klein_bottle_search.json')]
    t_start = time.time()
    BY = {r['name']: r for r in link_list()[1:]}
    pd = lambda k: json.loads(BY[k]['pd_notation'])
    items, reps = [], {}
    for f in files:
        src = os.path.basename(f)
        K = json.load(open(f))['knots']
        for name, r in K.items():
            rep = r['representative']
            for i, c in enumerate(r['certificates']):
                if c['contradiction']:
                    continue
                is_rep = bool(rep) and same_band(c, rep)
                items.append((src, name, i, c, pd(name), pd(c['J']), is_rep))
                if is_rep:
                    reps[(src, name)] = i
            if r.get('moebius_band_of_K'):
                items.append((src, name, 'K', r['moebius_band_of_K'], pd(name), None, True))
    res = defaultdict(dict)
    with Pool(workers) as p:
        for o in p.imap_unordered(job, items, chunksize=4):
            res[(o['file'], o['knot'])][o['certificate']] = o
    out, summary = {}, {}
    for f in files:
        src = os.path.basename(f)
        K = json.load(open(f))['knots']
        knots = {}
        for name, r in K.items():
            got = res.get((src, name), {})
            i = reps.get((src, name))
            rep = got.get(i) if i is not None else None
            mk = got.get('K')
            knots[name] = dict(
                group=r['group'], range=r.get('range'), certificates=len([x for x in got if x != 'K']),
                certificates_rebuilt=sum(o['ok'] for x, o in got.items() if x != 'K'),
                self_contained_rebuilt=sum(o['ok'] and o['self_contained'] for x, o in got.items() if x != 'K'),
                representative=rep, representative_ok=bool(rep and rep['ok'] and rep['self_contained']),
                representative_diagrammatic=bool(rep and rep['diagrammatic']),
                representative_b1=rep.get('b1') if rep else None,
                moebius_band_of_K=mk, moebius_band_of_K_ok=bool(mk and mk['ok']),
                moebius_band_of_K_b1=mk.get('b1') if mk else None)
        T = {k: v for k, v in knots.items() if v['group'] == 'target'}
        start = lambda o: ('KnotInfo PD' if o.get('K_is_knotinfo_pd') else
                           f"stored diagram ({o['diagram_source']}): isometric {o.get('K_isometric_to_knotinfo')}, "
                           f"simplifies to {o.get('K_simplifies_to')}")
        summary[src] = dict(
            knots=len(knots), targets=len(T),
            targets_representative_rebuilt=sum(v['representative_ok'] for v in T.values()),
            targets_representative_b1=dict(Counter(v['representative_b1'] for v in T.values())),
            targets_not_rebuilt=[k for k, v in T.items() if not v['representative_ok']],
            targets_representative_diagrammatic=sum(v['representative_diagrammatic'] for v in T.values()),
            targets_not_diagrammatic=[k for k, v in T.items() if not v['representative_diagrammatic']],
            representatives_not_self_contained={k: dict(J=v['representative']['J'], ok=v['representative']['ok'],
                                                        J_made_reaches_surface_diagram=v['representative'].get(
                                                            'J_made_reaches_surface_diagram'),
                                                        diagrammatic=v['representative']['diagrammatic'])
                                                for k, v in knots.items() if v['representative']
                                                and not v['representative']['self_contained']},
            representatives_from_knotinfo_pd=sum(bool(v['representative'] and v['representative'].get('K_is_knotinfo_pd'))
                                                 for v in knots.values()),
            starts_not_knotinfo_pd={f"{o['knot']} {o['certificate']}": start(o) for k, v in knots.items()
                                    for o in (v['representative'], v['moebius_band_of_K'])
                                    if o and not o.get('K_is_knotinfo_pd')},
            certificates=sum(v['certificates'] for v in knots.values()),
            certificates_rebuilt=sum(v['certificates_rebuilt'] for v in knots.values()),
            self_contained_certificates_rebuilt=sum(v['self_contained_rebuilt'] for v in knots.values()),
            certificates_failing=[[k, x] for k, v in knots.items() for x, o in res.get((src, k), {}).items()
                                  if not o['ok'] and x != 'K'],
            others={k: dict(group=v['group'], representative_ok=v['representative_ok'], b1=v['representative_b1'],
                            J=v['representative']['J'] if v['representative'] else None,
                            moebius_band_of_K_ok=v['moebius_band_of_K_ok'], moebius_band_of_K_b1=v['moebius_band_of_K_b1'])
                    for k, v in knots.items() if v['group'] != 'target' and (v['representative'] or v['moebius_band_of_K'])})
        tries = Counter()
        for k in knots:
            for o in res.get((src, k), {}).values():
                for ch in o['chains']:
                    if 'simplify_tries_spherogram' in ch:
                        tries[(ch['simplify_tries_spherogram'], ch['simplify_tries_rebuilt'])] += 1
        summary[src]['simplify_tries'] = {f'{a},{b}': n for (a, b), n in sorted(tries.items(), key=str)}
        if any(v['group'] == 'mobius' for v in knots.values()):              # output of mode mobius
            mk = {k: v for k, v in knots.items() if v['moebius_band_of_K']}
            kb = {k: v for k, v in knots.items() if v['representative']}
            ranges = ('[1,2]', '[1,3]')
            summary[src]['mobius'] = dict(
                moebius_bands_of_K_rebuilt={g: [sum(v['moebius_band_of_K_ok'] for v in mk.values() if v['range'] == g
                                                    and v['group'] == 'mobius'),
                                                sum(v['range'] == g and v['group'] == 'mobius' for v in mk.values())]
                                            for g in ranges},
                moebius_bands_of_K_b1=dict(Counter(v['moebius_band_of_K_b1'] for v in mk.values())),
                moebius_bands_of_K_diagrammatic=sum(v['moebius_band_of_K']['diagrammatic'] for v in mk.values()),
                moebius_bands_of_K_failing=[k for k, v in mk.items() if not v['moebius_band_of_K_ok']],
                moebius_bands_of_K_of_controls=[k for k, v in mk.items() if v['group'] == 'control'],
                klein_bottles=dict(representatives=len(kb),
                                   rebuilt=sum(v['representative_ok'] for v in kb.values()),
                                   b1=dict(Counter(v['representative_b1'] for v in kb.values())),
                                   diagrammatic=sum(v['representative_diagrammatic'] for v in kb.values()),
                                   not_rebuilt=[k for k, v in kb.items() if not v['representative_ok']]))
        out[src] = knots
    meta = dict(tries=TRIES, files=[os.path.basename(f) for f in files],
                versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                              database_knotinfo=database_knotinfo.version(), python=sys.version.split()[0]),
                seconds=round(time.time() - t_start))
    mobius = [os.path.basename(f) for f in files] == ['klein_bottle_search_mobius.json']
    path = os.environ.get('REBUILD_OUT', os.path.join(ROOT, 'data', 'klein_bottle_rebuild_mobius.json' if mobius
                                                      else 'klein_bottle_rebuild.json'))
    json.dump(dict(meta=meta, summary=summary, knots=out), open(path, 'w'), separators=(',', ':'), default=str)
    for src, s in summary.items():
        print(src, {k: v for k, v in s.items() if k not in ('others', 'mobius')})
        print('  mobius:', s['mobius']) if 'mobius' in s else print('  others:', s['others'])
    print('time %.0fs' % (time.time() - t_start))
