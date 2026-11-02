#!/usr/bin/env python
"""Donaldson's theorem for definite checkerboard surfaces ([C] Proposition 2.7): the search of [C] Proposition 3.2(c),
on the knots that Proposition 3.2(a), (b) and the second search for Mobius bands leave open, and its controls
([C] Section 3.5, Donaldson's theorem).

The obstruction.  Let s = sigma(K) + 4 Arf(K) mod 8 (KnotInfo Seifert matrix; field s of data/rank1_targets.json) and
eps = +1 if s = 2, eps = -1 if s = 6.  If K bounds a smooth Mobius band F, then Sign Sigma_2(B^4, F) = eps
(Gilmer-Livingston, Cor. 8).  If a diagram of K has a checkerboard surface V whose Gordon-Litherland form G_V is
(-eps)-definite, then Sigma_2(B^4, V) and -Sigma_2(B^4, F) glue along Sigma_K to a closed smooth definite 4-manifold of
rank b_1(V) + 1, so by Donaldson's theorem the positive definite lattice L = -eps G_V embeds in Z^{b_1(V) + 1}.  A
diagram for which L does not embed shows gamma_4(K) >= 2.  For s = 0 a Mobius band may have either sign and a
non-alternating knot has no definite surfaces of both signs, so only s = 2 and s = 6 are searched; s = 4 does not occur
among these knots (Yasuhara).
Conventions.  rank1_goeritz.goeritz_pair gives, for each colour, the matrix on the faces of that colour with one face
deleted; it is -G_V for the surface V of the other colour, G_V being the Goeritz matrix of Gordon and Litherland
(Thm. 1).  Check (convention_check): on the KnotInfo PD code of 3_1 (sigma = -2) the matrix of rank 2, which belongs
to the Seifert surface (orientable, so its form has signature sigma = -2), is positive definite.  Hence L is the matrix
M of goeritz_pair for s = 2 and -M for s = 6, as in jk_lattice.py, whose signed test has no contradiction on the
alternating controls.
Diagrams, in this order: the KnotInfo PD code; the diagrams stored by mode mobius of klein_bottle_search.py and by
mobius_open_search.py (all of them checked there to have the exterior of the KnotInfo knot); the certificate diagrams
kept in data/donaldson_pilot_diagrams.json (source 'pilot' in the records); and
the new diagrams of about NEW random walks of Reidemeister moves of spherogram from the KnotInfo PD code, as in
mobius_open_search.new_diagrams: kind walk (10 to 60 moves that do not decrease the crossing number, types II and III),
walk_level (20 to 100 moves, then simplify('level')) and walk_global (20 to 200 moves, then simplify('global')), kept
with c to c + X_MAX crossings and a new signed Gauss code.  Random moves of spherogram cannot be repeated from a seed,
so every certificate diagram is stored.
Search.  For each diagram and colour, L is tested for positive definiteness (Sylvester's criterion, exact) and, if
definite and new for the knot, for an embedding in Z^{rank L + 1} (jk_lattice.embeds, a complete search).  A target
stops at its first L that does not embed (the certificate); a control tests every definite L and any L that does not
embed is a contradiction.  Controls: the non-alternating knots that are not Montesinos with s = 2 or 6 that bound a
certified smooth Mobius band (mode mobius of klein_bottle_search.py, and mobius_open_search.py).
The certificates are checked by scripts/donaldson_check.py, which shares no code with this search.
A run appends each finished knot to data/donaldson_<mode>.partial (gitignored), from which a stopped run resumes, and
at the end writes its part of the output and removes that file.  Environment: DS_ONLY=knot,knot,... runs these knots
again and merges them into the existing output; DS_NEW replaces NEW (the output records it in meta).
Usage: python scripts/donaldson_surfaces.py targets | controls [workers]
Output: data/donaldson_surfaces.json (both parts; a mode replaces its own part)
"""
import gzip, json, os, random, sys, time, warnings, zlib
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'vendor', 'code', 'lower_bounds', 'owens'))
warnings.filterwarnings('ignore')
from rank1_goeritz import goeritz_pair, det
from jk_lattice import embeds
from owens_obstruction import positive_definite
from smooth_top_gap import gauss_key
from open_sets import upto

NEW = int(os.environ.get('DS_NEW', 1500))  # random walks per knot; repeated diagrams are scanned once
X_MAX = 10                                # new diagrams have at most c + X_MAX crossings
ATTEMPTS = 4                              # attempts per new diagram wanted
OUT = os.path.join(ROOT, 'data', 'donaldson_surfaces.json')


def data(name):
    return json.load(open(os.path.join(ROOT, 'data', name)))


def knots(mode):
    """The knots of a mode with their s.  Targets: the knots left open before [C] Proposition 3.2(c) (the set oleft3 of
    paper_counts.py) that are neither alternating nor Montesinos.  Controls: the knots of that kind with a Mobius band
    of mode mobius of klein_bottle_search.py or of mobius_open_search.py.  main() searches those with s = 2 or 6."""
    g = upto('oleft3')
    alt, mont = g['alt'], g['mont']
    R1 = data('rank1_targets.json')
    if mode == 'targets':
        pool = g['oleft3']
    else:
        old = data('klein_bottle_search_mobius.json')['knots']
        mos = json.load(gzip.open(os.path.join(ROOT, 'data', 'mobius_open_search.json.gz'), 'rt'))['summary']
        pool = {k for k, v in old.items() if v.get('moebius_band_of_K')} | set(mos['self_contained_up_to_reidemeister_moves_of_J'])
    names = sorted(k for k in pool if not alt(k) and not mont(k))
    assert all(k in R1 for k in names)
    return {k: R1[k]['s'] for k in names}


def stored_diagrams(names):
    """Stored diagrams of each knot: (source, pd), the KnotInfo diagram first."""
    old = data('klein_bottle_search_mobius.json')['knots']
    mos = json.load(gzip.open(os.path.join(ROOT, 'data', 'mobius_open_search.json.gz'), 'rt'))['knots']
    pilot = data('donaldson_pilot_diagrams.json')['diagrams']
    out = {}
    for k in names:
        ds = [('mobius search', d['pd']) for d in old.get(k, {}).get('diagrams', [])]
        ds += [('open search', d['pd']) for d in mos.get(k, {}).get('diagrams', [])]
        ds += [('pilot', pd) for pd in pilot.get(k, [])]
        out[k] = ds
    return out


def new_diagrams(pd0, seed, count):
    """count new diagrams of the knot of pd0 by random Reidemeister moves (see the docstring)."""
    import spherogram
    c = len(pd0)
    out, attempts = [], 0
    while len(out) < count and attempts < ATTEMPTS * count:
        random.seed(seed + attempts)
        kind = ('walk', 'walk_level', 'walk_global')[attempts % 3]
        attempts += 1
        L = spherogram.Link([list(x) for x in pd0])
        if kind == 'walk':
            st = random.randint(10, 60); L.backtrack(steps=st, prob_type_1=0, prob_type_2=random.randint(1, 5) / st)
        elif kind == 'walk_level':
            L.backtrack(steps=random.randint(20, 100), prob_type_1=0.05, prob_type_2=0.3); L.simplify('level')
        else:
            L.backtrack(steps=random.randint(20, 200), prob_type_1=0.1, prob_type_2=0.4); L.simplify('global')
        P = [tuple(x) for x in L.PD_code()]
        if len(L.link_components) == 1 and c <= len(P) <= c + X_MAX:
            out.append((kind, P))
    return out, attempts


def one(item):
    """Scan the diagrams of one knot.  item = (mode, name, s, det K, KnotInfo PD code, stored diagrams)."""
    mode, name, s, D, pd0, stored = item
    t0 = time.time()
    seen_d, seen_L = set(), set()
    rec = dict(knot=name, s=s, det=D, diagrams=Counter(), lattices=0, ranks=Counter(), embedded=0, certificate=None,
               contradictions=[])
    new, attempts = new_diagrams([tuple(x) for x in pd0], zlib.crc32(name.encode()), NEW)
    rec['new_attempts'] = attempts
    for src, pd in [('knotinfo', pd0)] + stored + new:
        P = [tuple(x) for x in pd]
        key = gauss_key(P)
        if key in seen_d:
            continue
        seen_d.add(key)
        rec['diagrams'][src] += 1
        for colour, M in enumerate(goeritz_pair(P)):
            if not M:
                continue
            L = M if s == 2 else [[-x for x in row] for row in M]
            k = str(L)
            if k in seen_L or not positive_definite(L):
                continue
            seen_L.add(k)
            rec['lattices'] += 1; rec['ranks'][len(L)] += 1
            assert abs(det(L)) == D, 'a Goeritz determinant differs from det K'
            if embeds(L):
                rec['embedded'] += 1
                continue
            cert = dict(source=src, crossings=len(P), colour=colour, rank=len(L), lattice=L, pd=[list(x) for x in P])
            if mode == 'targets':
                rec['certificate'] = cert
                break
            rec['contradictions'].append(cert)
        if rec['certificate']:
            break
    rec['diagrams'] = dict(rec['diagrams']); rec['ranks'] = dict(sorted(rec['ranks'].items()))
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def convention_check(BY):
    G = goeritz_pair([tuple(x) for x in json.loads(BY['3_1']['pd_notation'])])
    two = [M for M in G if len(M) == 2]
    return dict(sigma_3_1=int(BY['3_1']['signature']), matrices=G,
                ok=int(BY['3_1']['signature']) == -2 and len(two) == 1 and positive_definite(two[0]))


def main(mode, workers):
    import snappy, spherogram
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}
    conv = convention_check(BY)
    assert conv['ok'], conv
    S = knots(mode)
    stored = stored_diagrams(S)
    partial = os.path.join(ROOT, 'data', 'donaldson_%s.partial' % mode)
    done = {}
    if os.path.exists(partial):
        done = {r['knot']: r for r in map(json.loads, open(partial)) if r}
    only = os.environ.get('DS_ONLY')
    if only and os.path.exists(OUT):          # run some knots again and merge them into the existing output
        done = {k: r for k, r in json.load(open(OUT)).get(mode, {}).get('knots', {}).items() if k not in only.split(',')}
    items = [(mode, k, s, int(BY[k]['determinant']), json.loads(BY[k]['pd_notation']), stored[k]) for k, s in S.items()
             if s in (2, 6) and k not in done and (not only or k in only.split(','))]
    t0 = time.time()
    with Pool(workers) as P, open(partial, 'a') as f:
        for r in P.imap_unordered(one, items):
            f.write(json.dumps(r) + '\n'); f.flush()
            done[r['knot']] = r
            print(r['knot'], 's', r['s'], 'diagrams', sum(r['diagrams'].values()), 'lattices', r['lattices'],
                  'CERTIFICATE rank %d' % r['certificate']['rank'] if r['certificate'] else '',
                  'CONTRADICTION' if r['contradictions'] else '', r['seconds'], flush=True)
    recs = {k: done[k] for k in sorted(done)}
    if mode == 'targets':
        summary = dict(knots=len(S), s=dict(Counter(S.values())), searched=len(recs),
                       with_definite_surface=sum(r['lattices'] > 0 for r in recs.values()),
                       certified=sorted(k for k, r in recs.items() if r['certificate']),
                       by_s=dict(Counter(r['s'] for r in recs.values() if r['certificate'])),
                       certificate_sources=dict(Counter(r['certificate']['source'] for r in recs.values() if r['certificate'])),
                       certificate_ranks=dict(Counter(r['certificate']['rank'] for r in recs.values() if r['certificate'])),
                       diagrams=sum(sum(r['diagrams'].values()) for r in recs.values()),
                       lattices=sum(r['lattices'] for r in recs.values()))
    else:
        summary = dict(knots=len(S), s=dict(Counter(S.values())), searched=len(recs),
                       with_definite_surface=sum(r['lattices'] > 0 for r in recs.values()),
                       lattices=sum(r['lattices'] for r in recs.values()),
                       contradictions=sorted(k for k, r in recs.items() if r['contradictions']),
                       diagrams=sum(sum(r['diagrams'].values()) for r in recs.values()))
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    out.setdefault('meta', {})[mode] = dict(new=NEW, x_max=X_MAX, attempts=ATTEMPTS, convention_check=conv,
                                            seconds=round(time.time() - t0),
                                            versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__))
    out[mode] = dict(summary=summary, knots=recs)
    json.dump(out, open(OUT, 'w'), separators=(',', ':'))
    os.remove(partial)
    print('summary:', json.dumps(summary))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10)
