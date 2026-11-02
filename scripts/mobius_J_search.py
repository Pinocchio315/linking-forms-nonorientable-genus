#!/usr/bin/env python
"""Mobius bands for the knots of range [1,2] that [C] leaves open (the set oleft of scripts/paper_counts.py): a search
centred on the knots J, a band K -> J followed by a ribbon disc of J in a diagram of J that Reidemeister moves relate to
the diagram made by the band.

Why.  The second search (scripts/mobius_open_search.py) built the face bands K -> J in the diagrams of K that it stored
and looked for a ribbon disc of J by one band in the diagram of J made by the band.  It looked in other diagrams of J
(surface_of_J, J_VARIANTS = 100 attempts) only when SnapPy identified J with a knot that KnotInfo lists as slice; among
these knots that happened only for 13a_1412, with J = 12a_77 and no disc.  Yet 12a_77 has a ribbon disc by one band in
some diagrams reached from that diagram of J, in about 1 % of them.  Here every J of those bands is a candidate.

Stage 1, the knots J (one stored diagram of K per task): the bands of mobius_open_search.search_diagram, face bands
with |s| <= S_MAX built when det J = |alpha + n s| (arc_alphas) can be a nonzero square, and J a knot with square
det, Delta_J of Fox-Milnor form and signature 0.  The diagrams of K are those stored by the second search (field pd of
its output, all with the exterior of the KnotInfo diagram).
Stage 2, groups: the J of one K are grouped by (det, Delta_J, volume of the exterior to 1e-6); a J that is not
hyperbolic is a group of its own.  Each group keeps REPS representatives with fewest crossings (REPS_SLICE for a J
that KnotInfo lists as slice).  A group is dropped when SnapPy identifies J only with knots that KnotInfo lists as not
slice.  Necessary conditions for J to be smoothly slice order and prune the search, and no certificate rests on them:
the Levine-Tristram signatures at the roots of unity of order at most LT_ORDER that are not roots of Delta_J vanish
(numerical, from a Seifert matrix), and tau(J) = epsilon(J) = 0 (Szabo's program on a simplified diagram with at most
HFK_CROSSINGS crossings; a larger J passes).
Stage 3, discs: the groups that pass, KnotInfo-slice ones first, then by crossing number and size, each in turn: a
ribbon disc by one or two bands (mobius_open_search.ribbon) in the diagram of J made by the band, then by one band
(ribbon_one; every disc found by the exploration of 13a_1412 had one band) in up to N_VARIANTS diagrams made from it
as in mobius_open_search.surface_of_J (variant, and for odd attempts 1 to 8 further random moves of types II and III).
The first disc found ends the search for K, and so does the time limit KNOT_SECONDS, checked between diagrams.
Certificates have the format of the second search (kind 'band to a ribbon knot', field surface_of_J, the diagram of K
by its index there), with mobius_open_search.certificate_recheck, so that scripts/mobius_open_rebuild.py rebuilds them
(REBUILD_IN=data/mobius_J_search.json.gz).
Controls: the knots of mobius_open_search.SETTLED (gamma_4 = 2 by the exact correction terms of [C] Proposition 3.2(b))
and CONTROLS knots drawn by random.Random(0) from the controls of the second search (gamma_4 >= 2 by [C] Proposition
3.1(b) or (c)); a disc for a control would be a contradiction.
Checkpoints: stages 1 and 2 are saved to the .groups.json.gz file next to the output and reused by a later run with the
same parameters; stage 3 appends one line per knot to the .partial file, and a later run skips the knots in it.
Usage: python scripts/mobius_J_search.py [workers] | status
Environment: MJS_ONLY=knot,knot,... restricts the knots (pilots and tests); MJS_OUT replaces the output file.
Output: data/mobius_J_search.json.gz
"""
import cmath, gzip, json, math, os, random, sys, time, warnings, zlib
from collections import Counter, defaultdict
from multiprocessing import Pool
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import mobius_open_search as MOS
from bands import from_pd, faces, band, pd_code, components, graph_components
from smooth_top_gap import det_link, alexander_fox, gauss_key
from klein_bottle_search import arc_alphas, fox_milnor

S_MAX = MOS.S_MAX                         # half twists of the band K -> J, as in the second search
REPS, REPS_SLICE = 1, 2                   # representatives searched per group (per KnotInfo-slice group)
LT_ORDER = 12                             # Levine-Tristram signatures at roots of unity of order <= LT_ORDER
HFK_CROSSINGS = 20                        # tau and epsilon of J only for simplified diagrams this small
N_VARIANTS = 300                          # diagrams of J tried per representative (surface_of_J: 100)
KNOT_SECONDS = 600                        # stage 3 stops for a knot after this time
CONTROLS = 8                              # controls of the second search, besides SETTLED
OUT = os.environ.get('MJS_OUT', os.path.join(ROOT, 'data', 'mobius_J_search.json.gz'))
GROUPS = OUT.replace('.json.gz', '.groups.json.gz')
PARTIAL = OUT.replace('.json.gz', '.partial')
PARAMS = dict(s_max=S_MAX, reps=REPS, reps_slice=REPS_SLICE, lt_order=LT_ORDER, hfk_crossings=HFK_CROSSINGS,
              n_variants=N_VARIANTS, knot_seconds=KNOT_SECONDS, controls=CONTROLS)


def square(x):
    return x > 0 and isqrt(x) ** 2 == x


def knots():
    """Targets (oleft of paper_counts.py) and controls, with the diagrams stored by the second search."""
    from open_sets import upto
    g = upto('oleft3')
    dsf = json.load(open(os.path.join(ROOT, 'data', 'donaldson_surfaces.json')))
    oleft = g['oleft3'] - set(dsf['targets']['summary']['certified'])
    assert len(oleft) == 229
    S = json.load(gzip.open(os.path.join(ROOT, 'data', 'mobius_open_search.json.gz')))['knots']
    ctl = sorted(k for k, r in S.items() if r['group'] == 'control' and k not in MOS.SETTLED)
    ctl = set(random.Random(0).sample(ctl, CONTROLS)) | set(MOS.SETTLED)
    out = {}
    for k, r in S.items():
        group = 'target' if k in oleft else 'control' if k in ctl else None
        if group:
            assert all(d['isometric_to_knotinfo'] for d in r['diagrams'])
            out[k] = dict(group=group, det=r['det'], diagrams=[(d['index'], d['pd']) for d in r['diagrams']],
                          stored=r['diagrams'])
    assert sum(v['group'] == 'target' for v in out.values()) == 229
    only = os.environ.get('MJS_ONLY')
    return {k: v for k, v in out.items() if not only or k in only.split(',')}


# ---------- stage 1: the knots J of one diagram of K

def volume(P):
    import snappy, spherogram
    try:
        E = spherogram.Link([list(x) for x in P]).exterior()
        v = float(E.volume())
        return round(v, 6) if E.solution_type() == 'all tetrahedra positively oriented' and v > 0.5 else None
    except Exception:
        return None


def collect(task):
    name, idx, pd, n = task
    nbr, under, label = from_pd([tuple(x) for x in pd])
    arcs = MOS.face_arcs(nbr, label)
    alphas = [arc_alphas(nbr, under, f[i], f[j], n) for fi, f, i, j in arcs]
    st, out = Counter(arcs=len(arcs)), []
    for s in sorted(range(-S_MAX, S_MAX + 1), key=lambda x: (abs(x), x)):
        for k, (fi, f, i, j) in enumerate(arcs):
            if alphas[k] and not any(square(abs(a + n * s)) for a in alphas[k]):
                continue
            N, U = band(nbr, under, f[i], f[j], s)
            st['bands_built'] += 1
            if graph_components(N, len(U)) > 1 or len(components(N)) != 1:
                continue
            P = pd_code(N, U)
            d = det_link(P)
            if not square(d):
                continue
            st['J_square_det'] += 1
            A = alexander_fox(P)
            if not fox_milnor(A):
                continue
            st['J_fox_milnor'] += 1
            if MOS.sigma(P) != 0:
                continue
            st['J_signature_zero'] += 1
            out.append(dict(diagram=idx, arc=k, face=fi, darts=[i, j], s=s, det=d, A=A, vol=volume(P), c=len(P),
                            pd=[list(x) for x in P]))
    return name, idx, out, dict(st)


# ---------- stage 2: groups and necessary conditions for sliceness

def lt_signatures(P):
    """{'k/m': sigma_omega} at omega = e^(2 pi i k / m), m <= LT_ORDER, gcd(k, m) = 1, k <= m / 2, omega not a root of
    Delta (numerical)."""
    import numpy as np, spherogram
    V = np.array(spherogram.Link([list(x) for x in P]).seifert_matrix(), dtype=float)
    out = {}
    if not len(V):
        return out
    for m in range(3, LT_ORDER + 1):
        for k in range(1, m // 2 + 1):
            if math.gcd(k, m) != 1:
                continue
            w = cmath.exp(2j * math.pi * k / m)
            if abs(np.linalg.det(V - w * V.T)) < 1e-6:
                continue
            H = (1 - w) * V + (1 - w.conjugate()) * V.T
            ev = np.linalg.eigvalsh(H)
            out['%d/%d' % (k, m)] = int((ev > 1e-9).sum() - (ev < -1e-9).sum())
    return out


def hfk(P):
    import spherogram, knot_floer_homology as kfh
    L = spherogram.Link([list(x) for x in P])
    L.simplify('global')
    if len(L.crossings) > HFK_CROSSINGS:
        return dict(crossings=len(L.crossings))
    h = kfh.pd_to_hfk(L.PD_code())
    return dict(crossings=len(L.crossings), tau=h['tau'], epsilon=h['epsilon'])


def judge(task):
    """Identification, Levine-Tristram signatures and tau, epsilon of one group, on its first representative."""
    name, key, rep = task
    P = [tuple(x) for x in rep['pd']]
    ids, names = MOS.identify(P)
    ki = {x: MOS.G['slice'].get(x) for x in names}
    rec = dict(names=names, census=ids[:3], knotinfo_slice=sorted(x for x, v in ki.items() if v == '0'))
    if names and not rec['knotinfo_slice'] and not set(ki.values()) & {'0', '[0,1]', ''}:
        return name, key, dict(rec, verdict='KnotInfo: not slice')
    lt = lt_signatures(P)
    rec['lt_nonzero'] = {w: v for w, v in lt.items() if v}
    if rec['lt_nonzero']:
        return name, key, dict(rec, verdict='Levine-Tristram signature')
    rec['hfk'] = hfk(P)
    if rec['hfk'].get('tau') or rec['hfk'].get('epsilon'):
        return name, key, dict(rec, verdict='tau or epsilon')
    return name, key, dict(rec, verdict='candidate')


def stages_1_2(KN, pool):
    """Stages 1 and 2, or their checkpoint."""
    if os.path.exists(GROUPS):
        G = json.load(gzip.open(GROUPS, 'rt'))
        if G['params'] == PARAMS and set(KN) <= set(G['groups']):
            return {k: G['groups'][k] for k in KN}, {k: Counter(G['counts'][k]) for k in KN}
    t0 = time.time()
    tasks = [(k, i, pd, v['det']) for k, v in KN.items() for i, pd in v['diagrams']]
    Js, cnt = defaultdict(list), defaultdict(Counter)
    for name, idx, out, st in pool.imap_unordered(collect, tasks, chunksize=4):
        Js[name] += out
        cnt[name].update(st)
    print('stage 1: %d diagrams, %.0fs' % (len(tasks), time.time() - t0), flush=True)
    groups = {k: {} for k in KN}
    for name, L in Js.items():
        gr = defaultdict(list)
        for j in L:
            key = ('%d|%s|%s' % (j['det'], j['A'], j['vol']) if j['vol'] is not None else
                   'nonhyp|%d|%d|%d|%d' % (j['diagram'], j['arc'], j['s'], j['det']))
            gr[key].append(j)
        groups[name] = {key: dict(size=len(v), reps=sorted(v, key=lambda j: (j['c'], j['diagram']))[:REPS_SLICE])
                        for key, v in gr.items()}
    jt = [(name, key, g['reps'][0]) for name, gs in groups.items() for key, g in gs.items()]
    for name, key, rec in pool.imap_unordered(judge, jt, chunksize=8):
        groups[name][key]['judge'] = rec
    print('stage 2: %d groups, %.0fs' % (len(jt), time.time() - t0), flush=True)
    json.dump(dict(params=PARAMS, groups=groups, counts={k: dict(cnt[k]) for k in KN}), gzip.open(GROUPS, 'wt'),
              separators=(',', ':'))
    return groups, cnt


# ---------- stage 3: ribbon discs in diagrams of J

def disc_of_J(P, seed, st, deadline):
    """A ribbon disc of J: one or two bands in P, then one band in diagrams made from P as in
    mobius_open_search.surface_of_J, part (a), until the deadline."""
    import spherogram
    rb = MOS.ribbon(P, st)
    if rb:
        return None, rb
    keys = {gauss_key(P)}
    for i in range(N_VARIANTS):
        if time.time() > deadline:
            st['stopped_by_time'] += 1
            break
        random.seed(seed + i)
        V = MOS.variant(P, i)
        W, steps = V, 0
        if V and i % 2:
            L = spherogram.Link([list(x) for x in V])
            steps = 1 + (i // 2) % 8
            L.backtrack(steps=steps, prob_type_1=0, prob_type_2=0.5)
            W = [tuple(x) for x in L.PD_code()]
        if not W or gauss_key(W) in keys:
            continue
        keys.add(gauss_key(W))
        st['J_variants'] += 1
        rb = MOS.ribbon_one(W, st)
        if rb:
            return dict(diagram=MOS.REACHED, diagram_pd=[list(x) for x in W], variant=i, walk_steps=steps,
                        simplified_pd=[list(x) for x in V] if steps else None, ribbon_bands=[rb]), [rb]
    return None, None


def certificate(name, rep, key, g, s, rb, pd_of, src_of):
    import spherogram
    cert = dict(face=rep['face'], darts=rep['darts'], s=rep['s'], J_pd=rep['pd'], J_det=rep['det'], J_alexander=rep['A'],
                diagram=rep['diagram'], diagram_pd=pd_of[rep['diagram']], diagram_source=src_of[rep['diagram']],
                kind='band to a ribbon knot', group=key, J_names=g['judge']['names'])
    nbr, under, label = from_pd([tuple(x) for x in pd_of[rep['diagram']]])
    cert.update(MOS.band_record(MOS.face_arcs(nbr, label)[rep['arc']], label, s=rep['s']))
    if s:
        cert.update(surface_of_J=s, verification='self-contained up to Reidemeister moves of J')
    else:
        cert.update(ribbon_bands=rb, verification='self-contained')
    E_K = spherogram.Link([list(x) for x in pd_of[0]]).exterior()
    cert['recheck'] = MOS.certificate_recheck(pd_of[rep['diagram']], cert, E_K)
    return cert


def discs(task):
    """Stage 3 for one knot: its candidate groups in turn, until a disc or KNOT_SECONDS."""
    name, groups, pd_of, src_of = task
    t0, st, tried = time.time(), Counter(), 0
    deadline = t0 + KNOT_SECONDS
    for key, g in groups:
        if time.time() > deadline:
            break
        tried += 1
        for rep in g['reps'][:REPS_SLICE if g['judge']['knotinfo_slice'] else REPS]:
            P = [tuple(x) for x in rep['pd']]
            s, rb = disc_of_J(P, zlib.crc32(('%s|%s|%d' % (name, key, rep['diagram'])).encode()), st, deadline)
            if rb:
                return dict(knot=name, certificate=certificate(name, rep, key, g, s, rb, pd_of, src_of),
                            groups_tried=tried, groups=len(groups), counts=dict(st), seconds=round(time.time() - t0, 1))
    return dict(knot=name, certificate=None, groups_tried=tried, groups=len(groups), counts=dict(st),
                seconds=round(time.time() - t0, 1))


# ---------- the run

def init(G):
    MOS.worker_init(G)


def read_partial():
    if not os.path.exists(PARTIAL):
        return {}
    out = {}
    for line in open(PARTIAL):
        if line.strip():
            r = json.loads(line)
            out[r['knot']] = r
    return out


if __name__ == '__main__':
    if sys.argv[1:2] == ['status']:
        done = read_partial()
        print('knots done', len(done), '| discs', sorted(k for k, r in done.items() if r['certificate']))
        sys.exit()
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    t_start = time.time()
    KN = knots()
    G = MOS.knotinfo_data()
    print('knots', Counter(v['group'] for v in KN.values()), flush=True)
    with Pool(workers, initializer=init, initargs=(G,)) as pool:
        groups, cnt = stages_1_2(KN, pool)
        res = read_partial()
        order = lambda kv: (not kv[1]['judge']['knotinfo_slice'], kv[1]['reps'][0]['c'], -kv[1]['size'])
        dt = [(name, sorted(((k, g) for k, g in gs.items() if g['judge']['verdict'] == 'candidate'), key=order),
               dict(KN[name]['diagrams']), {d['index']: d['source'] for d in KN[name]['stored']})
              for name, gs in groups.items() if name not in res]
        dt.sort(key=lambda t: -len(t[1]))       # long tasks first
        print('stage 3: %d knots to do, %d done before' % (len(dt), len(res)), flush=True)
        with open(PARTIAL, 'a') as fp:
            for r in pool.imap_unordered(discs, dt):
                res[r['knot']] = r
                fp.write(json.dumps(r, separators=(',', ':'), default=str) + '\n'); fp.flush()
                print(r['knot'], KN[r['knot']]['group'], 'DISC' if r['certificate'] else '-', '%d/%d groups' %
                      (r['groups_tried'], r['groups']), r['seconds'], '(%d done, %.0fs)' % (len(res), time.time() - t_start),
                      flush=True)
    res = {k: r for k, r in res.items() if k in KN}
    has = lambda k: bool((res.get(k) or {}).get('certificate'))
    summary = dict(
        knots=dict(Counter(v['group'] for v in KN.values())),
        moebius_band=sorted(k for k in res if has(k) and KN[k]['group'] == 'target'),
        recheck_failing=sorted(k for k in res if has(k) and not res[k]['certificate']['recheck']['passes']),
        contradictions=sorted(k for k in res if has(k) and KN[k]['group'] == 'control'),
        stopped_by_time=sorted(k for k, r in res.items() if not has(k) and r['counts'].get('stopped_by_time')),
        groups=dict(Counter(g['judge']['verdict'] for gs in groups.values() for g in gs.values())),
        counts=dict(sum((cnt[k] for k in cnt), Counter())))
    out = dict(meta=dict(params=PARAMS, seconds=round(time.time() - t_start)), summary=summary,
               knots={k: dict(group=KN[k]['group'], det=KN[k]['det'], counts=dict(cnt[k]), result=res.get(k),
                              # the fields that scripts/mobius_open_rebuild.py reads (REBUILD_IN)
                              certificates=[res[k]['certificate']] if has(k) else [],
                              diagrams=KN[k]['stored'] if has(k) else [], knotinfo_slice=[],
                              groups={key: dict(size=g['size'], judge=g['judge'],
                                                reps=[{a: b for a, b in r.items() if a != 'pd'} for r in g['reps']])
                                      for key, g in groups.get(k, {}).items()})
                      for k in KN})
    json.dump(out, gzip.open(OUT, 'wt'), separators=(',', ':'), default=str)
    print(json.dumps(summary, indent=1, default=str))
    assert not summary['contradictions'] and not summary['recheck_failing'], summary
