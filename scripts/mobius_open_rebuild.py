#!/usr/bin/env python
"""Rebuild of the Mobius bands of the second search of [C] Section 3.3 (Proposition 3.3(a)), found by
scripts/mobius_open_search.py, and of its bands to knots that KnotInfo lists as slice.

Mobius bands.  Each is a chain of face bands from a stored diagram of K: a band K -> J (s half twists), then a band
J -> unknot, or a ribbon disc of J (one band to the two-component unlink, or two bands J -> L -> three-component unlink),
in the diagram of J made by the first band (verification 'self-contained') or in a diagram of J related to it by
Reidemeister moves ('self-contained up to Reidemeister moves of J'; the field surface_of_J of the record).  The first
band joins two knots, so the surface (bands at successive levels of a collar, a disc for each component of the final
unlink) is nonorientable with one boundary component, chi = discs - bands = 0 and b_1 = 1.
Rebuild.  No code of the searches (scripts/mobius_open_search.py, klein_bottle_search.py, bands.py,
smooth_top_gap.py) is used.  Bands are done by run_chain of scripts/klein_bottle_rebuild.py: the face from its number,
its edges and the two edges of the band compared with those stored, each band done twice (spherogram's add_one_band and
klein_bottle_rebuild.band_move), the two diagrams isomorphic to each other and to the stored diagram after the band,
the number of components checked, and the final diagram simplified to no crossings by Reidemeister moves alone with the
stated number of unlinked unknotted components.  The diagram of K: the KnotInfo PD code itself (database_knotinfo
pd_notation) or a stored diagram whose exterior is that of the KnotInfo diagram (SnapPy) and which spherogram's
Reidemeister moves relate to the KnotInfo diagram (klein_bottle_rebuild.reaches, with meets as fallback), directly or
through the stored diagrams it was made from (field parent).  For a surface of J in another diagram D of J, the
diagram of J made by the band is related to D in the same way (for a diagram W made from a simplified diagram V by
moves that do not decrease the crossing number: the diagram of J to V and W to V; for a diagram of the KnotInfo-slice
knot X, which the search related to J by a common diagram, meets first).  A certificate is diagrammatic when all these
relations are found: then it depends neither on SnapPy nor on any value of KnotInfo.
KnotInfo-slice bands (reported apart; no result uses them).  The band K -> J is rebuilt as above; the exterior of J is
compared with that of the KnotInfo diagram of the knot X named by the search (SnapPy is_isometric_to, the determinant
of the cusp map of each isometry that extends to the knots: +1 for J = X, -1 for the mirror image); the Jones
polynomial of a simplified diagram of J (at most JONES_MAX crossings) and of the KnotInfo PD code of X, by a state sum
over all 2^n smoothings written here, decides J = X or J = mirror image of X when the polynomial of X is not symmetric;
and X is smoothly slice in KnotInfo (database_knotinfo smooth_four_genus, the same as in the KnotInfo 2026.10.5 CSV).
Usage: python scripts/mobius_open_rebuild.py [workers]
Environment, for tests: REBUILD_IN and REBUILD_OUT replace the input and the output file.
Output: data/mobius_open_rebuild.json
"""
import gzip, json, os, sys, time, warnings
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
from klein_bottle_rebuild import run_chain, reaches, meets, isometric

JONES_MAX = 16                            # crossings of the simplified diagram of J for the state sum
SIMPLIFY_TRIES = 30
REACHED = 'diagram of J reached from J made by the band'


def jones_state_sum(pd):
    """Jones polynomial {exponent of t: coefficient} of a knot diagram (PD code, counterclockwise from the incoming
    under-strand) from the Kauffman bracket summed over all 2^n states: the A-smoothing of X[a, b, c, d] joins a to b
    and c to d, <D> = sum A^(#A - #B) (-A^2 - A^-2)^(loops - 1), V = (-A^3)^(-w) <D> with A = t^(-1/4), and w the
    writhe (orientation by traversal: the over-strand of a positive crossing enters at slot 3 when the under-strand
    enters at slot 0)."""
    n = len(pd)
    if n == 0:
        return {0: 1}
    labels = sorted({x for X in pd for x in X})
    pos = {x: i for i, x in enumerate(labels)}
    occ = {}
    for c, X in enumerate(pd):
        for s, x in enumerate(X):
            occ.setdefault(x, []).append((c, s))
    entering, (c, s) = set(), (0, 0)
    for _ in range(2 * n):                         # traverse the knot from slot 0 of crossing 0
        entering.add((c, s))
        out = pd[c][(s + 2) % 4]
        a, b = occ[out]
        c, s = b if a == (c, (s + 2) % 4) else a
    w = sum(1 if (c, 3) in entering else -1 for c in range(n))
    poly = Counter()
    for state in range(1 << n):
        parent = list(range(len(labels)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        na = 0
        for c, (a, b, cc, d) in enumerate(pd):
            if not (state >> c) & 1:
                na += 1
                parent[find(pos[a])] = find(pos[b]); parent[find(pos[cc])] = find(pos[d])
            else:
                parent[find(pos[a])] = find(pos[d]); parent[find(pos[b])] = find(pos[cc])
        loops = len({find(i) for i in range(len(labels))})
        term = Counter({na - (n - na): 1})
        for _ in range(loops - 1):
            t2 = Counter()
            for e, x in term.items():
                t2[e + 2] -= x; t2[e - 2] -= x
            term = t2
        poly.update(term)
    out = Counter()
    for e, x in poly.items():
        if x:
            assert (e - 3 * w) % 4 == 0
            out[-(e - 3 * w) // 4] += (-1) ** (w % 2) * x
    return {e: x for e, x in sorted(out.items()) if x}


def simplified(pd, limit):
    """A diagram of the same knot with at most limit crossings by spherogram's simplify('global'), or None."""
    import spherogram
    best = None
    for _ in range(SIMPLIFY_TRIES):
        L = spherogram.Link([list(x) for x in pd]); L.simplify('global')
        P = [tuple(x) for x in L.PD_code()]
        if best is None or len(P) < len(best):
            best = P
        if len(best) <= limit - 2:
            break
    return best if best is not None and len(best) <= limit else None


def first_step(c):
    """The band K -> J of a record, as a step of run_chain."""
    return dict(face=c['face'], darts=c['darts'], edges=c['edges'], face_edges=c.get('face_edges'), t=c['s'],
                after=c['J_pd'], components=1)


def ribbon_steps(bands):
    """The bands of a ribbon disc, as steps of run_chain; the last one ends in the unlink."""
    out = []
    for i, b in enumerate(bands):
        last = i == len(bands) - 1
        out.append(dict(face=b['face'], darts=b['darts'], edges=b['edges'], face_edges=b.get('face_edges'), t=b['r'],
                        after=None if last else b['after'], components=b['components'],
                        crossings=b.get('crossings') if last else None, final=last))
    return out


def relate_K(diags, idx, pd_K):
    """Reidemeister relation of the stored diagram idx of K to the KnotInfo diagram: directly, else through the stored
    diagrams it was made from (field parent)."""
    direct = reaches(diags[idx]['pd'], pd_K)
    if direct:
        return dict(ok=True, direct=direct)
    chain, i = [], idx
    while i != 0 and diags[i].get('parent') is not None:
        p = diags[i]['parent']
        r = reaches(diags[i]['pd'], diags[p]['pd'])
        chain.append([i, p, r])
        if not r:
            return dict(ok=False, direct=None, chain=chain)
        i = p
    if i != 0:
        r = reaches(diags[i]['pd'], pd_K)
        chain.append([i, 0, r])
    return dict(ok=bool(chain) and all(x[2] for x in chain), direct=None, chain=chain)


def job(item):
    """One certificate or KnotInfo-slice band: (knot, what, index, record, stored diagrams of K by index, KnotInfo PD
    of K, KnotInfo PD and smooth 4-genus of X or None)."""
    import snappy, spherogram
    name, what, idx, c, diags, pd_K, X = item
    t0 = time.time()
    out = dict(knot=name, what=what, index=idx, diagram=c['diagram'], diagram_source=c['diagram_source'])
    start = c['diagram_pd']
    assert diags[c['diagram']]['pd'] == start
    if c['diagram'] == 0:
        out['K_is_knotinfo_pd'] = [list(x) for x in start] == pd_K
    else:
        out['K_isometric_to_knotinfo'] = isometric(spherogram.Link(start).exterior(), spherogram.Link(pd_K).exterior())
        out['K_related_to_knotinfo'] = relate_K(diags, c['diagram'], pd_K)
    K_ok = bool(out.get('K_is_knotinfo_pd') or out.get('K_isometric_to_knotinfo'))
    K_diagrammatic = bool(out.get('K_is_knotinfo_pd') or (out.get('K_related_to_knotinfo') or {}).get('ok'))
    chains = []
    if what == 'moebius band':
        s = c.get('surface_of_J')
        if c['kind'] == 'band to the unknot':
            chains.append(('K', start, [dict(first_step(c), final=True, crossings=len(c['J_pd']))]))
        elif not s:
            chains.append(('K', start, [first_step(c)] + ribbon_steps(c['ribbon_bands'])))
        else:
            chains.append(('K', start, [first_step(c)]))
            chains.append(('J', s['diagram_pd'], ribbon_steps(s['ribbon_bands'])))
            if s.get('simplified_pd'):
                out['J_related'] = dict(J_to_simplified=reaches(c['J_pd'], s['simplified_pd']),
                                        diagram_to_simplified=reaches(s['diagram_pd'], s['simplified_pd']))
                J_rel = bool(out['J_related']['J_to_simplified'] and out['J_related']['diagram_to_simplified'])
            elif s['diagram'] == REACHED:
                out['J_related'] = dict(J_to_diagram=reaches(c['J_pd'], s['diagram_pd']))
                J_rel = bool(out['J_related']['J_to_diagram'])
            else:                                  # a diagram of X that the search related to J by a common diagram
                out['J_related'] = dict(J_to_diagram=meets(c['J_pd'], s['diagram_pd']) or
                                        reaches(c['J_pd'], s['diagram_pd']))
                J_rel = bool(out['J_related']['J_to_diagram'])
    else:                                          # KnotInfo-slice band: the first band only
        chains.append(('K', start, [first_step(c)]))
    out['chains'] = []
    for part, st, steps in chains:
        try:
            out['chains'].append(dict(part=part, **run_chain(st, steps)))
        except (AssertionError, IndexError, KeyError, ValueError) as e:
            out['chains'].append(dict(part=part, ok=False, error=f'{type(e).__name__}: {e}'))
    out['bands_ok'] = all(ch['ok'] for ch in out['chains'])
    if what == 'moebius band':
        final = out['chains'][-1]
        bands = sum(len(ch['steps']) for ch in out['chains'])
        out['b1'] = 1 - (final.get('discs', 0) - bands) if out['bands_ok'] else None
        out['ok'] = out['bands_ok'] and K_ok and out['b1'] == 1 and (not c.get('surface_of_J') or J_rel)
        out['diagrammatic'] = out['ok'] and K_diagrammatic
    else:
        pd_X, g4_X = X
        J = spherogram.Link([list(x) for x in c['J_pd']])
        E_J, E_X = J.exterior(), spherogram.Link(pd_X).exterior()
        signs = None
        for _ in range(8):
            try:
                signs = sorted({int(m.cusp_maps()[0].det()) for m in E_J.is_isometric_to(E_X, return_isometries=True)
                                if m.extends_to_link()})
                break
            except Exception:
                E_J.randomize(); E_X.randomize()
        out['isometry_signs'] = signs
        S = simplified(c['J_pd'], JONES_MAX)
        out['J_simplified_crossings'] = len(S) if S is not None else None
        VJ, VX = (jones_state_sum(S) if S is not None else None), jones_state_sum([tuple(x) for x in pd_X])
        mir = {-e: x for e, x in VX.items()}
        out['jones'] = (None if VJ is None else 'symmetric' if VJ == VX == mir else 'same' if VJ == VX else
                        'mirror' if VJ == mir else 'different')
        out['chirality'] = {'same': 'J = X', 'mirror': 'J = mirror image of X'}.get(out['jones'], out['jones'])
        out['smooth_four_genus_X'] = g4_X
        agree = (out['jones'] == 'same' and 1 in (signs or [])) or (out['jones'] == 'mirror' and -1 in (signs or [])) \
            or (out['jones'] == 'symmetric' and bool(signs))
        out['ok'] = out['bands_ok'] and K_ok and bool(signs) and agree and g4_X == '0'
        out['search_chirality'] = c['recheck']['chirality']
        out['chirality_agrees_with_search'] = out['chirality'] == c['recheck']['chirality'] or out['jones'] == 'symmetric'
    out['seconds'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    import snappy, spherogram, database_knotinfo
    from database_knotinfo import link_list
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    src = os.environ.get('REBUILD_IN', os.path.join(ROOT, 'data', 'mobius_open_search.json.gz'))
    dst = os.environ.get('REBUILD_OUT', os.path.join(ROOT, 'data', 'mobius_open_rebuild.json'))
    t_start = time.time()
    BY = {r['name']: r for r in link_list()[1:]}
    pd = lambda k: json.loads(BY[k]['pd_notation'])
    D = json.load(gzip.open(src, 'rt'))
    items = []
    for name, r in D['knots'].items():
        diags = {d['index']: d for d in r['diagrams']}
        for i, c in enumerate(r['certificates']):
            items.append((name, 'moebius band', i, c, diags, pd(name), None))
        for i, h in enumerate(r['knotinfo_slice']):
            items.append((name, 'KnotInfo-slice band', i, h, diags, pd(name),
                          (pd(h['J']), BY[h['J']]['smooth_four_genus'].strip())))
    res = []
    with Pool(workers) as p:
        for o in p.imap_unordered(job, items):
            res.append(o)
            print(o['knot'], o['what'], o['index'], 'ok', o['ok'], 'diagrammatic', o.get('diagrammatic'),
                  o.get('chirality', ''), f"({o['seconds']}s)", flush=True)
    by = {}
    for o in res:
        by.setdefault(o['knot'], []).append(o)
    T = {k for k, r in D['knots'].items() if r['group'] == 'target'}
    rep = {k: next((o for o in by.get(k, []) if o['what'] == 'moebius band' and o['index'] == 0), None) for k in T}
    mob = {k: o for k, o in rep.items() if o}
    ks = {k for k, r in D['knots'].items() if k in T and not r['certificates'] and r['knotinfo_slice']}
    ks_ok = {k for k in ks if any(o['what'] == 'KnotInfo-slice band' and o['ok'] for o in by.get(k, []))}
    summary = dict(
        moebius_bands=dict(targets=len(mob), rebuilt=sorted(k for k, o in mob.items() if o['ok']),
                           diagrammatic=sorted(k for k, o in mob.items() if o['diagrammatic']),
                           failing=sorted(k for k, o in mob.items() if not o['ok']),
                           b1=dict(Counter(o['b1'] for o in mob.values()))),
        all_certificates=[sum(o['ok'] for o in res if o['what'] == 'moebius band'),
                          sum(o['what'] == 'moebius band' for o in res)],
        knotinfo_slice_only=dict(targets=len(ks), rebuilt=sorted(ks_ok), failing=sorted(ks - ks_ok)),
        knotinfo_slice_bands=[sum(o['ok'] for o in res if o['what'] == 'KnotInfo-slice band'),
                              sum(o['what'] == 'KnotInfo-slice band' for o in res)],
        chirality=dict(Counter(o['chirality'] for o in res if o['what'] == 'KnotInfo-slice band')),
        chirality_disagreeing_with_search=[[o['knot'], o['index']] for o in res if o['what'] == 'KnotInfo-slice band'
                                           and not o['chirality_agrees_with_search']],
        controls={k: [o['ok'] for o in v] for k, v in by.items() if D['knots'][k]['group'] == 'control'})
    meta = dict(input=os.path.basename(src), jones_max=JONES_MAX, simplify_tries=SIMPLIFY_TRIES,
                versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__,
                              database_knotinfo=database_knotinfo.version(), python=sys.version.split()[0]),
                seconds=round(time.time() - t_start))
    json.dump(dict(meta=meta, summary=summary, records=sorted(res, key=lambda o: (o['knot'], o['what'], o['index']))),
              open(dst, 'w'), separators=(',', ':'), default=str)
    print('summary:', json.dumps(summary, default=str))
    print('time %.0fs' % (time.time() - t_start))
