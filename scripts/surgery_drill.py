#!/usr/bin/env python
"""[C] Section 3.2, with Proposition 3.2(c): surgery descriptions Sigma(K) = S^3_{p/q}(C) from short geodesics, and the
test of [C] Proposition 2.6, for the knots of Proposition 3.2(c) whose double branched cover is hyperbolic and bounds
no negative-definite plumbing tree and whose unknotting number is not one in KnotInfo.

The knots: those left open before Proposition 3.2(c) (the set oleft3 of paper_counts.py, from open_sets.upto('oleft3'))
that are neither alternating nor Montesinos, whose Sigma(K) is hyperbolic (data/plumbing_dinv.json, class hyperbolic),
and whose unknotting number is not 1 in KnotInfo (those are treated by surgery_dinv.py).
For each knot, Sigma(K) is built as in surgery_dinv.sigma_cover and its dual curves with at most
surgery_dinv.DRILL_SEGMENTS segments are drilled in turn, at most MAX_CURVES of them, by surgery_dinv.drill_route: if
the drilled manifold X has H_1 = Z and a filling of X that Regina recognises as S^3, then X is the exterior of a knot
kappa' in S^3 and Sigma(K) = S^3_{p/q}(kappa'), with the checks of drill_route (closed isometry, Casson-Walker
invariant from the signature and the Jones polynomial, the linking form against V + V^T, a diagram of kappa' by
exterior_to_link whose chirality is fixed by an isometry, V_j from the knot Floer homology of kappa' when it is that of
an L-space knot, of its mirror or thin, checked on the model complex).  The test of Proposition 2.6 is run with the
two implementations of surgery_dinv.py (rank1_dinv.test and montesinos_dinv_check.obstructed), and the d-functions of
two routes, when there are two, must agree.  At most MAX_FORMS routes per knot are carried to the end, and no route is
started after KNOT_SECONDS.
Usage: python scripts/surgery_drill.py [workers] [redo]
  redo: run the knots whose record in data/surgery_drill.json is an error again.
Output: data/surgery_drill.json (checkpoint data/surgery_drill.partial while it runs)
"""
import json, os, random, sys, time, traceback, warnings
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')

MAX_CURVES = 30                           # dual curves drilled per knot
MAX_FORMS = 2                             # routes carried to the end
KNOT_SECONDS = 1800                       # no new route after this time
OUT = os.path.join(ROOT, 'data', 'surgery_drill.json')
PARTIAL = os.path.join(ROOT, 'data', 'surgery_drill.partial')


def targets():
    from open_sets import upto
    g = upto('oleft3')
    pdv, su = g['pdv'], g['su']
    return sorted(k for k in g['oleft3'] if not g['alt'](k) and not g['mont'](k) and k not in su
                  and pdv['targets'][k]['cls'] == 'hyperbolic')


def one(name):
    import snappy
    import surgery_dinv as SD
    import montesinos_dinv_check as MC
    from database_knotinfo import link_list
    from cypari import pari
    pari.allocatemem(2 ** 28, 2 ** 31, silent=True)
    BY = {r['name']: r for r in link_list()[1:]}
    R1 = json.load(open(os.path.join(ROOT, 'data', 'rank1_targets.json')))
    snappy.set_rand_seed(SD.zlib_seed(name)); random.seed(SD.zlib_seed(name))
    t0 = time.time()
    pd = json.loads(BY[name]['pd_notation'])
    info = SD.knot_info(BY, name, R1[name]['s'])
    rec = dict(knot=name, det=info['det'], s=info['s'], routes=[])
    try:
        SigK = SD.sigma_cover(pd)
        rec['H1'] = str(SigK.homology())
        if rec['H1'] != 'Z/%d' % info['det']:          # not cyclic: Sigma(K) is no surgery on a knot
            rec['verdict'] = 'H_1 not cyclic'
            rec['seconds'] = round(time.time() - t0, 1)
            return rec
        R = SD.mullins_R(info['det'], info['sigma'], info['jones'])
        S = SD.improved(SigK, 100)
        rec['solution_type'] = S.solution_type()
        curves = S.dual_curves(max_segments=SD.DRILL_SEGMENTS)
        rec['curves'] = len(curves)
        forms = []
        for i in range(min(MAX_CURVES, len(curves))):
            if len(forms) == MAX_FORMS or time.time() - t0 > KNOT_SECONDS:
                break
            r, FZ = SD.drill_route(S, i, info, R)
            r.pop('traceback', None)
            rec['routes'].append(r)
            if FZ is not None:
                forms.append((i, FZ))
        if forms:
            i0, FZ = forms[0]
            rec['route_used'] = i0
            rec['routes_agree'] = [bool(SD.same_form(FZ, F)) for _, F in forms[1:]]
            assert all(rec['routes_agree']), 'the d-functions of two routes differ'
            rec['passing'], rec['passing_signed'] = SD.verdicts(FZ, info['s'])
            rec['obstructed'] = rec['passing'] is None
            ob2, _ = MC.obstructed(FZ.D, FZ.adj[0][0], {k: FZ.d[(k,)] for k in range(FZ.D)})
            rec['obstructed_independent_code'] = ob2
            assert ob2 == rec['obstructed'], 'the two implementations of the d-test differ'
            r0 = next(r for r in rec['routes'] if r.get('curve') == i0)
            rec.update(surgery=r0['surgery'], V_source=r0['V_source'], kappa_identify=r0['kappa'].get('identify'),
                       kappa_crossings=len(r0['kappa']['diagram']['pd']))
            rec['verdict'] = 'obstructed' if rec['obstructed'] else 'passes'
        else:
            s3 = [r for r in rec['routes'] if r.get('s3_slope')]
            rec['verdict'] = 'no surgery description' if not s3 else 'no d-function'
    except Exception as ex:
        rec.update(verdict='error', error=repr(ex)[:500], traceback=traceback.format_exc()[-1500:])
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def child(name, q):
    try:
        q.put(one(name))
    except BaseException as ex:
        q.put(dict(knot=name, routes=[], verdict='error', error=repr(ex)[:500]))


def main(workers, redo=False):
    import multiprocessing as mp
    import snappy, spherogram
    names = targets()
    if redo:                                  # run the knots whose record is an error again
        old = json.load(open(OUT))['knots']
        with open(PARTIAL, 'w') as f:
            for r in old.values():
                if r['verdict'] != 'error':
                    f.write(json.dumps(r, default=str) + '\n')
    done = {}
    if os.path.exists(PARTIAL):
        done = {r['knot']: r for r in map(json.loads, open(PARTIAL))}
    todo = [n for n in names if n not in done]
    ctx = mp.get_context('spawn')
    running, t_start = {}, time.time()
    with open(PARTIAL, 'a') as f:
        while todo or running:
            while todo and len(running) < workers:
                n = todo.pop(0); q = ctx.Queue(); p = ctx.Process(target=child, args=(n, q)); p.start()
                running[n] = (p, q, time.time())
            time.sleep(0.5)
            for n in list(running):
                p, q, t0 = running[n]
                rec = None
                try:
                    rec = q.get_nowait()
                except Exception:
                    if time.time() - t0 > 3 * KNOT_SECONDS:
                        p.terminate(); rec = dict(knot=n, routes=[], verdict='error', error='time limit')
                    elif not p.is_alive() and time.time() - t0 > 10:
                        rec = dict(knot=n, routes=[], verdict='error', error='process died')
                if rec is None:
                    continue
                p.join(5); del running[n]
                f.write(json.dumps(rec, default=str) + '\n'); f.flush()
                done[n] = rec
                print(n, rec['verdict'], rec.get('surgery', ''), rec.get('V_source', ''), rec['seconds'], flush=True)
    recs = {k: done[k] for k in names}
    summary = dict(knots=len(names), verdicts=dict(Counter(r['verdict'] for r in recs.values())),
                   with_surgery_description=sorted(k for k, r in recs.items() if any(x.get('s3_slope') for x in r['routes'])),
                   obstructed=sorted(k for k, r in recs.items() if r['verdict'] == 'obstructed'),
                   passes=sorted(k for k, r in recs.items() if r['verdict'] == 'passes'),
                   q=dict(Counter(x['slope_p_q'][1] for r in recs.values() for x in r['routes'] if x.get('slope_p_q'))),
                   h1_not_cyclic=sorted(k for k, r in recs.items() if r['verdict'] == 'H_1 not cyclic'),
                   errors=sorted(k for k, r in recs.items() if r['verdict'] == 'error'))
    import surgery_dinv as SD
    meta = dict(max_curves=MAX_CURVES, max_forms=MAX_FORMS, knot_seconds=KNOT_SECONDS, drill_segments=SD.DRILL_SEGMENTS,
                hfk_seconds=SD.HFK_SECONDS, seconds=round(time.time() - t_start),
                versions=dict(snappy=snappy.__version__, spherogram=spherogram.__version__))
    json.dump(dict(meta=meta, summary=summary, knots=recs), open(OUT, 'w'), separators=(',', ':'), default=str)
    os.remove(PARTIAL)
    print('summary:', json.dumps(summary))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6, sys.argv[2:3] == ['redo'])
