#!/usr/bin/env python
"""Exact V_j(C) for the surgery descriptions of scripts/surgery_drill.py whose d-function that script left open, and
the test of [C] Proposition 2.6 with them (the knots of range [1,2] left open by [C]).

Targets.  The knots of data/surgery_drill.json with verdict 'no d-function': drilling a dual curve of Sigma_K gave the
exterior of a knot kappa' in S^3, so Sigma_K = S^3_{p/q}(kappa'), but surgery_drill.py determines V_j(kappa') only when
the knot Floer homology of kappa' is that of an L-space knot, of its mirror, or thin.  Here V_j(C) comes from the UV = 0
knot Floer complex of C over every lift to F_2[U,V], by steps 1-3 of scripts/surgery_cfk.py, and the test is that of
surgery_cfk.one (steps 5 and 6: both implementations of the test, the checks of step 5, and the agreement of the
d-functions of two routes).  The chirality of kappa' is fixed by the isometry of exteriors that surgery_drill.py stored
(orientation 'diagram' or 'mirror').
Controls.  The eight knots of data/surgery_drill.json with a verdict (five pass, three are excluded), whose V_j were
determined there from staircase or thin complexes: the exact method must give the same verdict.  The tensor controls of
surgery_cfk.py (connected sums of torus knots) are run as well.
Usage: python scripts/drill_cfk.py [workers]
Output: data/drill_cfk.json (checkpoint data/drill_cfk.partial for the runs of Szabo's program)
"""
import json, os, sys, warnings
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import surgery_cfk as SC
import surgery_dinv as sd

OUT = os.path.join(ROOT, 'data', 'drill_cfk.json')
PARTIAL = os.path.join(ROOT, 'data', 'drill_cfk.partial')


def as_surgery_record(BY, rec):
    """A record of surgery_drill.py in the form that surgery_cfk.one reads: no crossing routes, its drilling routes."""
    r = BY[rec['knot']]
    R = sd.mullins_R(rec['det'], int(r['signature']), r['jones_polynomial'])
    return dict(det=rec['det'], s=rec['s'], R=str(R), routes=[], drill_routes=[x for x in rec['routes'] if x.get('kappa_hat')])


if __name__ == '__main__':
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    DR = json.load(open(os.path.join(ROOT, 'data', 'surgery_drill.json')))['knots']
    parts = dict(targets={n: as_surgery_record(BY, r) for n, r in DR.items() if r['verdict'] == 'no d-function'},
                 controls={n: as_surgery_record(BY, r) for n, r in DR.items() if r['verdict'] in ('passes', 'obstructed')})
    assert len(parts['targets']) == 13 and len(parts['controls']) == 8
    jobs = SC.jobs_of(parts, BY)
    cxs = SC.run_jobs(jobs, workers, PARTIAL)
    covers, res = {}, {}
    for part, recs in parts.items():
        res[part] = {}
        for n, rec in sorted(recs.items()):
            res[part][n] = SC.one(n, rec, SC.knot_inputs(BY, n, rec), cxs, covers)
            print(part, n, res[part][n]['verdict'], flush=True)
    ten = SC.tensor_controls(cxs, BY)
    expected = {n: {'passes': 'passes', 'obstructed': 'excluded'}[DR[n]['verdict']] for n in parts['controls']}
    T = res['targets']
    problems = dict(
        self_test=[k for k, v in SC.self_test().items() if not v],
        controls={n: v['verdict'] for n, v in res['controls'].items() if v['verdict'] != expected[n]},
        tensor=[k for k, t in ten.items() if not (t.get('agrees') and t.get('same_hat'))])
    summary = dict(targets=len(T), excluded=sorted(n for n, v in T.items() if v['verdict'] == 'excluded'),
                   passes=sorted(n for n, v in T.items() if v['verdict'] == 'passes'),
                   undetermined={n: SC.row(v) for n, v in T.items() if v['verdict'] == 'undetermined'},
                   table={n: SC.row(v) for n, v in sorted(T.items())},
                   kinds=dict(Counter(SC.row(v).get('kind') for v in T.values())),
                   controls={n: dict(verdict=v['verdict'], expected=expected[n]) for n, v in res['controls'].items()},
                   problems=problems)
    for part in res.values():
        for v in part.values():
            v.pop('form', None)
    meta = dict(versions=sd.versions(), parameters=dict(cfk_seconds=SC.CFK_SECONDS, max_slots=SC.MAX_SLOTS), jobs=len(jobs),
                jobs_not_ok={k: v['status'] for k, v in cxs.items() if k in jobs and v['status'] != 'ok'})
    out = dict(meta=meta, summary=summary, targets=T, controls=res['controls'], tensor=ten,
               complexes={k: dict(generators=v['generators'], differentials=[[x, y] for x, y, c in v['differentials']])
                          for k, v in cxs.items() if k in jobs and v['status'] == 'ok'})
    json.dump(out, open(OUT, 'w'), indent=0, default=str)
    print(json.dumps({k: summary[k] for k in ('excluded', 'passes', 'undetermined', 'kinds', 'controls', 'problems')},
                     indent=1, default=str))
    assert not any(v for v in problems.values() if v), problems
