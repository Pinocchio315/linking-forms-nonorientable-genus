#!/usr/bin/env python
"""The sets of scripts/paper_counts.py up to a stage, for the searches that take their knots from them.

paper_counts.py checks every count against all the outputs, while these searches take their knots from sets that it
computes before it reads their own outputs and those of the searches after them.  So upto(stop) runs paper_counts.py
only up to such a set: with the environment variable PAPER_COUNTS_UPTO=stop it stops there (SystemExit), and the
globals defined so far are returned.  The stops:
- newd: the knots that the first search and [C] Propositions 3.1 and 3.2(a), (b) leave open are left0 - newd.
  mobius_open_search.py ([C] Section 3.3) starts here.
- oleft3: the knots left open before [C] Proposition 3.2(c), left0 - newd - wide, which needs the outputs of
  mobius_open_search.py and mobius_open_rebuild.py.  donaldson_surfaces.py and surgery_drill.py ([C] Proposition 3.2(c)
  and Section 3.2) and mobius_J_search.py ([C] Section 3.3, after donaldson_surfaces.py) start here.
So each of these searches needs only the outputs of the stages before it.
Usage: from open_sets import upto (a module; it writes nothing)
"""
import contextlib, io, os

HERE = os.path.dirname(os.path.abspath(__file__))
STOPS = ('newd', 'oleft3')


def upto(stop):
    assert stop in STOPS, stop
    path = os.path.join(HERE, 'paper_counts.py')
    g = {'__name__': 'paper_counts', '__file__': path}
    os.environ['PAPER_COUNTS_UPTO'] = stop
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(open(path).read(), path, 'exec'), g)
    except SystemExit:
        pass
    finally:
        os.environ.pop('PAPER_COUNTS_UPTO', None)
    assert stop in g, f'paper_counts.py did not reach {stop}'
    return g
