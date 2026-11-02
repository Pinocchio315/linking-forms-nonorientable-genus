#!/usr/bin/env python
"""The KnotInfo values and ranges of gamma_4 ([C] Section 3.1), extracted from the data file of KnotInfo 2026.10.5.

The input is knotinfo_data_complete.csv of KnotInfo 2026.10.5 (not redistributed; see the README, Input data), with
columns separated by '|'; its second line holds column descriptions, not a knot.  The script writes
- gamma4_knotinfo_2026.10.5.json: the field smooth_4d_crosscap_number, without spaces, of every knot but 0_1 (a value
  such as '2', or a range such as '[1,3]');
- open_gamma4_knotinfo_2026.10.5.json: the fields KEEP of the knots for which that field is a range.
KnotInfo records the least b_1 of a nonorientable surface in B^4, which is 1 for slice knots ([C] Section 3.1).
The CSV is untrusted data: run this script with python -I.
Usage: python -I scripts/knotinfo_extract.py path/to/knotinfo_data_complete.csv [out_dir]
Output: the two files above, in out_dir (default data/)
"""
import csv, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
KEEP = ['name', 'crossing_number', 'alternating', 'signature', 'determinant', 'arf_invariant', 'unknotting_number',
        'smooth_four_genus', 'smooth_4d_crosscap_number', 'topological_4d_crosscap_number', 'crosscap_number',
        'montesinos_notation']

if __name__ == '__main__':
    path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'data')
    csv.field_size_limit(10 ** 9)
    rows = list(csv.DictReader(open(path), delimiter='|'))[1:]
    full = {r['name']: r['smooth_4d_crosscap_number'].replace(' ', '') for r in rows if r['name'] != '0_1'}
    json.dump(full, open(os.path.join(out_dir, 'gamma4_knotinfo_2026.10.5.json'), 'w'), indent=0, sort_keys=True)
    opened = {}
    for r in rows:
        if '[' in (r.get('smooth_4d_crosscap_number') or '').replace(' ', ''):
            opened[r['name']] = {k: r.get(k, '') for k in KEEP}
    json.dump(opened, open(os.path.join(out_dir, 'open_gamma4_knotinfo_2026.10.5.json'), 'w'), indent=1)
    print(len(full), 'knots,', len(opened), 'open')
