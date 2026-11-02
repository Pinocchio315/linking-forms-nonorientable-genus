#!/usr/bin/env python
"""Compare KnotInfo data files on the prime knots with at most 13 crossings (README, Input data).

Each argument is a file knotinfo_data_complete.csv of KnotInfo (columns separated by '|', the second line holding
column descriptions): for instance the data file of version 2026.10.5, the file csv_data/knotinfo_data_complete.csv of
the package database_knotinfo 2026.8.1, and that of the package of version 2026.11.2.  Every file is compared with the
last one: the script prints the knots and columns that only one of them has, and each column in which some knot differs,
ignoring spaces, with the number of such knots and the first few of them.  For the unknotting number it also prints
whether a value that differs is 1 in either file.
Usage: python -I scripts/knotinfo_compare.py OLD.csv [OLD2.csv ...] NEW.csv   (prints; writes nothing)
"""
import csv, re, sys

csv.field_size_limit(10 ** 9)


def read(path):
    with open(path, newline='', encoding='utf-8', errors='replace') as f:
        r = csv.reader(f, delimiter='|')
        head = next(r)
        next(r)
        rows = {}
        for row in r:
            if row:
                d = dict(zip(head, row))
                if d['name'] != '0_1' and int(re.match(r'\d+', d['name']).group()) <= 13:
                    rows[d['name']] = {k: re.sub(r'\s+', '', v) for k, v in d.items()}
    return head, rows


def main(paths):
    H, R = read(paths[-1])
    for path in paths[:-1]:
        h, r = read(path)
        print(f'{path}\n  against {paths[-1]}: {len(r)} and {len(R)} knots')
        print(f'  knots only in the first: {sorted(set(r) - set(R))}, only in the second: {sorted(set(R) - set(r))}')
        print(f'  columns only in the first: {[c for c in h if c not in H]}, only in the second: {[c for c in H if c not in h]}')
        for c in [c for c in h if c in H]:
            ks = sorted(k for k in r if k in R and r[k][c] != R[k][c])
            if ks:
                ones = sum('1' in (r[k][c], R[k][c]) for k in ks) if c == 'unknotting_number' else None
                print(f'  {c}: {len(ks)} knots, e.g. {[(k, r[k][c][:20], R[k][c][:20]) for k in ks[:3]]}'
                      + (f'; with the value 1 in either file: {ones}' if ones is not None else ''))


if __name__ == '__main__':
    main(sys.argv[1:])
