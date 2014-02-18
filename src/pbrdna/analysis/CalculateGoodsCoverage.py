#! /usr/bin/env python

import sys

file_list = sys.argv[1:-1]
max_dist = float(sys.argv[-1])

header = ['Filename']
rows = {}
for filename in file_list:
    row = [filename]
    with open( filename ) as handle:
        for line in handle:
            parts = line.strip().split()
            dist = parts[0]
            if dist == 'unique':
                dist = '0.00'   
            if float(dist) > max_dist:
                continue
            if dist not in header:
                header.append( dist )
            counts = [len(p.split(',')) for p in parts[2:]]
            singles = sum([c for c in counts if c == 1])
            G = round(1 - (singles / float(sum(counts))), 3)
            row.append(str(G))
    rows[filename] = row

print ','.join( header )
for row in rows:
    print ','.join(rows[row])
