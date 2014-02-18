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
            total = len(counts)
            singles = len([c for c in counts if c == 1])
            doubles = len([c for c in counts if c == 2])
            C = round( total + ((singles * (singles - 1))/ (2*(doubles+1))), 3)
            row.append(str(int(C)))
    rows[filename] = row

print ','.join( header )
for row in rows:
    print ','.join(rows[row])
