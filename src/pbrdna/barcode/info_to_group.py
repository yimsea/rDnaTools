#! /usr/bin/env python
import csv, sys

from collections import namedtuple

primer_info = namedtuple('primer_info', 'id strand seen5 seenA seen3 end5 endA end3 primer')

def info_to_group( filename, output=sys.stdout ):
    with open( filename ) as handle:
        for info in map(primer_info._make, csv.reader(handle, delimiter='\t')):
            if info.id == 'ID':
                continue
            if info.primer == 'NA':
                continue
            output.write( '{0}\tG{1}\n'.format(info.id, info.primer) )

if __name__ == '__main__':
    for filename in sys.argv[1:]:
        info_to_group( filename )
