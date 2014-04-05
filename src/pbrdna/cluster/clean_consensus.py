#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import os
from pbrdna.utils import write_dummy_file

def clean_consensus_outputs( dir_path, output_file ):
    for filename in os.listdir( dir_path ):
        file_path = os.path.join( dir_path, filename )
        if file_path.endswith('_input.fa'):
            os.remove( file_path )
        elif file_path.endswith('_input.fa.aln'):
            os.remove( file_path )
        elif file_path.endswith('_input.fa.aln_unsorted'):
            os.remove( file_path )
    write_dummy_file( output_file )

if __name__ == '__main__':
    import sys

    reseq_folder = sys.argv[1]
    output_file = sys.argv[2]

    clean_consensus_outputs( reseq_folder, output_file )