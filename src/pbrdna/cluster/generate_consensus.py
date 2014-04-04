#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbrdna.resequence.DagConTools import DagConRunner

def generate_consensus_files( cluster_list, consensus_tool, output_file ):
    consensus_files = []
    with open( cluster_list ) as handle:
        for line in handle:
            sequence_file, reference_file, count = line.strip().split()
            if reference_file.endswith('None'):
                consensus_files.append( (sequence_file, 'None') )
            else:
                consensus = consensus_tool( sequence_file, reference_file )
                consensus_files.append( (reference_file, consensus) )
    write_consensus_files( consensus_files, output_file )

def write_consensus_files( consensus_files, output_file ):
    with open( output_file, 'w' ) as handle:
        for filename_pair in consensus_files:
            handle.write('%s\t%s\n' % filename_pair)

if __name__ == '__main__':
    import sys

    cluster_list = sys.argv[1]
    consensus_tool = DagConRunner('gcon.py', 'r')
    output_file = sys.argv[2]

    generate_consensus_files( cluster_list, consensus_tool, output_file )