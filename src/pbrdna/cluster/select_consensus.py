#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbrdna.fasta.utils import fasta_count

def select_consensus_files( consensus_file, output_file ):
    selected_files = []
    with open( consensus_file ) as handle:
        for line in handle:
            sequence_file, reference_file, consensus_file = line.strip().split()
            if consensus_file.endswith('None'):
                pass
            elif fasta_count( consensus_file ) == 1:
                selected_files.append( consensus_file )
            else:
                selected_files.append( reference_file )
    output_selected_consensus( output_file, selected_files )

def select_combined_sequences( consensus_file, output_file ):
    selected_files = []
    with open( consensus_file ) as handle:
        for line in handle:
            sequence_file, reference_file, consensus_file = line.strip().split()
            if consensus_file.endswith('None'):
                selected_files.append( sequence_file )
            elif:
                selected_files.append( consensus_file )
    output_selected_consensus( output_file, selected_files )

def output_selected_consensus( output_file, selected_files ):
    with open( output_file, 'w' ) as handle:
        for filename in selected_files:
            handle.write(filename + '\n')

if __name__ == '__main__':
    import sys

    consensus_file = sys.argv[1]
    output_file = sys.argv[2]

    select_combined_sequences( consensus_file, output_file )