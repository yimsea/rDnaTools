#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

from pbrdna.fasta.utils import fasta_names

def create_name_file( consensus_file, selected_file, output_file ):
    selected = read_selected_files( selected_file )
    data_files = read_file_data( consensus_file, selected )
    sequence_names = read_sequence_names( data_files )
    output_names( sequence_names, output_file )

def read_selected_files( selected_file ):
    files = set([])
    with open( selected_file ) as handle:
        for line in handle:
            filename = line.strip()
            files.add( filename )
    return files

def read_file_data( consensus_file, selected ):
    file_data = {}
    with open( consensus_file ) as handle:
        for line in handle:
            source, ref, consensus = line.strip().split()
            if consensus in selected:
                file_data[consensus] = source
            elif ref in selected:
                file_data[ref] = source
            elif source in selected:
                file_data[source] = source
    return file_data

def read_sequence_names( data_files ):
    sequence_names = {}
    for reference, source in data_files.iteritems():
        ref_names = fasta_names( reference )
        if len( ref_names ) == 1:
            assert ref_names[0] not in sequence_names
            sequence_names[ref_names[0]] = fasta_names( source )
        elif reference == source:
            for name in ref_names:
                assert name not in sequence_names
                sequence_names[name] = [name]
        else:
            raise ValueError
    return sequence_names

def output_names( sequence_names, output_file ):
    with open( output_file, 'w') as handle:
        for reference, names in sequence_names.iteritems():
            handle.write( "%s\t%s\n" % (reference, ','.join(names)) )

if __name__ == '__main__':
    import sys

    consensus_file = sys.argv[1]
    selected_file = sys.argv[2]
    output_file = sys.argv[3]

    create_name_file( consensus_file, selected_file, output_file )