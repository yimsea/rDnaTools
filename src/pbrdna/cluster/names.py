#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

def create_name_file( consensus_file, selected_file, output_file ):
    selected = read_selected_files( selected_file )


def read_selected_files( selected_file ):
    return set([l.strip() for l in open(selected_file)])

if __name__ == '__main__':
    import sys

    consensus_file = sys.argv[1]
    selected_file = sys.argv[2]
    output_file = sys.argv[3]

    create_name_file( consensus_file, selected_file, output_file )