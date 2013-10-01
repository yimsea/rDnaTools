#! /usr/bin/env python
import csv
import sys
import logging

from collections import namedtuple
from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqReader, FastqWriter

barcode = namedtuple('barcode', 'id strand seen5 seenA seen3 end5 endA end3 primer')

log = logging.getLogger()

class SequenceSeparator( object ):

    def __init__( self, input_file, barcode_file, prefix=None, filetype=None ):
        self.input_file = input_file
        self.barcode_file = barcode_file
        self.prefix = prefix or get_prefix( input_file )
        self.filetype = filetype or get_filetype( input_file )
        self.writers = {}

    def run( self ):
        self.parse_barcode_data()
        self.parse_group_list()
        self.separate_sequences()

    def parse_barcode_data( self ):
        self.groups = {}
        with open( self.barcode_file ) as handle:
            for entry in map(barcode._make, csv.reader(handle, delimiter='\t')):
                if entry.id == 'ID':
                    continue
                self.groups[entry.id] = entry.primer

    def parse_group_list( self ):
        self.group_list = set()
        for read, group in self.groups.iteritems():
            self.group_list.add( group )

    def separate_sequences( self ):
        # Open the appropriate Sequence Reader
        if self.filetype == 'fasta':
            reader = FastaReader( self.input_file )
        elif self.filetype == 'fastq':
            reader = FastqReader( self.input_file )
        # Iterate through records, writing out the
        for record in reader:
            try:
                group = self.groups[record.name]
            except:
                continue
            self.write_record( record, group )

    def write_record( self, record, group ):
        if group not in self.writers:
            self.add_writer( group )
        self.writers[group].writeRecord( record )

    def add_writer( self, group ):
        if self.filetype == 'fasta':
            output_file = '%s.g%s.fasta' % (self.prefix, group)
            self.writers[group] = FastaWriter( output_file )
        if self.filetype == 'fastq':
            output_file = '%s.g%s.fastq' % (self.prefix, group)
            self.writers[group] = FastqWriter( output_file )

def get_prefix( filename ):
    return '.'.join( filename.split('.')[:-1] )

def get_filetype( filename ):
    if (filename.lower().endswith( '.fa' ) or
        filename.lower().endswith( '.fsa' ) or
        filename.lower().endswith( '.fasta' )):
        return 'fasta'
    elif (filename.lower().endswith( '.fq' ) or
          filename.lower().endswith( '.fastq' )):
        return 'fastq'
    else:
        msg = 'Input file is not a recognized filetype!'
        log.error( msg )
        raise TypeError( msg )

if __name__ == '__main__':
    log.basicConfig( level=log.INFO )

    sequence_file = sys.argv[1]
    barcode_file = sys.argv[2]

    SequenceSeparator( sequence_file, barcode_file ).run()
