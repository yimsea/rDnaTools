#! /usr/bin/env python
import csv, sys, log

from collections import namedtuple
from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord
from pbcore.io.FastqIO import FastqReader, FastqWriter, FastqRecord

barcode = namedtuple('barcode', 'id strand seen5 seenA seen3 end5 endA end3 primer')

log = log.getLogger()

class BarcodeTrimmer( object ):

    def __init__( self, input_file, barcode_file, prefix=None, filetype=None ):
        self.input_file = input_file
        self.barcode_file = barcode_file
        self.prefix = prefix or get_prefix( input_file )
        self.filetype = filetype or get_filetype( input_file )
        self.positions = {}

    def run( self ):
        self.parse_barcode_data()
        self.open_reader()
        self.open_writer()
        self.trim_sequences()

    def parse_barcode_data( self ):
        with open( self.barcode_file ) as handle:
            for entry in map(barcode._make, csv.reader(handle, delimiter='\t')):
                if entry.id == 'ID':
                    continue
                start = None if entry.end5 == 'NA' else int(entry.end5)
                end = None if entry.end3 == 'NA' else int(entry.end3)
                self.positions[entry.id] = (start, end)

    def open_reader( self ):
        if self.filetype == 'fasta':
            self.reader = FastaReader( self.input_file )
        elif self.filetype == 'fastq':
            self.reader = FastqReader( self.input_file )

    def open_writer( self ):
        if self.filetype == 'fasta':
            output_file = '%s.trim.fasta' % self.prefix
            self.writer = FastaWriter( output_file )
        elif self.filetype == 'fastq':
            output_file = '%s.trim.fastq' % self.prefix
            self.writer = FastqWriter( output_file )

    def trim_sequences( self ):
        for record in self.reader:
            try:
                start, end = self.positions[record.name]
            except:
                msg = 'Unknown sequence record "%s"!' % record.name
                log.error( msg )
                raise ValueError( msg )
            trimmed_record = trim_record( record, start, end )
            self.writer.writeRecord( trimmed_record )

def trim_record( record, start, end ):
    if isinstance(record, FastaRecord):
        return trim_fasta_record( record, start, end )
    elif isinstance(record, FastqRecord):
        return trim_fastq_record( record, start, end )
    else:
        msg = 'Unrecognized record type "%s"' % type(record)
        log.error( msg )
        raise TypeError( msg )

def trim_fasta_record( record, start, end ):
    if start is None and end is None:
        trimmed_sequence = record.sequence
    elif start is None:
        trimmed_sequence = record.sequence[:end]
    elif end is None:
        trimmed_sequence = record.sequence[start:]
    else:
        trimmed_sequence = record.sequence[start:end]
    return FastaRecord( record.name,
                        trimmed_sequence )

def trim_fastq_record( record, start, end ):
    if start is None and end is None:
        trimmed_sequence = record.sequence
        trimmed_quality = record.quality
    elif start is None:
        trimmed_sequence = record.sequence[:end]
        trimmed_quality = record.quality[:end]
    elif end is None:
        trimmed_sequence = record.sequence[start:]
        trimmed_quality = record.quality[start:]
    else:
        trimmed_sequence = record.sequence[start:end]
        trimmed_quality = record.quality[start:end]
    return FastqRecord( record.name,
                        trimmed_sequence,
                        trimmed_quality )

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

    BarcodeTrimmer( sequence_file, barcode_file ).run()
