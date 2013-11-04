#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import logging
from pbcore.io.BasH5IO import BasH5Collection
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbrdna.arguments import args, MIN_SNR
from pbrdna.log import initialize_logger

SNR = getattr(args, 'min_snr', MIN_SNR)

log = logging.getLogger(__name__)

def snr_filter(input_fastq, raw_data_file, output_fastq, min_snr=SNR):
    """
    Filter out sequences below a threshold of predicted accuracy
    """
    log.info("Filtering sequences below {0} Signal-To-Noise Ratio".format(min_snr))
    seq_count = 0
    pass_count = 0
    raw_data = BasH5Collection( raw_data_file )
    with FastqWriter( output_fastq ) as writer:
        for record in FastqReader( input_fastq ):
            seq_count += 1
            zmw_name = '/'.join( record.name.strip().split('/')[:2] )
            zmw = raw_data[zmw_name]
            zmw_snr = min( zmw.zmwMetric("HQRegionSNR") )
            print zmw_name, zmw_snr
            if zmw_snr >= min_snr:
                pass_count += 1
                writer.writeRecord( record )
    percentage = round(100.0*pass_count/seq_count)
    log.info("{0} sequences of {1} ({2}%) passed filtering".format(pass_count,
                                                                   seq_count,
                                                                   percentage))


if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]
    raw_data = sys.argv[2]
    output_file = sys.argv[3]
    min_snr = float(sys.argv[4])

    initialize_logger()
    snr_filter(input_file, raw_data, output_file, min_snr)