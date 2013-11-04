#! /usr/bin/env python

import logging
from numpy import mean
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbrdna.arguments import args, MIN_ACCURACY

ACCURACY = getattr(args, 'min_accuracy', MIN_ACCURACY)

log = logging.getLogger(__name__)

def quality_filter(input_fastq, output_fastq, min_accuracy=ACCURACY):
    """
    Filter out sequences below a threshold of predicted accuracy
    """
    log.info("Filtering sequences below {0}% predicted accuracy".format(100*min_accuracy))
    seq_count = 0
    pass_count = 0
    with FastqWriter( output_fastq ) as writer:
        for record in FastqReader( input_fastq ):
            seq_count += 1
            if predicted_accuracy(record) >= min_accuracy:
                pass_count += 1
                writer.writeRecord( record )
    percentage = round(100.0*pass_count/seq_count)
    log.info("{0} sequences of {1} ({2}%) passed filtering".format(pass_count,
                                                                   seq_count,
                                                                   percentage))

# Utility Functions
def predicted_accuracy(record):
    p_scores = convert_quality( record.quality )
    return round(mean(p_scores), 4)

def convert_quality(qv_scores):
    return [quality_to_p(qv) for qv in qv_scores]

def quality_to_p(qv):
    return 1-(10**(-1*qv/10.0))


if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_accuracy = float(sys.argv[3])

    quality_filter(input_file, output_file, min_accuracy)