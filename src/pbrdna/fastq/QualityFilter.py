#! /usr/bin/env python

from pbcore.io.FastqIO import FastqReader, FastqWriter
from numpy import mean, log10

MIN_ACCURACY = 0.99

class QualityFilter( object ):
    """
    A tool for filtering out low-quality fastq files
    """
    def __init__(self, input_fastq, output_fastq, min_accuracy=MIN_ACCURACY):
        self.input_reader = FastqReader(input_fastq)
        self.output_writer = FastqWriter(output_fastq)
        self.min_accuracy = min_accuracy

    def __call__(self):
        for fastq in self.input_reader:
            if predicted_accuracy(fastq) >= self.min_accuracy:
                self.output_writer.writeRecord( fastq )


# Utility Functions
def quality_to_p(qv):
    return 1-(10**(-1*qv/10.0))

def convert_quality(qv_scores):
    return [quality_to_p(qv) for qv in qv_scores]

def predicted_accuracy(record):
    p_scores = convert_quality( record.quality )
    return round(mean(p_scores), 4)


if __name__ == '__main__':
    import argparse, sys

    desc = "A tool for filtering Fastq sequences based on mean predicted accuracy"
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add("input_file",
        metavar="FASTQ",
        help="A Fastq file of PacBio CCS Sequence data to filter")
    add("-o", "--output_file",
        metavar="FILE",
        default=sys.stdout,
        help="Location to output the filtered sequence files to [STDOUT]")
    add("-m", "--min_accuracy",
        type=float,
        metavar="FLOAT",
        default=MIN_ACCURACY,
        help="Minimum average predicted accuracy to filter on [%s]" % MIN_ACCURACY)
    args = parser.parse_args()

    filter_tool = QualityFilter( args.input_file, 
                                 args.output_file,
                                 args.min_accuracy )
    filter_tool()
