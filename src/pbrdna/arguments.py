import os, argparse

from . import __VERSION__
from utils import validate_int, validate_float

NPROC = 1
MIN_DIST = 0.001
DIST = 0.03
MAX_DIST = 0.5
MIN_ACCURACY = 0.99
MIN_QV = 15
FRACTION = 0.8
MIN_LENGTH = 500
MIN_RATIO = 0.5
PRECLUSTER_DIFFS = 4
MIN_CLUSTER_SIZE = 3
MIN_SNR = 3
CLUSTER_METHODS = ('nearest', 'average', 'furthest')
DEFAULT_METHOD = 'average'

args = argparse.Namespace()

def parse_args():
    """
    Parse the options for running the HLA pipeline and
    """
    desc = 'A pipeline tool for analyzing PacBio sequenced rRNA amplicons'
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add('input_file',
        metavar='FILE',
        help="File of rRNA sequencing data to use")
    add('-r', '--raw_data',
        metavar='FILE',
        help='BasH5, BaxH5 or FOFN of raw H5-format sequence data')
    add('-a', '--min_accuracy',
        type=float, 
        metavar='FLOAT',
        default=MIN_ACCURACY,
        help='Minimum predicted sequence accuracy to allow (%s)' % MIN_ACCURACY)
    add('-l', '--min_length',
        type=int, 
        metavar='INT', 
        default=MIN_LENGTH,
        help='Minimum length sequence to allow (%s)' % MIN_LENGTH)
    add('-s', '--min_snr',
        type=float,
        metavar='FLOAT',
        default=MIN_SNR,
        help='Minimum Signal-to-Noise ratio to allow (%s)' % MIN_SNR)
    add('-d', '--distance', 
        type=float, 
        metavar='FLOAT', 
        default=DIST,
        help='Distance at which to cluster sequences (%s)' % DIST)
    add('-n', '--num_processes', 
        type=int,
        metavar='INT',
        default=NPROC, 
        dest='nproc', 
        help='Number of processors to use (%s)' % NPROC)
    add('-f', '--fraction', 
        type=float, 
        metavar='FLOAT', 
        default=FRACTION,
        help='Fraction of full-length to require of each read (%s)' % FRACTION)
    add('-o', '--output', 
        dest='output_dir', 
        metavar='DIR',
        default='rna_pipeline_run',
        help="Specify the output folder")
    add('-q', '--min_qv',
        type=int, 
        metavar='INT',
        default=MIN_QV,
        help='Minimum QV to allow after sequence masking (%s)' % MIN_QV)
    add('-c', '--min_cluster_size', 
        type=int, 
        metavar='INT',
        default=MIN_CLUSTER_SIZE,
        help='Minimum cluster to generate consensus sequences (%s)' % MIN_CLUSTER_SIZE)
    add('--clustering_method', 
        metavar='METHOD',
        dest='clusteringMethod', 
        default=DEFAULT_METHOD,
        choices=CLUSTER_METHODS,
        help="Distance algorithm to use in clustering (%s)" % DEFAULT_METHOD)
    add('--precluster_diffs', 
        type=int, 
        metavar='INT',
        default=PRECLUSTER_DIFFS,
        help='Maximum number of differences to allow in pre-clustering (%s)' % PRECLUSTER_DIFFS)
    add('-A', '--alignment_reference', 
        metavar='REF',
        default='silva.both.align',
        help="Reference MSA for aligning query sequences")
    add('-C', '--chimera_reference', 
        metavar='REF',
        default='silva.gold.align',
        help="Reference MSA for Chimera detection")
    add('--enable_masking', 
        action='store_true',
        help="Turn off the low-quality Masking step")
    add('--sub_cluster',
        action="store_true",
        help="Subcluster each OTU to separate individual rDNA alleles")
    add('--disable_clustering', 
        action='store_false',
        dest='enable_clustering',
        help="Turn off the Clustering and Resequencing steps")
    add('--disable_consensus', 
        action='store_false',
        dest='enable_consensus',
        help="Turn off the Consensus step")
    add('--blasr', 
        metavar='BLASR_PATH', 
        help="Specify the path to the Blasr executable")
    add('--mothur', 
        metavar='MOTHUR_PATH', 
        default='mothur',
        help="Specify the path to the Mothur executable")
    add('--debug', 
        action='store_true',
        help="Turn on DEBUG message logging")

    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "\tHLA Analysis Pipeline version: %s" % __VERSION__
            raise SystemExit

    add("--version",
        nargs=0,
        action=PrintVersionAction)

    parser.parse_args( namespace=args )

    # Validate numerical parameters
    validate_int( 'NumProc', args.nproc, minimum=0 )
    validate_float( 'Distance', args.distance, minimum=MIN_DIST, 
                                               maximum=MAX_DIST )
