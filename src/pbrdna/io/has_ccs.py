#! /usr/bin/env python

__author__ = 'bbowman@pacificbiosciences.com'

import logging
from pbcore.io.BasH5IO import BasH5Collection
from pbrdna.utils import is_fasta, is_fastq, is_bash5, is_fofn

log = logging.getLogger(__name__)

def file_has_ccs( filename ):
    if is_fasta( filename ) or is_fastq( filename ):
        return True
    elif is_bash5( filename ) or is_fofn( filename ):
        return h5_has_ccs( filename )
    else:
        msg = "Unrecognized file-type"
        log.error( msg )
        raise TypeError( msg )

def h5_has_ccs( filename ):
    """
    Check if a given input file has CCS Data
    """
    collection = BasH5Collection( filename )
    for movie in collection.movieNames:
        log.info('Testing movie "%s" for the presence of CCS data')
        # Make sure the movie has valid sequencing data
        if len(collection[movie].sequencingZmws) < 1:
            msg = 'Movie "%s" has no valid sequencing ZMWs' % movie
            log.error( msg )
            raise ValueError( msg )
        # Test a random ZMW to see if the movie has CCS data
        test_zmw = collection[movie].sequencingZmws[0]
        try:
            collection[movie][test_zmw].ccsRead
        except ValueError:
            log.info('Movie "%s" has no CCS data' % movie)
            return False
        log.info('Movie "%s" has valid CCS data' % movie)
    log.info('All supplied sequencing movies have valid CCS data\n')
    return True