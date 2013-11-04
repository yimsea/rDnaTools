#! /usr/bin/env python

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import os
import sys
import logging

from pbcore.io.BasH5IO import BasH5Collection
from pbcore.io.FastqIO import FastqRecord, FastqWriter 

from pbrdna.arguments import args, MIN_LENGTH, MIN_SNR
from pbrdna.utils import get_output_name

log = logging.getLogger(__name__)

LENGTH = args.min_length if hasattr(args, 'min_length') else MIN_LENGTH
SNR = args.min_snr if hasattr(args, 'min_snr') else MIN_SNR

def extract_ccs( input_file, output_file=None,
                             min_length=LENGTH,
                             min_snr=SNR):
    """
    Extract CCS reads from an input_file
    """
    output_file = output_file or get_output_name( input_file, 'fastq' )
    collection = BasH5Collection( input_file )
    extract_ccs_fastq( collection, output_file, min_length, min_snr )
    return output_file

def extract_ccs_fastq( collection, output_file, min_length, min_snr ):
    log.info('Extracting fastq CCS reads from input files')
    log.debug('    min_length: %s' % min_length)
    log.debug('    min_snr: %s' % min_snr)
    ccs_total = 0
    pass_total = 0
    with FastqWriter( output_file ) as writer:
        for movie in collection.movieNames:
            log.info('Extracting fastq CCS reads from %s' % os.path.basename(movie))
            ccs_count = 0
            pass_count = 0
            for well in collection[movie].sequencingZmws:
                zmw = collection[movie][well]

                # Skip non-CCS ZMWs
                if not zmw.ccsRead:
                    continue
                ccs_count += 1

                # Skip short and low-SNR sequences
                basecalls = zmw.ccsRead.basecalls()
                if len(basecalls) < min_length:
                    continue
                zmw_snr = min( zmw.zmwMetric("HQRegionSNR") )
                if zmw_snr < min_snr:
                    continue
                pass_count += 1

                # Finally write the CCS Fastq to file
                record = FastqRecord(zmw.ccsRead.readName,
                                     basecalls,
                                     zmw.ccsRead.QualityValue())
                writer.writeRecord( record )
            percentage = round(100.0*pass_count/ccs_count)
            log.info("Identified {0} CCS reads, of which {1} ({2}%) passed filter".format(ccs_count,
                                                                                           pass_count,
                                                                                           percentage))
            ccs_total += ccs_count
            pass_total += pass_count
    percentage = round(100.0*pass_total/ccs_total)
    log.info('Found a total of {0} CCS reads, of which {1} ({2}%) passed filter'.format(ccs_total,
                                                                                         pass_total,
                                                                                         percentage))

if __name__ == '__main__':
    import sys

    input_file = sys.argv[1]
    min_snr = float(sys.argv[2])

    logging.basicConfig( stream=sys.stdout, level=logging.DEBUG )
    extract_ccs( input_file, min_snr=min_snr )