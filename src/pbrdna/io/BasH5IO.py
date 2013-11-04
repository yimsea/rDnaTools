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

from pbcore.io.BasH5Reader import BasH5Reader
from pbcore.io.FastqIO import FastqRecord, FastqWriter

log = logging.getLogger(__name__)

class BasH5Extractor(object):
    """
    A tool for extracting sequence data from a BasH5 file    
    """

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, input_file=None, output_file=None):
        if input_file is None:
            self.initialize_from_args()
        else:
            self.initialize_from_call(input_file, output_file)
        self.validate_settings()
        self.initialize_readers()

    def initialize_from_args(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_file', metavar='FILE',
                            help="BasH5 or FOFN to extract from")
        parser.add_argument('-o', '--output', default=sys.stdout,
                            help="Specify a file to output the data to")
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--subreads',  action='store_true',
                            help="Output sequences from individual subreads")
        group.add_argument('--CCS', action='store_true',
                            help="Output sequences from CCS reads")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initialize_from_call(self, input_file, output_file):
        self.input_file = input_file
        if output_file is None:
            self.output = sys.stdout
        else:
            self.output = output_file

    def validate_settings(self):
        if self.input_file.endswith('.bas.h5') or \
           self.input_file.endswith('.fofn'):
            log.info('Creating a BasH5Extractor for "{0}"'.format(self.input_file))
            log.info('Outputing extracted reads to "{0}"'.format(self.output))
        else: 
            raise ValueError('Input files must be either FOFN or BasH5!')

    def initialize_readers(self):
        self.bash5_readers = []
        if self.input_file.endswith('.bas.h5'):
            self.add_bash5_reader( self.input_file )
        elif self.input_file.endswith('.fofn'):
            self.parse_fofn_file( self.input_file )

    def add_bash5_reader(self, bash5_file):
        filepath = os.path.realpath( bash5_file )
        path, filename = os.path.split( filepath )
        log.info('Creating a BasH5Reader for "{0}"'.format(filename))
        self.bash5_readers.append( BasH5Reader( filepath ) )

    def parse_fofn_file(self, fofn_file):
        with open(fofn_file, 'r') as handle:
            for line in handle:
                fofn_entry = line.strip()
                if not fofn_entry:
                    continue
                if fofn_entry.endswith('.bas.h5'):
                    self.add_bash5_reader( fofn_entry )
                elif fofn_entry.endswith('.bax.h5'):
                    self.add_bash5_reader( fofn_entry )
                else:
                    raise ValueError('FOFN must contain only BasH5 and BaxH5 files!')

    #################
    # Class Methods #
    #################

    @classmethod
    def writeCcsFastq(cls, basH5Reader, fastqWriter):
        log.info('Writing Fastq CCS reads from "%s"...' % basH5Reader.movieName)
        for zmw in basH5Reader:
            if zmw.ccsRead:
                fastqRecord = FastqRecord(zmw.ccsRead.readName,
                                          zmw.ccsRead.basecalls(),
                                          zmw.ccsRead.QualityValue())
                fastqWriter.writeRecord( fastqRecord )

    @classmethod
    def writeSubreadFastq(cls, basH5Reader, fastqWriter):
        log.info('Writing Fastq subreads from "%s"...' % basH5Reader.movieName)
        for zmw in basH5Reader:
            for subread in zmw.subreads():
                fastqRecord = FastqRecord(subread.readName,
                                          subread.basecalls(),
                                          subread.QualityValue())
                fastqWriter.writeRecord( fastqRecord )

    ####################
    # Instance Methods #
    ####################

    def outputCcsFastq(self):
        log.info('Parsing Fastq CCS reads from input BAS.H5 files')
        with FastqWriter(self.output) as writer:
            for reader in self.bash5_readers:
                self.writeCcsFastq( reader, writer )

    def outputSubreadFastq(self):
        log.info('Parsing Fastq subreads from input BAS.H5 files')
        with FastqWriter(self.output) as writer:
            for reader in self.bash5_readers:
                self.writeSubreadFastq( reader, writer )

    def __call__(self):
        if self.CCS:
            self.outputCcsFastq()
        elif self.subreads:
            self.outputSubreadFastq()

if __name__ == '__main__':
    extractor = BasH5Extractor()
    extractor()
