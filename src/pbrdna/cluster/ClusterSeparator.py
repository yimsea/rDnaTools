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
import subprocess

from pbcore.io.FastqIO import FastqRecord, FastqReader, FastqWriter
from pbcore.io.FastaIO import FastaRecord, FastaWriter  
from pbrdna.fastq.utils import meanPQv
from pbrdna.utils import get_zmw, create_directory, validate_input, validate_float

DEFAULT_DIST = 0.03
MIN_FULL_LENGTH = 1400
MIN_CLUSTER_SIZE = 3

class ClusterSeparator(object):
    """
    A tool for resequencing clusters of rDNA sequences    
    """

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, listFile=None, ccsFile=None, output=None, 
                                                    distance=None, 
                                                    min_cluster_size=None):
        if listFile is None or ccsFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(listFile, ccsFile, output, distance, min_cluster_size)
        self.validate_settings()

    def initializeFromArgs(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('listFile', metavar='FILE',
                            help="Mothur list file of cluster data")
        parser.add_argument('ccsFile', metavar='FILE',
                            help="Fasta or Fastq file of CCS sequences")
        parser.add_argument('-c', '--min_cluster_size', type=int,
                            default=MIN_CLUSTER_SIZE, metavar='INT',
                            help="Minimum number of sequences to require in a cluster")
        parser.add_argument('-d', '--distance', type=float,
                            default=DEFAULT_DIST, metavar='FLOAT',
                            help="Distance at which to cluster sequences")
        parser.add_argument('-m', '--minRefLength', type=int,
                            default=MIN_FULL_LENGTH, metavar='INT',
                            help="Minimum length to allow for reference sequences")
        parser.add_argument('--outputFastq', action='store_true',
                            help="Output FASTQ files for each cluster")
        parser.add_argument('--outputReference', action='store_true',
                            help="Output reference sequences for each cluster")
        parser.add_argument('-o', '--outputDir', default='reseq',
                            help="Specify a directory for output files")
        parser.add_argument('-s', '--summary', dest='output', default=None,
                            help="Output a summary of the individual clusters")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall(self, listFile, ccsFile, output, distance=None, 
                                                            min_cluster_size=None):
        self.listFile = listFile
        self.ccsFile = ccsFile
        self.output = output
        self.distance = DEFAULT_DIST if distance is None else distance
        self.min_cluster_size = MIN_CLUSTER_SIZE if distance is None else min_cluster_size
        self.outputDir = 'reseq'
        self.outputFastq = False
        self.outputReference = True
        self.minRefLength = MIN_FULL_LENGTH

    def validate_settings( self ):
        # Check the values of the supplied input files
        self.listFile = validate_input( self.listFile, ['.list'])
        self.ccsFile = validate_input( self.ccsFile, ['.fq', '.fastq'])
        # Check the value of the supplied distance
        validate_float( 'Distance', self.distance, minimum=0.001, 
                                                   maximum=0.5)

    def initializeOutputFolder( self ):
        # Create the output directory if needed and move into it
        self.origin=   os.getcwd()
        create_directory( self.outputDir )
        os.chdir( self.outputDir )

    #################
    # Class Methods #
    #################

    @classmethod
    def convertDistance(cls, distance):
        try:
            distance = 'unique' if distance == 'unique' else float(distance)
        except:
            raise ValueError('"%s" is not a valid distance!' % parts[0])
        return distance

    ####################
    # Instance Methods #
    ####################

    def parseSequenceData(self):
        self.sequenceData = {}
        self.qualityData = {}
        for record in FastqReader( self.ccsFile ):
            zmw = get_zmw( record.name )
            new_record = FastqRecord(zmw, record.sequence, record.quality)
            self.sequenceData[zmw] = new_record
            self.qualityData[zmw] = meanPQv(new_record)

    def parseDistances(self):
        distances = []
        with open( self.listFile, 'r' ) as handle:
            for line in handle:
                parts = line.split()
                distance = self.convertDistance( parts[0] )
                distances.append( distance )
        return distances

    def selectDistance(self, distances):
        # If our selected distance is present, simply return it
        print self.distance
        print distances
        if self.distance in distances:
            return self.distance
        # Otherwise find the largest clustering distance smaller than 
        #    the specified distance and return that
        possible = [d for d in distances if d != 'unique']
        smaller = [d for d in possible if d < self.distance]
        if not smaller:
            raise ValueError('No valid clustering distances found!')
        return max(smaller)

    def parseClusters( self, targetDist ):
        with open( self.listFile, 'r' ) as handle:
            for line in handle:
                # Skip lines until we find the target distance
                parts = line.split()
                currDist = self.convertDistance( parts[0] )
                if currDist != targetDist:
                    continue
                # Check that the number of clusters is concordant
                clusterCount = int(parts[1])
                clusters = parts[2:]
                assert len(clusters) == clusterCount
                # Convert the strings of clusters to Lists and return
                clusters = [c.split(',') for c in clusters]
                return clusters

    def trimClusterNames(self, clusters):
        trimmed = []
        for cluster in clusters:
            cluster = [get_zmw(c) for c in cluster]
            trimmed.append( frozenset(cluster) )
        return trimmed

    def getClusterReads(self, cluster):
        reads = []
        for ccsZmw in cluster:
            try:
                ccsRead = self.sequenceData[ccsZmw]
            except KeyError:
                #raise Warning("No CCS read found for '%s', skipping..." % ccsZmw)
                continue
            reads.append( ccsRead )
        return reads

    def outputClusterFastq(self, reads, count):
        fastqFile = 'cluster%s.fastq' % count
        if os.path.exists( fastqFile ):
            return fastqFile
        # Rename the "Reference" sequence to the cluster
        with FastqWriter( fastqFile ) as handle:
            for fastqRecord in reads:
                handle.writeRecord( fastqRecord )
        return fastqFile

    def outputClusterFasta( self, reads, count ):
        fastaFile = 'cluster%s.fasta' % count
        if os.path.exists( fastaFile ):
            return fastaFile
        # Rename the "Reference" sequence to the cluster
        with FastaWriter( fastaFile ) as handle:
            for fastqRecord in reads:
                fastaRecord = FastaRecord( fastqRecord.name,
                                           fastqRecord.sequence )
                handle.writeRecord( fastaRecord )
        return fastaFile

    def pickReference( self, reads ):
        longReads = [read for read in reads
                          if len(read.sequence) > self.minRefLength]
        if longReads:
            return self.findLowestErrorRead( reads )
        # If no 'full-length' reads are present, simply return the longest
        else:
            return self.findLongestRead( reads )

    def findLowestErrorRead( self, reads ):
        pQvs = [self.qualityData[read.name] for read in reads]
        maxQv = max(pQvs)
        lowestErrorReads = [read for read in reads
                                 if self.qualityData[read.name] == maxQv]
        return lowestErrorReads[0]

    def findLongestRead( self, reads ):
        lengths = [len(read.sequence) for read in reads]
        maxLength = max(lengths)
        longestReads = [read for read in reads
                             if len(read.sequence) == maxLength]
        return longestReads[0]

    def outputReferenceFasta( self, reference, count):
        print "Creating reference sequence for Cluster #%s" % count
        referenceFile = 'cluster%s_ref.fasta' % count
        reference_desc = 'cluster{0}_reference\t{1}'.format(count, reference.name)
        if os.path.exists( referenceFile ):
            return referenceFile
        with FastaWriter( referenceFile ) as handle:
            referenceFasta = FastaRecord( reference_desc,
                                          reference.sequence )
            handle.writeRecord( referenceFasta )
        return referenceFile

    def outputClusterFileList( self, clusterFiles ):
        print "Writing out the names of the individual cluster files"
        with open( self.output, 'w') as handle:
            for clusterFile, referenceFile, count in clusterFiles:
                clusterPath = os.path.join( self.outputDir, clusterFile )
                referencePath = os.path.join( self.outputDir, referenceFile )
                handle.write('{0}\t{1}\t{2}\n'.format(clusterPath, 
                                                      referencePath, 
                                                      count))
        return self.output

    def __call__( self ):
        self.initializeOutputFolder()
        self.parseSequenceData()
        # Select the appropriate distance, and parse the matching clusters
        distances = self.parseDistances()
        distance = self.selectDistance( distances )
        clusters = self.parseClusters( distance )
        # Trim the cluster neams and iterate, outputing each subset
        trimmedClusters = self.trimClusterNames( clusters )
        clusterFiles = []
        for count, cluster in enumerate( trimmedClusters ):
            count = str(count+1).zfill(4)
            print "Analyzing cluster #%s now..." % (count)
            reads = self.getClusterReads( cluster )
            clusterFile = self.outputClusterFasta( reads, count )
            if self.outputFastq:
                self.outputClusterFastq( reads, count )
            if len(reads) >= self.min_cluster_size and self.outputReference:
                reference = self.pickReference( reads )
                referenceFile = self.outputReferenceFasta( reference, count )
                clusterFiles.append( (clusterFile, referenceFile, len(reads)) )
            else:
                clusterFiles.append( (clusterFile, 'None', len(reads)) )
        # Return to the origin directory and output the results summary
        os.chdir( self.origin )
        if self.output:
            return self.outputClusterFileList( clusterFiles )
        return

if __name__ == '__main__':
    separator = ClusterSeparator()
    separator()
