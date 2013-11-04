#! /usr/bin/env python

import os
import logging

from pbcore.io.FastaIO import FastaWriter

from pbrdna.log import initialize_logger
from pbrdna.arguments import args, parse_args
from pbrdna.io.has_ccs import file_has_ccs
from pbrdna.io.extract_ccs import extract_ccs
from pbrdna.io.MothurIO import SummaryReader
from pbrdna.fasta.utils import fasta_count, copy_fasta_sequences
from pbrdna.fastq.quality_filter import quality_filter
from pbrdna.fastq.QualityAligner import QualityAligner
from pbrdna.fastq.QualityMasker import QualityMasker
from pbrdna.mothur.MothurTools import MothurRunner
from pbrdna.cluster.ClusterSeparator import ClusterSeparator
from pbrdna.resequence.DagConTools import DagConRunner
from pbrdna.utils import (validate_executable,
                          create_directory,
                          split_root_from_ext,
                          get_output_name,
                          file_exists,
                          all_files_exist)

log = logging.getLogger(__name__)

class rDnaPipeline( object ):
    """
    A tool for running a community analysis pipeline on PacBioData
    """

    def __init__(self):
        parse_args()
        self.__dict__.update( vars(args) )
        self.validate_settings()
        self.initialize_output()
        if self.debug:
            initialize_logger( self.log_file, logging.DEBUG )
        else:
            initialize_logger( self.log_file, logging.INFO )

    def validate_settings(self):
        # Validate the input file
        root, ext = split_root_from_ext( self.input_file )
        if ext in ['.bas.h5', '.fofn']:
            self.data_type = 'bash5'
        elif ext in ['.fq', '.fastq']:
            self.data_type = 'fastq'
        elif ext in ['.fa', '.fsa', '.fasta']:
            self.data_type = 'fasta'
            self.enable_masking = False
            self.enable_consensus = False
        else:
            raise TypeError('Sequence file must be a bas.h5 file, a ' + \
                            'fasta file, or a fofn of multiple such files')
        # If Clustering was disabled, also disable the consensus process
        if not self.enable_clustering:
            self.enable_consensus = False
        # If Consensus is enabled, initialize the appropriate tool
        if self.enable_consensus:
            self.consensusTool = DagConRunner('gcon.py', 'r')
        # Searching for Mothur executable, and set the Mothur Process counter
        self.mothur = validate_executable( self.mothur )
        self.processCount = 0

    def initialize_output(self):
        # Create the Output directory
        create_directory( self.output_dir )
        # Create a symbolic link from the data file to the output dir
        baseName = os.path.basename( self.input_file )
        symlinkPath = os.path.join( self.output_dir, baseName )
        if os.path.exists( symlinkPath ):
            pass
        else:
            absPath = os.path.abspath( self.input_file )
            os.symlink( absPath, symlinkPath )
        self.sequenceFile = baseName
        # Move into the Output directory and create Log directory and files
        os.chdir( self.output_dir )
        create_directory( 'log' )
        stdoutLog = os.path.join('log', 'mothur_stdout.log')
        stderrLog = os.path.join('log', 'mothur_stderr.log')
        self.log_file = os.path.join('log', 'rna_pipeline.log')
        # Instantiate the MothurRunner object
        self.factory = MothurRunner( self.mothur, 
                                     self.nproc, 
                                     stdoutLog, 
                                     stderrLog)

    def getProcessLogFile(self, process, isMothurProcess=False):
        if isMothurProcess:
            logFile = 'process%02d.mothur.%s.logfile' % (self.processCount, 
                                                         process)
        else:
            logFile = 'process%02d.%s.logfile' % (self.processCount, process)
        return os.path.join('log', logFile)

    def process_setup(self, inputFile, processName, suffix=None, suffixList=None):
        """ 
        Return a tuple containing the output file and a boolean flag describing
        whether the output file already exists
        """
        log.info('Preparing to run %s on "%s"' % (processName, inputFile))
        self.processCount += 1
        if suffix:
            outputFile = get_output_name(inputFile, suffix)
            return outputFile
        elif suffixList:
            outputFiles = []
            for suffix in suffixList:
                outputFile = get_output_name( inputFile, suffix )
                outputFiles.append( outputFile )
            return outputFiles

    def output_files_exist( self, outputFile=None, outputList=None ):
        if outputFile:
            if file_exists( outputFile ):
                log.info('Output files detected, skipping process...\n')
                return True
            else:
                log.info('Output files not found, running process...')
                return False
        elif outputList:
            if all_files_exist( outputList ):
                log.info('Output files detected, skipping process...\n')
                return True
            else:
                log.info('Output files not found, running process...')
                return False

    def checkOutputFile( self, outputFile ):
        if file_exists( outputFile ):
            log.info('Expected output "%s" found' % outputFile)
        else:
            msg = 'Expected output "%s" not found!' % outputFile
            log.info( msg )
            raise IOError( msg )

    def process_cleanup(self, outputFile=None, outputList=None):
        """
        Log if the process successfully created it's output, and raise an
        error message if not
        """
        if outputFile:
            self.checkOutputFile( outputFile )
        elif outputList:
            for outputFile in outputList:
                self.checkOutputFile( outputFile )
        log.info('All expected output files found - process successful!\n')

    def write_dummy_file(self, dummyFile):
        with open(dummyFile, 'w') as handle:
            handle.write('DONE')
        return dummyFile

    def extract_raw_ccs(self, inputFile):
        outputFile = self.process_setup( inputFile, 
                                         'extractCcsFromBasH5',
                                         suffix='fastq' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        elif file_has_ccs( inputFile ):
            extract_ccs(inputFile, outputFile, self.raw_data)
        else:
            msg = 'Raw data file has no CCS data!'
            log.error( msg )
            raise ValueError( msg )
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def filter_fastq(self, fastqFile):
        outputFile = self.process_setup( fastqFile, 
                                         'FilterQuality',
                                         suffix='filter.fastq' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        quality_filter( fastqFile, outputFile )
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def separate_fastq(self, fastqFile):
        outputList = self.process_setup( fastqFile, 
                                        'Fastq.Info', 
                                        suffixList=['fasta', 'qual'] )
        if self.output_files_exist( outputList=outputList ):
            return outputList
        mothurArgs = {'fastq':fastqFile, 'fasta':'T', 'qfile':'T'}
        logFile = self.getProcessLogFile('fastq.info', True)
        self.factory.runJob('fastq.info', mothurArgs, logFile)
        self.process_cleanup( outputList=outputList )
        return outputList

    def align_sequences(self, fastaFile):
        outputFile = self.process_setup( fastaFile, 
                                        'Align.Seqs', 
                                        suffix='align' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':fastaFile,
                      'reference':self.alignment_reference,
                      'flip':'t'}
        logFile = self.getProcessLogFile('align.seqs', True)
        self.factory.runJob('align.seqs', mothurArgs, logFile)
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def screen_sequences(self, alignFile, start=None, end=None, min_length=None):
        if alignFile.endswith('.align'):
            outputExt = 'good.align'
        elif alignFile.endswith('.fasta'):
            outputExt = 'good.fasta'
        outputFile = self.process_setup( alignFile, 
                                         'Screen.Seqs', 
                                         suffix=outputExt )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'start':start,
                      'end':end,
                      'minlength':min_length}
        logFile = self.getProcessLogFile('screen.seqs', True)
        self.factory.runJob('screen.seqs', mothurArgs, logFile)
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def summarize_sequences(self, fastaFile):
        outputFile = self.process_setup( fastaFile, 
                                        'Summary.Seqs', 
                                        suffix='summary' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':fastaFile}
        logFile = self.getProcessLogFile('summary.seqs', True)
        self.factory.runJob('summary.seqs', mothurArgs, logFile)
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def parse_summary_file(self, summaryFile):
        log.info('Preparing to run SummaryReader...')
        parser = SummaryReader(summaryFile, self.fraction)
        log.info('Identifying full-length alignment positions...')
        start, end = parser.getFullLengthPositions()
        log.info('Full-length start is NAST Alignment position %s' % start)
        log.info('Full-length end is NAST Alignment position %s' % end)
        log.info('Calculating minimum allowed alignment positions...')
        maxStart, minEnd = parser.getAllowedPositions()
        log.info('Maximum allowed start is NAST Alignment position %s' % maxStart)
        log.info('Minimum allowed end is NAST Alignment position %s\n' % minEnd)
        return maxStart, minEnd

    def find_chimeras(self, alignFile):
        outputFile = self.process_setup( alignFile, 
                                        'UCHIME', 
                                        suffix='uchime.accnos' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'reference':self.chimera_reference}
        logFile = self.getProcessLogFile('chimera.uchime', True)
        self.factory.runJob('chimera.uchime', mothurArgs, logFile)
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def remove_sequences(self, alignFile, idFile):
        outputFile = self.process_setup( alignFile, 
                                        'Remove.Seqs', 
                                        suffix='pick.align' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'accnos':idFile}
        logFile = self.getProcessLogFile('remove.seqs', True)
        self.factory.runJob('remove.seqs', mothurArgs, logFile)
        self.process_cleanup( outputFile=outputFile )
        return outputFile
  
    def filter_sequences(self, alignFile):
        outputFile = self.process_setup( alignFile, 
                                        'Filter.Seqs', 
                                        suffix='filter.fasta' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'vertical':'T'}
        logFile = self.getProcessLogFile( 'filter.seqs', True )
        self.factory.runJob( 'filter.seqs', mothurArgs, logFile )
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def add_quality_to_alignment(self, fastqFile, alignFile):
        outputFile = self.process_setup( alignFile, 
                                        'QualityAligner', 
                                        suffix='fastq' )
        if self.output_files_exist( outputFile=output ):
            return output
        aligner = QualityAligner( fastqFile, alignFile, outputFile )
        aligner.run()
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def mask_fastq_sequences(self, fastqFile):
        outputFile = self.process_setup( fastqFile, 
                                        'QualityMasker', 
                                        suffix='masked.fastq' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        masker = QualityMasker(fastqFile, outputFile, self.minQv)
        masker.run()
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def unique_sequences( self, alignFile ):
        if alignFile.endswith('.align'):
            outputSuffixes = ['unique.align', 'names']
        elif alignFile.endswith('.fasta'):
            outputSuffixes = ['unique.fasta', 'names']
        outputList = self.process_setup( alignFile,
                                        'Unique.Seqs',
                                        suffixList=outputSuffixes )
        if self.output_files_exist( outputList=outputList ):
            return outputList
        mothurArgs = {'fasta':alignFile}
        logFile = self.getProcessLogFile('unique.seqs', True)
        self.factory.runJob('unique.seqs', mothurArgs, logFile)
        self.process_cleanup( outputList=outputList )
        return outputList

    def precluster_sequences( self, alignFile, nameFile ):
        if alignFile.endswith('.align'):
            outputSuffixes = ['precluster.align', 'precluster.names']
        elif alignFile.endswith('.fasta'):
            outputSuffixes = ['precluster.fasta', 'precluster.names']
        outputList = self.process_setup( alignFile,
                                        'Pre.Cluster',
                                        suffixList=outputSuffixes )
        if self.output_files_exist( outputList=outputList ):
            return outputList
        mothurArgs = { 'fasta':alignFile,
                       'name': nameFile,
                       'diffs':self.precluster_diffs }
        logFile = self.getProcessLogFile('pre.cluster', True)
        self.factory.runJob('pre.cluster', mothurArgs, logFile)
        self.process_cleanup( outputList=outputList )
        return outputList

    def calculate_distance_matrix( self, alignFile ):
        outputFile = self.process_setup( alignFile, 
                                        'Dist.Seqs', 
                                        suffix='phylip.dist' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = { 'fasta':alignFile,
                       'calc':'onegap',
                       'countends':'F',
                       'output':'lt' }
        logFile = self.getProcessLogFile('dist.seqs', True)
        self.factory.runJob('dist.seqs', mothurArgs, logFile)
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def cluster_sequences(self, distanceMatrix):
        if self.clusteringMethod == 'nearest':
            outputSuffix = 'nn.list'
        elif self.clusteringMethod == 'average':
            outputSuffix = 'an.list'
        elif self.clusteringMethod == 'furthest':
            outputSuffix = 'fn.list'
        outputFile = self.process_setup( distanceMatrix, 
                                        'Cluster', 
                                        suffix=outputSuffix )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'phylip':distanceMatrix,
                      'method':self.clusteringMethod}
        logFile = self.getProcessLogFile( 'cluster', True )
        self.factory.runJob( 'cluster', mothurArgs, logFile )
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def separate_cluster_sequences(self, listFile, sequenceFile):
        outputFile = self.process_setup( listFile, 
                                        'ClusterSeparator', 
                                        suffix='list.clusters')
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        separator = ClusterSeparator( listFile, 
                                      sequenceFile,
                                      outputFile,
                                      self.distance, 
                                      self.min_cluster_size )
        separator()
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def generate_consensus_sequences(self, clusterListFile):
        outputFile = self.process_setup( clusterListFile, 
                                        'ClusterResequencer', 
                                        suffix='consensus')
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        consensusFiles = []
        with open( clusterListFile ) as handle:
            for line in handle:
                sequenceFile, referenceFile, count = line.strip().split()
                if referenceFile.endswith('None'):
                    consensusFiles.append( (sequenceFile, 'None') )
                else:
                    root_name = os.path.basename( sequenceFile )
                    consensus = self.consensusTool( sequenceFile, 
                                                    referenceFile )
                    consensusFiles.append( (referenceFile, consensus) )
        with open( outputFile, 'w' ) as handle:
            for filenamePair in consensusFiles:
                handle.write('%s\t%s\n' % filenamePair)
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def cleanup_consensus_folder( self, consensusFile ):
        outputFile = self.process_setup( consensusFile, 
                                        'ConsensusCleanup', 
                                        suffix='consensus.cleanup' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        reseqPath = os.path.join( os.getcwd(), 'reseq' )
        for filename in os.listdir( reseqPath ):
            filePath = os.path.join( reseqPath, filename )
            if filePath.endswith('_input.fa'):
                os.remove( filePath )
            elif filePath.endswith('_input.fa.aln'):
                os.remove( filePath )
            elif filePath.endswith('_input.fa.aln_unsorted'):
                os.remove( filePath )
        self.write_dummy_file( outputFile )
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def select_final_sequences( self, consensusFile ):
        outputFile = self.process_setup( consensusFile, 
                                        'SequenceSelector', 
                                        suffix='consensus.selected' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        selectedFiles = []
        with open( consensusFile ) as handle:
            for line in handle:
                referenceFile, consensusFile = line.strip().split()
                if consensusFile.endswith('None'):
                    pass
                elif fasta_count( consensusFile ) == 1:
                    selectedFiles.append( consensusFile )
                else:
                    selectedFiles.append( referenceFile )
        with open( outputFile, 'w' ) as handle:
            for filename in selectedFiles:
                handle.write(filename + '\n')
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def output_final_sequences( self, finalSequenceList ):
        outputFile = self.process_setup( finalSequenceList, 
                                        'SequenceWriter',
                                        suffix='fasta' )
        if self.output_files_exist( outputFile=outputFile ):
            return outputFile
        with FastaWriter( outputFile ) as writer:
            with open( finalSequenceList ) as handle:
                for line in handle:
                    sequenceFile = line.strip()
                    copy_fasta_sequences( sequenceFile, writer )
        self.process_cleanup( outputFile=outputFile )
        return outputFile

    def run(self):
        if self.data_type == 'bash5':
            fastqFile = self.extract_raw_ccs( self.sequenceFile )
        elif self.data_type == 'fastq':
            fastqFile = self.sequenceFile
        elif self.data_type == 'fasta':
            fastqFile = None
            fastaFile = self.sequenceFile
        # If we have a Fastq, filter low-quality reads and convert to FASTA
        if fastqFile:
            filteredFastq = self.filter_fastq( fastqFile )
            fastaFile, qualFile = self.separate_fastq( filteredFastq )
        # Align the Fasta sequences and remove partial reads
        alignedFile = self.align_sequences( fastaFile )
        summaryFile = self.summarize_sequences( alignedFile )
        maxStart, minEnd = self.parse_summary_file( summaryFile )
        screenedFile = self.screen_sequences(alignedFile, 
                                            start=maxStart,
                                            end=minEnd)
        # Identify and remove chimeric reads
        chimeraIds = self.find_chimeras( screenedFile )
        noChimeraFile = self.remove_sequences( screenedFile, chimeraIds )
        # Filter out un-used columns to speed up re-alignment and clustering
        filteredFile = self.filter_sequences( noChimeraFile )
        # If masking is enabled, create an aligned FASTQ, mask the 
        # low-quality bases and remove over-masked reads
        if self.enable_masking:
            alignedFastqFile = self.add_quality_to_alignment( fastqFile, filteredFile )
            maskedFastq = self.mask_fastq_sequences( alignedFastqFile )
            maskedFasta = self.convert_fastq_to_fasta( maskedFastq )
            screenedFasta = self.screen_sequences( maskedFasta,
                                                  min_length=self.min_length)
            fileForClustering = screenedFasta
        # Otherwise if masking is disabled, we'll use unique-ify and 
        #    pre-cluster our sequences
        else:
            uniqueFile, nameFile = self.unique_sequences( filteredFile )
            preclusteredFile, nameFile = self.precluster_sequences( uniqueFile, nameFile )
            fileForClustering = preclusteredFile
        # If enabled, calculate sequence distances and cluster
        if self.enable_clustering:
            distanceMatrix = self.calculate_distance_matrix( fileForClustering )
            listFile = self.cluster_sequences( distanceMatrix )
        # If enabled, generate a consensus for each cluster from above
        if self.enable_consensus:
            clusterListFile = self.separate_cluster_sequences( listFile, fastqFile )
            consensusFile = self.generate_consensus_sequences( clusterListFile )
            self.cleanup_consensus_folder( consensusFile )
            selectedFile = self.select_final_sequences( consensusFile )
            finalFile = self.output_final_sequences( selectedFile )

if __name__ == '__main__':
    rDnaPipeline().run()
