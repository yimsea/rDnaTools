#! /usr/bin/env python

import os, sys, logging

from pbrdna import __VERSION__
from pbrdna.arguments import args, parse_args
from pbrdna.io.FastaIO import FastaWriter
from pbrdna.io.BasH5IO import BasH5Extractor
from pbrdna.io.MothurIO import SummaryReader
from pbrdna.fasta.utils import fasta_count, copy_fasta_sequences
from pbrdna.fastq.QualityFilter import QualityFilter
from pbrdna.fastq.QualityAligner import QualityAligner
from pbrdna.fastq.QualityMasker import QualityMasker
from pbrdna.mothur.MothurTools import MothurRunner
from pbrdna.cluster.ClusterSeparator import ClusterSeparator
from pbrdna.resequence.DagConTools import DagConRunner
from pbrdna._utils import *

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
        self.initialize_logger()

    def validate_settings(self):
        # Validate the input file
        root, ext = splitRootFromExt( self.sequenceFile )
        if ext in ['.bas.h5', '.fofn']:
            self.dataType = 'bash5'
        elif ext in ['.fq', '.fastq']:
            self.dataType = 'fastq'
        elif ext in ['.fa', '.fsa', '.fasta']:
            self.dataType = 'fasta'
            self.enableMasking = False
            self.enableConsensus = False
        else:
            raise TypeError('Sequence file must be a bas.h5 file, a ' + \
                            'fasta file, or a fofn of multiple such files')
        # If Clustering was disabled, also disable the consensus process
        if not self.enableClustering:
            self.enableConsensus = False
        # If Consensus is enabled, initialize the appropriate tool
        if self.enableConsensus:
            self.consensusTool = DagConRunner('gcon.py', 'r')
        # Searching for Mothur executable, and set the Mothur Process counter
        self.mothur = validateExecutable( self.mothur )
        self.processCount = 0

    def initialize_output(self):
        # Create the Output directory
        createDirectory( self.output )
        # Create a symbolic link from the data file to the output dir
        baseName = os.path.basename( self.sequenceFile )
        symlinkPath = os.path.join( self.output, baseName )
        if os.path.exists( symlinkPath ):
            pass
        else:
            absPath = os.path.abspath( self.sequenceFile )
            os.symlink( absPath, symlinkPath )
        self.sequenceFile = baseName
        # Move into the Output directory and create Log directory and files
        os.chdir( self.output )
        createDirectory( 'log' )
        stdoutLog = os.path.join('log', 'mothur_stdout.log')
        stderrLog = os.path.join('log', 'mothur_stderr.log')
        self.log_file = os.path.join('log', 'rna_pipeline.log')
        # Instantiate the MothurRunner object
        self.factory = MothurRunner( self.mothur, 
                                     self.nproc, 
                                     stdoutLog, 
                                     stderrLog)

    def initialize_logger(self):
        if self.debug:
            log_level = logging.DEBUG
        else:
            log_level = logging.INFO 
        # Initialize the LogHandler for the master log file
        logging.basicConfig( level=log_level, 
                             stream=sys.stdout )
        logging.basicConfig( level=log_level, 
                             filename=self.log_file )
        # Record the initialization of the pipeline
        log.info("INFO logger initialized")
        log.debug("DEBUG logger initialized")
        log.info("Initializing rDnaPipeline v%s" % __VERSION__)
        log.debug("Using the following parameters:")
        for param, value in self.__dict__.iteritems():
            log.debug("\t%s = %s" % (param, value))
        log.info("Initialization of rDnaPipeline completed\n")

    def getProcessLogFile(self, process, isMothurProcess=False):
        if isMothurProcess:
            logFile = 'process%02d.mothur.%s.logfile' % (self.processCount, 
                                                         process)
        else:
            logFile = 'process%02d.%s.logfile' % (self.processCount, process)
        return os.path.join('log', logFile)

    def processSetup(self, inputFile, processName, suffix=None, suffixList=None):
        """ 
        Return a tuple containing the output file and a boolean flag describing
        whether the output file already exists
        """
        log.info('Preparing to run %s on "%s"' % (processName, inputFile))
        self.processCount += 1
        if suffix:
            outputFile = predictOutputFile(inputFile, suffix)
            return outputFile
        elif suffixList:
            outputFiles = []
            for suffix in suffixList:
                outputFile = predictOutputFile( inputFile, suffix )
                outputFiles.append( outputFile )
            return outputFiles

    def outputFilesExist( self, outputFile=None, outputList=None ):
        if outputFile:
            if fileExists( outputFile ):
                log.info('Output files detected, skipping process...\n')
                return True
            else:
                log.info('Output files not found, running process...')
                return False
        elif outputList:
            if allFilesExist( outputList ):
                log.info('Output files detected, skipping process...\n')
                return True
            else:
                log.info('Output files not found, running process...')
                return False

    def checkOutputFile( self, outputFile ):
        if fileExists( outputFile ):
            log.info('Expected output "%s" found' % outputFile)
        else:
            msg = 'Expected output "%s" not found!' % outputFile
            log.info( msg )
            raise IOError( msg )

    def processCleanup(self, outputFile=None, outputList=None):
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

    def writeDummyFile(self, dummyFile):
        with open(dummyFile, 'w') as handle:
            handle.write('DONE')
        return dummyFile

    def extractCcsFromBasH5(self, inputFile):
        outputFile = self.processSetup( inputFile, 
                                        'extractCcsFromBasH5', 
                                        suffix='fastq' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        extractor = BasH5Extractor( inputFile, outputFile )
        extractor.outputCcsFastq()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def filterFastqFile(self, fastqFile):
        outputFile = self.processSetup( fastqFile, 
                                        'FilterQuality', 
                                        suffix='filter.fastq' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        filter_tool = QualityFilter( fastqFile, outputFile, self.min_accuracy )
        filter_tool()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def separateFastqFile(self, fastqFile):
        outputList = self.processSetup( fastqFile, 
                                        'Fastq.Info', 
                                        suffixList=['fasta', 'qual'] )
        if self.outputFilesExist( outputList=outputList ):
            return outputList
        mothurArgs = {'fastq':fastqFile, 'fasta':'T', 'qfile':'T'}
        logFile = self.getProcessLogFile('fastq.info', True)
        self.factory.runJob('fastq.info', mothurArgs, logFile)
        self.processCleanup( outputList=outputList )
        return outputList

    def alignSequences(self, fastaFile):
        outputFile = self.processSetup( fastaFile, 
                                        'Align.Seqs', 
                                        suffix='align' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':fastaFile,
                      'reference':self.alignmentRef,
                      'flip':'t'}
        logFile = self.getProcessLogFile('align.seqs', True)
        self.factory.runJob('align.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def screenSequences(self, alignFile, start=None, end=None, min_length=None):
        if alignFile.endswith('.align'):
            outputExt = 'good.align'
        elif alignFile.endswith('.fasta'):
            outputExt = 'good.fasta'
        outputFile = self.processSetup( alignFile, 
                                         'Screen.Seqs', 
                                         suffix=outputExt )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'start':start,
                      'end':end,
                      'minlength':min_length}
        logFile = self.getProcessLogFile('screen.seqs', True)
        self.factory.runJob('screen.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def summarizeSequences(self, fastaFile):
        outputFile = self.processSetup( fastaFile, 
                                        'Summary.Seqs', 
                                        suffix='summary' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':fastaFile}
        logFile = self.getProcessLogFile('summary.seqs', True)
        self.factory.runJob('summary.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def parseSummaryFile(self, summaryFile):
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

    def findChimeras(self, alignFile):
        outputFile = self.processSetup( alignFile, 
                                        'UCHIME', 
                                        suffix='uchime.accnos' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'reference':self.chimeraRef}
        logFile = self.getProcessLogFile('chimera.uchime', True)
        self.factory.runJob('chimera.uchime', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def removeSequences(self, alignFile, idFile):
        outputFile = self.processSetup( alignFile, 
                                        'Remove.Seqs', 
                                        suffix='pick.align' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'accnos':idFile}
        logFile = self.getProcessLogFile('remove.seqs', True)
        self.factory.runJob('remove.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile
  
    def filterSequences(self, alignFile):
        outputFile = self.processSetup( alignFile, 
                                        'Filter.Seqs', 
                                        suffix='filter.fasta' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'vertical':'T'}
        logFile = self.getProcessLogFile( 'filter.seqs', True )
        self.factory.runJob( 'filter.seqs', mothurArgs, logFile )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def addQualityToAlignment(self, fastqFile, alignFile):
        outputFile = self.processSetup( alignFile, 
                                        'QualityAligner', 
                                        suffix='fastq' )
        if self.outputFilesExist( outputFile=output ):
            return output
        aligner = QualityAligner( fastqFile, alignFile, outputFile )
        aligner.run()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def maskFastqSequences(self, fastqFile):
        outputFile = self.processSetup( fastqFile, 
                                        'QualityMasker', 
                                        suffix='masked.fastq' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        masker = QualityMasker(fastqFile, outputFile, self.minQv)
        masker.run()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def uniqueSequences( self, alignFile ):
        if alignFile.endswith('.align'):
            outputSuffixes = ['unique.align', 'names']
        elif alignFile.endswith('.fasta'):
            outputSuffixes = ['unique.fasta', 'names']
        outputList = self.processSetup( alignFile,
                                        'Unique.Seqs',
                                        suffixList=outputSuffixes )
        if self.outputFilesExist( outputList=outputList ):
            return outputList
        mothurArgs = {'fasta':alignFile}
        logFile = self.getProcessLogFile('unique.seqs', True)
        self.factory.runJob('unique.seqs', mothurArgs, logFile)
        self.processCleanup( outputList=outputList )
        return outputList

    def preclusterSequences( self, alignFile, nameFile ):
        if alignFile.endswith('.align'):
            outputSuffixes = ['precluster.align', 'precluster.names']
        elif alignFile.endswith('.fasta'):
            outputSuffixes = ['precluster.fasta', 'precluster.names']
        outputList = self.processSetup( alignFile,
                                        'Pre.Cluster',
                                        suffixList=outputSuffixes )
        if self.outputFilesExist( outputList=outputList ):
            return outputList
        mothurArgs = { 'fasta':alignFile,
                       'name': nameFile,
                       'diffs':self.preclusterDiffs }
        logFile = self.getProcessLogFile('pre.cluster', True)
        self.factory.runJob('pre.cluster', mothurArgs, logFile)
        self.processCleanup( outputList=outputList )
        return outputList

    def calculateDistanceMatrix( self, alignFile ):
        outputFile = self.processSetup( alignFile, 
                                        'Dist.Seqs', 
                                        suffix='phylip.dist' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = { 'fasta':alignFile,
                       'calc':'onegap',
                       'countends':'F',
                       'output':'lt' }
        logFile = self.getProcessLogFile('dist.seqs', True)
        self.factory.runJob('dist.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def clusterSequences(self, distanceMatrix):
        if self.clusteringMethod == 'nearest':
            outputSuffix = 'nn.list'
        elif self.clusteringMethod == 'average':
            outputSuffix = 'an.list'
        elif self.clusteringMethod == 'furthest':
            outputSuffix = 'fn.list'
        outputFile = self.processSetup( distanceMatrix, 
                                        'Cluster', 
                                        suffix=outputSuffix )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'phylip':distanceMatrix,
                      'method':self.clusteringMethod}
        logFile = self.getProcessLogFile( 'cluster', True )
        self.factory.runJob( 'cluster', mothurArgs, logFile )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def separateClusterSequences(self, listFile, sequenceFile):
        outputFile = self.processSetup( listFile, 
                                        'ClusterSeparator', 
                                        suffix='list.clusters')
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        separator = ClusterSeparator( listFile, 
                                      sequenceFile,
                                      outputFile,
                                      self.distance, 
                                      self.min_cluster_size )
        separator()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def generateConsensusSequences(self, clusterListFile):
        outputFile = self.processSetup( clusterListFile, 
                                        'ClusterResequencer', 
                                        suffix='consensus')
        if self.outputFilesExist( outputFile=outputFile ):
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
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def cleanupConsensusFolder( self, consensusFile ):
        outputFile = self.processSetup( consensusFile, 
                                        'ConsensusCleanup', 
                                        suffix='consensus.cleanup' )
        if self.outputFilesExist( outputFile=outputFile ):
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
        self.writeDummyFile( outputFile )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def selectFinalSequences( self, consensusFile ):
        outputFile = self.processSetup( consensusFile, 
                                        'SequenceSelector', 
                                        suffix='consensus.selected' )
        if self.outputFilesExist( outputFile=outputFile ):
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
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def outputFinalSequences( self, finalSequenceList ):
        outputFile = self.processSetup( finalSequenceList, 
                                        'SequenceWriter',
                                        suffix='fasta' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        with FastaWriter( outputFile ) as writer:
            with open( finalSequenceList ) as handle:
                for line in handle:
                    sequenceFile = line.strip()
                    copy_fasta_sequences( sequenceFile, writer )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def __call__(self):
        if self.dataType == 'bash5':
            fastqFile = self.extractCcsFromBasH5( self.sequenceFile )
        elif self.dataType == 'fastq':
            fastqFile = self.sequenceFile
        elif self.dataType == 'fasta':
            fastqFile = None
            fastaFile = self.sequenceFile
        # If we have a Fastq, filter low-quality reads and convert to FASTA
        if fastqFile:
            filteredFastq = self.filterFastqFile( fastqFile )
            fastaFile, qualFile = self.separateFastqFile( filteredFastq )
        # Align the Fasta sequences and remove partial reads
        alignedFile = self.alignSequences( fastaFile )
        summaryFile = self.summarizeSequences( alignedFile )
        maxStart, minEnd = self.parseSummaryFile( summaryFile )
        screenedFile = self.screenSequences(alignedFile, 
                                            start=maxStart,
                                            end=minEnd)
        # Identify and remove chimeric reads
        chimeraIds = self.findChimeras( screenedFile )
        noChimeraFile = self.removeSequences( screenedFile, chimeraIds )
        # Filter out un-used columns to speed up re-alignment and clustering
        filteredFile = self.filterSequences( noChimeraFile )
        # If masking is enabled, create an aligned FASTQ, mask the 
        # low-quality bases and remove over-masked reads
        if self.enableMasking:
            alignedFastqFile = self.addQualityToAlignment( fastqFile, filteredFile )
            maskedFastq = self.maskFastqSequences( alignedFastqFile )
            maskedFasta = self.convertFastqToFasta( maskedFastq )
            screenedFasta = self.screenSequences( maskedFasta,
                                                  min_length=self.min_length)
            fileForClustering = screenedFasta
        # Otherwise if masking is disabled, we'll use unique-ify and 
        #    pre-cluster our sequences
        else:
            uniqueFile, nameFile = self.uniqueSequences( filteredFile )
            preclusteredFile, nameFile = self.preclusterSequences( uniqueFile, nameFile )
            fileForClustering = preclusteredFile
        # If enabled, calculate sequence distances and cluster
        if self.enableClustering:
            distanceMatrix = self.calculateDistanceMatrix( fileForClustering )
            listFile = self.clusterSequences( distanceMatrix )
        # If enabled, generate a consensus for each cluster from above
        if self.enableConsensus:
            clusterListFile = self.separateClusterSequences( listFile, fastqFile )
            consensusFile = self.generateConsensusSequences( clusterListFile )
            self.cleanupConsensusFolder( consensusFile )
            selectedFile = self.selectFinalSequences( consensusFile )
            finalFile = self.outputFinalSequences( selectedFile )

if __name__ == '__main__':
    pipeline = rDnaPipeline()
    pipeline()
