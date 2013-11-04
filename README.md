# rDnaTools #

rDnaTools is a python package of tools and pipelines for working with
ribosomal DNA sequence data generated with the PacBio(R) SMRT sequencing.
rDnaTools works by wrapping existing tools from microbial ecology,
primarily the Mothur suite of utilities.  Currently rDnaTools implements
a single pipeline for the export, filtering, and cluster of 16S sequences.

Though primarily intended for use in analyzing 16S rDNA sequences, the
same tools and approaches should apply equally well to 18S, 23S, or ITS
sequences, provided that suitable reference sequences are supplied.


## Requirements ##

The core functionality of rDnaTools is built upon Python2.7 using the
pbcore framework for accessing PacBio data files.  In addition rDnaTools
wraps the functionality from a number of stand-alone commandline tools
that must available for the package to function
* Python 2.7 (www.python.org)
* pbcore (www.github.com/PacificBiosciences/pbcore)
* pbdagcon (www.github.com/PacificBiosciences/pbdagcon)
* Blasr (www.github.com/PacificBiosciences/blasr)
* Mothur (www.mothur.org)


## rDnaPipeline ##

The primary tool for analyzing rDNA sequence data is a script called
"rDnaPipeline", which takes as an input PacBio sequence data from ribosomal
DNA amplicons.  The pipeline accepts data in either FOFN, BAS.H5, FASTA, or
FASTQ format, and runs a sequential series of analyses, similar to Mothur`s
Batch Mode.  The analysis is based on Mothur's recommended SOP for analyzing
454 rDNA sequence tags, with some modifications to account for the unique
nature of PacBio's data.

Since rDnaPipeline analyzes PacBio Circular-Consensus Sequence (CCS) data,
The basic call to the rDnaPipeline will vary slightly depending on what version
of SMRT Analysis is being used, as the methods by which CCS data is generated
and presented has changed.

If the rDNA sequence generated on a PacBio RS running SMRT Analysis v2.0, then
the basic call to rDnaPipeline.py will look as follows:
rDnaPipeline.py FOFN -n PROCS -A ALIGNMENT_REF -C CHIMERA_REF

If the rDNA sequence generated on a PacBio RS running SMRT Analysis v2.1, then
running rDnaPipeline is a two-step process.  First the Reads_of_Insert protocol
must be run as normal on the sequencing data to generate CCS data for the data-set.
Then rDnaPipeline.py is called as follows:
rDnaPipeline.py READS_OF_INSERT -r FOFN -n PROCS -A ALIGNMENT_REF -C CHIMERA_REF

For the reference files, we recommend using the curated SILVA alignments provided
on the Mothur website (http://www.mothur.org/wiki/Silva_reference_files).


## Citation ##

rDnaTools would not have been possible were it not for the hard work of the
existing Microbial Ecology community, and their existing tools for analyzing
ribosomal DNA sequence data.  Since the core of the analyses wrapped by 
rDnaTools come from the Mothur suite, please cite their publication if you 
use rDnaTools in your work:

Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, 
community-supported software for describing and comparing microbial communities. 
Appl Environ Microbiol, 2009. 75(23):7537-41.