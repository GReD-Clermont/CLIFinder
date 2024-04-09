.. image:: https://github.com/GReD-Clermont/CLIFinder/actions/workflows/planemo.yml/badge.svg
    :target: https://github.com/GReD-Clermont/CLIFinder/actions/workflows/planemo.yml

CLIFinder v0.5.2
================


Description
-----------

L1 Chimeric Transcripts (LCTs)  are transcribed from LINE 1 antisense promoter and include the L1 5’UTR sequence in antisense orientation followed by the adjacent genomic region.  
CLIFinder v0.5.2 is a Galaxy tool, specifically designed to identify  potential LCTs from one or several oriented RNA-seq paired-end reads in the human genome.
CLIFinder v0.5.2 is customizable to detect transcripts initiated by different types of repeat elements.



Installation
------------

Some tools are used by CLIFinder that must be installed and added to the PATH (listed in CLIFinder.xml). This is easily done through conda with the command:
::

    conda create -n clifinder samtools=1.9 bedtools=2.26.0gx repeatmasker=4.0.9_p2 bwa=0.7.17 fastx_toolkit=0.0.14 perl=5.26.2 perl-getopt-long=2.50 perl-file-copy-recursive=0.45 perl-parallel-forkmanager=2.02 perl-statistics-r=0.34 r-base=3.5.1 r-plyr=1.8.5 bioconductor-genomicranges=1.34.0 wget=1.20.1

You should then be able to use CLIFinder by activating the conda environment and running the script with:
::

    conda activate clifinder
    perl script/CLIFinder.pl

Galaxy uses conda to solve tool requirements starting with Galaxy release 16.01 so that is the minimum version required to install this tool (which can be done through the toolshed: https://toolshed.g2.bx.psu.edu/repository?repository_id=5c73d1cf20ab37c3).



Usage
-----

The command you need to use to run the script is as follows:
::

    CLIFinder.pl --first <first fastq of paired-end set 1> --name <name 1> --second <second fastq of paired-end set 1> [--first <first fastq of paired-end set 2> --name <name 2> --second <second fastq of paired-end set 2> ...] --ref <reference genome> [--build_ref] --TE <transposable elements> [--build_TE] --html <results.html> --html-path <results directory> [options]


    Arguments:
        --first <fastq file>    First fastq file to process from paired-end set
        --name <name>           Name of the content to process
        --second <fastq file>   Second fastq file to process from paired-end set
        --ref <reference>       Fasta file containing the reference genome
        --TE <TE>               Fasta file containing the transposable elements
        --rmsk <text file>      Tab-delimited text file (with headers) containing reference repeat sequences (e.g. rmsk track from UCSC)
        --refseq <text file>    Tab-delimited file (with headers) containing reference genes (e.g. RefGene.txt from UCSC)
        --html                  Main HTML file where results will be displayed
        --html-path             Folder where results will be stored

    For any fasta file, if a bwa index is not provided, you should build it through the corresponding '--build_[element]' argument

    Options:
        --rnadb <RNA db>        Blast database containing RNA sequences (default: empty)
        --estdb <EST db>        Blast database containing EST sequences (default: empty)
        --size_read <INT>       Size of reads (default: 100)
        --BDir <0|1|2>          Orientation of reads (0: undirectional libraries, 1: TEs sequences in first read in pair, 2: TEs sequences in second read in pair) (default: 0)
        --size_insert <INT>     Maximum size of insert tolerated between R1 and R2 for alignment on the reference genome (default: 250)
        --min_L1 <INT>          Minimum number of bp matching for L1 mapping (default: 50)
        --mis_L1 <INT>          Maximum number of mismatches tolerated for L1 mapping (default: 1)
        --min_unique <INT>      Minimum number of consecutive bp not annotated by RepeatMasker (default: 33)
        --species <STRING>      Species to use in RepeatMasker (default: human)
        --threads <INT>         Number of threads (default: 1)

    For Blast database files, if a fasta is provided, the database can be built with '--build_[db]'. Otherwise, provide a path or URL. "tar(.gz)" files are acceptable, as well as wild card (rna*) URLs.

