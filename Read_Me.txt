####################################################################################
##############################CLiFinder v0.4.1 README###############################
####################################################################################

#############
#Description#
#############
L1 Chimeric Transcripts (LCTs)  are transcribed from LINE 1 antisense promoter and include the L1 5’UTR sequence in antisense orientation followed by the adjacent genomic region.  
CLIFinder v0.4.1 is a Galaxy tool, specifically designed to identify  potential LCTs from one or several oriented RNA-seq paired-end reads in the human genome.  
CLIFinder v0.4.1 is customizable to detect transcripts initiated by different types of repeat elements.


###############
#Prerequisites#
###############

1.Unix system with A Galaxy server (release jully 2014 or later installed)

2. Some tools are used by CLiFinder and must be installed and added to the Path. They are generally already present when you install Galaxy server.

	- Bwa aligner: you can obtain it here: https://sourceforge.net/projects/bio-bwa/files/ . Please download version  0.7.12-r1039 or higher
	- BLAST software suite: Installers and source code are available from  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ . we use version 2.2.28+
	- BedTools powerful toolset for genome arithmetic is also needed. It should be found here: http://bedtools.readthedocs.io/en/latest/ . We recommend to use v2.17.0 or higher.
	- FastX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing. Here is the installation address: http://hannonlab.cshl.edu/fastx_toolkit/download.html/ Version higher yhan v.0.0.6 is needed.

3. RepeatMasker: is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences. You can obtain a free copy here: http://www.repeatmasker.org/RMDownload.html/ We use version 4.0.6 to develop CLifinder. It must be added to the PATH.

4. R project version higher than 3.1 is needed with libraries "GenomicRanges" and "plyr" installed. You can find respectively these two libraries here: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html and https://cran.r-project.org/web/packages/plyr/index.html

##############
#Installation#
##############

The process has to be completed by an administrator of your Galaxy server to install CLiFinder.

1. Download CLiFinder
You can find Clinfinder here: https://github.com/GReD-Clermont/CLIFinder/

2. Put the tool into Galaxy's tools directory
You need to add files into tools/ directory , where all tool-related files are stored, within your Galaxy installation.

3. Make Galaxy aware of the new tool CLiFinder
Now that the tool and its definition file are ready, the final step is to make Galaxy aware of the new files. 
Galaxy recognizes installed tools by reading the tool_conf.xml tool configuration file. Thus, letting Galaxy know about the new tool is as easy as adding a few lines to the tool_conf.xml file located in the config/ directory of the Galaxy installation. New tools can either be added to existing sections or added to new sections defined in the following way:

 <section name="NewTools" id="mTools">
    <tool file="LINE_html_4_1.xml" />
 </section>
 
 4. Start or Restart Galaxy to use it.
 
 
 


