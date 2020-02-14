#!/usr/bin/env perl

################################################
#Declaration of necessary libraries#############
################################################

use strict;
use warnings;
use Parallel::ForkManager;
use POSIX;
use Statistics::R;
use Getopt::Long qw(HelpMessage VersionMessage);
use File::Basename;
use File::Copy::Recursive;
use FindBin qw($Bin);
use Archive::Tar;

our $VERSION = '0.5.0';


#####################################################################
#Definition options of execution according to the previous variables#
#####################################################################

GetOptions(
  "first|1=s"     => \my @fastq1,
  "second|2=s"    => \my @fastq2,
  "name=s"        => \my @name,
  "html=s"        => \my $html,
  "html_path=s"   => \my $html_repertory,
  "TE=s"          => \my $TE,
  "ref=s"         => \my $ref,
  "rnadb:s"       => \my $rna_source,
  "estdb:s"       => \my $est_source,
  "build_TE"      => \my $build_TE,
  "build_ref"     => \my $build_ref,
  "build_rnadb"   => \my $build_rnadb,
  "build_estdb"   => \my $build_estdb,
  "rmsk=s"        => \my $rmsk_source,
  "refseq=s"      => \my $refseq,
  "min_unique:i"  => \(my $prct = 33),
  "size_insert:i" => \(my $maxInsertSize = 250),
  "size_read:i"   => \(my $size_reads = 100),
  "BDir:i"        => \(my $Bdir = 0),
  "min_L1:i"      => \(my $min_L1 = 50),
  "mis_L1:i"      => \(my $mis_L1 = 2),
  "threads:i"     => \(my $threads = 1),
  "help"          => sub { HelpMessage(0); },
  "version"       => sub { VersionMessage(0); },
) or HelpMessage(1);

HelpMessage(1) unless @fastq1 && @fastq2 && @name && defined($TE) && defined($ref) && defined($rmsk_source) && defined($refseq) && defined($html) && defined($html_repertory);

my $iprct = 100 - (($prct / $size_reads)*100) ;
my $mis_auth = $size_reads - $min_L1 + $mis_L1 ;
my $eprct = ($iprct * $size_reads) /100;
my $dprct = ((100-$iprct) * $size_reads) / 100;

################################################
#Construct index of ref and TE if doesn't exist#
################################################

`(bwa index $ref)` if ($build_ref);
`(bwa index $TE)` if ($build_TE);

############################################
#Create repository to store resulting files#
############################################

mkdir $html_repertory;

##########################################
#Define hash                             #
##########################################

my %frag_exp_id;

##########################
#Data file we have to use#
##########################

print STDOUT "Extracting data from rmsk file\n";
my $line_only=$html_repertory.'/'.'line_only.txt';
my $rmsk = $html_repertory.'/rmsk.bed'; 
filter_convert_rmsk($rmsk_source, $rmsk, $line_only);

##############################
# Analyse of each fastq file #
##############################

my @garbage; my $num = 0;
foreach my $tabR (0..$#fastq1)
{
  ###################################################
  # Paired end mapping against L1 promoter sequences#
  ###################################################
  
  print STDOUT "Alignment of $name[$tabR] to L1\n";
  my $sam = $html_repertory.'/'.$name[$tabR]."_L1.sam"; push(@garbage, $sam);
  halfmap_paired($TE, $fastq1[$tabR], $fastq2[$tabR], $sam, $threads, $mis_auth);
  print STDOUT "Alignment done\n";
  
  ##################################################
  # Creation of two fastq for paired halfed mapped:#
  # - _1 correspond to sequences mapped to L1      #
  # - _2 correspond to sequences unmapped to L1    #
  ##################################################

  my $hm_reads_1 = $html_repertory.'/'.$name[$tabR]."_halfmapped_1.fastq"; push(@garbage, $hm_reads_1);
  my $hm_reads_2 = $html_repertory.'/'.$name[$tabR]."_halfmapped_2.fastq"; push(@garbage, $hm_reads_2);
  
  print STDOUT "Getting pairs with one mate matched to L1 and the other mate undetected by repeatmasker as a repeat sequence\n";
  my $out_ASP_1 = $html_repertory.'/'.$name[$tabR]."_1.fastq"; push(@garbage, $out_ASP_1);
  my $out_ASP_2 = $html_repertory.'/'.$name[$tabR]."_2.fastq"; push(@garbage, $out_ASP_2);
  
  ##split mate that matched to L1 and others##
  my $half_num_out = get_half($sam, $mis_L1, $min_L1, $Bdir, $hm_reads_1, $hm_reads_2);
  print STDOUT "Number of half mapped pairs: $half_num_out\n";
  
  ##pairs obtained after repeatmasker on the other mate##
  my $left = sort_out($threads, $hm_reads_1, $hm_reads_2, $out_ASP_1, $out_ASP_2, $dprct, $eprct, $html_repertory);
  print STDOUT "Number of pairs after repeatmasker: $left\n";
  
  ##################################################
  # Alignment of halfed mapped pairs on genome     #
  ##################################################
  print STDOUT "Alignment of potential chimeric sequences to the genome\n";
  $sam = $html_repertory.'/'.$name[$tabR]."_genome.sam";
  push(@garbage, $sam);
  align_genome($ref, $out_ASP_1, $out_ASP_2, $sam, $maxInsertSize, $threads);
  print STDOUT "Alignment done\n";
  
  ##compute the number of sequences obtained after alignment ##
  
  $left = `samtools view -@ $threads -Shc $sam`;
  chomp $left; $left = $left/2;
  print STDOUT "Number of sequences: $left\n";

  ##################################################
  # Create bedfiles of potential chimerae          #
  # and Know repeat sequences removed              #
  ##################################################
  
  print STDOUT "Looking for chimerae\n";
  results($html_repertory, $sam, $name[$tabR], $iprct, \%frag_exp_id, $rmsk, $num, \@garbage);
  $num++;
}

##define files variables ##

my $repfirst = $html_repertory.'/first.bed'; push(@garbage,$repfirst);
my $repsecond = $html_repertory.'/second.bed'; push(@garbage,$repsecond);
my $repMfirst = $html_repertory.'/firstM.bed'; push(@garbage,$repMfirst);
my $repMsecond = $html_repertory.'/secondM.bed'; push(@garbage,$repMsecond);
#my $covRepMsecond = $html_repertory.'/covSecondM.bed'; push(@garbage,$covRepMsecond);

##Concate all files for first and second mate results ##

`cat $html_repertory/*-first.bed > $repfirst`; #*/
`cat $html_repertory/*-second.bed > $repsecond`; #*/

## Sort Files and generate files that merge reads in the same locus ##
print STDOUT "Sort files and merge reads in the same locus\n";
`bedtools sort -i $repfirst | bedtools merge -c 4,5 -o collapse,max -d 100 -s > $repMfirst `;
`bedtools sort -i $repsecond | bedtools merge -c 4,5 -o collapse,max -d 100 -s > $repMsecond `;

my (%frag_uni, @second_R, @second_exp, @results);
my $merge_target = $html_repertory.'/target_merged.bed'; push(@garbage, $merge_target);
my $merge = $html_repertory.'/merged.bed'; push(@garbage, $merge);

open (my $mT, ">".$merge_target) || die "Cannot open $merge_target\n";
open (my $m, ">".$merge) || die "Cannot open $merge\n";
open (my $in, $repMsecond) || die "Cannot open $repMsecond\n";
my $cmp = 0;
while (<$in>)
{
  chomp $_;
  my @tmp = (0) x scalar(@fastq1);
  my @line = split /\t/, $_;
  my @names =split /,/, $line[4];
  foreach my $n (@names){$n =~/(.*?)\/[12]/; $frag_uni{$1} = $cmp; $tmp[$frag_exp_id{$1}]++; }
  $second_exp[$cmp] = \@tmp;
  $cmp++;
  push @second_R, [$line[0],$line[1],$line[2],$line[3]];
}

$cmp = 0;
open ($in, $repMfirst) || die "Cannot open $repMfirst\n";
while (<$in>)
{
  chomp $_;
  my %sec;
  my @line = split /\t/, $_;
  my @names =split /,/, $line[4];
  my @tmp = (0) x scalar(@fastq1);
  foreach my $n (@names){$n =~/(.*?)\/[12]/; $tmp[$frag_exp_id{$1}]++; }
  foreach my $n (@names)
  {
    $n =~/(.*?)\/[12]/;
    unless (exists ($sec{$frag_uni{$1}}) )
    {
      my @lmp = ($line[0], $line[1], $line[2], $line[3]);
      foreach my $exp_N (@tmp){ push @lmp, $exp_N;}
      push (@lmp, $second_R[$frag_uni{$1}]->[0], $second_R[$frag_uni{$1}]->[1], $second_R[$frag_uni{$1}]->[2], $second_R[$frag_uni{$1}]->[3]);
      foreach my $exp_N (@{$second_exp[$frag_uni{$1}]}){ push @lmp, $exp_N;}
      
      my $name = $cmp.'-'.$second_R[$frag_uni{$1}]->[0].'-'.$second_R[$frag_uni{$1}]->[1].'-'.$second_R[$frag_uni{$1}]->[2];
      print $mT $second_R[$frag_uni{$1}]->[0]."\t".$second_R[$frag_uni{$1}]->[1]."\t".$second_R[$frag_uni{$1}]->[2]."\t$name\t29\t".$second_R[$frag_uni{$1}]->[3]."\n";
      
      my ($b, $e) = (0,0);
      if ($line[1] < $second_R[$frag_uni{$1}]->[1])
      {
        $b = $line[1] - 1000; $e = $second_R[$frag_uni{$1}]->[2] + 1000;
        $name = $cmp.'-'.$line[0].'-'.$b.'-'.$e;
        print $m $line[0]."\t".$b."\t".$e."\t$name\t29\t".$second_R[$frag_uni{$1}]->[3]."\n";
      }
      else
      {
        $b = $second_R[$frag_uni{$1}]->[1] - 1000; $e = $line[2] + 1000;
        $name = $cmp.'-'.$line[0].'-'.$b.'-'.$e;
        print $m $line[0]."\t".$b."\t".$e."\t$name\t29\t".$second_R[$frag_uni{$1}]->[3]."\n";
      }
      $results[$cmp] = \@lmp;
      $cmp++;
    }
    $sec{$frag_uni{$1}} = undef;
  }
}
close $mT; close $m;

my $fasta = $html_repertory.'/target_merged.fasta'; push(@garbage, $fasta);
my $extend = $html_repertory.'/extend.fasta'; push(@garbage, $extend);

`bedtools getfasta -name -fi $ref -bed $merge -fo $extend`;
`bedtools getfasta -name -fi $ref -bed $merge_target -fo $fasta`;

################################################
#Blast against human rna and est, if provided  #
################################################

my $rna;
my $est;
if(defined($rna_source))
{
  ##get databases for est and rna
  print STDOUT "Getting blast databases for rna\n";
  my $rna_db = get_blastdb_from_source($rna_source, $build_rnadb, 'rna', $html_repertory);

  print STDOUT "Blast against human rna\n";
  my $tabular = $html_repertory."/chimerae_rna.tab"; push(@garbage, $tabular);
  blast($rna_db, $fasta, $tabular, $threads);
  $rna = extract_blast($tabular);

  # Clean RNA blast database if in html dir
  if(rindex($rna_db, $html_repertory, 0) == 0)
  {
    my $toErase = $rna_db.'.*';
    unlink glob "$toErase";
  }
}
if(defined($est_source))
{
  print STDOUT "Getting blast databases for est\n";
  my $est_db = get_blastdb_from_source($est_source, $build_estdb, 'est', $html_repertory);

  print STDOUT "Blast against human est\n";
  my $tabular2 = $html_repertory."/chimerae_est.tab"; push(@garbage, $tabular2);
  blast($est_db, $fasta, $tabular2, $threads);
  $est = extract_blast($tabular2);

  # Clean EST blast database if in html dir
  if(rindex($est_db, $html_repertory, 0) == 0)
  {
    my $toErase = $est_db.'.*';
    unlink glob "$toErase";
  }
}

################################################
#Create Results html file                      #
################################################
print STDOUT "Save results in file\n";
save_csv(\@fastq1, \@name, \@results, $line_only, $refseq, $html_repertory);

print STDOUT "Create HTML\n";
html_tab(\@fastq1, \@name, \@results, $rna, $est, $html, $html_repertory);

$extend = $extend.'*';
push(@garbage, glob($extend));
push(@garbage, $line_only);
push(@garbage, $rmsk);
unlink @garbage;

print STDOUT "Job done!\n";
  
########################################### END MAIN ##########################################################


##############################################################################################################
############################################     FUNCTIONS   #################################################
##############################################################################################################


############################################################
## Function to convert rmsk table to bed and line_only #####
############################################################
## @param:                                                 #
##       $source: rmsk table file                          #
##       $bed: rmsk bed file                               #
##       $line_only: rmsk table file with only LINE        #
############################################################

sub filter_convert_rmsk
{
  my ($source, $bed, $line_only) = @_;
  open(my $input, $source) || die "Cannot open rmsk file! $!\n"; ## Open source file
  open(my $bedfile, ">".$bed) || die "Cannot open output bed file for rmsk! $!\n"; ## Open bed file
  open(my $linefile, ">".$line_only) || die "Cannot open output LINE-only file for rmsk! $!\n"; ## Open line_only file
  my @headers;
  my %indices;

  print $linefile "#filter: rmsk.repClass = 'LINE'\n";

  while(<$input>)
  {
    chomp $_;
    if($. == 1)
    {
      if(substr($_, 0, 1) ne "#")
      {
        die "rmsk file does not have header starting with #\n";
      }
      else
      {
        print $linefile "$_\n";
        my $firstline = substr($_, 1);
        @headers = split(/\t/, $firstline);
        @indices{@headers} = 0..$#headers;
      }
    }
    else
    {
      my @line = split(/\t/,$_);
      if($line[$indices{"repClass"}] eq "LINE")
      {
        print $linefile "$_\n";
      }
      print $bedfile "$line[$indices{'genoName'}]\t$line[$indices{'genoStart'}]\t$line[$indices{'genoEnd'}]\t$line[$indices{'repName'}]\t$line[$indices{'swScore'}]\t$line[$indices{'strand'}]\n";
    }
  }
  close $input;
  close $bedfile;
  close $linefile;
}


############################################################
## Function to get blast db from the specified source ######
############################################################
## @param:                                                 #
##       $source: db source (URL or path)                  #
##       $build_db: whether the db should be created       #
##       $name: name of the db that could be created       #
##       $dest_dir: where the data can be placed           #
## @return:                                                #
##       $path: blast db path                              #
############################################################

sub get_blastdb_from_source
{
  my ($source, $build_db, $name, $dest_dir) = @_;
  # Assume source is just db path
  my $path = $source;
  my ($file) = $path =~ m~([^/\\]*)$~;
  my $dbname = $file;
  my @garbage;

  if($build_db)
  {
    $dbname = $name;
    $path = $dest_dir.'/'.$name;
    print STDOUT "Making $dbname blast database\n";
    `makeblastdb -in $source -dbtype nucl -out $path`;
  }
  else
  {
    # Check if source is URL
    if(index($source, ":/") != -1)
    {
      my $url = $source;
      if($file =~ /\*/)
      {
        $url =~ s/\Q$file//;
        print STDOUT "Downloading blast database from $url\n";
        `wget -q -N -r -nH -nd -np --accept=$file $url -P $dest_dir`;

        # Assume regexp matches db name
        $dbname =~ s/\*$//;
      }
      else
      {
        print STDOUT "Downloading blast database from $url\n";
        `wget -q -N $source -P $dest_dir`;
        push(@garbage, $dest_dir.'/'.$file);
      }
      if($? == 0)
      {
        print "Downloaded database\n";
      }
      else
      {
        print "Error while downloading database\n";
      }
      $path = $dest_dir.'/'.$dbname;
    }
    if(index($file, ".") != -1)
    {
      if(index($file, ".tar") != -1)
      {
        ## Extract tar files
        print STDOUT "Extracting blast database from $file\n";
        my @properties = ('name');
        my $tar=Archive::Tar->new();
        $tar->setcwd($dest_dir);
        $tar->read($path);
        my @files = $tar->list_files(\@properties);
        $tar->extract();
        $tar->clear();
        unlink @garbage;

        ## Get dbname from filenames
        my @parts = split(/\./, $files[0]);
        $dbname = $parts[0];
        $path = $dest_dir.'/'.$dbname;
        print STDOUT "Extracted database\n";
      }
      else
      {
        print STDOUT "Unexpected file format for database"
      }
    }
  }
  print "Using $dbname database\n";
  return $path;
}

############################################################
##Function that aligned paired-end reads on a genome########
############################################################
## @param:                                                 #
##       $index: referential genome                        #
##       $fasq1: first paired end file                     #
##       $fasq2: second paired end file                    #
##       $sam: alignment output file                       #
##       $maxInsertSize: maximum size of insert            #
##       $threads: number of threads used                  #
############################################################

sub align_genome
{
  my ($index, $fastq1, $fastq2, $sam, $maxInsertSize, $threads) = @_ ;
  my @L_garbage =();
  my $sai1 = $sam."_temporary.sai1"; push @L_garbage,$sai1;
  my $sai2 = $sam."_temporary.sai2"; push @L_garbage,$sai2;
  `bwa aln -o4 -e1000 -t $threads $index $fastq1 > $sai1`;
  `bwa aln -o4 -e1000 -t $threads $index $fastq2 > $sai2`;
  ## -A force the insertion size
  `bwa sampe -s -A -a $maxInsertSize $index $sai1 $sai2 $fastq1 $fastq2 | samtools view -@ $threads -F4 -f 2 -Sh /dev/stdin -o $sam`;
  unlink @L_garbage;
}


############################################################
##Function get_half get alignement on TE                 ###
############################################################
## @param:                                                 #
##       $sam: name of alignement file                     #
##       $mis_L1: maximum number of mismatches             #
##       $min_L1: minimum number of bp matching            #
##       $Bdir: reads orientation                          #
##       $fastq1: fastq file with matching reads (1)       #
##       $fastq2: fastq file with matching reads (2)       #
##                                                         #
## @return:                                                #
##       $half_num_out: number of alignment saved          #
############################################################

sub get_half
{
  ## store name of file
  my $sam = shift;
  my $mis_L1 = shift;
  my $min_L1 = shift;
  my $Bdir = shift;
  my $fastq1 = shift;
  my $fastq2 = shift;
  my $filtered_sam = $sam.'_filtered.sam'; push(@garbage, $filtered_sam);
  my $cmp = 0;
  my $sequence = '';
  my $score = '';
  open(my $in, $sam) || die "Cannot open $sam file! ($!)\n";
  open(my $out, '>'.$filtered_sam) || die "Cannot open $filtered_sam file! ($!)\n";
  
  ##read file##
  while(<$in>)
  {
    chomp $_;
    ##We copy the headers##
    if($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
    {
      print $out "$_\n";
      next;
    }
    ##We split line in several part##
    my @line = split (/\t/,$_);
   
    ##Find if alignemets have min_L1 consecutives bases mapped on R1 ##
    if ($_ =~/NM:i:(\d+)\t.*MD:Z:(.*)/)
    {
       my $misT = $1; my $MD = $2; my $maxT = 0;
       $MD = $1 if ($MD =~ /(.*)\t/);
       $MD =~ s/\^//g;
       my @tab = split /[ATCG]/,$MD;
       my $tot = 0;
       my $accept = 0;
       if ($misT <= $mis_L1){$accept = 1;}
       else
       {
         if ( $mis_L1 > scalar(@tab) ) { $maxT = scalar(@tab); }
         else{ $maxT = $mis_L1; }
         for (my $t = 0; $t < $maxT ; $t++)
         {
           $tot += $tab[$t] if $tab[$t] ne '';
         }
         $accept = 1 if $tot >= $min_L1;
       }
       ## if sequence is not accepted we go to the next sequence ##
       next if $accept == 0;
    }
    
    ##looking for flag of the alignment and keep only good reads##
    ##Find if it aligns on R1 or on R2##
    unless($line[1] & 2 || ($line[1] & 4 && $line[1] & 8))
    {
      if ( $Bdir == 0
              || ($Bdir == 1 && (($line[1] & 64 && $line[1] & 8) || ($line[1] & 128 && $line[1] & 4)))
              || ($Bdir == 2 && (($line[1] & 128 && $line[1] & 8) || ($line[1] & 64 && $line[1] & 4))) )
      {
        $cmp++;
      }
      # If Bdir is 0, modify flags so that reads that map are first in pair
      if($Bdir == 0)
      {
        if($line[1] & 128 && $line[1] & 8) { $line[1] = $line[1] - 64; } # Set second mapping read as first in pair
        else { if($line[1] & 64 && $line[1] & 4) { $line[1] = $line[1] + 64; }} # Set unmapping read with mapping mate as second in pair
      }
      print $out join("\t", @line), "\n";
    }
  }
  close $in;
  close $out;

  # Generate fastq files
  my $report;
  if($Bdir == 2)
  {
    `samtools view -h -G 72 $filtered_sam | samtools fastq -G 132 -1 $fastq2 -2 $fastq1 -s /dev/null /dev/stdin`;
  }
  else
  {
    `samtools view -h -G 136 $filtered_sam | samtools fastq -G 68 -1 $fastq1 -2 $fastq2 -s /dev/null /dev/stdin`;
  }
  unlink $filtered_sam;
  return $cmp;
}

############################################################
##Function sort_out: extract paired end reads            ###
############################################################
## @param:                                                 #
##       $threads: number of threads used                  #
##       $in1: input file halfmapped 1                     #
##       $in2: input file halfmapped 2                     #
##       $out1: output file accepted 1                     #
##       $out2: output file accepted 2                     #
##       $dprct: number of bp not annotated by RepeatMasker#
##       $eprct: number of repeated bases tolerated        #
##       $html_repertory: folder for html files            #
############################################################

sub sort_out
{
  my ($threads, $in1, $in2, $out1, $out2, $dprct, $eprct, $html_repertory) = @_;
  my ($name,$path) = fileparse($out2,'.fastq');
  my %repeat;
  my @garbage = (); my $cmp = 0;
  my $repout = $html_repertory.'/'.$name."_repout/";
  my $list = $html_repertory.'/'.$name.".list"; push(@garbage, $list);
  my $fa = $html_repertory.'/'.$name.".fa"; push(@garbage, $fa);
  mkdir $repout;
  my %notLine;
  
  ## Transform fastq file to fasta
  `fastq_to_fasta -i $in2 -o $fa -Q33`;
  
  ##Launch RepeatMasker on fasta file
  `RepeatMasker -s -pa $threads -dir $repout -engine hmmer -species human $fa`;
  my $repfile = $repout.$name.".fa.out";
  open (my $rep, $repfile) || die "Cannot open $repfile ($!)\n";
  while(<$rep>)
  {
    chomp;
    ## test the percent of repeats ##
    my $string = $_;
    $string =~ s/^\s+//;
    next unless ($string =~ /^\d/);
    $string =~ s/\s+/ /g;
    my @line = split (' ',$string);
    if ( exists($repeat{$line[4]}) )
    {
      $repeat{$line[4]}->[0] = $line[5] if $repeat{$line[4]}->[0] > $line[5];
      $repeat{$line[4]}->[1] = $line[6] if $repeat{$line[4]}->[1] < $line[6];
    }
    else{ $repeat{$line[4]} = [$line[5],$line[6]];}
  }
  close $rep;
  
  ## store in list if pair passed the repeat test ##
  open (my $fq, $in2) || die "Cannot open $in2 ($!)\n";
  open (my $lst, '>'.$list) || die "Cannot open $list ($!)\n";
  while(<$fq>)
  {
    chomp $_;
    next unless(index($_, '@') == 0 && $. % 4 == 1);
    my $seqname = substr($_, 1);
    my $repseq = $repeat{$seqname};
    unless(defined($repseq) && ($repseq->[0] <= $dprct && $repseq->[1] >= $eprct))
    {
      print $lst "$seqname\n";
      $cmp++;
    }
  }
  close $fq;
  close $lst;

  ##write resulting reads in both files for paired ##
  `seqtk subseq $in1 $list > $out1`;
  `seqtk subseq $in2 $list > $out2`;
  
  ##drop files and directories generated by repeatmasker##
  my $toErase = $repout.'*';
  unlink glob "$toErase";
  unlink @garbage; rmdir $repout;
  return $cmp;
}

############################################################
##Function to get half-mapped paired-end reads on a ref    #
############################################################
## @param:                                                 #
##       $index: referential file                          #
##       $fasq1: first file paired end reads               #
##       $fasq2: second file paired end reads              #
##       $sam: output alignment file                       #
##       $threads: number of threads used                  #
##       $mis: tolerated mismatches                        #
############################################################
sub halfmap_paired
{
  my ($index, $fastq1, $fastq2, $sam, $threads, $mis) = @_ ;
  my @garbage = ();
  my $sai1 = $sam."_temporary.sai1"; push @garbage,$sai1;
  my $sai2 = $sam."_temporary.sai2"; push @garbage,$sai2;
  
  ##alignement with bwa
  
  `bwa aln -n $mis -t $threads $index $fastq1 > $sai1`;
  `bwa aln -n $mis -t $threads $index $fastq2 > $sai2`;
  `bwa sampe $index $sai1 $sai2 $fastq1 $fastq2 | samtools view -@ $threads -h -F 2 -G 12 -o $sam`;
  
  ## delete temporary single aligned files
  unlink @garbage;
}

############################################################
##Function that aligned reads on a referential           ###
############################################################
## @param:                                                 #
##       $index: referential file                          #
##       $fasq: reads file                                 #
##       $sam: output alignment file                       #
##       $maxInsertSize: maximum size of insert            #
##       $threads: number of threads used                  #
############################################################

sub align
{
  my ($index, $fastq, $sam, $maxInsertSize, $threads ) = @_ ;
  `bwa aln -o4 -e$maxInsertSize -t $threads $index $fastq | bwa samse $index /dev/stdin $fastq > $sam `;
}

############################################################
##Function results computes bed files for result         ###
############################################################
## @param:                                                 #
##       $out_repository: repository to store results      #
##       $file: sam file resulting of alignement           #
##       $name: name of paireds end reads file             #
##       $iprct: percentage of repeats tolerated           #
##       $hashRef: store number for each first read value  #
##       $rmsk: UCSC repeat sequences                      #
##       $ps: number of the paired end file                #
##       $garbage_ref: reference to garbage array          #
############################################################

sub results
{
  my ($out_repertory, $file, $name, $iprct, $hashRef, $rmsk, $ps, $garbage_ref) = @_;
  my $namefirst = $out_repertory.'/'.$name.'-first.bed'; push(@$garbage_ref, $namefirst);
  my $namesecond = $out_repertory.'/'.$name.'-second.bed'; push(@$garbage_ref, $namesecond);
  
  ## store reads mapped in proper pair respectively first and second in pair in bam files and transform in bed files##
  `samtools view -Sb -f66 $file | bedtools bamtobed -i /dev/stdin > temp_name_first`;
  `samtools view -Sb -f130 $file | bedtools bamtobed -i /dev/stdin > temp_name_second`;
  
  ##compute converage of second mate on rmsk##
  my $baseCov = 0;
  my %IdCov = ();
  my @coverage = `bedtools coverage -a temp_name_second -b $rmsk`;
  
  
  ## store coverage fraction ##
  foreach my $covRmsk (@coverage)
  {
    chomp $covRmsk;
    my @split_cov = split /\t/, $covRmsk;
    ##generate identifier for IdCov ##
    $split_cov[3] =~ /(.*?)\/[12]/;
    ##store value in IdCov ##
    if (!exists($IdCov{$1}))
    {
      $IdCov{$1} = $split_cov[-1];
    }
    else
    {
      $IdCov{$1} = $split_cov[-1] if $split_cov[-1] > $IdCov{$1};
    }
  }
  
  ## get only first mate that have less tant $iprct repeats ##
  open (my $tmp_fi, 'temp_name_first') || die "Cannot open $namefirst!\n";
  open (my $nam_fi, ">".$namefirst) || die "Cannot open $namefirst!\n";
  while (<$tmp_fi>)
  {
    my @line = split /\t/, $_;
    $line[3] =~ /(.*?)\/[12]/;
    
    if ($IdCov{$1} <= $iprct/100)
    {
      print $nam_fi $_;
      
      ${$hashRef}{$1}= $ps;
    }
  }
  close $tmp_fi; close $nam_fi;
  
  
  ## get only  second mate that have less than $iprct repeats ##

  open (my $tmp_sec, 'temp_name_second') || die "Cannot open $namesecond!\n";
  open (my $nam_sec, ">".$namesecond) || die "Cannot open $namesecond!\n";
  while (<$tmp_sec>)
  {
    my @line = split /\t/, $_;
    $line[3] =~ /(.*?)\/[12]/;
    if ($IdCov{$1} <= $iprct/100)
    {
      print $nam_sec $_;
    }
  }
  close $tmp_sec; close $nam_sec;
}


############################################################
##Function blast: blast nucleotide sequences on ref      ###
############################################################
## @param:                                                 #
##       $db: database where to search                     #
##       $fasta: file containing nucleotide sequences      #
##       $tabular: out file name                           #
##       $threads: number of threads used                  #
############################################################



sub blast
{
  my ($db, $fasta, $tabular, $threads) = @_;
  `blastn -db $db -query $fasta -num_threads $threads -out $tabular -outfmt 6 -evalue 10e-10`;
}

############################################################
##Function extract_blast: extract result from blast      ###
############################################################
## @param:                                                 #
##       $file: name of sequences file                     #
## @return: hash that contains sequences                   #
############################################################


sub extract_blast
{
  my $file = shift;
  my %hash = ();
  open (my $f, $file) || die "Cannot open $file\n";
  while (<$f>)
  {
    chomp $_;
    my ($seq,$id) = split /\t/,$_;
    $seq = $1 if ($seq =~ /(\d+)-(.*?)-(\d+)-(\d+)/);
    $hash{$seq} = [] unless exists $hash{$seq};
    push @{$hash{$seq}}, $id;
  }
  close $f;
  return \%hash;
}
  
############################################################
##Function print_header: header of html file             ###
############################################################
## @param:                                                 #
############################################################

sub print_header
{
  my $fileR = shift; my $title = shift;
  print $fileR "<!DOCTYPE html>\n<html>\n<head>\n\t<title>$title</title>\n";
  print $fileR "\t<style type=\"text/css\">\n";
  print $fileR "\t\tbody { font-family:Arial, Helvetica, Sans-Serif; font-size:0.8em;}\n";
  print $fileR "\t\t#report { border-collapse:collapse;}\n";
  print $fileR "\t\t#report h4 { margin:0px; padding:0px;}\n";
  print $fileR "\t\t#report img { float:right;}\n";
  print $fileR "\t\t#report ul { margin:10px 0 10px 40px; padding:0px;}\n";
  print $fileR "\t\t#report th { background:#7CB8E2 url(static/images/header_bkg.png) repeat-x scroll center left; color:#fff; padding:7px 15px; text-align:left;}\n";
  print $fileR "\t\t#report td { background:#C7DDEE none repeat-x scroll center left; color:#000; padding:7px 15px; }\n";
  print $fileR "\t\t#report tr.odd td { background:#fff url(static/images/row_bkg.png) repeat-x scroll center left; cursor:pointer; }\n";
  print $fileR "\t\t#report div.arrow { background:transparent url(static/images/arrows.png) no-repeat scroll 0px -16px; width:16px; height:16px; display:block;}\n";
  print $fileR "\t\t#report div.up { background-position:0px 0px;}\n";
  print $fileR "\t</style>\n";
  print $fileR "\t<script src=\"./js/jquery.min.js\" type=\"text/javascript\"></script>\n";
  print $fileR "\t<script type=\"text/javascript\">\n";
  print $fileR "\t\t\$(document).ready(function(){\n";
  print $fileR "\t\t\t\$(\"#report tr:odd\").addClass(\"odd\");\n";
  print $fileR "\t\t\t\$(\"#report tr:not(.odd)\").hide();\n";
  print $fileR "\t\t\t\$(\"#report tr:first-child\").show();\n";
  print $fileR "\t\t\t\$(\"#report tr.odd\").click(function(){\n";
  print $fileR "\t\t\t\t\$(this).next(\"tr\").toggle();\n";
  print $fileR "\t\t\t\t\$(this).find(\".arrow\").toggleClass(\"up\");\n";
  print $fileR "\t\t\t});\n";
  print $fileR "\t\t\t//\$(\"#report\").jExpand();\n";
  print $fileR "\t\t});\n\t</script>\n";
  print $fileR "</head>\n<body>\n\t<table id=\"report\">\n";
}
  
############################################################
##Function html_tab: definition of html file             ###
############################################################
## @param:                                                 #
##       $fastq1_ref: reference to first paired end files  #
##       $name_ref: reference to names of reads files      #
##       $results_ref: reference to results files          #
##       $rna: results for known RNA                       #
##       $est: results for known EST                       #
##       $html: html results file                          #
##       $html_repertory: repository to store results      #
############################################################

sub html_tab
{
  my ($fastq1_ref, $name_ref, $results_ref, $rna, $est, $html, $html_repertory) = @_;
  my $out = $html_repertory;
  my @fastq1 = @{$fastq1_ref};
  my @name = @{$name_ref};
  my @results = @{$results_ref};
  
  # Copy HTML resources to results folder
  File::Copy::Recursive::dircopy "$Bin/js/", "$out/js" or die "Copy failed: $!";
  File::Copy::Recursive::dircopy "$Bin/static/", "$out/static" or die "Copy failed: $!";

  my $chimOut = $html;
  
  open(my $tab, ">".$chimOut) || die "Cannot open $chimOut";
  print_header($tab,"Chimerae");
  print $tab "\t\t<tr>\n\t\t\t<th>L1 chromosome</th>\n\t\t\t<th>L1 start</th>\n\t\t\t<th>L1 end</th>\n\t\t\t<th>L1 strand</th>\n";
  for my $i (0..$#fastq1)
  {
    print $tab "\t\t\t<th>$name[$i] read #</th>\n";
  }
  print $tab "\t\t\t<th>Chimera chromosome</th>\n\t\t\t<th>Chimera start</th>\n\t\t\t<th>Chimera end</th>\n\t\t\t<th>Chimera strand</th>\n";
  for my $i (0..$#fastq1)
  {
    print $tab "\t\t\t<th>$name[$i] read #</th>\n";
  }
  if(defined($rna))
  {
    print $tab "\t\t\t<th>Known RNA</th>\n";
  }
  else
  {
    print $tab "\t\t\t<th></th>\n";
  }
  if(defined($est))
  {
    print $tab "\t\t\t<th>Known EST</th>\n";
  }
  else
  {
    print $tab "\t\t\t<th></th>\n";
  }
  print $tab "\t\t\t<th></th>\n\t\t</tr>\n";
  
  for my $i (0..$#results)
  {
    print $tab "\t\t<tr>\n";
    foreach my $j (@{$results[$i]})
    {
      print $tab "\t\t\t<td>$j</td>\n";
    }
    my ($Hrna, $Hest) = ('','');
    $Hrna = ${$rna}{$i}[0] if exists(${$rna}{$i});
    $Hest = ${$est}{$i}[0] if exists(${$est}{$i});
    chomp $Hrna; chomp $Hest;
    if($Hrna)
    {
      print $tab "\t\t\t<td><a target=\"_blank\" rel=\"noopener noreferrer\" href=\"https://www.ncbi.nlm.nih.gov/nuccore/$Hrna\">$Hrna</a></td>\n";
    }
    else
    {
      print $tab "\t\t\t<td></td>\n";
    }
    if($Hest)
    {
      print $tab "\t\t\t<td><a target=\"_blank\" rel=\"noopener noreferrer\" href=\"https://www.ncbi.nlm.nih.gov/nuccore/$Hest\">$Hest</a></td>\n";
    }
    else
    {
      print $tab "\t\t\t<td></td>\n";
    }
    print $tab "\t\t\t<td><div class=\"arrow\"></div></td>\n\t\t</tr>\n";
    my $colspan = scalar(@fastq1) * 2 + 8 ;
    print $tab "\t\t<tr>\n\t\t\t<td valign=top colspan=$colspan></td>\n\t\t\t<td valign=top>\n";
    if (exists(${$rna}{$i}))
    {
      for (my $w = 1; $w <= $#{${$rna}{$i}}; $w++)
      {
        $Hrna = '';
        $Hrna = ${$rna}{$i}[$w];
        chomp $Hrna;
        print $tab "\t\t\t\t<a target=\"_blank\" rel=\"noopener noreferrer\" href=\"https://www.ncbi.nlm.nih.gov/nuccore/$Hrna\">$Hrna</a><br>\n";
      }
      delete ${$rna}{$i};
    }
    print $tab "\t\t\t</td>\n\t\t\t<td valign=top>\n";
    if (exists (${$est}{$i}))
    {
      for (my $w = 1; $w <= $#{${$est}{$i}}; $w++)
      {
        $Hest = '';
        $Hest = ${$est}{$i}[$w];
        chomp $Hest;
        print $tab "\t\t\t\t<a target=\"_blank\" rel=\"noopener noreferrer\" href=\"https://www.ncbi.nlm.nih.gov/nuccore/$Hest\">$Hest</a><br>\n";
      }
      delete ${$est}{$i};
    }
    print $tab "\t\t\t</td>\n\t\t\t<td></td>\n\t\t</tr>\n";
  }
  print $tab "\t</table>\n</body>\n</html>\n";
  close $tab;
}
  
############################################################
##Function save_csv: save results in different formats   ###
############################################################
## @param:                                                 #
##       $fastq1_ref: reference to first paired end files  #
##       $name_ref: reference to names of reads files      #
##       $results_ref: reference to results files          #
##       $line_only: Line only database                    #
##       $refseq: refseq text file                         #
##       $out: repository to store results                 #
############################################################
sub save_csv{
  my ($fastq1_ref, $name_ref, $results_ref, $line_only, $refseq, $out) = @_;
  my @fastq1 = @{$fastq1_ref};
  my @name = @{$name_ref};
  my @results = @{$results_ref};
  my $out1= $out.'/results.txt';
  my $out2= $out.'/first_results.txt';
  my $out3= $out.'/final_result_chimerae.txt';

  # save result in csv file ##
  
  my $filed = $out1;
  open(my $tab, ">".$filed) || die "Cannot open $filed";
  print $tab "L1 chromosome \t L1 start \t L1 end \t L1 strand";;
  for my $i (0..$#fastq1)
  {
    print $tab "\t $name[$i] read #";
  }
  print $tab "\t Chimera chromosome\t Chimera start \t Chimera end \t Chimera strand";
  for my $i (0..$#fastq1)
  {
    print $tab "\t $name[$i] read #";
  }
  print $tab "\n";
  for my $i ( 0 .. $#results )
  {
    my $rowref = $results[$i];
    my $n = @$rowref - 1;
    for my $j ( 0 .. $n-1 )
    {
      print $tab "$results[$i][$j]\t";
    }
    print $tab "$results[$i][$n]\n";
  }
  close $tab;
  
  ##Add some information via R Scripts##
  
  # Create bridge between Perl and R
  my $R = Statistics::R->new();
  $R->start();
  $R->set('out1', $out1);
  $R->set('out2', $out2);
  $R->set('out3', $out3);
  $R->set('nfastq', $#fastq1);
  $R->set('line_only', $line_only);
  $R->set('refseq', $refseq);
  my $R_out = $R->run_from_file("$Bin/CLIFinder_results.R");
  $R->stop();
  print STDOUT "$R_out\n";
}

__END__

=head1 NAME

        CLIFinder - Identification of L1 Chimeric Transcripts in RNA-seq data

=head1 SYNOPSIS

        CLIFinder.pl --first <first fastq of paired-end set 1> --name <name 1> --second <second fastq of paired-end set 1> [--first <first fastq of paired-end set 2> --name <name 2> --second <second fastq of paired-end set 2> ...] --ref <reference genome> [--build_ref] --TE <transposable elements> [--build_TE] --html <results.html> --html-path <results directory> [options]

        Arguments:
                --first <fastq file>    First fastq file to process from paired-end set
                --name <name>           Name of the content to process
                --second <fastq file>   Second fastq file to process from paired-end set
                --ref <reference>       Fasta file containing the reference genome
                --TE <TE>               Fasta file containing the transposable elements
                --rmsk <text file>      Tab-delimited text file (with headers) containing reference repeat sequences (e.g. rmsk track from UCSC)
                --refseq <text file>    Tab-delimited file (with headers) containing reference genes (e.g. RefGene.txt from UCSC)
                --html <html file>      Main HTML file where results will be displayed
                --html_path <dir>       Folder where results will be stored

        For any fasta file, if a bwa index is not provided, you should build it through the corresponding '--build_[element]' argument

        Options:
                --rnadb <RNA db>        Blast database containing RNA sequences (default: empty)
                --estdb <EST db>        Blast database containing EST sequences (default: empty)
                --size_read <INT>       Size of reads (default: 100)
                --BDir <0|1|2>          Orientation of reads (0: undirectional libraries, 1: TEs sequences in first read in pair, 2: TEs sequences in second read in pair) (default: 0)
                --size_insert <INT>     Maximum size of insert tolerated between R1 and R2 for alignment on the reference genome (default: 250)
                --min_L1 <INT>          Minimum number of bp matching for L1 mapping (default: 50)
                --mis_L1 <INT>          Maximum number of mismatches tolerated for L1 mapping (default: 2)
                --min_unique <INT>      Minimum number of consecutive bp not annotated by RepeatMasker (default: 33)
                --threads <INT>         Number of threads (default: 1)

        For Blast database files, if a fasta is provided, the database can be built with '--build_[db]'. Otherwise, provide a path or URL. \"tar(.gz)\" files are acceptable, as well as wild card (rna*) URLs.

