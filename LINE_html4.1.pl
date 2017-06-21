 #/usr/bin/perl

################################################
#Declaration of necessary libraries#############
################################################ 

use strict;
use warnings;
use Sys::CPU;
use Parallel::ForkManager;
use POSIX;
use Statistics::R ;
use Getopt::Long;
use File::Basename;

################################################
#Declaration of necessary global variables######
################################################

my (@fastq1, @fastq2, @name, $html, $size_reads, $ref, $TE, $build_index, $build_TE, $html_repertory, $maxInsertSize, $prct, $help, $image, $Bdir, $minL1, $mis_L1);

#####################################################################
#Definition options of execution according to the previous variables#
#####################################################################

GetOptions (
"first:s"  =>  \@fastq1,
"second:s"  =>  \@fastq2,
"name:s"  =>  \@name,
"html:s"  => \$html,
"TE:s" => \$TE,
"ref:s" => \$ref,
"build_TE" => \$build_TE,
"build_index" => \$build_index,
"pourcentage:i" => \$prct,
"size_insert:i" => \$maxInsertSize,
"size_read:i" => \$size_reads,
"html_path:s" => \$html_repertory,
"image:s" => \$image,
"BDir:i" => \$Bdir,
"minL1:i" => \$minL1,
"mis_L1:i" => \$mis_L1,
);
my $iprct = 100 - (($prct / $size_reads)*100) ;
my $mis_auth = $size_reads - $minL1 + $mis_L1 ;
my $eprct = ($iprct * $size_reads) /100; 
my $dprct = ((100-$iprct) * $size_reads) / 100;

################################################
#Construct index of ref and TE if doesn't exist#
################################################

`(bwa index $ref) 2> /dev/null ` if ($build_index);
`(bwa index $TE) 2> /dev/null ` if ($build_TE);

############################################
#Create repository to store resulting files#
############################################

mkdir $html_repertory;

######################################################
#Define  number of Cpu process we can use = total /5 #
######################################################

my $cpu = int(Sys::CPU::cpu_count() /5) ;

##########################################
#Define  hash                            #
##########################################

my %frag_exp_id;

##########################
#Data file we have to use#
##########################

my $NCBI_est = "https://galaxy.gred-clermont.fr/clifinder/est_human"; # NCBI Human est
my $NCBI_rna = "https://galaxy.gred-clermont.fr/clifinder/rna"; # NCBI Human rna
my $rmsk = "https://galaxy.gred-clermont.fr/clifinder/rmsk.bed"; # UCSC repeat sequences
####################################
# Redirection standart error output#
####################################

my $file = $html_repertory.'/report.txt';
open STDERR, '>', $file or die "Can't redirect STDERR: $!";

##############################
# Analyse of each fastq file #
##############################


my @garbage;  my $num = 0;
foreach my $tabR (0..$#fastq1)
{
  ###################################################
  # Paired end mapping against L1 promoter sequences#
  ###################################################
  
  print STDERR "alignement of $name[$tabR] to L1\n";
  my $sam = $html_repertory.'/'.$name[$tabR]."_L1.sam"; push(@garbage, $sam);
  align_paired( $TE, $fastq1[$tabR], $fastq2[$tabR], $sam, $cpu, $mis_auth);
  print STDERR "alignement done\n";
  
  ##################################################
  # Creation of two fastq for paired halfed mapped:#
  # - _1 correspond to sequences mapped to L1      #
  # - _2 correspond to sequences unmapped to L1    #
  ##################################################
  print STDERR "getting pairs with one mate matched to L1 and the other mate undetected by repeatmasker as a repeat sequence\n";
  
  my $out_ASP_1 = $html_repertory.'/'.$name[$tabR]."_1.fastq"; push(@garbage, $out_ASP_1);
  my $out_ASP_2 = $html_repertory.'/'.$name[$tabR]."_2.fastq"; push(@garbage, $out_ASP_2);
  
  ##split mate that matched to L1 and others##
  my ($ASP_readsHashR, $half_num_out) = get_half($sam); 
  # $ASP_reads{$line[0]}[0] mapped - $ASP_reads{$line[0]}[1] unmapped
  
  ##pairs obtained after repeatmasker on the other mate## 
  my $left = sort_out($cpu, $out_ASP_1, $out_ASP_2, $ASP_readsHashR);

  print STDERR "number of half mapped pairs : $half_num_out\n";
  print STDERR "number of pairs after repeatmasker: $left\n";
  
  ##################################################
  # Alignment of halfed mapped pairs on genome     #
  ##################################################
  print STDERR "alignement of potential chimeric sequences to the genome\n";
  $sam = $html_repertory.'/'.$name[$tabR]."_genome.sam";  
  push(@garbage, $sam);
  align_genome( $ref, $out_ASP_1, $out_ASP_2, $sam, $cpu);
  print STDERR "alignement done\n";
  
  ##compute the number of sequences obtained after alignement ##
  
  $left = `samtools view -Shc $sam 2> /dev/null`;
  chomp $left; $left = $left/2;
  print STDERR "number of sequences....: $left\n";

  ##################################################
  # Create bedfiles of potential chimerae          #
  # and Know repeat sequences removed              #
  ##################################################
  
  print STDERR "looking for chimerae\n";
  results($html_repertory, $sam, $name[$tabR], \%frag_exp_id, $num);
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
`bedtools sort -i $repfirst 2> /dev/null | bedtools merge -scores max -d 100 -s -nms -i /dev/stdin > $repMfirst  2> /dev/null`;
`bedtools sort -i $repsecond  2> /dev/null | bedtools merge  -scores max -d 100 -s -nms -i /dev/stdin > $repMsecond  2> /dev/null`;

my (%Gviz, %frag_uni, @second_R, @second_exp, @results);
my $merge_target = $html_repertory.'/target_merged.bed'; push(@garbage, $merge_target);
my $merge = $html_repertory.'/merged.bed'; push(@garbage, $merge);

open (my $mT, ">".$merge_target) || die "cannot open $merge_target\n";
open (my $m, ">".$merge) || die "cannot open $merge\n";
open (my $in, $repMsecond) || die "cannot open secondM\n";
my $cmp = 0;
while (<$in>)
{
  chomp $_;
  my @tmp = (0) x scalar(@fastq1);
  my @line = split /\t/, $_;
  my @names =split /;/, $line[3];
  foreach my $n (@names){$n =~/(.*?)\/[12]/; $frag_uni{$1} = $cmp; $tmp[$frag_exp_id{$1}]++; }
  $second_exp[$cmp] = \@tmp;
  $cmp++;
  push @second_R, [$line[0],$line[1],$line[2],$line[5]];
}

$cmp = 0;
open ($in, $repMfirst) || die "cannot open firstM\n";
while (<$in>)
{
  chomp $_;
  my %sec;
  my @line = split /\t/, $_;
  my @names =split /;/, $line[3];
  my @tmp = (0) x scalar(@fastq1);
  foreach my $n (@names){$n =~/(.*?)\/[12]/; $tmp[$frag_exp_id{$1}]++; }
  foreach my $n (@names)
  {
    $n =~/(.*?)\/[12]/;
    unless (exists ($sec{$frag_uni{$1}}) )
    {
      my @lmp = ($line[0], $line[1], $line[2], $line[5]);
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
#Blast against human rna and est               #
################################################

print STDERR "blast against human rna\n";
my $tabular = $html_repertory."/chimerae_rna.tab"; push(@garbage, $tabular);
blast($NCBI_rna, $fasta, $tabular, $cpu);
my $rna = extract_blast($tabular);

print STDERR "blast against human est\n";
my $tabular2 = $html_repertory."/chimerae_est.tab";push(@garbage, $tabular2);
blast($NCBI_est, $fasta, $tabular, $cpu);
my $est = extract_blast($tabular);

################################################
#Create Results html file                      #
################################################
print STDERR "save result in file\n";
save_csv();
print STDERR "create HTML\n";
html_tab($rna,$est,$html_repertory);
$extend = $extend.'*';
push(@garbage,glob($extend));
unlink @garbage;
print STDERR "Job done!\n";

########################################### END MAIN ##########################################################


##############################################################################################################
############################################     FUNCTIONS   #################################################   
##############################################################################################################



############################################################
##Function that aligned paired-end reads on a genome########
############################################################
## @param:												   #
##       $index: referential genome						   #
##       $fasq1: first paired end file					   #
##       $fasq2: second paired end file				       #
##       $sam:	Alignment output file 				       #
##       $number_of_cpus: Number of Cpu used			   #
############################################################

sub align_genome
{
  my ($index, $fastq1, $fastq2, $sam, $number_of_cpus) = @_ ;
  my @L_garbage =();
  my $sai1 = $sam."_temporary.sai1"; push @L_garbage,$sai1;
  my $sai2 = $sam."_temporary.sai2"; push @L_garbage,$sai2;
  `bwa aln -o4 -e1000  -t $number_of_cpus $index $fastq1 > $sai1 2> /dev/null`;
  `bwa aln -o4 -e1000 -t $number_of_cpus $index $fastq2 > $sai2 2> /dev/null`;
  ## -A force the insertion size
  `bwa sampe -s -A -a $maxInsertSize $index  $sai1 $sai2 $fastq1 $fastq2 2> /dev/null | samtools view -F4 -f 2  -Sh /dev/stdin > $sam`;
  unlink @L_garbage;
}


############################################################
##Function get_half get alignement on TE                 ###
############################################################
## @param:												   #
##       $sam: Name of alignement file   				   #
##														   #	
## @return:												   #
##		 $ASP_readsHashR: table to store sequences         #
##		 $half_num_out:	 number of alignment saved         #
############################################################

sub get_half
{
  ## store name of file	
  my $sam = shift; 
  open(my $fic, $sam) || die "cannot open sam file! $!\n"; ## Open file
  my (%ASP_reads); my $cmp = 0; ## Declare variables for
  my $sequence = '';
  
  ##read file##
  while(<$fic>)
  {
    chomp $_;
    ##We don't consider lines of sam files that are in header##
    next if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ );
    ##We split line in several part##
    my @line = split (/\t/,$_);
   
	##Find if alignemets have min_L1 consecutives bases mapped on R1 ##
	##Corriger enlever la condition mis_L1 = 0 ##   
#   if ($mis_L1 != 0)
#   {
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
          $accept = 1 if $tot >= $minL1;
        }
        ## if sequence is not accepted we go to the next sequence ##
        next if $accept == 0;
     }
#   }
    
    ##looking for flag of the alignment and keep only good reads##
    ##Find if it aligns on R1 or on R2##
    
    if ($line[1] == 73 || $line[1] == 89 || $line[1] == 117 || $line[1] == 69  || $line[1] == 133 || $line[1] == 181 || $line[1] == 153|| $line[1] == 137)
    {
      if ( $Bdir == 0 || ($Bdir == 1 && $line[1] & 64) || ($Bdir = 2 && $line[1] & 128))
      {
        $cmp++;
        $sequence = $line[9];
        ## if sequence is reversed aligned then reverse sequence ## 
        if ($line[1] & 16)
        {
          $sequence =reverse($sequence);
          $sequence =~ tr/atgcuATGCU/tacgaTACGA/;
        }
        ## define table contains ## 
        $ASP_reads{$line[0]} = [undef,undef] unless exists( $ASP_reads{$line[0]} );
        
        ##split if first mate (R1) is mapped on L1 or not (R2) ## 
        if ($line[1] & 8)
        {
          $ASP_reads{$line[0]}[0] = "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
        }
        else
        {
          $ASP_reads{$line[0]}[1] = "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
        }
      }
    }
  }
  close $fic;
  return ( \%ASP_reads, $cmp);
}

############################################################
##Function sort_out: extract paired end reads            ###
############################################################
## @param:												   #
##       $cpus:	number of Cpu used						   #
##       $out1:	output file accepted 1				       #
##       $out2:	output file accepted 2			           #
##       $readsHashTabR: reads to consider 				   #
############################################################

sub sort_out
{
  my ($cpus, $out1, $out2, $readsHashTabR) = @_;
  my ($name,$path) = fileparse($out2,'.fastq');
  my %repeat;
  my @garbage = (); my $cmp = 0;
  my $repout = $html_repertory.'/'.$name."_repout/";
  my $fa = $html_repertory.'/'.$name.".fa"; push (@garbage,$fa );
  my $second = $html_repertory.'/'.$name."_temporary.fastq";  push (@garbage,$second);
  mkdir $repout;
  my %notLine;
  
  ##Write on file containing of readssHashTabR
  
  open(my $tmp, ">".$second) || die "cannot open temp file $second\n";
  while ( my ($k,$v) = each %{$readsHashTabR} )
  {
    print $tmp ${$v}[1] if defined(${$v}[1]);
  }
  close $tmp;
  
  ## Transform fastq file to fasta
  
  `fastq_to_fasta -i $second -o $fa -Q33`;
  
  ##Launch RepeatMasker on fasta file
  
  `/data/oldgalaxy/Repeat/RepeatMasker/RepeatMasker -s -pa $cpus -dir $repout -species human $fa`;
  my $repfile = $repout.$name.".fa.out";
  open (my $rep, $repfile) || die "cannot open $repfile $!\n";
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
  
  ## store in table if pair passed the repeat test ##
  while (my ($k, $v) = each %repeat)
  {
    $notLine{$k} = 1 unless ($v->[0] > $dprct || $v->[1] < $eprct);
  }
  
  ##write resulting reads in  both files for paired ##
  open(my $accepted_1, ">".$out1 ) || die "cannot open $out1 file $!\n";
  open(my $accepted_2, ">".$out2 ) || die "cannot open $out2 file $!\n";
  while ( my ($k,$v) = each %{$readsHashTabR} )
  {
    if ( defined (${$v}[0]) && defined (${$v}[1]) )
    {
      unless (defined ($notLine{$k}) && $notLine{$k} == 1)
      {
        $cmp++;
        print $accepted_1 ${$v}[0];
        print $accepted_2 ${$v}[1];
      }
    }
  }
  close $accepted_1; close $accepted_2;
  
  ##drop files and directories generated by repeatmasker##
  my $toErase = $repout.'*';
  unlink glob "$toErase";
  unlink @garbage; rmdir $repout;
  return $cmp;
}

############################################################
##Function that aligned paired-end reads on a referential###
############################################################
## @param:												   #
##       $index: referential file						   #
##       $fasq1: first file paired end reads			   # 
##       $fasq2: second file paired end reads			   #
##       $sam: output alignment file				       #
##       $number_of_cpus: number of CPU used			   #
##		 $mis: tolerated mismatches					       #
############################################################
sub align_paired
{
  my ($index, $fastq1, $fastq2, $sam, $number_of_cpus, $mis) = @_ ;
  my @garbage = ();
  my $sai1 = $sam."_temporary.sai1"; push @garbage,$sai1;
  my $sai2 = $sam."_temporary.sai2"; push @garbage,$sai2;
  
  ##alignement with bwa
  
  `bwa aln -n $mis -t $number_of_cpus $index $fastq1 > $sai1 2> /dev/null`;
  `bwa aln -n $mis -t $number_of_cpus $index $fastq2 > $sai2 2> /dev/null`;
  `bwa sampe $index  $sai1 $sai2 $fastq1 $fastq2  > $sam 2> /dev/null`;
  
  ## delete temporary single aligned files
  unlink @garbage; 
}

############################################################
##Function that aligned reads on a referential           ###
############################################################
## @param:												   #
##       $index: referential file						   #
##       $fasq: reads file								   #
##       $sam:	output alignment file				       #
##       $number_of_cpus: number of CPU used			   #
############################################################

sub align
{
  my ($index, $fastq, $sam, $number_of_cpus ) = @_ ;
  `bwa aln -o4 -e$maxInsertSize -t $number_of_cpus $index $fastq  2> /dev/null | bwa samse $index /dev/stdin $fastq > $sam  2> /dev/null`;
}

############################################################
##Function results computes bed files for result         ###
############################################################
## @param:												   #
##       $out_repository: repository to store results	   #
##       $file:	sam file resulting of alignement		   #
##       $name:	name of paireds end reads file             #
##       $hashRef: store number for each first read value  #
##       $ps: number of the paired end file		           #
############################################################

sub results
{
  my ($out_repertory, $file, $name, $hashRef,$ps) = @_;
  my $namefirst = $out_repertory.'/'.$name.'-first.bed'; push(@garbage, $namefirst);
  my $namesecond = $out_repertory.'/'.$name.'-second.bed'; push(@garbage, $namesecond);
  
  ## store reads mapped in proper pair respectively  first and second in pair in bam files and transform in bed files## 
  `samtools view -Sb -f66 $file 2> /dev/null | bedtools bamtobed -i /dev/stdin > temp_name_first 2> /dev/null`;
  `samtools view -Sb -f130 $file 2> /dev/null | bedtools bamtobed -i /dev/stdin > temp_name_second 2> /dev/null`;
  
  ##compute converage of second mate on rmsk##
  my $baseCov = 0;
  my %IdCov = ();
  my @coverage = `bedtools coverage -b temp_name_second -a $rmsk`;
  
  
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
      IdCov{$1} = $split_cov[-1] if $split_cov[-1] > IdCov{$1};
    }
  }
  
  ## get only  first mate that have less tant $iprct repeats ##
  open (my $tmp_fi, 'temp_name_first') || die "cannot open $namefirst!\n";
  open (my $nam_fi, ">".$namefirst) || die "cannot open $namefirst!\n";
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

  open (my $tmp_sec, 'temp_name_second') || die "cannot open $namesecond!\n";
  open (my $nam_sec, ">".$namesecond) || die "cannot open $namesecond!\n";
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

#sub results
#{
#  my ($out_repertory, $file, $name, $hashRef,$ps) = @_;
#  my $namefirst = $out_repertory.'/'.$name.'-first.bed'; push(@garbage, $namefirst);
#  my $namesecond = $out_repertory.'/'.$name.'-second.bed'; push(@garbage, $namesecond);
#  `samtools view -Sb -f66 $file 2> /dev/null | bedtools bamtobed -i /dev/stdin > $namefirst 2> /dev/null`;
#  `samtools view -Sb -f130 $file 2> /dev/null | bedtools bamtobed -i /dev/stdin > $namesecond 2> /dev/null`;
#  open( my $in, $out_repertory.'/'.$name.'-first.bed') || die "cannot open first read bed\n";
#  while (<$in>)
#  {
#    my @line = split /\t/, $_;
#    $line[3] =~ /(.*?)\/1/;
#    ${$hashRef}{$1}= $ps;
#  }
#}


############################################################
##Function blast: blast nucleotide sequences on ref      ###
############################################################
## @param:												   #
##       $db:database where to search					   #
##       $fasta: file containing nucleotide sequences	   #
##       $tabular: out file name						   #
##       $cpu:	Number of Cpu used					       #
############################################################



sub blast
{
  my ($db, $fasta, $tabular, $cpus) = @_;
  `blastn -db $db -query $fasta -num_threads $cpus -out $tabular -outfmt 6 -evalue 10e-10`;
}

############################################################
##Function extract_blast: extract result from blast      ###
############################################################
## @param:  											   #
##       $file: Name of sequences file   				   #
## @return: hash that contains sequences                   #
############################################################


sub extract_blast
{
  my $file = shift;
  my %hash = ();
  open (my $f, $file) || die "cannot open $file\n";
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
## @param:												   #
############################################################



sub print_header
{
  my $fileR = shift; my $title = shift;
  print $fileR "<!DOCTYPE html> <html> <head> <title>$title</title>";
  print $fileR "<style type=\"text/css\">
  body { font-family:Arial, Helvetica, Sans-Serif; font-size:0.8em;}
  #report { border-collapse:collapse;}
  #report h4 { margin:0px; padding:0px;}
  #report img { float:right;}
  #report ul { margin:10px 0 10px 40px; padding:0px;}
  #report th { background:#7CB8E2 url(header_bkg.png) repeat-x scroll center left; color:#fff; padding:7px 15px; text-align:left;}
  #report td { background:#C7DDEE none repeat-x scroll center left; color:#000; padding:7px 15px; }
  #report tr.odd td { background:#fff url(row_bkg.png) repeat-x scroll center left; cursor:pointer; }
  #report div.arrow { background:transparent url(arrows.png) no-repeat scroll 0px -16px; width:16px; height:16px; display:block;}
  #report div.up { background-position:0px 0px;}
  </style>\n";
  print $fileR " <script src=\"./jquery.min.js\" type=\"text/javascript\"></script>\n";
  print $fileR "<script type=\"text/javascript\">
  \$(document).ready(function(){
    \$(\"#report tr:odd\").addClass(\"odd\");
    \$(\"#report tr:not(.odd)\").hide();
    \$(\"#report tr:first-child\").show();
    
    \$(\"#report tr.odd\").click(function(){
    \$(this).next(\"tr\").toggle();
    \$(this).find(\".arrow\").toggleClass(\"up\");
  });
  //\$(\"#report\").jExpand();
});
</script>";
print $fileR "</head> <body> <table id=\"report\" >\n";
}

############################################################
##Function html_tab: definition of html file             ###
############################################################
## @param:												   #	
##		 $html_repository: repository to store results	   #
##       $rna: results for know RNA						   #
##       $est: results for known EST				       #
############################################################

sub html_tab
{
  my ($rna,$est) = @_;
  my $out = $html_repertory;
  
  `cp https://galaxy.gred-clermont.fr/clifinder/arrows.png $out && cp https://galaxy.gred-clermont.fr/clifinder/row_bkg.png $out && cp https://galaxy.gred-clermont.fr/clifinder/jquery.min.js $out`;
  my $chimOut = $html;
  
  open(my $tab, ">".$chimOut) || die "cannot open $chimOut";
  print_header($tab,"Chimerae");
  print $tab "<tr>
  <th>L1 chromosome</th>
  <th>L1 start</th>
  <th>L1 end</th>
  <th>L1 strand</th>";
  for my $i (0..$#fastq1)
  {
    print $tab "\t<th>$name[$i] read #</th>\n";
  }
  print $tab "
  <th>Chimera chromosome</th>
  <th>Chimera start</th>
  <th>Chimera end</th>
  <th>Chimera strand</th>";
  for my $i (0..$#fastq1)
  {
    print $tab " <th>$name[$i] read #</th>\n";
  }
  print $tab "\t<th>Known RNA</th>
  \t<th>Known EST</th>\n\t<th></th>\n</tr>";  
  
  for my $i (0..$#results)
  {
    print $tab "<tr>";
    foreach my $j (@{$results[$i]})
    {
      print $tab " <td>$j</td>";
    }
    my ($Hrna, $Hest) = ('','');
    $Hrna = ${$rna}{$i}[0]  if exists(${$rna}{$i});
    $Hest = ${$est}{$i}[0] if exists(${$est}{$i});
    my $Lrna ='link break'; my $Lest = 'link break';
    chomp $Hrna; chomp $Hest;
    $Lrna = $3 if $Hrna =~/gi\|(.*?)\|(.*?)\|(.*)\|$/;
    $Lest = $3 if $Hest =~/gi\|(.*?)\|(.*?)\|(.*)\|$/;
    print $tab "\t<td><A HREF=\"http://www.ncbi.nlm.nih.gov/nuccore/$Lrna\">$Hrna</A></td>\n";
    print $tab "\t<td><A HREF=\"http://www.ncbi.nlm.nih.gov/nuccore/$Lest\">$Hest</A></td>\n";
    print $tab "\t<td><div class=\"arrow\"></div></td>\n</tr>\n";
    my $img = 'link break';
    $img = $i.'.svg';
    my $colspan = scalar(@fastq1) * 2 + 8 ;
    print $tab "<tr>\n\t<td valign=top  colspan = $colspan><img src=\"$img\"/></td>\n\t<td valign=top>";
    if (exists(${$rna}{$i}))
    {
      for (my $w = 1; $w <= $#{${$rna}{$i}}; $w++)
      {
        $Hrna = '';
        $Hrna = ${$rna}{$i}[$w];
        $Lrna ='link break';
        chomp $Hrna;
        $Lrna = $3 if $Hrna =~/gi\|(.*?)\|(.*?)\|(.*)\|$/;
        print $tab "<A HREF=\"http://www.ncbi.nlm.nih.gov/nuccore/$Lrna\">$Hrna</A><br>\n";
      }
      delete ${$rna}{$i};
    }
    print $tab "</td>\n\t<td valign=top>";
    if (exists (${$est}{$i}))
    {
      for (my $w = 1; $w <= $#{${$est}{$i}}; $w++)
      {
        $Hest = '';
        $Hest = ${$est}{$i}[$w];
        chomp $Hest;
        $Lest ='link break';
        $Lest = $3 if $Hest =~/gi\|(.*?)\|(.*?)\|(.*)\|$/;
        print $tab "\t<A HREF=\"http://www.ncbi.nlm.nih.gov/nuccore/$Lest\">$Hest</A><br>\n";
      }
      delete ${$est}{$i};
    }
    print $tab "</td>\n\t<td></td>\n</tr>\n";
  }
  print $tab qw{
    </table>
  };
  print $tab "<a href=\"report.txt\">Report</a>";
  print $tab qw{
    </body>
    </html>
  };
  close $tab;
}

############################################################
##Function save_csv: save results in different formats   ###
############################################################
## @param:												   #	
############################################################
sub save_csv{	
	my $Line_only="https://galaxy.gred-clermont.fr/clifinder/Line_only_hg19.txt.gz"; #Line Only H19 database
	my $Hg19_refseq="https://galaxy.gred-clermont.fr/clifinder/hg19_refseq.bed"; #h19 refseq bed file
	my $out1= $html_repertory.'/results.txt';
	my $out2= $html_repertory.'/first_results.txt';
	my $out3 =$html_repertory.'/final_result_chimerae.txt';
	
	# save result in csv file ##
	
    my $filed = $out1;
	open(my $tab, ">".$filed) || die "cannot open $file";
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
	
	use Statistics::R;
	# Create bridge between Perl and R
	my $R = Statistics::R->new();
	$R->startR;
	eval{
		$R->send(
		qq '
		rm(list=ls())
		library(GenomicRanges)
		library(plyr)
		chim<-read.delim("$out1")
		
		chim<-chim[order(chim[,$#fastq1+7],decreasing=F),]
		chim<-chim[order(chim[,2],decreasing=F),]
		chr<-sub("chr","",as.character(chim[,1]))
		chim<-chim[order(as.numeric(chr)),]
		
		grchim <- GRanges(seqnames = chim[,1],
		              IRanges(start = chim[,2],
		                      end = chim[,3]),strand=chim[,4])
		
		
		grchim\$ID<-paste("Id_",1:length(chim[,1]),sep="")
		
		mcols(grchim)<-cbind(mcols(grchim),chim[,5:($#fastq1+9)])
		
		
		grfusR<- union(grchim,grchim)
		
		position<-as.data.frame(findOverlaps(grfusR,grchim))
		
		grfusR\$dup<- table(position[,1])
		
		position2<-as.data.frame(findOverlaps(grfusR[grfusR\$dup>1],grchim))
		
		
		grfusR2<-grfusR
		strand(grfusR2)<- "+"
		gr3<-union(grfusR2,grfusR2)
		position3<- as.data.frame(findOverlaps(gr3,grfusR2))
		gr3\$dup<-table(position3[,1])
		
		position3<-as.data.frame(findOverlaps(gr3[gr3\$dup>1],grfusR2))
		grfusR\$info<-"no"
		grfusR\$info [position3[,2]]<-"overlap sens opposé"
		
		
				
		grfusR\$ID<-"Id"
		grfusR\$ID[position[!duplicated(position[,1]),1]]<-grchim\$ID[position[!duplicated(position[,1]),2]]
		
		if(nrow(position2)!=0)
		{
			result <- aggregate(position2[,2] ~ position2[,1], data = position2, paste, collapse = "_")
			grfusR\$ID[grfusR\$dup>1]<-paste("ID",result[,2],sep="_")
		}
		
		mcols(grfusR)<-cbind(mcols(grfusR), mcols(grchim[position[!duplicated(position[,1]),2]]))
		
		
		
		min<-ddply(as.data.frame(grchim), .(seqnames,end,strand), function(x)x[x\$Chimera.start==min(x\$Chimera.start),])
		min<-ddply(as.data.frame(min), .(seqnames,start,strand), function(x)x[x\$Chimera.start==min(x\$Chimera.start),])
		max<-ddply(as.data.frame(grchim), .(seqnames,end,strand), function(x)x[x\$Chimera.end==max(x\$Chimera.end),])
		max<-ddply(as.data.frame(max), .(seqnames,start,strand), function(x)x[x\$Chimera.end==max(x\$Chimera.end),])

		grfusR<-as.data.frame(grfusR)
		grfusR<-grfusR[order(grfusR[,1],grfusR[,2],grfusR[,3],grfusR[,4],decreasing=F),]

		grfusR\$Chimera.start<- min\$Chimera.start
		grfusR\$Chimera.end<-max\$Chimera.end

	
		datax<-as.data.frame(grfusR)
		colnames(datax)[1:3]<-colnames(chim)[1:3]
		colnames(datax)[5]<-colnames(chim)[4]
		
		
		grchim2 <- GRanges(seqnames = datax[,$#fastq1+12],
		              IRanges(start = datax[,$#fastq1+13],
		                      end = datax[,$#fastq1+14]),strand=datax[,$#fastq1+15])
		mcols(grchim2)<-datax[,-c(4,$#fastq1+12:$#fastq1+15)]
		
		
		
		grfus<- union(grchim2,grchim2)
		
		position<-as.data.frame(findOverlaps(grfus,grchim2))
		
		grfus\$dup<- table(position[,1])
		
		position2<-as.data.frame(findOverlaps(grfus[grfus\$dup>1],grchim2))
		
		grchim2[ position2[,2] ]
		
		
		
		
		grfus\$ID_final<-"Id"
		grfus\$ID_final[position[!duplicated(position[,1]),1]]<-grchim2\$ID[position[!duplicated(position[,1]),2]]
		
		if(nrow(position2)!=0)
		{
			result <- aggregate(position2[,2] ~ position2[,1], data = position2, paste, collapse = "_")
			grfus\$ID_final[grfus\$dup>1]<-paste("Id",result[,2],sep="_")
		}
		
		
		
		mcols(grfus)<-cbind(mcols(grfus), mcols(grchim2[position[!duplicated(position[,1]),2]]))
		
		 for (i in 0:$#fastq1)
		 {
		  mcols(grfus)[grfus\$dup>1,12+i] <- mcols(grfus)[grfus\$dup>1,12+i] +  mcols(grchim2)[position[duplicated(position[,1]),2],10+i]
		 }
		
		
		grfus2<-grfus
		strand(grfus2)<-"+"
		gr3<-union(grfus2,grfus2)
		
		position3<- as.data.frame(findOverlaps(gr3,grfus2))
		gr3\$dup<-table(position3[,1])
		
		position3<-as.data.frame(findOverlaps(gr3[gr3\$dup>1],grfus2))
		
		grfus\$info [position3[,2]]<-"overlap sens opposé"

		min<-ddply(as.data.frame(grchim2), .(seqnames,end,strand), function(x)x[x\$L1.start==min(x\$L1.start),])
		min<-ddply(data.frame(min), .(seqnames,start,strand), function(x)x[x\$L1.start==min(x\$L1.start),])
		max<-ddply(as.data.frame(grchim2), .(seqnames,end,strand), function(x)x[x\$L1.end==max(x\$L1.end),])
		max<-ddply(data.frame(max), .(seqnames,start,strand), function(x)x[x\$L1.end==max(x\$L1.end),])
		
		grfus1<-as.data.frame(grfus)
		grfus1<-grfus1[order(grfus1[,1],grfus1[,2],grfus1[,3],grfus1[,4],decreasing=F),]
		
		grfus1\$L1.start<- min\$L1.start
		grfus1\$L1.end<-max\$L1.end
		
		dataf<-as.data.frame(grfus1)		
		
		result<-( data.frame("Chimera.Chr"= dataf\$L1.chromosome, "Chimera.Start"=apply(data.frame( dataf\$start,dataf\$end,dataf\$L1.start,dataf\$L1.end  ),1,min) , "Chimera.End"= apply(data.frame( dataf\$start,dataf\$end,dataf\$L1.start,dataf\$L1.end  ),1,max) ,"Chimera.Strand"=dataf\$L1.strand ,"L1.Chr"= dataf\$L1.chromosome, "L1.Start"=dataf\$L1.start ,"L1.End"= dataf\$L1.end , "L1.Strand"=dataf\$L1.strand , "Unique.Chr"= dataf\$seqnames, "Unique.Start"=dataf\$start , "Unique.End"= dataf\$end , "Unique.Strand"=dataf\$strand , "ID_final"=dataf\$ID_final,"info"=dataf\$info, dataf[,18:($#fastq1+18)]   )  )
		
		result<-result[order(result[,2],decreasing=F),]
		chr<-sub("chr","",as.character(result[,1]))
		result<-result[order(as.numeric(chr)),]
		options(scipen=10) 
		write.table(result,"$out2",sep="\t",row.names = F,quote = F) 
	    grchim <- GRanges(seqnames = result\$L1.Chr,
	              IRanges(start = result\$L1.Start,
	              		  end = result\$L1.End),strand=result\$L1.Strand)
	    mcols(grchim)<-result
	    
		Rep<-read.delim("/data/oldgalaxy/database/Human/Line_only_hg19.txt.gz",skip=1)
		
		Gene<-read.delim("/data/oldgalaxy/database/Human/hg19_refseq.bed") 
	    grLINE <- GRanges(seqnames = Rep\$genoName,
	              IRanges(start = Rep\$genoStart,
	                      end = Rep\$genoEnd),
	              repStrand=as.character(Rep\$strand),
	                     repName =as.character(Rep\$repName))
	   
	    grGene <- GRanges(seqnames = Gene\$chrom,
	                      IRanges(start = Gene\$txStart,
	                              end = Gene\$txEnd),
	                      geneStrand=as.character(Gene\$strand),
	                              geneName = as.character(Gene\$name2))
	                              
		position<-as.data.frame(findOverlaps(grchim,grLINE))
	
		position2<-as.data.frame(findOverlaps(grchim,grGene))
	
		table(grLINE\$repName[position[,2]])
		write.table(table(grLINE\$repName[position[,2]]))
		
		grchim\$GeneName<-"no_gene"
	
		grchim\$GeneName[position2[,1]]<- grGene\$geneName[position2[,2]]
		
		grchim\$GeneStrand<-"*"
		
		grchim\$GeneStrand[position2[,1]]<- grLINE\$repStrand[position2[,2]]
		
		grchim\$repName<-"no"
		
		grchim\$repName[position[,1]]<- grLINE\$repName[position[,2]]
		
		grchim\$repStart<-0
	
		grchim\$repStart[position[,1]]<-start(grLINE[position[,2]])
		
		grchim\$repEnd<-0
		
		grchim\$repEnd[position[,1]]<-end(grLINE[position[,2]])
		
		grchim\$repWidth<-0
		
		grchim\$repWidth[position[,1]]<-width(grLINE[position[,2]])
		
		grchim\$repStrand<-"*"
		
		grchim\$repStrand[position[,1]]<- grLINE\$repStrand[position[,2]]
		
		dup<-position[duplicated(position[,1]),1]
		if(length(dup != 0))
		{
			for (i in 1:length(dup))
			{
			  
				grchim\$repName[dup[i]] <-paste(grLINE\$repName[position[position[,1]==dup[i],2]],collapse="/")
				grchim\$repStart[dup[i]] <-paste(start(grLINE[position[position[,1]==dup[i],2]]),collapse="/")
				grchim\$repEnd[dup[i]] <-paste(end(grLINE[position[position[,1]==dup[i],2]]),collapse="/")
				grchim\$repWidth[dup[i]] <-paste(width(grLINE[position[position[,1]==dup[i],2]]),collapse="/")
				grchim\$repStrand[dup[i]] <-paste(grLINE\$repStrand[position[position[,1]==dup[i],2]],collapse="/")	
				
			}
		}
		
		final_result<-as.data.frame(grchim)
		options(scipen=10)
		write.table(final_result[,-c(1:5)],"$out3",sep="\t",row.names = F,quote = F)
		'
		);
		};
	
		$R->stop();
		
		
		
		
			
}	


