#!/usr/bin/perl

=head1 Name

pepMfa_to_cdsMfa.pl  --  convert protein alignment to cds alignment

=head1 Description

generate cds multiple alignment based on protein multiple alignment

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  perl pepMfa_to_cdsMfa.pl [option] <pep.mfa> <cds.fa>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/pepMfa_to_cdsMfa.pl ./leptin.cds.fa.pep.fa.muscle ../input/leptin.cds.fa > ./leptin.cds.fa.muscle

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);


my $aa_align_file = shift;
my $cds_file = shift;

my %aa;

##read the protein alignment information
open IN, $aa_align_file || die "fail open $aa_align_file\n";
$/=">"; <IN>; $/="\n";
while (<IN>) {
	my $name = $1 if(/^(\S+)/);
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
	$aa{$name} = $seq;
}
close IN;

##convert protein alignment to cds alignment
my $out;
open IN, $cds_file || die "fail open $cds_file\n";
$/=">"; <IN>; $/="\n";
while (<IN>) {
	$out .= ">".$_;
	my $name = $1 if(/^(\S+)/);
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\s//g;
	$seq =~ s/---//g;
	$seq = uc($seq);
	$/="\n";
	
	my $cds;
	my $prot = $aa{$name};
	my $len_prot = length($prot);
	my $j = 0;
        my (@prot,@CDS);
	for (my $i=0; $i<$len_prot; $i++) {
		my $aa = substr($prot,$i,1);
		if ($aa ne '-') {
                        my $coden=substr($seq,$j,3);
                        if($aa ne &translate($coden)){print STDERR "$aa != $coden at $name, use original seq\n"}
                        $cds .=$coden; 
                        push @prot," ".$aa." ";
                        push @CDS,$coden;
			$j += 3;
		}else{
			$cds .= '---';
                        push @prot," - ";
                        push @CDS,"---";
		}
	}
	Display_seq(\$cds);
	$out .= $cds;
        if(scalar @prot == scalar @CDS){print STDERR "@CDS","\n","@prot","\n"}
}
close IN;


print  $out;

	
#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################

sub translate(){
#standard code table:http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
#    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#  Starts = ---M---------------M---------------M----------------------------
#  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
   my $orf=shift;
   $orf=uc($orf);
   my %code=(
             'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S',
             'TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
             'TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*',
             'TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W',
             'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
             'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
             'CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q',
             'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
             'ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M',
             'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
             'AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
             'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R',
             'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
             'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
             'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E',
             'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G'
            );

   die"les or more code nucleotide\n" if(length $orf!=3);
   if($orf=~/N|Y|M/i){return "X"}else{return ($code{$orf})}
}
