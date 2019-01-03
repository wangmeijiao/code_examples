use warnings;
use strict;

my %align;
$/=">";
while(<stdin>){
  chomp;
  next if ($_ eq "");
  my @box=split/\n/,$_;
  my $head=shift @box;
  my $seq=join "",@box;
  if(!exists $align{$head}){$align{$head}=$seq}else{die "dup $head\n"}

}
$/="\n";

my @list=sort keys %align;
my (@first,@second,@first_degap,@second_degap);
my ($identity,$cnt_total,$cnt_similar);
if(length $align{$list[0]} == length $align{$list[1]}){
   @first=split//,$align{$list[0]};
   @second=split//,$align{$list[1]};
   for(0..$#first){
      if($first[$_] eq "-" || $second[$_] eq "-"){}else{
         push @first_degap,$first[$_]; 
         push @second_degap,$second[$_];
         $cnt_total++;
         if($first[$_] eq $second[$_]){$cnt_similar++}}
   }
}
$identity=100*$cnt_similar/$cnt_total;

my $identity_aa;
($cnt_total,$cnt_similar)=(0,0);
if(scalar @first_degap == scalar @second_degap && scalar(@second_degap)%3==0){
   for (my $i=0;$i<=$#first_degap-2;$i+=3){
      my $codon_first=$first_degap[$i].$first_degap[$i+1].$first_degap[$i+2] ;
      my $codon_second=$second_degap[$i].$second_degap[$i+1].$second_degap[$i+2] ;
      $cnt_total++;   
      if( &translate($codon_first) eq &translate($codon_second) ){ $cnt_similar++}
   }
}else{die "cds length diff or length % 3 != 0\n"}

$identity_aa=100*$cnt_similar/$cnt_total;


print "$list[0]&$list[1]$identity%$identity_aa%\n",join("",@first_degap),"\n",join("",@second_degap),"\n";


#####sub#######
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

