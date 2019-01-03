use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#input: gff+fa
#output:  -gene 
#         -mrna(longest mrna, exons joined) 
#         -cds (CDS of longest mrna, CDSs joined)
#         -pep (translated pep of longest mrna CDSs, considering frameshift and check premature stop codon)
# 
#        gene(total region) -> mrna(exons) -> cds(CDSs) -> pep(cds_segs2pep joined and checked)


#usage: perl gffExtractAny_v1.0.pl -gff Sbicolor_255_v2.1.gene_exons.gff3 -fa ~/data/all_genomes/sorgJGI2.0/assembly/Sbicolor_255_v2.0.softmasked.fa -prefix sorg2 -gene -mrna -cds -pep -outdir sorg2_out -width 100 

#note: modify gff_read part first for compatiblity /ken'paetibl/
#
#

my($gene,$mrna,$cds,$pep,$prefix,$outdir,$gff,$fa,$width);
GetOptions("gene!",\$gene,"mrna!",\$mrna,"cds!",\$cds,"pep!",\$pep,"gff=s",\$gff,"fa=s",\$fa,"prefix:s",\$prefix,"outdir:s",\$outdir,"width:i",\$width);
$width||=50;

my %gff;#main structure table:  chr->gene->mRNAs->(elements:CDS,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)
my %gene;#gene info
my %mrna;#mrna info
#1,read&fill without order consideration, line by line
open GFF, $gff or die "$!";
while(<GFF>){
  chomp;
  next if($_ eq "" || $_=~/^#/);
  my ($chr,undef,$ftype,$start,$end,undef,$strand,$frameshift,$attr)=split/\t/,$_;
  die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
  next if($chr eq "chrUn" || $chr eq "chrSy"); #get rid of other chrs
  next if($ftype eq "chromosome" || $ftype eq "Chromosome");
  ($start,$end)=($start<$end)?($start,$end):($end,$start);
  my ($ID,$pID,$name,$note,$longest);
  if($ftype eq "gene"){
    $attr=~/ID=([^;]+);Name=([^;]+)/;
    ($ID,$name)=($1,$2);
    #$alias=$ID;
    #$name=~s/%20|%2C|%2D|%2E|%2F|%5F|%3A|%3B|%2A|%2B/_/g;
    if(!exists $gff{$chr}->{$ID}){
      my %temp;
      $gff{$chr}->{$ID}=\%temp;
      $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
    }else{die "dup gene name $ID\n"}
   }elsif($ftype eq "mRNA"){
      $attr=~/ID=([^;]+);Name=([^;]+);pacid=([^;]+);longest=([^;]+);Parent=([^;]+)/;
      ($ID,$name,$longest,$pID)=($1,$2,$4,$5);
      #$alias=$ID;
      if(!exists $gff{$chr}->{$pID}->{$ID}){
        my %temp;
        $gff{$chr}->{$pID}->{$ID}=\%temp;
        $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID\t$longest";
      }else{die "dup mRNA ID at $ID\n"}
     }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"){
         $attr=~/ID=([^;]+);Parent=([^;]+);pacid=([^;]+)/;
         $pID=$2;
         $name=$1;
         my $geneID;
         if(exists $mrna{$pID}){
          my @box=split/\t/,$mrna{$pID};
          if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
         }else{die "elements comes first, I can't assign it to gene\n"}
         if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
            my @temp;
            if($ftype eq "CDS"){push @temp,"$start-$end-$frameshift"}else{push @temp,"$start-$end"}
            $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
         }else{($ftype eq "CDS")?(push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$frameshift"):(push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end")}
       }else{die "unknow feature type\n"}
}#while end
close GFF;
print STDERR "read gff done\n";

#print Dumper \%gff;
#print Dumper \%gene;
#print Dumper \%mrna;
#exit;

#2, read in fa

my $fa_index=&faRead($fa);
print STDERR "read fa done\n";

#3, start to extract and output
$outdir||=".";
if(-d $outdir){}else{mkdir $outdir}
$prefix||="out";
print STDERR "start to extract from fa: $fa according to gff $gff :\n";
if($gene){open GENE, ">$outdir/${prefix}_genes.fa" or die "$!"}
if($mrna){open MRNA, ">$outdir/${prefix}_mrnas.longest.fa" or die "$!";}
if($cds){open CDS, ">$outdir/${prefix}_CDSs.longest.fa" or die "$!";}
if($pep){open PEP, ">$outdir/${prefix}_protein.longest.fa" or die "$!"; open PEP_LOST, ">$outdir/${prefix}_protein.lost.fa" or die "$!"}
foreach my$chr(sort keys %gff){
 print STDERR "$chr ..";
 foreach my $gene(sort keys %{$gff{$chr}}){
    #1,output gene region
    my ($chr,$start,$end,$ID,$strand,$note)=split/\t/,$gene{$gene};
    if($gene){
     my $seq = &faFormat(&faExtract($fa_index,"$chr:$start-$end:$strand"),$width);
     print GENE ">$ID\n$seq";
    }
    #my $max_len=0;
    my $max_len_mrna;
    foreach my $mrna(keys %{$gff{$chr}->{$gene}}){
      my @box=split/\t/,$mrna{$mrna};
      #if($max_len<($box[2]-$box[1])){$max_len=$box[2]-$box[1];$max_len_mrna=$mrna;} 
      if($box[6] == 1){$max_len_mrna=$mrna}
    }
    #mrna (the longest one, which all of the following caculation are based on)
    my $index=$gff{$chr}->{$gene}->{$max_len_mrna};
    if(!defined $index || $index eq ""){print STDERR "the longest mrna of $gene is lost:$max_len_mrna\n"; next}
    my ($chr_mrna,$start_mrna,$end_mrna,$id_mrna,$strand_mrna,$pID)=split/\t/,$mrna{$max_len_mrna};
    $id_mrna=~s/\.v2\.1$//;
    #2,output mrna region longest
    #if($mrna){
    # my $seq = &faFormat(&faExtract($fa_index,"$chr_mrna:$start_mrna-$end_mrna:$strand_mrna"),50);
    # print MRNA ">$alias_mrna\n$seq";
    #}
    #3,output mrna transcritps (exons with utrs)
    if(exists $index->{"exon"}){
      my @exon=@{&segSort($index->{"exon"})};
      my @exon_seq;
      foreach(@exon){
        my($exon_start,$exon_end)=split/-/,$_;
        my $exon_seq=&faExtract($fa_index,"$chr_mrna:$exon_start-$exon_end:$strand_mrna");
        ($strand_mrna eq "+")?(push @exon_seq, $exon_seq):(unshift @exon_seq, $exon_seq);
      }
      if($mrna){
         my $seq = &faFormat(join("",@exon_seq),$width);
         print MRNA ">$id_mrna\n$seq";
      }
    }else{print STDERR "no exon records for $id_mrna, mrna transcripts output nothing\n"}
    #4,output CDS and/or pep
    if(exists $index->{"CDS"}){
       my @CDS=@{&segSort($index->{"CDS"})};
       #my $i;
       my @CDS_seq;
       my @PEP_seq;
       foreach(@CDS){
         #$i++;
         my($cds_start,$cds_end,$frameshift)=split/-/,$_;
         my $cds_seq=&faExtract($fa_index,"$chr_mrna:$cds_start-$cds_end:$strand_mrna");
         ($strand_mrna eq "+")?(push @CDS_seq, $cds_seq):(unshift @CDS_seq, $cds_seq);
         my $pep_seq=&translate_CDS($cds_seq,$frameshift);
         ($strand_mrna eq "+")?(push @PEP_seq, $pep_seq):(unshift @PEP_seq, $pep_seq);
       }
       if($cds){
         my $seq=&faFormat(join("",@CDS_seq),$width);
         print CDS ">$id_mrna\n$seq";
       }
       if($pep){
         #check PEP validity
         my $pep=join("",@PEP_seq);
         my $flag=&checkPEP($pep);
         if($flag eq "OK"){ 
           my $seq=&faFormat($pep,$width);
           print PEP ">$id_mrna\n$seq";
         }else{print PEP_LOST "$id_mrna PEP illegal: $pep, $flag\n"}
       }
    }else{die "no CDS at $ID\n"}
 }#foreach gene
  print STDERR " done\n";
}#foreach chr end
if($gene){close GENE}
if($mrna){close MRNA}
if($cds){close CDS}
if($pep){close PEP; close PEP_LOST}


###subs###

sub faRead(){
  my $file = shift;
  my %fa;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my $seq = join "", @box;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return \%fa;
}

sub segSort(){
 my @segs=@{shift @_};
 my @sorted=sort{
      my ($s1,$e1,$s2,$e2);#$s1<=$e1 && $s2<=$e2 and deOverlaped beforehand by default
      ($s1,$e1)=split/-/,$a;
      ($s2,$e2)=split/-/,$b;
      if($e2<=$s1){1}elsif($e1<=$s2){-1}else{print STDERR "unknown circumstance at $s1,$e1:$s2,$e2\n"}
    } @segs;
# print "@CDS\n@sorted\n";
 return \@sorted;
}


sub faExtract(){
    my ($index,$coord)=@_;
    my ($chr,$s,$e,$strand)=split/:|-/,$coord;
    die "not exists $chr in your fa file\n"if(!exists $index->{$chr});
    my $seq_extract=($strand eq "+")?(substr($index->{$chr},$s-1,$e-$s+1)):(&RC(substr($index->{$chr},$s-1,$e-$s+1)));
    #return uc($seq_extract);
    return $seq_extract;
}

sub faFormat(){
  my ($seq,$num)=@_;
  my @seq=split//,$seq;
  my ($string,$cnt);
  foreach(@seq){
   $cnt++;
   if($cnt%$num==0){$string.=$_."\n"}else{$string.=$_}
  }
  if($cnt%$num!=0){return $string."\n"}else{return $string}
}

sub RC(){
  my $seq=shift;
  $seq=reverse $seq;
  $seq=~tr/atcgATCG/tagcTAGC/;
  return $seq;
}


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
   my $aa="";
   my $codon;
   if($orf=~/^ATG/ && $orf=~/(TAA|TAG|TGA)$/ && (length $orf)%3==0){
       for(my$i=1;$i<=((length $orf)-2);$i+=3){
        $codon=substr($orf,$i-1,3);
        if($i!=(length $orf)-2 && $codon=~/TAA|TAG|TGA/){return "##earlier_stop"}
        if(exists $code{$codon}){$aa.=$code{$codon}}else{return "##wrong_codon:$codon"}
       }
       return $aa;
     }else{return "##truncted_orf"}
}



sub translate_CDS(){
#standard code table:http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
#    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#  Starts = ---M---------------M---------------M----------------------------
#  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
   my $orf=shift;
   my $frameshift=shift; #0,1,2: delete 0,1,2 nucleiotide before translating
   if($frameshift eq "." || $frameshift eq "" || !defined $frameshift){$frameshift=0}
   $orf=uc($orf);
   $orf=substr($orf,$frameshift);
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
   my $aa="";
   my $codon;
   for(my$i=1;$i<=((length $orf)-2);$i+=3){
     $codon=substr($orf,$i-1,3);
     if($i!=(length $orf)-2 && $codon=~/TAA|TAG|TGA/){return "##earlier_stop"}
     if(exists $code{$codon}){$aa.=$code{$codon}}else{return "##wrong_codon"}
   }
   return $aa;
}

sub checkPEP(){
   my $pep=shift;
   my @comments;
   if($pep!~/^M/){push @comments,"##start codon not M"}
   if($pep!~/\*$/){push @comments,"##stop codon not exists"}
   if(length $pep <= 10){push @comments,"##short than 50 aa"}
   chop($pep);
   if($pep=~/\*/){push @comments,"##stop codon in pep"}
   if(scalar @comments == 0){ return "OK"}else{return join",",@comments}

}


