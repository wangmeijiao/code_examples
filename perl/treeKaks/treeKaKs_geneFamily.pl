use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#1, extract pep and cds 
#2, align peps with muscle
#3, transform pep_alignment to CDS_alignment
#4, build phylogenetic trees with treebest and draw tree svg
#5, calculate Ka/Ks with KaKs_calculator

#usage: cat pairs.list.cluster.homoTab | perl treeKaKs.pl -cds cds.fa -pep pep.fa -outdir outdir2 2> err > log &


my($pep,$cds,$outdir,$tmpdir);
GetOptions("pep=s",\$pep,"cds=s",\$cds,"outdir:s",\$outdir,"tmp:s",\$tmpdir);
$outdir||="out";
#$tmpdir||="tmp";
my $pep_index=&readFa($pep);
my $cds_index=&readFa($cds);
print STDERR "read in fa done \n";

#print Dumper $pep_index;
#print STDERR Dumper $cds_index;

#read in cluster tab
my $cnt=0;
while(<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_ =~/^$/);
  $cnt++;
  my @ids=split/[ \t]+/,$_;
#1, extract prot and cds 
  my @seqs_prot;
  my @seqs_cds;
  foreach my $id(@ids){
       if($id ne "-"){
         if(exists $pep_index->{$id}){ push @seqs_prot, ">$id\n$pep_index->{$id}\n" }else{die "not exists $id in  pep file $pep\n"} 
         if(exists $cds_index->{$id}){ push @seqs_cds, ">$id\n$cds_index->{$id}\n" }else{die "not exists $id in cds file $cds\n"} 
       }
  }
  print "\n";
  my $prot=join"",@seqs_prot;
  my $cds=join"",@seqs_cds;
  #print "$prot";
  #print "$cds";

#2, align aa by muscle and check quality
my $aln=`echo -ne "$prot"|muscle -quiet `;
#`clustalw2 -INFILE=temp.fa -TYPE=DNA -OUTFILE=temp.alnFa -OUTPUT=FASTA `;
#print $aln;


#stat pepAlign 
my ($coverage,$identity)=split/[ |\t]+/,`echo -ne "$aln" |perl bin/mfaStats_v1.2.pl -quiet`;
#print "$coverage,$identity\n";


#3, transform prot_alignment to cds_aligment, remove aligment gaps and transform to axt
#`perl ../bin/pepMfa_to_cdsMfa.pl temp_prot.muscle temp_cds.fa 2>temp_prot2cds.codeDetail|perl ../bin/degap2Axt.pl > temp_prot2cds.degap.axt ` ;
#if(! -d $tmpdir ){mkdir $tmpdir};
if(! -d $outdir ){if(! mkdir $outdir){die "can't mkdir $outdir: dir or file already exists\n"}};
open PEPALN, ">$outdir/orth${cnt}.pepAlign.mfa" or die "$!"; 
open ALN, ">$outdir/orth${cnt}.cdsAlign.mfa" or die "$!"; 
open ALN_degap, ">$outdir/orth${cnt}.cdsAlign.degap.mfa" or die "$!"; # not used in the following steps 

print PEPALN "$aln\n";
print ALN &pepAlign2cdsAlign($aln,$cds);
print ALN_degap &pepAlign2cdsAlign_degap($aln,$cds);

close PEPALN;
close ALN;
close ALN_degap;

#4, build tree with treebest and draw tree svg file
#`treebest nj $outdir/orth${cnt}.cdsAlign.mfa > $outdir/orth${cnt}.cdsAlign.mfa.nhx 2> /dev/null`;
`treebest filter -n $outdir/orth${cnt}.cdsAlign.mfa > $outdir/orth${cnt}.cdsAlign.filter.mfa`;
`treebest nj $outdir/orth${cnt}.cdsAlign.filter.mfa > $outdir/orth${cnt}.cdsAlign.filter.mfa.nhx 2>/dev/null`;
`nw_display -s -w 1000 -v 100 -S -b 'opacity:0' -i 'font-size:10' -d 'stroke-width:2;stroke:blue' $outdir/orth${cnt}.cdsAlign.filter.mfa.nhx > $outdir/orth${cnt}.cdsAlign.filter.mfa.nhx.svg`;


#5, calculate kaks by KaKs_Calculator and check Ka/Ks
#`cat $outdir/orth${cnt}.cdsAlign.mfa |perl bin/mfa2axt.pl > $outdir/orth${cnt}.cdsAlign.mfa.axt`;
#`KaKs_Calculator -i $outdir/orth${cnt}.cdsAlign.mfa.axt -o $outdir/orth${cnt}.cdsAlign.mfa.axt.kaks -m NG `;
`cat $outdir/orth${cnt}.cdsAlign.filter.mfa |perl bin/mfa2axt.pl -degap > $outdir/orth${cnt}.cdsAlign.filter.mfa.axt`;
`KaKs_Calculator -i $outdir/orth${cnt}.cdsAlign.filter.mfa.axt -o $outdir/orth${cnt}.cdsAlign.filter.mfa.axt.kaks -m NG `;
my @kaks=`awk -vOFS='\t' '(NR>1){print \$1,\$5}' $outdir/orth${cnt}.cdsAlign.filter.mfa.axt.kaks`;


#6, combine pepalign quality and kaks values to report into a single file 
print join ("|",@ids),"\n","$coverage,$identity",@kaks,"\n";

=pod
my @out=`awk -vOFS='\t' '(NR!=1){print \$3,\$4,\$5,\$2,\$1}' temp_prot2cds.degap.axt.kaks`;
print "$seqQ\t$seqS\n@out";
chomp @out;
my @data;
foreach(@out){
  my ($ka,$ks,$kas,$method,undef)=split/\t/,$_;
  push @data, $ks;

}
my $sd=&SD(\@data);
print "SD: $sd\n"

#my $ks=`awk -vOFS='\t' '(NR==2){print \$2,\$3,\$4,\$5,\$1}' temp_prot2cds.degap.axt.kaks`;
#chomp $ks;
#if(defined $ks && $ks ne "" && $ks<0.7){print "$seqQ\t$seqS\t$ks\n"}else{print STDERR "$seqQ\t$seqS\t$ks\n"}
=cut
#  if($cnt % 1000 ==0){
#     print STDERR "#";
#     `mkdir $outdir/index$cnt`;    
#     `mv $outdir/orth*  $outdir/index$cnt`;
#  }
  #last if($cnt>5);
}#read stdin while end here
print STDERR "\ntotal $cnt orthologous groups done\n";

#####sub########

sub readFa(){
  my $file=shift;
  my %fa;
  $/=">";
  open FA, $file or die "$!";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my $seq=join"",@box;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id in file $file\n"}
  }
  close FA;
  $/="\n";
  return \%fa;
}


sub pepAlign2cdsAlign(){
    my ($pepAlign,$cds)=@_;
    my (%pep,%cds);
    my @cdsAlign;
    my @pepAlign=split/>/,$pepAlign;
    foreach (@pepAlign){
      next if($_ eq "" || ! defined $_);
      my @box=split/\n+/,$_;
      my $id=shift @box;
      my $seq=join"",@box;
      if(!exists $pep{$id}){$pep{$id}=$seq}else{die "dup $id in pepAlign $pepAlign\n"};
    }
    my @cds=split/>/,$cds; 
    foreach (@cds){
      next if($_ eq "" || ! defined $_);
      my @box=split/\n+/,$_;
      my $id=shift @box;
      my $seq=join"",@box;
      if(!exists $cds{$id}){$cds{$id}=$seq}else{die "dup $id in cds $cds\n"};
    }
    foreach my $id(sort keys %pep){ # pepAlign -> cdsAlign
     my $cdsAlign=""; 
     if(exists $cds{$id}){
        my $pepAlign=$pep{$id};
        my $cds=$cds{$id};
        my $j=0;
        #my (@prot,@CDS);
        for (my $i=0; $i< length $pepAlign; $i++) {
                my $aa = substr($pepAlign,$i,1);
                if ($aa ne '-') {
                        my $codon=substr($cds,$j,3);
                        if($aa ne &translate($codon)){print STDERR "$aa != $codon at $id, use original seq\n"}
                        $cdsAlign.=$codon;
                        #push @prot," ".$aa." ";
                        #push @CDS,$coden;
                        $j += 3;
                }else{
                        $cdsAlign.= '---';
                        #push @prot," - ";
                        #push @CDS,"---";
                     }    
        }#foreach pepAlign
      }else{die "pep $id not exists in cds\n"}
      my $string=&faFormat($cdsAlign,50);
      push @cdsAlign,">$id\n$string";   
   }#pep id foreach
   my $result=join("",@cdsAlign);
   return $result;
}#sub end


sub pepAlign2cdsAlign_degap(){
    my ($pepAlign,$cds)=@_;
    my (%pep,%cds);
    my @order_ids;
    my %count; #use for aa -> cds count coord
    my %cdsAlign; #use for aa -> cds container
    my @pepAlign=split/>/,$pepAlign;
    foreach (@pepAlign){
      next if($_ eq "" || ! defined $_ || $_=~/^\s+$/);
      my @box=split/\n+/,$_;
      my $id=shift @box;
      my $seq=join"",@box;
      my @seq=split//,$seq;
      if(!exists $pep{$id}){$pep{$id}=\@seq; $count{$id}=0; $cdsAlign{$id}=""}else{die "dup $id in pepAlign $pepAlign\n"};
      push @order_ids,$id;
    }
    my @cds=split/>/,$cds; 
    foreach (@cds){
      next if($_ eq "" || ! defined $_ || $_=~/^\s+$/);
      my @box=split/\n+/,$_;
      my $id=shift @box;
      my $seq=join"",@box;
      if(!exists $cds{$id}){$cds{$id}=$seq}else{die "dup $id in cds $cds\n"};
    }

    my ($flag,$len)=&checkAlignLenAll(\%pep);
    die "pepAlign not equal: $pepAlign" if($flag == 0);

    for my $i(0..$len-1){ # pepAlign -> cdsAlign, process all alignen seqs at the same time, with degap

      my $flag_gap=0; #has gap or not ?
      foreach my $id(@order_ids){
        if(exists $cds{$id}){
          my $aa=$pep{$id}->[$i];
          if($aa eq "-"){$flag_gap=1;last}
        }else{die"pep $id not exists in cds\n"}
      }

      foreach my $id(@order_ids){
          my $aa=$pep{$id}->[$i];
          if($aa ne "-"){
            my $codon=substr($cds{$id},$count{$id},3);
            if($aa ne &translate($codon)){print STDERR "$aa != $codon at $id, use original seq\n"}
            if($flag_gap==0){$cdsAlign{$id}.=$codon}
            $count{$id}+=3;                        
          }
      }
    }#for 0-len end 
 
    my @cdsAlign;
    foreach my $id(@order_ids){
      my $string=&faFormat($cdsAlign{$id},50);
      push @cdsAlign,">$id\n$string";
    }

=pod
    foreach my $id(sort keys %pep){ # pepAlign -> cdsAlign
     my $cdsAlign=""; 
     if(exists $cds{$id}){
        my $pepAlign=$pep{$id};
        my $cds=$cds{$id};
        my $j=0;
        #my (@prot,@CDS);
        for (my $i=0; $i< length $pepAlign; $i++) {
                my $aa = substr($pepAlign,$i,1);
                if ($aa ne '-') {
                        my $codon=substr($cds,$j,3);
                        if($aa ne &translate($codon)){print STDERR "$aa != $codon at $id, use original seq\n"}
                        $cdsAlign.=$codon;
                        #push @prot," ".$aa." ";
                        #push @CDS,$coden;
                        $j += 3;
                }else{
                        $cdsAlign.= '---';
                        #push @prot," - ";
                        #push @CDS,"---";
                     }    
        }#foreach pepAlign
      }else{die "pep $id not exists in cds\n"}
      my $string=&faFormat($cdsAlign,50);
      push @cdsAlign,">$id\n$string";   
   }#pep id foreach
=cut

   my $result=join("",@cdsAlign);
   return $result;
}#sub end

sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=scalar(@{$index->{$id}})}else{if($len != scalar(@{$index->{$id}}) ){$flag=0}}
  }
  return ($flag,$len);
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


sub translate(){
#standard code table:http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
#    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#  Starts = ---M---------------M---------------M----------------------------
#  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
   my $codon=shift;
   $codon=uc($codon);
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

   if($codon=~/N/){return "X"}  #gap N tolerant, better if degap first
   return $code{$codon};

}



sub SD(){
  my $index=shift;
  my $total;
  my $cnt;
  foreach(@{$index}){
    return (1000) if($_ eq "NA");
    $total+=$_;
    $cnt++;
  }
  my $mean=$total/$cnt;
  #print "\n$mean\n";
  $total=0;
  foreach(@{$index}){
    $total+=($mean-$_)*($mean-$_);
  }
  #print "\n$total\n";
  return sqrt($total/($cnt-1));
}




