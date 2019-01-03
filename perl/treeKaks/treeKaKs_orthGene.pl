use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#1, extract pep and cds 
#2, align peps with muscle
#3, transform pep_alignment to CDS_alignment
#4, build phylogenetic trees with treebest and draw tree svg
#5, calculate Ka/Ks with KaKs_calculator

#usage: cat pairs.list.cluster | perl treeKaKs.pl -cds cds.list -pep pep.list -outdir outdir2 2> err > log &


my($list_pep,$list_cds,$outdir,$tmpdir);
GetOptions("pep=s",\$list_pep,"cds=s",\$list_cds,"outdir:s",\$outdir,"tmp:s",\$tmpdir);
$outdir||="out";
#$tmpdir||="tmp";
my $pep_index=&readFa_list(&readList($list_pep));
my $cds_index=&readFa_list(&readList($list_cds));
print STDERR "read in fa done \n";

#print Dumper $pep_index;
#print STDERR Dumper $cds_index;

#read in cluster tab
my @order=("japo","glab","punc","brac","lper","sorg");
my $cnt=0;
while(<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/);
  $cnt++;
  my @groups=split/\t/,$_;
  die "truncted line $_\n" if(scalar @groups != 6);    
  my @ids=split/\t|,/,$_;
  if(scalar @ids > 6){next}#print "big cluster at $_\n";next};
#1, extract prot and cds 
  my @seqs_prot;
  my @seqs_cds;
  for my $i(0..5){
    my @members=split/,/,$groups[$i];
    #print "@members|";
    foreach my $id(@members){
       if($id ne "-"){
         if(exists $pep_index->{$order[$i]}->{$id}){ push @seqs_prot, ">$id\n$pep_index->{$order[$i]}->{$id}\n" }else{die "not exists $id in $order[$i] pep file\n"} 
         if(exists $cds_index->{$order[$i]}->{$id}){ push @seqs_cds, ">$id\n$cds_index->{$order[$i]}->{$id}\n" }else{die "not exists $id in $order[$i] cds file\n"} 
       }
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
my ($coverage,$identity)=split/[ |\t]+/,`echo -ne "$aln" |perl bin/mfaStats_v1.2.pl -quiet`;
#print "$coverage,$identity\n";

#3, transform prot_alignment to cds_aligment, remove aligment gaps and transform to axt
#`perl ../bin/pepMfa_to_cdsMfa.pl temp_prot.muscle temp_cds.fa 2>temp_prot2cds.codeDetail|perl ../bin/degap2Axt.pl > temp_prot2cds.degap.axt ` ;
#if(! -d $tmpdir ){mkdir $tmpdir};
if(! -d $outdir ){if(! mkdir $outdir){die "can't mkdir $outdir: dir or file already exists\n"}};
open ALN, ">$outdir/orth${cnt}.cdsAlign.mfa" or die "$!"; 
#print  "$aln\n";
print ALN &pepAlign2cdsAlign($aln,$cds);
close ALN;

#4, build tree with treebest and draw tree svg file
#`treebest nj $outdir/orth${cnt}.cdsAlign.mfa > $outdir/orth${cnt}.cdsAlign.mfa.nhx 2> /dev/null`;
`treebest filter -n $outdir/orth${cnt}.cdsAlign.mfa > $outdir/orth${cnt}.cdsAlign.filter.mfa`;
`treebest nj $outdir/orth${cnt}.cdsAlign.filter.mfa > $outdir/orth${cnt}.cdsAlign.filter.mfa.nhx 2>/dev/null`;
`nw_display -s -w 1000 -v 100 -S -b 'opacity:0' -i 'font-size:10' -d 'stroke-width:2;stroke:blue' $outdir/orth${cnt}.cdsAlign.filter.mfa.nhx > $outdir/orth${cnt}.cdsAlign.filter.mfa.nhx.svg`;


#5, calculate kaks by KaKs_Calculator and check Ka/Ks
#`cat $outdir/orth${cnt}.cdsAlign.mfa |perl bin/mfa2axt.pl > $outdir/orth${cnt}.cdsAlign.mfa.axt`;
#`KaKs_Calculator -i $outdir/orth${cnt}.cdsAlign.mfa.axt -o $outdir/orth${cnt}.cdsAlign.mfa.axt.kaks -m NG `;
`cat $outdir/orth${cnt}.cdsAlign.filter.mfa |perl bin/mfa2axt.pl > $outdir/orth${cnt}.cdsAlign.filter.mfa.axt`;
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
  if($cnt % 1000 ==0){
     print STDERR "#";
     `mkdir $outdir/index$cnt`;    
     `mv $outdir/orth*  $outdir/index$cnt`;
  }
  #last if($cnt>5);
}#read stdin while end here
print STDERR "\ntotal $cnt orthologous groups done\n";

#####sub########

sub readList(){
  my $file=shift;
  my %list;
  open LIST,$file or die"$!";
  while(<LIST>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my($spec,$file)=split/[\t| ]+/,$_;
    if(!exists $list{$spec}){$list{$spec}=$file}else{die "dup species $spec\n"}
 }
  close LIST;
  return \%list;
}

sub readFa_list(){
    my $list=shift;
    my %fa;
    $/=">";
    foreach my $spec(sort keys %{$list}){
      open FA, $list->{$spec} or die "$!";
      while(<FA>){
        chomp;
        next if($_ eq "" || $_=~/^#/);
        my @box=split/\n+/,$_;
        my $id=shift @box;
        my $seq=join"",@box;
        if(!exists $fa{$spec}->{$id}){$fa{$spec}->{$id}=$seq}else{die "dup $id in file $list->{$spec}\n"}
      }
      close FA;
    }
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




