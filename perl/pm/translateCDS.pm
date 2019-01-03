














sub translate_CDS_simple(){
#standard code table:http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
#    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#  Starts = ---M---------------M---------------M----------------------------
#  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
   my $orf=shift;
   #my $frameshift=shift; #0,1,2: delete 0,1,2 nucleiotide before translating
   print STDERR "Warning! orf undefined\n" if(!defined $orf);
   #if($frameshift eq "." || $frameshift eq "" || !defined $frameshift){$frameshift=0}
   $orf=uc($orf);
   my $len = length $orf;
   if($len % 3 != 0){my $num = $len % 3; $orf=~s/.{$num}$//;    print STDERR "\nWarning! ##not3x at orf $orf, chomp tails to fit 3x\n"}
   if($len < 3){die "err orf <3 :$orf"}
   #$orf=substr($orf,$frameshift);
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
   for(my$i=1;$i<=$len-2;$i+=3){
     $codon=substr($orf,$i-1,3);
     if($codon=~/n+/i){$aa.="X";next} #X for gap
     if($i!=(length $orf)-2 && $codon=~/TAA|TAG|TGA/){$aa.="*";print STDERR "\nWarning! ##earlier_stop at orf $orf\n";next}
     if(exists $code{$codon}){$aa.=$code{$codon}}else{$aa.="X";print STDERR "\nWarning! ##not_in_codon_tab $codon at orf $orf\n";next} #X for unknown
     #if(exists $code{$codon}){$aa.=$code{$codon}}else{return "##wrong_codon"}
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


