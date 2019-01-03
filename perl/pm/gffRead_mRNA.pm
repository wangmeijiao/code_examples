
sub gffRead_mRNA(){
  ##standard nest gff3 parsing code v1.0,  @ Dec. 26, 2015
  my $file=shift;
  my %gff;#main structure table:  chr->gene->mRNAs->(elements:CDS,exon,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)
  my %gene;#gene info
  my %mrna;#mrna info
  #1,read&fill without order consideration
  open GFF, $file or die "$!";
  while(<GFF>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
    die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
    next if($chr eq "chrUn" || $chr eq "chrSy" || $chr eq "Pt" || $chr eq "Mt"); #get rid of other chrs
    next if($ftype eq "chromosome" || $ftype eq "Chromosome");
    ($start,$end)=($start<$end)?($start,$end):($end,$start);
    my ($ID,$pID,$name);
    if($ftype eq "mRNA"){
      $attr=~/ID=([^;]+);Note=([^;]+)/;
      ($ID,$name)=($1,$2);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
      }else{die "dup gene name $ID\n"}
     #}elsif($ftype eq "mRNA"){
      #  $attr=~/Parent=([^;]+)/;
      #  $pID = $1;
      #  if(!exists $gff{$chr}->{$pID}->{$ID}){
      #    my %temp;
      #    $gff{$chr}->{$pID}->{$ID}=\%temp;
      #    $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID";
      #  }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"){
           $attr=~/Parent=([^;]+)/;
           $pID=$1;

           #my $geneID;
           #if(exists $mrna{$pID}){
           # my @box=split/\t/,$mrna{$pID};
           # if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
           #}else{die "elements comes first, I can't assign it to gene\n"}
           
           if(!exists $gff{$chr}->{$pID}->{$ftype}){
              my @temp;
              if($ftype eq "CDS"){push @temp,"$start-$end-$fshift"}else{push @temp,"$start-$end"};
              $gff{$chr}->{$pID}->{$ftype}=\@temp;
           }else{ if($ftype eq "CDS"){push @{$gff{$chr}->{$pID}->{$ftype}},"$start-$end-$fshift"}else{push  @{$gff{$chr}->{$pID}->{$ftype}},"$start-$end"} }
         }else{die "unknow feature type\n"}
  }#while end
  close GFF;
  return (\%gff,\%mrna);

}


