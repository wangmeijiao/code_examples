






 my ($gff,$gene,$mrna) = &gffRead($ARGV[0]);

 print Dumper $gff,$gene,$mrna;










###sub

sub gffRead_processed(){
  ## gff parse code for gffSimplyv1.3 processed 
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
    my ($ID,$pID,$name,$note);
    if($ftype eq "gene"){
      $attr=~/ID=([^;]+);Name=([^;]+)/;
      ($ID,$name)=($1,$2);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
      }else{die "dup gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+);Note=([^;]+)/;
        ($ID,$pID,$note)=($1,$2,$3);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID\t$note";
        }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR" || $ftype eq "tss" || $ftype eq "tts"){
           $attr=~/Parent=([^;]+);Note=([^;]+)/;
           ($pID,$note)=($1,$2);
           my $geneID;
           if(exists $mrna{$pID}){
            my @box=split/\t/,$mrna{$pID};
            if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
           }else{die "elements comes first, I can't assign it to gene\n"}
           if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
              my @temp;
              if($ftype eq "CDS"){push @temp,"$start-$end-$fshift-$note"}else{push @temp,"$start-$end-$note"};
              $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
           }else{ if($ftype eq "CDS"){push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$fshift-$note"}else{push  @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$note"} }
         }else{die "unknow feature type\n"}
  }#while end
  close GFF;
  return (\%gff,\%gene,\%mrna);

}


