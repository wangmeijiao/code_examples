
my %gtf;#main structure table:  chr->gene->mRNAs->(elements:CDS,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)
my %gene;#gene info
my %mrna;#mrna info
#1,read&fill without order consideration
while(<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/);
  my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
  die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
  #$attr=~s/"//g;
  next if($chr eq "chrUn" || $chr eq "chrSy" || $chr eq "Pt" || $chr eq "Mt" ); #get rid of other chrs
  next if($ftype eq "chromosome" || $ftype eq "Chromosome");
  ($start,$end)=($start<$end)?($start,$end):($end,$start);
  my ($ID,$tID,$name,$tname);
  if($ftype eq "gene"){
    $attr=~/gene_id "([^;]+)"; gene_symbol "([^;]+)";/;
    ($ID,$name)=($1,$2);
    if(!exists $gtf{$chr}->{$ID}){
      my %temp;
      $gtf{$chr}->{$ID}=\%temp;
      $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
    }else{die "dup gene name $ID\n"}
   }elsif($ftype eq "mRNA" || $ftype eq "miRNA" || $ftype eq "ncRNA" || $ftype eq "pre_miRNA" || $ftype eq "snoRNA"|| $ftype eq "snRNA"|| $ftype eq "tRNA"|| $ftype eq "rRNA"|| $ftype eq "pseudogene" ){
      $attr=~/gene_id "([^;]+)"; gene_symbol "([^;]+)"; transcript_id "([^;]+)"; transcript_symbol "([^;]+)";/;
      ($ID,$name,$tID,$tname)=($1,$2,$3,$4);
      if(!exists $gtf{$chr}->{$ID}->{$tID}){
        my %temp;
        $gtf{$chr}->{$ID}->{$tID}=\%temp;
        $mrna{$tID}="$chr\t$start\t$end\t$tID\t$strand\t$ID\t$name";
      }else{die "dup mRNA ID at $ID\n"}
     }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "5UTR" ||  $ftype eq "3UTR" || $ftype eq "start_codon" || $ftype eq "stop_codon"){
         $attr=~/gene_id "([^;]+)"; gene_symbol "([^;]+)"; transcript_id "([^;]+)"; transcript_symbol "([^;]+)";/;
         ($ID,$name,$tID,$tname)=($1,$2,$3,$4);
         if(!exists $gtf{$chr}->{$ID}->{$tID}->{$ftype}){
            my @temp;
            if($ftype eq "CDS"){push @temp,"$start-$end-$fshift"}else{push @temp,"$start-$end"};
            $gtf{$chr}->{$ID}->{$tID}->{$ftype}=\@temp;
         }else{ if($ftype eq "CDS"){push @{$gtf{$chr}->{$ID}->{$tID}->{$ftype}},"$start-$end-$fshift"}else{push  @{$gtf{$chr}->{$ID}->{$tID}->{$ftype}},"$start-$end"} }
       }else{die "unknow feature type\n"}
}#while end


#print Dumper \%gtf;
#print Dumper \%gene;
#print Dumper \%mrna;
#exit;


