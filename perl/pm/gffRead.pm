use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;







my ($gff_idx,$gene_idx,$mrna_idx)=&gffRead($ARGV[0]);


print Dumper $gff_idx,$gene_idx,$mrna_idx;





###sub##




sub gffRead(){
  ##standard nest gff3 parsing code v1.0,  @ Dec. 26, 2015
  ## modified 29 Nov, 2016
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
    if($ftype eq "gene"){
      $attr=~/ID=([^;]+);Name=([^;]+)/;
      ($ID,$name)=($1,$2);
      die "not match at gene $_" if(!defined $ID || !defined $name);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
      }else{die "dup gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+)/;
        ($ID,$pID)=($1,$2);
        die "not match at mRNA $_" if(!defined $ID || !defined $pID);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID";
        }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR" || $ftype eq "tss" || $ftype eq "tts"){
           $attr=~/Parent=([^;]+)/;
           $pID=$1;
           die "not match at feature $ftype: $_" if(!defined $pID);
           my $geneID;
           if(exists $mrna{$pID}){
            my @box=split/\t/,$mrna{$pID};
            if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
           }else{die "elements comes first, I can't assign it to gene\n"}
           if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
              my @temp;
              if($ftype eq "CDS"){push @temp,"$start-$end-$fshift"}else{push @temp,"$start-$end"};
              $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
           }else{ if($ftype eq "CDS"){push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$fshift"}else{push  @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end"} }
         }else{die "unknow feature type\n"}
  }#while end
  close GFF;
  return (\%gff,\%gene,\%mrna);

}



sub gffRead_old(){
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
    if($ftype eq "gene"){
      $attr=~/ID=([^;]+);Name=([^;]+)/;
      ($ID,$name)=($1,$2);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
      }else{die "dup gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+)/;
        ($ID,$pID)=($1,$2);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID\t$name";
        }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"){
           $attr=~/Parent=([^;]+)/;
           $pID=$1;
           my $geneID;
           if(exists $mrna{$pID}){
            my @box=split/\t/,$mrna{$pID};
            if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
           }else{die "elements comes first, I can't assign it to gene\n"}
           if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
              my @temp;
              if($ftype eq "CDS"){push @temp,"$start-$end-$fshift"}else{push @temp,"$start-$end"};
              $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
           }else{ if($ftype eq "CDS"){push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$fshift"}else{push  @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end"} }
         }else{die "unknow feature type\n"}
  }#while end
  close GFF;
  return (\%gff,\%gene,\%mrna);

}

sub gffRead_oldold(){
   #may change read format for each species gff
   my ($file,$region,$outdir,$prefix)=@_;
   die "file $file not exists\n" if(! -f $file);
   my ($chr,$s,$e)=split/:|-/,$region;
   my %gff_region;
   my $outfile="$outdir/${prefix}_${chr}_${s}_${e}.gff";
   open GFF, $file or die "$!";
   open GFFO, ">$outfile" or die "$!";
   while(<GFF>){  #only capture mrna and CDS lines
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my ($chr_in,$foo1,$ftype,$start,$end,$foo2,$strand,$fshift,$attr)=split/\t/,$_;
     die"empty value at $_\n" if(!defined $chr_in || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
     ($start,$end)=($start<$end)?($start,$end):($end,$start);
     next if($chr ne $chr_in); #not the same chr
     next if( $start >= $e || $end <= $s); #not overlap at all
     if($start < $s ){$start = $s}; # cut pending ends
     if($end > $e){$end = $e};
     my ($id,$pid);
     if($ftype eq "mRNA"){
        $attr=~/ID=([^;]+)/;
        $id=$1;
        if(!exists $gff_region{$id}){
          $gff_region{$id}->{"range"}="$start-$end";
          $gff_region{$id}->{"strand"}=$strand;
          my ($start_real,$end_real);
          $start_real=$start-$s+1; #change to real coord for gff+fa -> embl 
          $end_real=$end-$s;
          print GFFO "$chr_in\t$foo1\t$ftype\t$start_real\t$end_real\t$foo2\t$strand\t$fshift\tID=$id\n";
        }else{die "CDS comes before mrna, quit\n"}
     }elsif($ftype eq "CDS"){
        $attr=~/Parent=([^;]+)/;
        $pid=$1;
        if(exists $gff_region{$pid}){
          push @{$gff_region{$pid}->{"CDS"}},"$start-$end";
          my ($start_real,$end_real);
          $start_real=$start-$s+1; #change to real coord for gff+fa -> embl 
          $end_real=$end-$s;
          print GFFO "$chr_in\t$foo1\t$ftype\t$start_real\t$end_real\t$foo2\t$strand\t$fshift\tParent=$pid\n";
        }else{die "CDS comes before mrna, quit\n"}
      }
   }#while end
   close GFF;
   close GFFO;
   #print Dumper \%gff_region;  
   return (\%gff_region,$outfile);

}




