use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#embl format:
#
#
#FH       key         Location/Qualifiers
#
#FT       feature_key   regions
#FT                     /qualifier1="  "
#FT                     /qualifier2="  "
#   ...
#
#SQ     sequence 1000 Bp ; ...
#       CGTAGCTACGAT   GTAGTTTCGATG  ACGTAGTCACAG  36
#          ...

#//

# valid feature_keys/qualifiers groups see  ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html


# usage: cat test.gene.gff |perl thisfile.pl  > test.gene.gff.embl


my ($gff_idx,$mrna_idx)=&gffRead_mRNA($ARGV[0]);

#print Dumper $gff_idx,$mrna_idx;
#exit;


printf("%-5s%-18s%s\n", "FH","Key","Location/Qualifier");
printf("%-23s\n", "FH");

foreach my $chr(keys %{$gff_idx}){
  #foreach my $geneID(sort keys %{$gff_idx->{$chr}}){
     foreach my $mrnaID(keys %{$gff_idx->{$chr}}){
         #my ($chr_gene,$start_gene,$end_gene,$ID_gene,$strand_gene,$name_gene)=split/\t/,$gene_idx->{$geneID};
         my ($chr_mrna,$start_mrna,$end_mrna,$ID_mrna,$strand_mrna,$pID)=split/\t/,$mrna_idx->{$mrnaID};
         #die "one of them not eq: $chr_gene ne $chr_mrna || $strand_gene ne $strand_mrna || $ID_gene ne $pID" if($chr_gene ne $chr_mrna || $strand_gene ne $strand_mrna || $ID_gene ne $pID);
         my $element=$gff_idx->{$chr}->{$mrnaID};
         #embl mRNA feature and related qualifiers    
         if(exists $element->{"exon"}){
            my @regions=@{$element->{"exon"}};
            if(exists $element->{"five_prime_UTR"}){foreach(@{$element->{"five_prime_UTR"}}){push @regions,$_}}
            if(exists $element->{"three_prime_UTR"}){foreach(@{$element->{"three_prime_UTR"}}){push @regions,$_}}
            my $coord=&getCoord(\@regions,$strand_mrna);
            printf("%-5s%-18s%s\n","FT","mRNA",$coord);         
            printf("%-23s%s\n","FT","/gene=\"$mrnaID\"");         
            printf("%-23s%s\n","FT","/note=\"$mrnaID\"");         
         }
         #embl CDS feature and related qualifiers
         if(exists $element->{"CDS"}){
             my $coord=&getCoord($element->{"CDS"},$strand_mrna);
             if(!exists $element->{"exon"}){
               printf("%-5s%-18s%s\n","FT","mRNA",$coord);
               printf("%-23s%s\n","FT","/gene=\"$mrnaID\"");
               printf("%-23s%s\n","FT","/note=\"$mrnaID\"");
             }
             printf("%-5s%-18s%s\n","FT","CDS",$coord);        
             printf("%-23s%s\n","FT","/gene=\"$mrnaID\""); 
             printf("%-23s%s\n","FT","/colour=3"); 
         }

         die "exon or CDS must have one at $mrnaID" if(!exists $element->{"exon"} && !exists $element->{"CDS"});

     }#foreach mrna
  #}#foreach gene
}#foreach chr







###sub##



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
      $attr=~/ID=([^;]+);/;
      $ID=$1;
      $name = $ID;
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
           $attr=~/Parent=([^;]+);/;
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


sub gffRead(){
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
      $attr=~/ID=([^;]+);Description=([^;]+);/;
      ($ID,$name)=($1,$2);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
      }else{die "dup gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+);/;
        ($ID,$pID)=($1,$2);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID";
        }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"){
           $attr=~/Parent=([^;]+);/;
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


sub getCoord(){
   my ($region,$strand)=@_;
   my $string;
   my @regions_sorted=@{&segSort_start($region)};
   if($strand eq "+"){
     if(scalar @regions_sorted == 1){
       my ($s,$e)=split/-/,$regions_sorted[0];
       $string="$s..$e";
     }elsif(scalar @regions_sorted > 1){
       my @str;
       foreach my $seg(@regions_sorted){
         my ($s,$e)=split/-/,$seg;
         push @str,"$s..$e";
       }
       $string=join ",",@str;
       $string="join($string)";
     }else{die "epmty segs at @regions_sorted"}
   }elsif($strand eq "-"){
       if(scalar @regions_sorted == 1){
         my ($s,$e)=split/-/,$regions_sorted[0];
         $string="complement($s..$e)";
       }elsif(scalar @regions_sorted > 1){
          my @str;
          foreach my $seg(@regions_sorted){
            my ($s,$e)=split/-/,$seg;
            push @str,"$s..$e";
          }
          $string=join ",",@str;
          $string="complement(join($string))";
        }else{die "epmty segs at @regions_sorted"}
     }else{die "unknow strand $strand at @{$region}"}
   return $string;
}


sub segSort_start(){ #simple sort by start
     #start < end
     my @segs=@{shift @_};
     my @segs_sorted=sort{
       my ($s1,$e1)=split/-/,$a;
       my ($s2,$e2)=split/-/,$b;
       die "s > e at @segs" if($s1 > $e1 || $s2 > $e2);
       if($s1<$s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
     } @segs;
     return \@segs_sorted;
}




