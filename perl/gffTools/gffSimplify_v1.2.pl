use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

# the same method in gff2beds; 
# filter raw gff as following:
#       filter short gene with len<cutoff (200)
#       remove TE-related
#       output longest mrna isoforms
#       calculate and add the exon and/or intron if not exist in original file
#       haven't deal with overlap genes and split gene fragments
#       gff attribute string formatted as "ID=...;Note=...  Parent=...;"
# usage: cat all.gff3 |perl thisfile.pl -cutoff 200 -exon -intron >all_tigr6.simplified.gff



my ($cutoff,$exon,$intron);
GetOptions("cutoff=i",\$cutoff,"exon!",\$exon,"intron!",\$intron);


my %gff;#main structure table:  chr->gene->mRNAs->(elements:CDS,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)
my %gene;#gene info
my %mrna;#mrna info
#1,read&fill without order consideration
while(<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
  my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
  die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
  next if($chr eq "ChrUn" || $chr eq "ChrSy"); #get rid of other chrs
  #next if($ftype eq "chromosome" || $ftype eq "Chromosome");
  ($start,$end)=($start<$end)?($start,$end):($end,$start);
  my ($ID,$pID,$note,$name);
  if($ftype eq "gene"){
    $attr=~/ID=([^;]+);Name=([^;]+);Note=([^;]+)/;
    ($ID,$name,$note)=($1,$2,$3);
    #$alias=$ID;
    $note=~s/(%20|%2C|%2D|%2E|%2F|%5F|%3A|%3B|%2A|%2B)+/_/g;
    if(!exists $gff{$chr}->{$ID}){
      my %temp;
      $gff{$chr}->{$ID}=\%temp;
      $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$note\t$name";
    }else{die "dup gene name $ID\n"}
   }elsif($ftype eq "mRNA"){
      $attr=~/ID=([^;]+);Name=([^;]+);Parent=([^;]+)/;
      ($ID,$name,$pID)=($1,$2,$3);
      if(!exists $gff{$chr}->{$pID}->{$ID}){
        my %temp;
        $gff{$chr}->{$pID}->{$ID}=\%temp;
        $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID\t$name";
      }else{die "dup mRNA ID at $ID\n"}
     }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"){
         $attr=~/ID=([^;]+);Parent=([^;]+)/;
         ($ID,$pID)=($1,$2);
         my (undef,$order) = split/:/,$ID;
         die "can't fetch feature order at $pID" if($order eq "" || !defined $order);
         my $geneID;
         if(exists $mrna{$pID}){
          my @box=split/\t/,$mrna{$pID};
          if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
         }else{die "elements comes first, I can't assign it to gene\n"}
         if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
            my @temp;
            if($ftype eq "CDS"){push @temp,"$order:$start-$end-$fshift"}else{push @temp,"$order:$start-$end"};
            $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
         }else{ if($ftype eq "CDS"){push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$order:$start-$end-$fshift"}else{push  @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$order:$start-$end"} }
       }else{die "unknow feature type\n"}
}#while end

#print Dumper \%gff;
#print Dumper \%gene;
#print Dumper \%mrna;
#exit;


#2,filter to output a simplified gff
open TE_G,">te_related.gff" or die "$!";
open SHORT,">short.gff" or die "$!";
$cutoff||=200; #gene length biger than this
#print STDERR "$cutoff\n";

foreach my$chr(sort keys %gff){
 #gene filter the short ones and te-related
 foreach my $gene(sort keys %{$gff{$chr}}){
   my ($chr,$start,$end,$id,$strand,$note,$name)=split/\t/,$gene{$gene};
   if($end-$start+1<$cutoff){print SHORT "$chr\tgff_simplify\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$note\n";next}
   if($note!~/retrotransposon|transposon/){
   #if($note!~/retrotransposon|transposon|rRNA|tRNA|mitochondrion|chloroplast/){
    print "$chr\tgff_simplify\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Note=$note\n";
    my $max_len=0;
    my $max_len_mrna;
    foreach my $mrna(keys %{$gff{$chr}->{$gene}}){
      my @box=split/\t/,$mrna{$mrna};
      if($max_len<($box[2]-$box[1])){$max_len=$box[2]-$box[1];$max_len_mrna=$mrna;}
    }

    #mrna (the longest one, which all of the following caculation are based on)
    my $index=$gff{$chr}->{$gene}->{$max_len_mrna};
    my ($chr_mrna,$start_mrna,$end_mrna,$id_mrna,$strand_mrna,$pID,$name_mrna)=split/\t/,$mrna{$max_len_mrna};
    print "$chr_mrna\tgff_simplify\tmRNA\t$start_mrna\t$end_mrna\t.\t$strand_mrna\t.\tID=$id_mrna;Parent=$id\n";

    #CDS
    my @CDS_store;# for exon and intron calculation
    if(exists $index->{"CDS"}){
       my @CDS=@{&segSort($index->{"CDS"})};
       my $i;
       foreach(@CDS){
         $i++;
         my($start,$end,$fshift)=split/-/,$_;
         print "$chr_mrna\tgff_simplify\tCDS\t$start\t$end\t.\t$strand_mrna\t$fshift\tParent=$id_mrna\n";
         push @CDS_store,"$start-$end";
       } 
    }else{die "no CDS at $id\n"} #CDS is a must

    #exon 
    if(exists $index->{"exon"}){
       my @exon=@{&segSort($index->{"exon"})};
       my $i;
       foreach(@exon){
         $i++;
         my($start,$end)=split/-/,$_;
         print "$chr_mrna\tgff_simplify\texon\t$start\t$end\t.\t$strand_mrna\t.\tParent=$id_mrna\n";
       } 
    }
            
    #intron
    if(exists $index->{"intron"}){
       my @intron=@{&segSort($index->{"intron"})};
       my $i;
       foreach(@intron){
         $i++;
         my($start,$end)=split/-/,$_;
         print "$chr_mrna\tgff_simplify\tintron\t$start\t$end\t.\t$strand_mrna\t.\tParent=$id_mrna\n";
       } 
    }

    #UTRs (utrs are start or end exons)
    my @utr5;# for exon and intron calculation
    if(exists $index->{"five_prime_UTR"}){
       my @box=@{&segSort($index->{"five_prime_UTR"})};
       my $i;
       foreach(@box){
         $i++;
         my($start,$end,$fshift)=split/-/,$_;
         print "$chr_mrna\tgff_simplify\tfive_prime_UTR\t$start\t$end\t.\t$strand_mrna\t.\tParent=$id_mrna\n";
         push @utr5,"$start-$end";
       }
    }
    my @utr3;# for exon and intron calculation
    if(exists $index->{"three_prime_UTR"}){
       my @box=@{&segSort($index->{"three_prime_UTR"})};
       my $i;
       foreach(@box){
         $i++;
         my($start,$end,$fshift)=split/-/,$_;
         print "$chr_mrna\tgff_simplify\tthree_prime_UTR\t$start\t$end\t.\t$strand_mrna\t.\tParent=$id_mrna\n";
         push @utr3,"$start-$end";
       }
    }

    #calculate exon and/or intron from mrna, CDS and UTRs if not exist
    if(! exists $index->{"exon"} && $exon ){
       my $exon;
       die "empty CDS in $id_mrna" if(scalar @CDS_store == 0);
       $exon=&getExon("$start_mrna-$end_mrna",\@CDS_store,\@utr5,\@utr3,$strand_mrna);
       die "empty exon in $id_mrna" if(scalar @$exon == 0);
       foreach(@$exon){
          my($start,$end)=split/-/,$_;
          print "$chr_mrna\tgff_simplify\texon\t$start\t$end\t.\t$strand_mrna\t.\tParent=$id_mrna\n";
       }
       #calculate intron if need
       if(scalar @{$exon} >1 && ! exists $index->{"intron"} && $intron){
         my $intron=&getIntron($exon);
         foreach(@$intron){
             my($start,$end)=split/-/,$_;
             print "$chr_mrna\tgff_simplify\tintron\t$start\t$end\t.\t$strand_mrna\t.\tParent=$id_mrna\n";
         }
       }elsif(scalar @{$exon} == 1 && ! exists $index->{"intron"} && $intron){
              print STDERR "$id_mrna has one exon :@$exon and can't calculate intron\n"
           }
    }#calculate end

  }else{print TE_G "$chr\tgff_simplifyv1\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Note=$note\n";next} #if gene te end
  print "##\n";
 }#foreach gene end
}#foreach chr end

close TE_G;
close SHORT;



#####subs####
sub segSort(){
 #sort with s < e and nonOverlap check
 my @segs=@{shift @_}; #s-e-frameshift or s-e
 my @sorted=sort{
      my ($s1,$e1,$s2,$e2);#$s1<=$e1 && $s2<=$e2 and deOverlaped beforehand by default
      my ($order1,$string1) = split/:/,$a;
      ($s1,$e1)=split/-/,$string1;
      my ($order2,$string2) = split/:/,$b;
      ($s2,$e2)=split/-/,$string2;
      if($s1 > $e1 || $s2 > $e2){die "start > end at $s1 > $e1 || $s2 > $e2"}
      if($s1 < $s2 && $e1 > $s2 || $s2 < $s1 && $e2 > $s1){die "overlap found at $a,$b"}  #unOverlap
      if($e2<=$s1){1}elsif($e1<=$s2){-1}else{die "unknown circumstance at $s1,$e1:$s2,$e2\n"}
    } @segs;
# print "@CDS\n@sorted\n";
 return \@sorted;
}


sub segSort_order(){
 #sort with order tag "exon_1" 
 my @segs=@{shift @_}; #order:s-e-frameshift or order:s-e
 my @sorted=sort{
      my ($order1, $order2);# order1 order2 .. 
      ($order1,undef)=split/:/,$a;
      ($order2,undef)=split/:/,$b;
      if($order1 gt $order2){1}elsif($order1 le $order2 ){-1}else{die "unknown circumstance at $order1,$order2\n"}
    } @segs;
# print "@CDS\n@sorted\n";
 return \@sorted;
}



sub getExon(){
    #calculate sorted exon regions from CDS, mrna and 5UTR, 3UTR (UTRs are parts of exons)
    #  mrna_range       ------------------------------------------------------------
    #     CDSs                     ------     --------------  ---    ----
    #     UTRs          ----   ----                                      ------   --
    #                  5utr_1   2                                        3utr_1    2
    #     exons         ----   ----------     --------------  ---    ----------   --   
    #splicing_out_introns   \ /          \   /              \/   \  /          \ /
    #                     1         2                3         4          5        6
    #                                 
    my ($mrna_range,$cds,$utr5,$utr3,$strand)=@_;
    my @exon;
    my ($mrna_s, $mrna_e)=split/-/,$mrna_range;
    die "mrna start > end in getExon function\n" if($mrna_s >= $mrna_e);
    die "empty cds" if(scalar @{$cds} == 0);
    #utr5
    foreach (@{$utr5}){
      push @exon,$_;

    } 
    #cds
    foreach (@{$cds}){
      push @exon,$_;

    } 
    #utr3
    foreach (@{$utr3}){
      push @exon,$_;

    } 
   my @exon_sorted_linked=@{&segLink(\@exon,1)};

   my ($exon_s,undef)=split/-/,$exon_sorted_linked[0];
   my (undef,$exon_e)=split/-/,$exon_sorted_linked[-1];
   #if($exon_s != $mrna_s || $exon_e != $mrna_e){ die "exon start or end != mrna start or end: $exon_s, $exon_e != $mrna_s, $mrna_e"}

   return \@exon_sorted_linked;
}

sub getIntron(){
  my @box=@{shift @_}; 
  my @intron;
  for(0..$#box-1){
    my ($s1,$e1)=split/-/,$box[$_];
    my ($s2,$e2)=split/-/,$box[$_+1];
    push @intron,($e1+1)."-".($s2-1);
  }
  return \@intron;
}



sub segLink(){
 #link small segments to longer ones with gap<=$gapSize
    my ($index,$gapSize)=@_;
    my @segments=@{&segSort($index)};#sorted
    my @segments_linked;
    for (my $i=1;$i<=$#segments;$i++){
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        if(($tmp2[0]-$tmp1[1]-1)<=$gapSize){
          my $start=$tmp1[0];
          my $end=$tmp2[1];
          $segments[$i]=$start."-".$end;
        }else{push @segments_linked, $segments[$i-1];}
        if($i==$#segments){push @segments_linked,$segments[$i];}
    }
    if($#segments==0){push @segments_linked,$segments[0];}
    return(\@segments_linked);
}



