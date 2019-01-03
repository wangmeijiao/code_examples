use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;


#####minor modification: add intron calculation from already exist exons  Jun 6, 2017
##### modify getIntron_strand() to skip intron calculation if two adjacent exon too near ( dist < -2)

#######version 1.3 update @ Jul 30, 2016
#
#add mrna attribute tag "longest_AS",  add exon cds intron attribute tag "note=calculated_exon$i" or "note=cds$i" if not calculated, order by transcription direction(strand)
#add tss and tts
#report TE_related gene details (for TE structure examination)
#do exon annotation:  exon1=5utr1 exon2=5utr2+cds1  exon3=cds2 ... based on two annotation system: exon-intron,  5utr-cds-3utr
#usage: cat all.gff3 | perl gffSimplify_v1.3.pl -cutoff 200 -exon -intron -tss -tts -anno > all.simple.gff3  2> nointron.err &        & grep -v "^#" all.simple.gff3 |sort -k1,1 -k4,4n > all.simple.sorted.gff3
#note: can't deal with overlap gene and potential gene fragments



# version 1.2
# the same method in gff2beds; 
# filter raw gff as following:
#       omit short gene with len<cutoff (200)
#       remove TE-related
#       choose the longest mrna isoforms
#       calculate and add the exon and/or intron if not exist in original file, append to the records
# usage: cat all.gff3 |perl thisfile.pl -cutoff 200 -exon -intron >all_tigr6.simplified.gff



my ($cutoff,$exon,$intron,$tss,$tts,$exon_anno);
GetOptions("cutoff=i",\$cutoff,"exon!",\$exon,"intron!",\$intron,"tss!",\$tss,"tts!",\$tts,"anno!",\$exon_anno);


my %gff;#main structure table:  chr->gene->mRNAs->(elements:CDS,intron,UTR)
my %gene;#gene info
my %mrna;#mrna info
#1,read&fill without order consideration
while(<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
  my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
  die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
  next if($chr eq "chrUn" || $chr eq "chrSy"); #get rid of other chrs
  #next if($ftype eq "chromosome" || $ftype eq "Chromosome");
  die "start > end :$start > $end at $_" if($start > $end);
  #($start,$end)=($start<$end)?($start,$end):($end,$start);
  my ($ID,$pID,$note,$alias);
  if($ftype eq "gene"){
    $attr=~/ID=([^;]+);Name=([^;]+);Alias=([^;]+)/;
    ($ID,$note,$alias)=($1,$2,$3);
    #$alias=$ID;
    $note=~s/(%20|%2C|%2D|%2E|%2F|%5F|%3A|%3B|%2A|%2B)+/_/g;
    if(!exists $gff{$chr}->{$ID}){
      my %temp;
      $gff{$chr}->{$ID}=\%temp;
      $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$note\t$alias";
    }else{die "dup gene name $ID\n"}
   }elsif($ftype eq "mRNA"){
      $attr=~/ID=([^;]+);Parent=([^;]+);Alias=([^;]+)/;
      ($ID,$pID,$alias)=($1,$2,$3);
      if(!exists $gff{$chr}->{$pID}->{$ID}){
        my %temp;
        $gff{$chr}->{$pID}->{$ID}=\%temp;
        $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID\t$alias";
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
   my ($chr,$start,$end,$id,$strand,$note,$alias)=split/\t/,$gene{$gene};
   if($end-$start+1<$cutoff){print SHORT "$chr\tgff_simplifyv1.3\tgene\t$start\t$end\t.\t$strand\t.\tID=$alias;Name=$note\n";next}
   if($note!~/retrotransposon|transposon/){
   #if($note!~/retrotransposon|transposon|rRNA|tRNA|mitochondrion|chloroplast/){
    print "$chr\tgff_simplifyv1.3\tgene\t$start\t$end\t.\t$strand\t.\tID=$alias;Name=$note\n";
    my $max_len=0;
    my $max_len_mrna;
    my @group;
    foreach my $mrna(keys %{$gff{$chr}->{$gene}}){
      my @box=split/\t/,$mrna{$mrna};
      if($max_len<($box[2]-$box[1])){$max_len=$box[2]-$box[1];$max_len_mrna=$mrna} #choose the first if equal
      push @group,$mrna;
    }
    my $group = join"-",@group;
    #mrna (the longest one, which all of the following caculation are based on)
    my $index=$gff{$chr}->{$gene}->{$max_len_mrna};
    my ($chr_mrna,$start_mrna,$end_mrna,$id_mrna,$strand_mrna,$pID,$alias_mrna)=split/\t/,$mrna{$max_len_mrna};
    die "$id_mrna ne $max_len_mrna at $pID" if($id_mrna ne $max_len_mrna);
    print "$chr_mrna\tgff_simplifyv1.3\tmRNA\t$start_mrna\t$end_mrna\t.\t$strand_mrna\t.\tID=$alias_mrna;Parent=$alias;Note=longest_AS;Group=$group;\n";

    #the 5utr-cds-3utr system
    #CDS
    my @CDS_store;# for exon and intron calculation
    if(exists $index->{"CDS"}){
       my @CDS=@{&segSort_strand($index->{"CDS"},$strand_mrna)};
       my $i;
       foreach(@CDS){
         $i++;
         my($start,$end,$fshift)=split/-/,$_;
         print "$chr_mrna\tgff_simplifyv1.3\tCDS\t$start\t$end\t.\t$strand_mrna\t$fshift\tParent=$alias_mrna;Note=CDS$i\n";
         push @CDS_store,"$start-$end";
       } 
    }else{die "no CDS at $id\n"} #CDS is a must

    #UTRs (utrs are start or end exons)
    my @utr5;# for exon and intron calculation
    if(exists $index->{"five_prime_UTR"}){
       my @box=@{&segSort_strand($index->{"five_prime_UTR"},$strand_mrna)};
       my $i;
       foreach(@box){
         $i++;
         my($start,$end)=split/-/,$_;
         print "$chr_mrna\tgff_simplifyv1.3\tfive_prime_UTR\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=utr5_$i\n";
         push @utr5,"$start-$end";
       }
    }
    my @utr3;# for exon and intron calculation
    if(exists $index->{"three_prime_UTR"}){
       my @box=@{&segSort_strand($index->{"three_prime_UTR"},$strand_mrna)};
       my $i;
       foreach(@box){
         $i++;
         my($start,$end)=split/-/,$_;
         print "$chr_mrna\tgff_simplifyv1.3\tthree_prime_UTR\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=utr3_$i\n";
         push @utr3,"$start-$end";
       }
    }

    #the exon-intron system
    #exon 
    if(exists $index->{"exon"}){
       my @exon=@{&segSort_strand($index->{"exon"},$strand_mrna)};
       my $i;
       foreach(@exon){
         $i++;
         my($start,$end)=split/-/,$_;
         print "$chr_mrna\tgff_simplifyv1.3\texon\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=exon$i\n";
       } 
    }
            
    #intron
    if(exists $index->{"intron"}){
       my @intron=@{&segSort_strand($index->{"intron"},$strand_mrna)};
       my $i;
       foreach(@intron){
         $i++;
         my($start,$end)=split/-/,$_;
         print "$chr_mrna\tgff_simplifyv1.3\tintron\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=intron$i\n";
       } 
    }

    ####do calculation if required and not exist

    #calculate TSS and TTS from longest mrna AS (not gene)   
    if(! exists $index->{"tss"} && $tss ){
       if($strand_mrna eq "+"){ print "$chr_mrna\tgff_simplifyv1.3\ttss\t$start\t$start\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated\n"  }elsif($strand_mrna eq "-"){print "$chr_mrna\tgff_simplify\ttss\t$end\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated\n"}else{die "unknow mrna strand $strand_mrna"}
    }
    if(! exists $index->{"tts"} && $tts ){
       if($strand_mrna eq "+"){ print "$chr_mrna\tgff_simplifyv1.3\ttts\t$end\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated\n"  }elsif($strand_mrna eq "-"){print "$chr_mrna\tgff_simplify\ttts\t$start\t$start\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated\n"}else{die "unknow mrna strand $strand_mrna"}
    }#calculate end


    #calculate exon and/or intron from mrna, CDS and UTRs if not exist
    if(! exists $index->{"exon"} && $exon ){
       my $exon;
       die "empty CDS in $alias_mrna" if(scalar @CDS_store == 0);
       $exon=&getExon("$start_mrna-$end_mrna",\@CDS_store,\@utr5,\@utr3,$strand_mrna);
       die "empty exon in $alias_mrna" if(scalar @$exon == 0);
       my $i;
       foreach(@$exon){
          my($start,$end,$anno)=split/-|\|/,$_;
          $i++;
          print "$chr_mrna\tgff_simplifyv1.3\texon\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated_exon$i;Anno=$anno\n";
       }
       #calculate intron if need
       if(scalar @{$exon} >1 && ! exists $index->{"intron"} && $intron){
         my $intron=&getIntron_strand($exon,$strand_mrna);
         my $i;
         foreach(@$intron){
             my($start,$end)=split/-/,$_;
             $i++;
             print "$chr_mrna\tgff_simplifyv1.3\tintron\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated_intron$i\n";
         }
       }elsif(scalar @{$exon} == 1 && ! exists $index->{"intron"} && $intron){
              print STDERR "$alias_mrna has one exon :@$exon and can't calculate intron\n"
           }
    }#calculate end

    #calculate intron from already exists exons
    if(!exists $index->{"intron"} && exists $index->{"exon"} && $intron){
       my $exon = $index->{"exon"};
       my $exon_sort = &segSort_strand($exon,$strand_mrna);
       my $intron_calc = &getIntron_strand($exon_sort,$strand_mrna);
       my $cnt;
       foreach my $seg(@{$intron_calc}){
          my ($s,$e) = split/-/,$seg;
          $cnt++;
          print "$chr_mrna\tgff_simplifyv1.3\tintron\t$s\t$e\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=calculated_intron$cnt\n";
      }
    }#calculate end 

  }else{
         #TE related gff, output as is; useful for DNA TE detail annotation, i.e CACTA transpoase 
         print TE_G "$chr\tgff_simplify\tgene\t$start\t$end\t.\t$strand\t.\tID=$alias;Name=$note\n";
         foreach my $mrna(sort keys %{$gff{$chr}->{$gene}}){
           my $index=$gff{$chr}->{$gene}->{$mrna};
           my ($chr_mrna,$start_mrna,$end_mrna,$id_mrna,$strand_mrna,$pID,$alias_mrna)=split/\t/,$mrna{$mrna};
           print TE_G "$chr_mrna\tgff_simplifyv1.3\tmRNA\t$start_mrna\t$end_mrna\t.\t$strand_mrna\t.\tID=$alias_mrna;Parent=$alias\n";

           #CDS
           if(exists $index->{"CDS"}){
              my @CDS=@{&segSort_strand($index->{"CDS"},$strand_mrna)};
              my $i;
              foreach(@CDS){
                $i++;
                my($start,$end,$fshift)=split/-/,$_;
                print TE_G "$chr_mrna\tgff_simplifyv1.3\tCDS\t$start\t$end\t.\t$strand_mrna\t$fshift\tParent=$alias_mrna;Note=cds$i\n";
              }
           }else{die "no CDS at $id\n"} #CDS is a must
           #exon 
           if(exists $index->{"exon"}){
              my @exon=@{&segSort_strand($index->{"exon"},$strand_mrna)};
              my $i;
              foreach(@exon){
                $i++;
                my($start,$end)=split/-/,$_;
                print TE_G "$chr_mrna\tgff_simplifyv1.3\texon\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=exon$i\n";
              }
           }
           #intron
           if(exists $index->{"intron"}){
              my @intron=@{&segSort_strand($index->{"intron"},$strand_mrna)};
              my $i;
              foreach(@intron){
                $i++;
                my($start,$end)=split/-/,$_;
                print TE_G "$chr_mrna\tgff_simplifyv1.3\tintron\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=intron$i\n";
              }
           }
          #5UTR 3UTR
          if(exists $index->{"five_prime_UTR"}){
             my @box=@{&segSort_strand($index->{"five_prime_UTR"},$strand_mrna)};
             my $i;
             foreach(@box){
               $i++;
               my($start,$end)=split/-/,$_;
               print TE_G "$chr_mrna\tgff_simplifyv1.3\tfive_prime_UTR\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=utr5_$i\n";
             }
          }
          if(exists $index->{"three_prime_UTR"}){
             my @box=@{&segSort_strand($index->{"three_prime_UTR"},$strand_mrna)};
             my $i;
             foreach(@box){
               $i++;
               my($start,$end)=split/-/,$_;
               print TE_G "$chr_mrna\tgff_simplifyv1.3\tthree_prime_UTR\t$start\t$end\t.\t$strand_mrna\t.\tParent=$alias_mrna;Note=utr3_$i\n";
             }
          }
       }#foreach TE mrna end
       next;        
     } #if gene te end
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
      ($s1,$e1)=split/-/,$a;
      ($s2,$e2)=split/-/,$b;
      if($s1 > $e1 || $s2 > $e2){die "start > end at $s1 > $e1 || $s2 > $e2"}
      if($s1 < $s2 && $e1 > $s2 || $s2 < $s1 && $e2 > $s1){die "overlap found at $a,$b"}  #unOverlap
      if($e2<=$s1){1}elsif($e1<=$s2){-1}else{die "unknown circumstance at $s1,$e1:$s2,$e2\n"}
    } @segs;
# print "@CDS\n@sorted\n";
 return \@sorted;
}


sub segSort_strand(){
 #sort with s < e and nonOverlap check
 my @segs=@{shift @_}; #s-e-frameshift or s-e
 my $strand = shift;
 my @sorted=sort{
      my ($s1,$e1,$s2,$e2);#$s1<=$e1 && $s2<=$e2 and deOverlaped beforehand by default
      ($s1,$e1)=split/-/,$a;
      ($s2,$e2)=split/-/,$b;
      if($s1 > $e1 || $s2 > $e2){die "start > end at $s1 > $e1 || $s2 > $e2"}
      if($s1 < $s2 && $e1 > $s2 || $s2 < $s1 && $e2 > $s1){die "overlap found at $a,$b"}  #unOverlap
      if($strand eq "+"){
         if($e2<=$s1){1}elsif($e1<=$s2){-1}else{die "unknown circumstance at $s1,$e1:$s2,$e2\n"}
      }elsif($strand eq "-"){ if($e2<=$s1){-1}elsif($e1<=$s2){1}else{die "unknown circumstance at $s1,$e1:$s2,$e2\n"} }else{die"unknow strand:$strand"}
    } @segs;
# print "@CDS\n@sorted\n";
 return \@sorted;
}


sub getExon(){
    #calculate exon regions from strand sorted CDS, 5UTR, 3UTR of a given mrna (exon=utr+cds)
    #check mrna_range and exon range
    #check +/- strand, if + strand, seg order incremental, decremental on the contrary
    #output exon order idx
    #output exon annotation details

    ###+ strand
    #  mrna_range       ------------------------------------------------------------
    #     CDSs                     ------     --------------  ---    ----
    #     UTRs          ----   ----                                      ------   --
    #                  5utr_1   2                                        3utr_1    2
    #     exons         ----   ----------     --------------  ---    ----------   --   
    #splicing_out_introns   \ /          \   /              \/   \  /          \ /
    #                     1         2                3         4          5        6

    ###- strand          
    #  mrna_range       ------------------------------------------------------------
    #     CDSs                     ------     --------------  ---    ----
    #     UTRs          ----   ----                                      ------   --
    #                  3utr_2   1                                        5utr_2    1
    #     exons         ----   ----------     --------------  ---    ----------   --   
    #splicing_out_introns   \ /          \   /              \/   \  /          \ /
    #                     6         5                4         3          2        1
    #                        
    my ($mrna_range,$cds,$utr5,$utr3,$strand)=@_;
    my ($mrna_s, $mrna_e)=split/-/,$mrna_range;
    die "mrna start > end in getExon function\n" if($mrna_s >= $mrna_e);
    die "empty cds" if(scalar @{$cds} == 0);
 
    ##1 collect segs
    my @exon;
    if($strand eq "+"){
      #utr5, ok if has no utr5
      if(scalar @$utr5 > 1){
        die "Incremental check failed in utr5" if(!&isIncremental($utr5));
        for my$i(0..$#{$utr5}){
          my ($s,$e) = split/-/, $utr5->[$i];
          my $j=$i+1;
          push @exon,"$s-$e|5utr$j";
        }#for end 
      }elsif(scalar @$utr5 == 1){
           my ($s,$e) = split/-/, $utr5->[0]; 
           push @exon,"$s-$e|5utr1";
         }

      #cds, exit if no cds found
      if(scalar @$cds > 1){
        die "Incremental check failed in cds" if(!&isIncremental($cds));
        for my$i(0..$#{$cds}){
           my ($s,$e) = split/-/, $cds->[$i];
           my $j=$i+1;
           push @exon,"$s-$e|cds$j";
        }#for end 
      }elsif(scalar @$cds == 1){
           my ($s,$e) = split/-/, $cds->[0];
           push @exon,"$s-$e|cds1";
         }else{die"empty cds:@$cds"}

      #utr3, ok if no utr3 found
      if(scalar @$utr3 > 1){
        die "Incremental check failed in utr3" if(!&isIncremental($utr3));
        for my$i(0..$#{$utr3}){
           my ($s,$e) = split/-/, $utr3->[$i];
           my $j=$i+1;
           push @exon,"$s-$e|3utr$j";
        }#for end 
      }elsif(scalar @$utr3 == 1){
           my ($s,$e) = split/-/, $utr3->[0];
           push @exon,"$s-$e|3utr1";
         }
    }elsif($strand eq "-"){
      #utr5, ok if has no utr5
      if(scalar @$utr5 > 1){
        die "Decremental check failed in utr5" if(!&isDecremental($utr5));
        for my$i(0..$#{$utr5}){
          my ($s,$e) = split/-/, $utr5->[$i];
          my $j=$i+1;
          push @exon,"$s-$e|5utr$j";
        }#for end 
      }elsif(scalar @$utr5 == 1){
           my ($s,$e) = split/-/, $utr5->[0]; 
           push @exon,"$s-$e|5utr1";
         }

      #cds, exit if no cds found
      if(scalar @$cds > 1){
        die "Decremental check failed in cds" if(!&isDecremental($cds));
        for my$i(0..$#{$cds}){
           my ($s,$e) = split/-/, $cds->[$i];
           my $j=$i+1;
           push @exon,"$s-$e|cds$j";
        }#for end 
      }elsif(scalar @$cds == 1){
           my ($s,$e) = split/-/, $cds->[0];
           push @exon,"$s-$e|cds1";
         }else{die"empty cds:@$cds"}

      #utr3, ok if no utr3 found
      if(scalar @$utr3 > 1){
        die "Decremental check failed in utr3" if(!&isDecremental($utr3));
        for my$i(0..$#{$utr3}){
           my ($s,$e) = split/-/, $utr3->[$i];
           my $j=$i+1;
           push @exon,"$s-$e|3utr$j";
        }#for end 
      }elsif(scalar @$utr3 == 1){
           my ($s,$e) = split/-/, $utr3->[0];
           push @exon,"$s-$e|3utr1";
         }

       }else{die "unknow strand: $strand"}


   ##2 join segs into  exon by 1bp overlap and  record structure detail
   my $exon_linked=&segLink_plus(\@exon,1,$strand);


   ##3 check mrna range and exon range
   if($strand eq "+"){
     my ($exon_s,undef,undef)=split/-|\|/,$exon_linked->[0];
     my (undef,$exon_e,undef)=split/-|\|/,$exon_linked->[-1];
     if($exon_s != $mrna_s || $exon_e != $mrna_e){ die "exon start or end != mrna start or end: $exon_s, $exon_e != $mrna_s, $mrna_e"}
   }elsif($strand eq "-"){
      my (undef,$exon_e,undef)=split/-|\|/,$exon_linked->[0];
      my ($exon_s,undef,undef)=split/-|\|/,$exon_linked->[-1];
      if($exon_s != $mrna_s || $exon_e != $mrna_e){ die "exon start or end != mrna start or end: $exon_s, $exon_e != $mrna_s, $mrna_e"}   
    }else{die "unknow strand:$strand"}

   return $exon_linked;
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

sub getIntron_strand_old(){
  #input exon segs already sorted by strand
  my @box=@{shift @_};
  my $strand = shift;
  my @intron;
  for(0..$#box-1){
    my ($s1,$e1)=split/-|\|/,$box[$_];
    my ($s2,$e2)=split/-|\|/,$box[$_+1];
    if($strand eq "+"){push @intron,($e1+1)."-".($s2-1)}elsif($strand eq "-"){push @intron,($e2+1)."-".($s1-1)}else{die "unknow strand:$strand"}
  }
  return \@intron;
}

sub getIntron_strand(){
  #input exon segs already sorted by strand
  my @box=@{shift @_};
  my $strand = shift;
  my @intron;
  for(0..$#box-1){
    my ($s1,$e1)=split/-|\|/,$box[$_];
    my ($s2,$e2)=split/-|\|/,$box[$_+1];
    if($strand eq "+"){
      if($e1-$s2 >= -2){print STDERR "exon too near or overlap at $s1,$e1 and $s2,$e2\n";next}
      push @intron,($e1+1)."-".($s2-1)
    }elsif($strand eq "-"){
       if($e2-$s1 >= -2){print STDERR "exon too near or overlap at $s2,$e2 and $s1,$e1\n";next}
       push @intron,($e2+1)."-".($s1-1)
     }else{die "unknow strand:$strand"}
    #if($strand eq "+"){push @intron,($e1+1)."-".($s2-1)}elsif($strand eq "-"){push @intron,($e2+1)."-".($s1-1)}else{die "unknow strand:$strand"}
  }
  return \@intron;
}


sub isIncremental(){
   #non overlap and incremental order   
   #must has at least two segs
   my $segs = shift;
   my $flag = 1;
   for my$i(0..$#{$segs}-1){
     my ($s1,$e1) = split/-/, $segs->[$i];
     my ($s2,$e2) = split/-/, $segs->[$i+1];
     if($e1 < $s2){}else{$flag = 0};
   }
   return $flag;
}


sub isDecremental(){
   #non overlap and decremental order   
   #must has at least two segs
   my $segs = shift;
   my $flag = 1;
   for my$i(0..$#{$segs}-1){
     my ($s1,$e1) = split/-/, $segs->[$i];
     my ($s2,$e2) = split/-/, $segs->[$i+1];
     if($e2 < $s1){}else{$flag = 0};
   }
   return $flag;
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


sub segLink_plus(){
 #link small segments to longer ones with gap<=$gapSize, consider strand, record link details
    my ($index,$gapSize,$strand)=@_;
    #my @segments=@{&segSort($index)};#sorted
    my @segments = @{$index};
    my @segments_linked;
    for (my $i=1;$i<=$#segments;$i++){
        my ($s1,$e1,$anno1)=split /-|\|/,$segments[$i-1];
        my ($s2,$e2,$anno2)=split /-|\|/,$segments[$i];
        if($strand eq "+"){
          if(($s2-$e1-1)<=$gapSize){
            my $start=$s1;
            my $end=$e2;
            my $anno = $anno1.":".$anno2;
            $segments[$i]=$start."-".$end."|".$anno;
          }else{push @segments_linked, $segments[$i-1]}
        }elsif($strand eq "-"){
              if(($s1-$e2-1)<=$gapSize){
                my $start=$s2;
                my $end=$e1;
                my $anno = $anno1.":".$anno2;
                $segments[$i]=$start."-".$end."|".$anno;
              }else{push @segments_linked, $segments[$i-1]}

          }else{ die "unknow strand in segLink_plus function:$strand"}

        if($i==$#segments){push @segments_linked,$segments[$i];}
    }
    if($#segments==0){push @segments_linked,$segments[0];}
    return(\@segments_linked);
}




