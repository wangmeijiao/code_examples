use strict;
use warnings;
use Data::Dumper;


# mfa aln call SNP (mismatch) and InDel (alignment gap -)
# only two species

#need debug for snp detect on chr03:22545355-22545376, which align gap in Chr03



  my ($file,$order) = ($ARGV[0],$ARGV[1]);
  $order ||= "japo-brac";
  my @order = split/-/,$order;
  die "order empty: $order" if(scalar @order == 0);
  #my %maf;
  my $cnt=0;
  my $cnt_snp = 0;
  my $cnt_indel = 0;
  open MFA, $file or die "$!";
  open SNP, ">mut.$file.snp" or die "$!"; #mismatch S/V type
  open INDEL1, ">mut.$file.InDel1" or die "$!"; #deletion in ref
  open INDEL2, ">mut.$file.InDel2" or die "$!"; #deletion in que

  $/ = "\n\n";
  while(<MFA>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
    my @temp = split/>/,$_;
    shift @temp;
    die "not two seqs at $_" if(scalar @temp != 2);
    my %align;
    #my @order;
    foreach my $ctrl(@temp){
      my @box=split/\n+/,$ctrl;
      my $head=shift @box;
      my @tmp = split/[\t ]+/,$head;
      my $id;
      if($tmp[0] =~/^chr/){$id = "japo"}elsif($tmp[0] =~/^Chr/){$id = "brac"}else{die "unknow spec at $_"}
      #push @order,$id;
      my $align=join"",@box;
      $align{$id}->{'aln'} = $align;
      $align{$id}->{'alnLen'} = $tmp[1];
      $align{$id}->{'gene'} = $tmp[2];
      $align{$id}->{'peak'} = $tmp[3];
      #link multi-orth-regions
      if($tmp[4]=~/\|/){
        print STDERR "link multi-orth-regions when readin mfa at $tmp[2] with orth $tmp[4]\n";
        my @orth = split/\|/,$tmp[4];
        my $result_link = &segLink_all(\@orth);print STDERR "@$result_link\n";
        if(scalar @$result_link == 1){
           my ($chr,$range,$strand) = split/:/,$result_link->[0];
           my ($start,$end) = split/-/,$range;
           my @orth = ("$chr","$start","$end","$strand");
           $align{$id}->{'orth'} = \@orth;
        }else{die "link failed at $tmp[2] with orth $tmp[4]"}
      }else{
             my ($chr,$range,$strand) = split/:/,$tmp[4];
             my ($start,$end) = split/-/,$range;
             my @orth = ("$chr","$start","$end","$strand");
             $align{$id}->{'orth'} = \@orth;
           }


    }
  
    $cnt++;
    #last if($cnt == 10000);
    #print STDERR "#" if($cnt % 1 ==0);
    
    my %aln;
    foreach my $spec(@order){$aln{$spec}=$align{$spec}->{'aln'}}
    
    my ($flag,$len)=&checkAlignLenAll(\%aln); 
    if($flag != 1){die "len diffs"}

    my ($densi_score, $cons_seq, $snp, $flagN, $flagGap)=&getSimilarDensi(\%aln,\@order,$len,"score");
    my $stars=&densi2star($densi_score,0.75*2);
    #print Dumper $flag_gaps, $densi_score, $cons_seq,$stars;
#    foreach my $id(@order){printf("%-20s  %s\n",$id,$align{$id})}
#    printf "%-20s  %s\n"," ",join("",@{$stars});
#    printf "%-20s  %s\n"," ",join("",@{$densi_score});
#    printf "%-20s  %s\n"," ",join("",@{$cons_seq});
#    printf "%-20s  %s\n"," ",join("",@{$snp}); #snp
#    printf "%-20s  %s\n"," ",join("",@{$flagGap}); #InDel
#    printf "%-20s  %s\n"," ",join("",@{$flagN});
#    print "#all equal#\n" if($flag == 1);


    #annotation and output as bed

    #snp 
    #$name,$start,$end,$strand,$seq_ali,$seq_ori
    my $chr_ref = $align{'japo'}->{'orth'}->[0];
    my $start_ref = $align{'japo'}->{'orth'}->[1];
    my $strand_ref = $align{'japo'}->{'orth'}->[3];
    my $chr_que = $align{'brac'}->{'orth'}->[0];
    my $start_que = $align{'brac'}->{'orth'}->[1];
    my $strand_que = $align{'brac'}->{'orth'}->[3];
    my $cnt_base_ref = 0;
    my $cnt_base_que = 0;
    my $align_ref = $aln{'japo'};
    my $align_que = $aln{'brac'};

    print SNP ">block$cnt\n";
    print SNP "#",join("|",@{$align{'japo'}->{'orth'}})," $align{'japo'}->{'gene'}"," $align{'japo'}->{'peak'}\n";
    print SNP "#",join("|",@{$align{'brac'}->{'orth'}})," $align{'brac'}->{'gene'}"," $align{'brac'}->{'peak'}\n";
    for my $i(0..$#{$snp}){
       my $base = substr($align_ref,$i,1);
       if($base=~/[ATGCN]/){$cnt_base_ref++}
       $base = substr($align_que,$i,1);
       if($base=~/[ATGCN]/){$cnt_base_que++}
       #if( ($flagN->[$i] != 1 && $flagGap->[$i] != 1 && $flagGap->[$i] != 3)  || ($flagN->[$i] == 1 && $flagGap->[$i] != 1 && $flagGap->[$i] != 3) ){$cnt_base_ref++}  
       #if($flagGap->[$i] != 2 ){$cnt_base_que++}  

       if($snp->[$i] eq "S" || $snp->[$i] eq "V"){
          $cnt_snp++;
          my $start_ref_real = $start_ref + $cnt_base_ref;
          my $end_ref_real = $start_ref_real + 1;
          my $start_que_real = $start_que + $cnt_base_que - 1;
          my $end_que_real = $start_que_real + 1;
          #my $end_que_real = $start_que_real + 1; #comment out 
          print SNP "$chr_ref\t$start_ref_real\t$end_ref_real\t$strand_ref\t$snp->[$i]\tblock${cnt}_snp$cnt_snp\t$chr_que\t$start_que_real\t$end_que_real\t$strand_que\t$snp->[$i]\tblock${cnt}_snp$cnt_snp\n";
       }
    }    
    print SNP "\n\n";

    #InDel
    print INDEL1 ">block$cnt\n";
    print INDEL1 "#",join("|",@{$align{'japo'}->{'orth'}})," $align{'japo'}->{'gene'}"," $align{'japo'}->{'peak'}\n";
    print INDEL1 "#",join("|",@{$align{'brac'}->{'orth'}})," $align{'brac'}->{'gene'}"," $align{'brac'}->{'peak'}\n";
    my($coverage,$segs1,$segs2)=&statBlocks($flagGap);  
    my $cnt_indel1 = 0;
    foreach my $seg(@$segs1){#transform align coordinates to real 
        my ($align_s,$align_e) = split/-/,$seg;
        my $len_indel = $align_e-$align_s+1;
        $cnt_indel1++;
        my $cnt_base_ref_in = 0;
        my $cnt_base_que_in = 0;
        for my $i(0..$align_s){ #get start real
          my $base = substr($align_ref,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_ref_in++}
          $base = substr($align_que,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_que_in++}
        }
        my $start_ref_real = $start_ref + $cnt_base_ref_in;
        my $start_que_real = $start_que + $cnt_base_que_in;


        for my $i($align_s..$align_e){#get end real
          my $base = substr($align_ref,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_ref_in++}
          $base = substr($align_que,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_que_in++}
        }
        my $end_ref_real = $start_ref + $cnt_base_ref_in;
        my $end_que_real = $start_que + $cnt_base_que_in;

        print INDEL1 "$chr_ref\t$start_ref_real\t$end_ref_real\t$strand_ref\t$seg:$len_indel\tblock${cnt}_InDel$cnt_indel1\t$chr_que\t$start_que_real\t$end_que_real\t$strand_que\t$seg:$len_indel\tblock${cnt}_InDel$cnt_indel1\n";

    }
    print INDEL1 "\n\n";

    print INDEL2 ">block$cnt\n";
    print INDEL2 "#",join("|",@{$align{'japo'}->{'orth'}})," $align{'japo'}->{'gene'}"," $align{'japo'}->{'peak'}\n";
    print INDEL2 "#",join("|",@{$align{'brac'}->{'orth'}})," $align{'brac'}->{'gene'}"," $align{'brac'}->{'peak'}\n";
    my $cnt_indel2 = 0;
    foreach my $seg(@$segs2){#transform align coordinates to real 
        my ($align_s,$align_e,$seg2_cut) = split/-|:/,$seg;
        #print "$seg\t$seg2_cut\n";
        my $len_indel = $align_e-$align_s+1;
        $cnt_indel2++;
        my $cnt_base_ref_in = 0;
        my $cnt_base_que_in = 0;
        for my $i(0..$align_s){ #get start real
          my $base = substr($align_ref,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_ref_in++}
          $base = substr($align_que,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_que_in++}
        }
        my $start_ref_real = $start_ref + $cnt_base_ref_in;
        my $start_que_real = $start_que + $cnt_base_que_in;


        for my $i($align_s..$align_e){#get end real
          my $base = substr($align_ref,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_ref_in++}
          $base = substr($align_que,$i,1);
          if($base=~/[ATGCN]/){$cnt_base_que_in++}
        }
        my $end_ref_real = $start_ref + $cnt_base_ref_in;
        my $end_que_real = $start_que + $cnt_base_que_in;

        print INDEL2 "$chr_ref\t$start_ref_real\t$end_ref_real\t$strand_ref\t$align_s-$align_e:$len_indel:$seg2_cut\tblock${cnt}_InDel$cnt_indel2\t$chr_que\t$start_que_real\t$end_que_real\t$strand_que\t$align_s-$align_e:$len_indel:$seg2_cut\tblock${cnt}_InDel$cnt_indel2\n";

    }
     print  INDEL2 "\n\n";


    #last;
  }#while end
  close MFA;
  close SNP;
  close INDEL1;
  close INDEL2;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks for $order\n";




##sub


sub segLink_all(){
  #link and delete  original overlap segs
  my $idx = shift;
  my $gap ||= 5;
  my ($chr,$strand);
  die "segs empty when link orth" if(scalar @$idx == 0);
  die "segs == 1 when link orth" if(scalar @$idx == 1);
  my @segs_sort = sort { #sort and check
                         my ($chr1,$range1,$strand1) = split/:/,$a;
                         my ($s1,$e1) = split/-/,$range1;
                         my ($chr2,$range2,$strand2) = split/:/,$b;
                         my ($s2,$e2) = split/-/,$range2;
                         die "chr diffs when link orth $a and $b" if($chr1 ne $chr2);
                         die "strand diffs when link orth $a and $b" if($strand1 ne $strand2);
                         die "s1 > e1 when link orth $a and $b" if($s1 > $e1);
                         die "s2 > e2 when link orth $a and $b" if($s2 > $e2);
                         if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                         #$chr = $chr1;print STDERR "see??? $chr\n"; #not valid in sort {}
                         #$strand = $strand1;
                       } @$idx;
  for(my $i = 0;$i<= $#segs_sort-1;$i++){
     my ($chr1,$range1,$strand1) = split/:/,$segs_sort[$i];
     my ($s1,$e1) = split/-/,$range1;
     my ($chr2,$range2,$strand2) = split/:/,$segs_sort[$i+1];
     my ($s2,$e2) = split/-/,$range2;
     if(!defined $chr && !defined $strand){$chr = $chr1; $strand= $strand1}else{if($chr1 eq $chr && $strand1 eq $strand){}else{die "chr or strand diffs when link orth $segs_sort[$i] and $segs_sort[$i+1]"}    }
     if($s2 - $e1 <= $gap){$segs_sort[$i] = "NA";$segs_sort[$i+1] = "$chr:$s1-$e2:$strand"}else{die "broken seg when link orth at @segs_sort"}

  }

  my @seg_link;
  foreach my $ctrl(@segs_sort){if($ctrl ne "NA"){push @seg_link,$ctrl}}
  return \@seg_link;

}



sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=length($index->{$id})}else{if($len != length($index->{$id}) ){$flag=0}}
  }
  return ($flag,$len);
}


sub getSimilarDensi(){
  my ($index,$order,$len,$type)=@_;
  my @densi; #use score or perc;
  my @cons;
  my @Tsnp;
  my @flag_gaps; #len == aligned_len ;
  my @flag_N; #len == aligned_len ;
  my %align;
  my @ids;
  foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align{$id}=\@box;
      push @ids,$id;
  }
  for my $i (0.. $len-1){
     my @elements;
     foreach my $id(@$order){
       push @elements, $align{$id}->[$i];
     }
     my($score,$cons, $Tsnp,$flag_N, $flag_gap)=&checkIdentity_pair(\@elements);
     #my($flag, $score, $perc, $cons)=&checkIdentity_densi(\@elements);
     #if($type eq "score"){push @densi, $score}elsif($type eq "perc"){push @densi, $perc}else{die "don't understand densi type $type"}
     push @densi, $score;
     push @cons, $cons;
     push @Tsnp,$Tsnp;
     push @flag_gaps,$flag_gap;
     push @flag_N, $flag_N;
  }
  return(\@densi, \@cons, \@Tsnp, \@flag_N, \@flag_gaps ) ;
}

sub densi2star(){
  my ($densi,$cutoff)=@_;
  my @stars;
  my $len=scalar @{$densi};
  for(0..$len-1){
    die "undefed value in @{$densi}" if(!defined $densi->[$_]);
    if($densi->[$_] >= $cutoff){push @stars, "*"}else{push @stars, " "}
  }
  return \@stars;
}

sub checkIdentity_pair(){
    #only compare two species
    #upcase and check all then compare    
    #use bracRS2 as ref (has very few N, only ATGC-)
    #return 2(all identical), 0 (not identical)
    #return flag_gap flag_N cons_base score
    my @elements=@{shift @_};
    die "more than 2 or empty at @elements" if(scalar @elements != 2);
    my ($ref,$que) = @elements;
    $ref = uc($ref);
    $que = uc($que);
    if($ref=~/[^-ATGCN]/i){die "illegal character found in ref $ref"} #nucleotide only
    if($que=~/[^-ATGCN]/i){die "illegal character found in que $que"} 

    die "both gap -" if($ref eq "-" && $que eq "-");
    die "both N " if($ref eq "N" && $que eq "N");

    my ($score,$cons,$Tsnp,$flag_N,$flag_gap);
    # score: 0(not match),2(match)
    # Tsnp : S(A<->G or C<->T transition,Ts),V(other,transversion,Tv), relate to ref
    # flag_N: 0(none) , 1 (ref), 2(que)
    # flag_gap: 0(none) , 1 (ref), 2(que), 3(both N or -)
    $Tsnp = " ";
    if($ref=~/[ATGC]/){ 
       if($ref eq $que){ $score = 2; $cons=$ref; $flag_N = 0; $flag_gap = 0 } #identical and not N or -
       elsif($ref ne $que && $que eq "-"){  # not identical and is  -, #InDel, N similar to -
         $score = 0; $cons=" "; $flag_N = 0; $flag_gap = 2             
       }elsif( $ref ne $que && $que eq "N"){ # not identical and is  N
          $score = 0; $cons=" "; $flag_N = 2; $flag_gap = 2
        }elsif($ref ne $que && $que=~/[ATGC]/){ # not identical and is ATCG, SNP
           $score = 0; $cons=" "; if($ref eq "A" && $que eq "G" || $ref eq "G" && $que eq "A" || $ref eq "C" && $que eq "T" || $ref eq "T" && $que eq "C"){$Tsnp = "S"}else{$Tsnp = "V"}; $flag_N = 0; $flag_gap = 0 
         }else{die "unknow situation at ref:$ref and que:$que"}

    }elsif($ref eq "-"){ 
       if($ref eq $que){ die "both gap -"}
       elsif($ref ne $que && $que=~/[ATGC]/){  # not identical and is  ATCG
         $score = 0; $cons=" "; $flag_N = 0; $flag_gap = 1
       }elsif( $ref ne $que && $que eq "N"){ # not identical and is  N
          $score = 0; $cons=" "; $flag_N = 2; $flag_gap = 3
         }else{die "unknow situation at ref:$ref and que:$que"}

     }elsif($ref eq "N"){ #very few
       if($ref eq $que){ die "both gap N"}
       elsif($ref ne $que && $que=~/[ATGC]/){  # not identical and is  ATCG
         $score = 0; $cons=" "; $flag_N = 1; $flag_gap = 1
       }elsif( $ref ne $que && $que eq "-"){ # not identical and is  -
          $score = 0; $cons=" "; $flag_N = 1; $flag_gap = 3
         }else{die "unknow situation at ref:$ref and que:$que"}

      }else{die "ref err: $ref:$que"} #should not happen

    return ($score,$cons,$Tsnp,$flag_N,$flag_gap);
}



sub statBlocks(){
    # 0,1,2,3
    # link 1; link 2/3
    my $star=shift;
    my $len=scalar @{$star};
    my ($coverage,$coverage_total)=(0,0);
    my $gap=4; #for link nearby similar blocks

    #find continual segments for ref gap (value = 1, most are -, N very few)
    my @segment1;
    my ($cnt_0,$cnt_1,$cnt_2,$cnt_3)=(0,0,0,0);
    my ($start,$end,$cnt)=(0,0,0);
    for my $i(0..$len-1){
      if($star->[$i] == 1){$cnt_1++}elsif($star->[$i] == 0){$cnt_0++}elsif($star->[$i] == 2){$cnt_2++}elsif($star->[$i] == 3){$cnt_3++}     
      if($star->[$i] == 1){
         $cnt++;
         if($i==($len-1)){ #last element
           $start=$i-$cnt+1;
           $end=$i;
           push @segment1,$start."-".$end;
           $cnt=0;
         }
      }else{
            $start=$i-$cnt;
            $end=$i-1;
            if($start<=$end){push @segment1,$start."-".$end}
            $cnt=0;
           }
    }#for end
    my @seg1=@segment1; #@segments changed in next steps
    #print "@segs\n";

    #find continual segments for ref and que gap (value = 2,3)
    my @segment2;
    ($cnt_0,$cnt_1,$cnt_2,$cnt_3)=(0,0,0,0);
    my $cnt_3_inner = 0;
    ($start,$end,$cnt)=(0,0,0);
    for my $i(0..$len-1){
      if($star->[$i] == 1){$cnt_1++}elsif($star->[$i] == 0){$cnt_0++}elsif($star->[$i] == 2){$cnt_2++}elsif($star->[$i] == 3){$cnt_3++}
      if($star->[$i] == 2 || $star->[$i] == 3){
         $cnt++;
         if($star->[$i] == 3){$cnt_3_inner++}
         if($i==($len-1)){ #last element
           $start=$i-$cnt+1;
           $end=$i;
           push @segment2,$start."-".$end.":$cnt_3_inner";
           $cnt=0;
           $cnt_3_inner = 0
         }
      }else{
            $start=$i-$cnt;
            $end=$i-1;
            if($start<=$end){push @segment2,$start."-".$end.":$cnt_3_inner"}
            $cnt=0;
            $cnt_3_inner = 0
           }
    }#for end
    my @seg2=@segment2; #@segments changed in next steps
    #print "@segs\n";

=pod
    #link small segments to long ones with gap<=4
    my @segments_linked;
    for (my $i=1;$i<=$#segments;$i++){
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        if(($tmp2[0]-$tmp1[1]-1)<= $gap){
        $start=$tmp1[0];
        $end=$tmp2[1];
        $segments[$i]=$start."-".$end;
        }else{push @segments_linked, $segments[$i-1];}
        if($i==$#segments){push @segments_linked,$segments[$i];}
    }
    if($#segments==0){push @segments_linked,$segments[0];}
    #print "@segments_linked\n";
=cut

    #start to calculate percentage etc
    $coverage=sprintf "%.3f",($cnt_1+$cnt_2+$cnt_3)/$len; 


    #return ($coverage,$identity,\@segs,\@segments_linked);
    return ($coverage,\@seg1,\@seg2);
}




