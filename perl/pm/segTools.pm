
sub segSortCheck_fast(){
  #code block to sort by start and check overlap at the same time
  my $segs = shift;
  my @segs_sorted = sort{
                          my ($s1,$e1) = split/-/,$a;
                          die "s1 > e1 in segs sort" if($s1 > $e1);
                          my ($s2,$e2) = split/-/,$b;
                          die "s2 > e2 in segs sort" if($s2 > $e2);
                          if($e1 <= $s2 || $e2<=$s1){}else{die "overlap segs $a,$b"}
                          if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                        }@$segs;


  return \@segs_sorted;
}



sub segCheckOverlap(){
   #Jan. 5, 2017; need relate2seg_detailed and segSort_start
   #function to only check segments for overlap, which belong to a single sequence; return 1 if overlap found else 0
   #use 3seg method if (seg_num >=3) or if (seg_num ==2 ) use 2seg method . do nothing if (seg_num <2)
   #input format @segs=("s1..e1","s2..e2","x..x")
   my @segs=@{shift @_};
   die "emtpy segs:@segs" if(scalar @segs == 0);
   if (scalar @segs<2){
     my ($s,$e) = split/-/,$segs[0];
     die "segs has one segment but s >= e: $s >= $e" if($s >= $e);
     return 0;
    }elsif(scalar @segs ==2){
       my @segs_sorted=@{&segSort_start(\@segs)};
       my $result = &relate2seg_detailed($segs_sorted[0],$segs_sorted[1]);
       if($result == -1){return 0}else{return 1}
      }elsif(scalar @segs >2 ){
         my @segs_sorted=@{&segSort_start(\@segs)};
         #print "sorted: @segs_sorted\n";
         my $flag = 0;
         for(1..$#segs_sorted-1){
             my $result1=&relate2seg_detailed($segs_sorted[$_-1],$segs_sorted[$_]);
             my $result2=&relate2seg_detailed($segs_sorted[$_],$segs_sorted[$_+1]);
             if($result1 == -1 && $result2 == -1){}else{$flag = 1;die  "overlap at $segs_sorted[$_-1],$segs_sorted[$_],$segs_sorted[$_+1]\n"}
         }
         if($flag == 1){return 1}else{return 0}
       }
}

sub segDeDup_checkOverlap(){
  #input raw bed3 segs, dedup, sort & checkOverlap, then output single, nonoverlap bed3 segs
  #input format @segs=("s1-e1","s2-e2","..-..")
  my @segs=@{shift @_};
  #dedup and SortByStart
  my %temp;
  foreach (@segs){$temp{$_}++}
  my @segs_deDup_sortByS=sort {
                               my($s1,$e1)=split/-/,$a;
                               my($s2,$e2)=split/-/,$b;
                               if($s1<$s2){return -1}elsif($s1>$s2){return 1}elsif($s1==$s2){return 0}
                         } keys %temp;
  #check nonOverlap
  for(0..$#segs_deDup_sortByS-1){
    my($s1,$e1)=split/-/,$segs_deDup_sortByS[$_];
    my($s2,$e2)=split/-/,$segs_deDup_sortByS[$_+1];
    if($s1<=$s2){if($e1>=$s2){die"Overlap segs found at $s1-$e1|$s2-$e2\n"}}else{die "S1 > S2 at $s1-$e1|$s2-$e2\n"}
  }
  #output 3 columns segs
  foreach(@segs_deDup_sortByS){my($s,$e)=split/-/,$_;push @result, "$s-$e"}
  return \@result;   
}


sub segDeOverlap(){
   #standard function to sort&link(deOverlap) segments, which belong to a single sequence
   #use 3seg method if (seg_num >=3) or if (seg_num ==2 ) use 2seg method . do nothing if (seg_num <2)
   #input format @segs=("s1-e1","s2-e2","..-..")
   my @segs=@{shift @_};
   if (scalar @segs<2){
     return \@segs
    }elsif(scalar @segs ==2){
       my @segs_sorted=@{&segSort(\@segs)};
       my @segs_sorted_linked=&relate2segs($segs_sorted[0],$segs_sorted[1]);
       return \@segs_sorted_linked;
      }elsif(scalar @segs >2 ){
         my @segs_sorted=@{&segSort(\@segs)};
         #print "sorted: @segs_sorted\n";
         my @result=shift @segs_sorted;
         for(0..$#segs_sorted){
             my @relate;
             if($_ != $#segs_sorted){
               @relate=&relate2segs($segs_sorted[$_],$segs_sorted[$_+1]);
               @relate=&relate2segs($result[-1],$relate[0]);
               pop @result;
               push @result,@relate;
             }elsif($_ == $#segs_sorted){
                   @relate=&relate2segs($result[-1],$segs_sorted[$_]);
                   pop @result;
                   push @result,@relate;
                 }
         }
         return \@result;
       }
}




sub segSort_start(){ #simple sort by start
     #start < end
     my @segs=@{shift @_};
     my @segs_sorted=sort{
       my ($s1,$e1)=split/-/,$a;
       my ($s2,$e2)=split/-/,$b;
       if($s1<$s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
     } @segs;
     return \@segs_sorted;
}


sub segSort_withCheck(){
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





sub relate2segs(){
    #input paras ("s1-e1","s2-e2") must be sorted by start
    my ($s1,$e1)=split/-/,(shift @_);
    my ($s2,$e2)=split/-/,(shift @_);
    die "start1>start2 at $s1-$e1|$s2-$e2\n" if ($s1>$s2);
    #three types of relationship between two segments
    #1,separated
    if($s2>$e1){
       return ("$s1-$e1","$s2-$e2");
     #2,overlap
     }elsif($s2<=$e1 && $e1<=$e2){
        return ("$s1-$e2")
      #3,included
      }elsif($e2<$e1 && $e2<=$e1){
         return ("$s1-$e1");
        }else{print STDERR "unknown relationship at $s1-$e1|$s2-$e2, return original\n"; return ("$s1-$e1","$s2-$e2")}
}


sub relate2seg_detailed(){
  #input A and B seg as "s1-e1","s2-e2";
  #start must < end
  #A as target
  # return -1 if seperated
  #         2 if one end overlap 
  #         1 if A includes  B
  #         0 if A covered by B
  my ($s1,$e1)=split/-/,(shift @_);
  my ($s2,$e2)=split/-/,(shift @_);
  die "start > end at $s1-$e1|$s2-$e2\n" if ($s1 >= $e1 || $s2 >= $e2);
  if($e1 <= $s2 && $e1 <= $e2 || $e2 <= $s1 && $e2 <= $e1){ # seperated
     return -1
   }elsif($s2 > $s1 && $s2 < $e1 && $e1 > $s2 && $e1 > $s2 && $e1 < $e2 || $s1 > $s2 && $s1 < $e2 && $e2 > $s1 && $e2 < $e1){
      return 2 #one end overlap
    }elsif($s2 >= $s1 && $e2 <= $e1){ #A includes B
      return 1
     }elsif( $s1 >= $s2 && $e1 <= $e2){ #A covered by B
        return 0
       }else{die "unknown circurrence! $s1-$e1,$s2-$e2"}
}










sub linkCheckSegs(){
#link small segments to longer ones with gap<=$gapSize and control percH >= $percLH_TH
    my ($index,$gapSize)=@_;
    my $percH_before;
    my $percH_swallowed;
    my @segments=@{$index};
    my @segments_linked;
    my @percH=("","");#two elements only
    for (my $i=1;$i<=$#segments;$i++){#more than one segments
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        die "not ordered at @tmp1 and  @tmp2\n" if($tmp1[1] >= $tmp2[0]);

        if($i == 1){my $len=$tmp1[1]-$tmp1[0]+1;$percH[0]="$len|$len" } #percH_init
        my ($totalH0,$total0)=split/\|/,$percH[0]; #percH_before
        my $totalH1=$totalH0+($tmp2[1]-$tmp2[0]+1);
        my $total1=$total0+($tmp2[1]-$tmp1[1]+1);
        $percH[1]="$totalH1|$total1";  #percH_swallowed
        my $percLH=$totalH1/$total1;

        if(($tmp2[0]-$tmp1[1]-1)<=$gapSize && $percLH >= $percLH_TH){
          my $start=$tmp1[0];
          my $end=$tmp2[1];
          $segments[$i]=$start."-".$end;
          $percH[0]=$percH[1];$percH[1]=""; #slide percH
        }else{
              push @segments_linked, $segments[$i-1];
              my $len=$tmp2[1]-$tmp2[0]+1;
              $percH[0]="$len|$len";$percH[1]=""; #init @percH with element i
             }
        if($i==$#segments){push @segments_linked,$segments[$i]} #the last seg
    }
    if($#segments==0){push @segments_linked,$segments[0]} #only one segment
   return(\@segments_linked);
}





sub segLink(){
 #link small segments to longer ones with gap<=$gapSize
 #must make sure segs are sorted and unOverlap
    my ($index,$gapSize)=@_;
    my @segments=@{$index};
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



sub binarize(){#turn densi to 0/1 value ###very slow when apply to genomewide practice
   my ($index,$ave_g,$given)=@_;#dynamic thresholds: local lambda or local ave???  local=10k up/down region or 1/10 data points
   my @binDen;
   my $last=$#{$index};
   for my $i(0..$last){
      if($i%1000==0){print STDERR "#"}
      my @localAve;
      if(defined $index->[$i]){
         my @localAve;
         foreach my $flank_half(500,2500,5000){
           my $regionS=$i-$flank_half; # define local region as 1k -500/+500 at present data point, shift if reach ends
           my $regionE=$i+$flank_half;
           my @localData;
           my $localAve;
           if($regionS<=0){
             @localData=@{$index}[1..2*$flank_half];
             $localAve=&calcMedian(\@localData);
           }elsif($regionE>=$last){
                @localData=@{$index}[$last-2*$flank_half..$last];
                $localAve=&calcMedian(\@localData);
               }else{
                      @localData=@{$index}[$i-$flank_half..$i+$flank_half];
                      $localAve=&calcMedian(\@localData);
                    }
           push @localAve,$localAve;
        }#local region 1k, 5k, 10k ...
        my $threshold=&max(@localAve,$ave_g,$given);
        #$threshold=$localAve;
        if($index->[$i] <=$threshold){$binDen[$i]=0}else{$binDen[$i]=1}
      }else{die "err\n"}
   }
   print STDERR "\n";
   return \@binDen;
}


sub findLinkBlocks(){
    my $star=shift;
    my $gapSize=shift;
    my $win=shift;
    my $len=$#{$star};
    my @segments;
    my $cnt=0;
    #find continual segments( hump data points)
    my ($start,$end);
    for my $i(0..$len){
       if($$star[$i]==1){
          $cnt++;
          if($i==$len){
            $start=$i-$cnt+1;
            $end=$i;
            push @segments,(($start-1)*$win)."-".(($end+1)*$win);
            $cnt=0;
          }
       }else{
          $start=$i-$cnt;
          $end=$i-1;
          if($start<=$end){push @segments,(($start-1)*$win)."-".(($end+1)*$win)}
          $cnt=0;
        }
    }#for end 

    #link small segments to long ones with gap<=4
    my @segments_linked;
    for (my $i=1;$i<=$#segments;$i++){
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        if(($tmp2[0]-$tmp1[1]-1)<=$gapSize){
        $start=$tmp1[0];
        $end=$tmp2[1];
        $segments[$i]=$start."-".$end;
        }else{push @segments_linked, $segments[$i-1];}
        if($i==$#segments){push @segments_linked,$segments[$i];}
    }
    if($#segments==0){push @segments_linked,$segments[0];}

    return(\@segments,\@segments_linked);
}#sub end



sub segJoin(){
    #filter -> join neighbor  small segments to longer ones with link_gap<=$gapSize -> filter again
    #must be sorted and deoverlaped
    my ($index,$gapSize,$filterSize1,$filterSize2)=@_;
    my @segments_ori=@{$index};
    #1,filter short ones
    my @segments_filtered;
    foreach(@segments_ori){
      my($s,$e)=split/-/,$_;
      die "start >= end at $_" if($s >= $e);
      if(($e-$s+1)> $filterSize1){push @segments_filtered, "$s-$e"}
    }
    my @segments=@segments_filtered;
    #2, link by simple_neighbor_join
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

    #3, filter short ones again
    my @segments_linked_filtered;
    foreach(@segments_linked){
      my($s,$e)=split/-/,$_;
      die "start >= end at $_" if($s >= $e);
      if(($e-$s+1)> $filterSize2){push @segments_linked_filtered, "$s-$e"}
    }

    return(\@segments_linked_filtered);
}




sub segComplement(){
   #must be sorted and deoverlaped
   my ($index,$len)=@_;
   die "too few segs @{$index}" if(scalar @{$index} < 2);
   my @segs_comp;
   foreach(0..$#{$index}-1){
     my($s1,$e1)=split/-/,$index->[$_];
     my($s2,$e2)=split/-/,$index->[$_+1];
     $e1+=1;
     $s2-=1;
     push @segs_comp,"$e1-$s2";
   }
   my ($s_first,$e_first)=split/-/,$index->[0];
   my ($s_last,$e_last)=split/-/,$index->[-1];
   if($s_first != 1){$s_first-=1; unshift @segs_comp, "1-$s_first" }
   if($e_last != $len ){$e_last+=1; push @segs_comp, "$e_last-$len" }
   return \@segs_comp;
}

sub segComplement_with_startend(){
   #must be sorted and deoverlaped
   my ($index,$range)=@_;
   die "too few segs @{$index}" if(scalar @{$index} < 2);
   my ($chr,$region) = split/_/,$range;
   my ($start,$end) = split/-/,$region;
   my @segs_comp;
   foreach(0..$#{$index}-1){
     my($s1,$e1)=split/-/,$index->[$_];
     my($s2,$e2)=split/-/,$index->[$_+1];
     $e1+=1;
     $s2-=1;
     push @segs_comp,"$e1-$s2";
   }
   my ($s_first,$e_first)=split/-/,$index->[0];
   my ($s_last,$e_last)=split/-/,$index->[-1];
   if($s_first != $start){$s_first-=1; unshift @segs_comp, "$start-$s_first" }
   if($e_last != $end ){$e_last+=1; push @segs_comp, "$e_last-$end" }
   return \@segs_comp;
}




sub segCheck_overlap_quick(){
   #must be ascendly sorted by start beforehand, 
   #then check nonOverlap or not
   my @segs=@{shift @_};
   my @points;
   my $flag=0;
   foreach (@segs){my ($s,$e)=split/-/,$_;push @points,($s,$e)}
   return("err:points <=1") if(scalar @points <= 1);
   for my $i(0 .. $#points-1){
      if($points[$i] <= $points[$i+1]){}else{$flag = 1}
   }
   return $flag;
}


sub segOverlap(){
   my ($exon_coord,$region)=@_;
   my @exon_coord_new;
   my ($chr,$s,$e)=split/:|-/,$region;
   die "$s >= $e at $region" if($s >= $e);
   die "exon coord not sorted and deOverlap at @{$exon_coord}" if(&segCheck($exon_coord));
   foreach(@{$exon_coord}){
     my ($start,$end)=split/-/,$_;
     if($start >= $s && $end <= $e){
       push @exon_coord_new, "$start-$end"
      }elsif($start < $s && $end > $s && $end <= $e){
        push @exon_coord_new, "$s-$end"
       }elsif($end > $e && $start > $s && $start < $e){
          push @exon_coord_new, "$start-$e"
         }elsif($start <= $s && $end >= $e){push @exon_coord_new, "$s-$e"}
   }
   die "empty region overlap at @{$exon_coord}"if (scalar @exon_coord_new == 0);
   return \@exon_coord_new;
}

sub segSumLen(){
  my @segs=@{shift @_};
  my $sum_len=0;
  foreach(@segs){
    my ($s,$e)=split/-/,$_;
    die "s > e" if($s > $e);
    $sum_len+=($e-$s+1);
  }
  return $sum_len;
}


sub minMax_seg(){
  my $index=shift;
  die "seg empty" if(scalar @$index == 0);
  my ($min,$max)=split/-|\|/,$index->[0];
  foreach (@{$index}){
    my ($s,$e) = split/-|\|/,$_;
    die "s > e: $s > $e" if($s > $e);
    if($min>$s){$min=$s}
    if($max<$e){$max=$e}
  }
  #print ("$min,$max\n");
  return ($min,$max);
}#end of minMax sub

sub segSort_ascend(){
    #sort by start in ascend order
    my $idx = shift;    
    my @sorted = sort {
       my ($s1,$e1) = split/-|\|/,$a;
       my ($s2,$e2) = split/-|\|/,$b;
       if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
    } @$idx;
    return \@sorted;
}

sub segSort_descend(){
    #sort by start in descend order
    my $idx = shift;    
    my @sorted = sort {
       my ($s1,$e1) = split/-|\|/,$a;
       my ($s2,$e2) = split/-|\|/,$b;
       if($s1 < $s2){return 1}elsif($s1 > $s2){return -1}else{return 0}
    } @$idx;
    return \@sorted;
}

sub segCheck_overlap(){
    #no need to sort, check all segment at one time
    my $idx = shift;
    my ($min,$max) = &minMax_seg($idx);
    my @box;
    for (0..$max-$min){$box[$_] = 0}  
    foreach my $ctrl(@$idx){
      my ($s,$e) = split/-|\|/,$ctrl;
      for ($s-$min..$e-$min){$box[$_]++}
    }
    my $flag = 1;
    foreach (@box){if($_ > 1){$flag = 0}}
    return $flag;
}



sub segCheck_sort(){
     #to check: start < end; sorted (ascendly)
     my $idx = shift;
     die "seg empty" if(scalar @$idx == 0);
     if(scalar @$idx == 1){
        my ($s1,$e1) = split/-|\|/,$idx->[0];
        die "s1 > e1:$s1 > $e1" if($s1 > $e1);
        return 1
     }
     my $flag = 1;
     for my $i(0..$#$idx-1){
        my ($s1,$e1) = split/-|\|/,$idx->[$i];
        die "s1 > e1:$s1 > $e1" if($s1 > $e1);
        my ($s2,$e2) = split/-|\|/,$idx->[$i+1];
        die "s2 > e2:$s2 > $e2" if($s2 > $e2);
        if($e1 <= $s2){}else{
          #if($s1 >= $s2 || $e1 >= $e2){print STDERR "not sorted at $s1-$e1, $s2-$e2 \n"};
          #if($s1<$s2 && $e1 >= $s2){print STDERR "overlap at $s1-$e1, $s2-$e2 \n"};
          $flag = 0
        }
    }
    return $flag;
}



sub segCheck_sort_type(){
     #sorted ascend or descend or unsort
     #must nonOverlap
     my $idx = shift;
     die "seg empty" if(scalar @$idx == 0);
     if(scalar @$idx == 1){return "single"}     
     my $flag_descend = 0;
     my $flag_ascend = 0;
     for my $i(0..$#$idx-1){
        my ($s1,$e1) = split/-|\|/,$idx->[$i];
        die "s1 > e1:$s1 > $e1" if($s1 > $e1);
        my ($s2,$e2) = split/-|\|/,$idx->[$i+1];
        die "s2 > e2:$s2 > $e2" if($s2 > $e2);
        if($e1 <= $s2){$flag_ascend = 1};
        if($e2 <= $s1){$flag_descend = 1};
     }
     if($flag_ascend == 1 && $flag_descend == 0){return "ascend"}elsif($flag_ascend == 0 && $flag_descend == 1){return "descend"}elsif($flag_ascend == 1 && $flag_descend == 1){return "unsort"}else{return "unknow"}
}

sub segShareMerge(){
   #must clustered as an overlap group first
   #fill in overlap array and scan to report
   my $segs = shift;
   die "empty segs :@$segs;" if(scalar @$segs == 0);
   #get range
   my ($min,$max);
   my $depth=0;
   foreach my $seg(@{$segs}){
      my ($chr,$start,$end,$id)= split/\t/,$seg;
      die "start >= end: $start >= $end" if($start >= $end);
      $depth++;
      if(!defined $min){$min = $start}else{if($min > $start){$min = $start} }
      if(!defined $max){$max = $end}else{if($max < $end){$max = $end} }
   }
   #start to fill
   my @overlap;
   for (0..$max-$min){$overlap[$_]=0}
   foreach my $seg(@{$segs}){
      my ($chr,$start,$end,$id)= split/\t/,$seg;
      for my $i(($start-$min)..($end-$min)){$overlap[$i]++}
   }
   #start to scan and report full depth overlap regions, use the same method in statBlock()
   my @segments;
   my ($start,$end,$cnt)=(0,0,0);
   for my $i(0..$max-$min){
     die "not filled at $i from min $min, max $max" if(!defined $overlap[$i]);
     if($overlap[$i] == $depth){
        $cnt++;
        if($i==($max-$min)){ #last element
          $start=$i-$cnt+1;
          $end=$i;
          $start+=$min;
          $end+=$min;
          push @segments,$start."-".$end;
          $cnt=0;
        }
     }else{
           $start=$i-$cnt;
           $end=$i-1;
           if($start<=$end){$start+=$min;$end+=$min;push @segments,$start."-".$end}
           $cnt=0;
          }
   }#for end

   return (\@segments,"$min-$max");
}

