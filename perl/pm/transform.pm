###
sub readTransform(){
  my $file = shift @_;
  my %tab;
  open TAB, $file or die "$!";
  while(<TAB>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    $_=~s/^\s+//;
    my($chr1,$start1,$end1,$chr2,$start2,$end2)=split/[ \t]+/,$_;
    if(!exists $tab{$chr1} && $start1<$end1 && $start2<$end2){
      $tab{$chr1}="$chr1:$start1-$end1"."|"."$chr2:$start2-$end2";
      #$tab{$chr2}="$chr1:$start1-$end1"."|"."$chr2:$start2-$end2";
    }else{die "more than one transform lines for $chr1\n"}
  }
  close TAB;
  #print Dumper \%tab;
  return \%tab;
}


####
sub transform_dropdown(){
    my ($region,$tab)=@_;
    my %tab=%{$tab};
    #get real_region from $region
    my($chr,$start,$end)=split/:|-/,$region;
    my($chr_real,$start_real,$end_real);
    if(exists $tab{$chr}){
      my($region1,$region2)=split/\|/,$tab{$chr};
      my($chr1,$start1,$end1)=split/:|-/,$region1; #local region
      my($chr2,$start2,$end2)=split/:|-/,$region2; #genome region
      if($chr eq $chr2){$chr_real=$chr1}else{die};
      $start_real=$start-$start2;
      if($start_real < $start1){ $start_real=$start1}
      $end_real=$end-$start2;
      if($end_real > $end1){ $end_real = $end1}
    }else{die "can't found $chr in transform tab\n"}
    my $region_real="$chr_real:$start_real-$end_real";
    #return $region_real;
    return ($chr_real,$start_real,$end_real);
}


sub transform_liftup(){
    my ($region,$tab)=@_;
    my %tab=%{$tab};
    #get real_region from $region
    my($chr,$start,$end)=split/:|-/,$region;
    my($chr_real,$start_real,$end_real);
    if(exists $tab{$chr}){
      my($region1,$region2)=split/\|/,$tab{$chr}; 
      my($chr1,$start1,$end1)=split/:|-/,$region1;
      my($chr2,$start2,$end2)=split/:|-/,$region2;
      if($chr eq $chr1){$chr_real=$chr2}else{die};
      if($start<=$end1){$start_real=$start2+($start-$start1)}else{die}
      if($end<=$end1){$end_real=$start2+($end-$start1)}else{die}
    }else{die "can't found $chr in transform tab\n"}
    my $region_real="$chr_real:$start_real-$end_real";
    #return $region_real;
    return ($chr_real,$start_real,$end_real);
}


sub coordTransform_reverse(){
# transform coordinate 180 degree, from minus to plus : first len-pos then exchange s and e
# for display RCed sequence, not for sequence extraction
# gff - strand features always numbered from the + strand
  my ($segs,$len) = @_; #input chr lenght
  my @segs_transform;
  foreach my $seg(@{$segs}){
    my ($s,$e) = split/-/,$seg;
    push @segs_transform,($len-$e)."-".($len-$s);
  }
  return \@segs_transform;

}

sub coordTransform_reverse_point(){
# transform coordinate 180 degree from minus to plus : first len-pos then exchange s and e
# for display RCed sequence, not for sequence extraction
# gff - strand features always numbered from the + strand
  my ($pos,$len) = @_; #input chr lenght
  return $len-$pos;

}


  
sub alignTransform2real(){
  #transform aligment coordinates to real sequence coordinates
  #considering CDS coordinates jumpping 
  my ($align,$align_header,$region) = @_;
  my ($start,$end) = split/-/,$region; 
  my %region_real;
  foreach my $spec(keys %{$align}){
     my @seq_align = split//,$align->{$spec};
     my $cnt=0;
     my $cnt_real=0;
     my $start_real = $start;
     my $end_real = $end;
     foreach my$base(@seq_align){
       if($base ne "-"){$cnt_real++;$cnt++}else{$cnt++}
       if($cnt == $start){$start_real = $cnt_real}
       if($cnt == $end){$end_real = $cnt_real}
     }
     
     my ($chr,$s,$e,$id,$strand,$segs) = split/\|/,$align_header->{$spec};
     my $coord_tab = &realTransform2Align($segs,$strand);
     foreach my $ctrl(keys %{$coord_tab}){
        my ($s,$e) = split/-/,$ctrl;
        my ($s_bp_real,$e_bp_real) = split/-/,$coord_tab->{$ctrl};
        if($start_real + 1 >= $s && $start_real <= $e){my $d = $start_real - $s +1;die "d < 0" if($d < 0);$start_real = $s_bp_real + $d  }
        if($end_real >= $s && $end_real <= $e){my $d = $end_real - $s +1;die "d < 0" if($d < 0);$end_real = $s_bp_real + $d  }
     }


     $region_real{$spec} = "$start_real-$end_real";
  }
  return \%region_real;
}


sub realTransform2Align(){
   #for CDS joined seqence aligment
   # if strand is +, check for overlap and sorted, if strand is -, do not sort or reverse
   my ($regions,$strand) = @_;
   my @segs = split/:/,$regions;
   my @segs_sorted = sort{
                           my ($s1,$e1) = split/-/,$a;
                           die "s1 > e1 in segs sort" if($s1 > $e1);
                           my ($s2,$e2) = split/-/,$b;
                           die "s2 > e2 in segs sort" if($s2 > $e2);
                           if($e1 <= $s2 || $e2<=$s1){}else{die "overlap segs $a,$b "}
                           if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                          }@segs;
   my %segs_transform;
   my $len=0;
   my $len_sum=0;
   foreach my $seg(@segs){
     my ($s,$e) = split/-/,$seg;
     $len=($e-$s+1);
     my $start_align = $len_sum+1;
     my $end_align = $len_sum+$len;
     $segs_transform{"$start_align-$end_align"} = "$s-$e";
     $len_sum+=$len;
   }
   return \%segs_transform;
}
