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


