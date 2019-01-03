
sub segsOverlap_multi(){
  #compare multisegments and output perct%
  #input s-e-other;
  #no need to sort
  #1,brief overlap: left1-right1 vs left2-right2, get star range
  #2overlap in details: fill star range
  my ($idx1,$idx2) = @_;
  my ($perc1,$perc2);
  my ($sum1,$sum2);
  my $cnt_1 = 0;
  my $cnt_2 = 0;
  my ($left1,$right1) = &getRange($idx1);  
  my ($left2,$right2) = &getRange($idx2);  
  die "left1 >= right1 at $left1 >= $right1" if($left1 >= $right1);
  die "left2 >= right2 at $left2 >= $right2" if($left2 >= $right2);
  if($right1 < $left2 || $right2 < $left1){
    print STDERR "  step 1, brief overlap:not overlap\n";
    return (0,0,0,0,0);
  }else{
     print STDERR "  step 1, brief overlap:overlap detected";
     my ($min,$max);
     if($left1 < $left2){$min = $left1}else{$min = $left2}
     if($right1 > $right2){$max = $right1}else{$max = $right2}
     my $length = $max-$min+1;
     print STDERR " with length $length\n";
     #start to fill star
     print STDERR "  step 2, overlap in details\n";
     my @star;
     foreach my $seg(@{$idx1}){
        my ($s,$e) = split/-/,$seg;
        die "s>e at $s > $e" if($s > $e); 
        $sum1+=($e-$s+1);
        die "coordinate err at $s-$min" if($s-$min < 0);
        die "coordinate err at $e-$min" if($e-$min < 0);        
        for my $i(($s-$min)..($e-$min)){$star[$i]++ }
     }
     for my $i(0..$#star){if(!defined $star[$i]){}elsif($star[$i] >= 2){print STDERR "dup segs at @$idx1\n"}}

     foreach my $seg(@{$idx2}){
        my ($s,$e) = split/-/,$seg;
        die "s>e at $s > $e" if($s > $e); 
        $sum2+=($e-$s+1);
        die "coordinate err at $s-$min" if($s-$min < 0);
        die "coordinate err at $e-$min" if($e-$min < 0);        
        for my $i(($s-$min)..($e-$min)){$star[$i]++ }
     }

     #start to summary star
     for my $i(0..$#star){
       if(!defined $star[$i]){}elsif($star[$i] == 1){$cnt_1++}elsif($star[$i] == 2){ $cnt_2++}elsif($star[$i] >= 2){die "unknow situation at $star[$i]"}
     }
     $perc1 = sprintf ("%.3f",$cnt_2/$sum1);
     $perc2 = sprintf ("%.3f",$cnt_2/$sum2);

   }
   return ($perc1,$perc2,$sum1,$sum2,$cnt_2);

}


sub getRange(){
  my $idx = shift;
  my ($left,$right);
  foreach my $i(0..$#$idx){
    my ($s,$e) = split/-/,$idx->[$i];
    die "s>e" if($s > $e);
    if(!defined $left && !defined $right){($left,$right) = ($s,$e)}else{
       if($s < $left){$left = $s}
       if($e > $right){$right = $e}
    }
  }
  return ($left,$right);
}


