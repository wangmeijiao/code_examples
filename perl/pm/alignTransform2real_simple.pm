

sub alignTransform2real_simple(){
  #transform aligment coordinates to real sequence coordinates
  my ($align,$align_header,$region) = @_;
  my ($start,$end) = split/-/,$region;#align start-end after cut ends  
  my %region_real;
  foreach my $spec(keys %{$align}){
     my @seq_align = split//,$align->{$spec};
     my ($left,$right) = split/-/,$align_header->{$spec};
     my $cnt=0;
     my $cnt_real=0;
     my $start_real = $left;
     my $end_real = $right;
     foreach my$base(@seq_align){
       if($base ne "-"){$cnt_real++;$cnt++}else{$cnt++}
       if($cnt == $start){$start_real = $left + $cnt_real}
       if($cnt == $end){$end_real = $left + $cnt_real}
     }
     $region_real{$spec} = "$start_real-$end_real";
  }
  return \%region_real;
}



