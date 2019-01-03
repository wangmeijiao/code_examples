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

