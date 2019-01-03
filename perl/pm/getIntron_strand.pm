

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



