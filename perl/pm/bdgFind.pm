sub bdgFind(){
  #find position in bdg segs, return bdg index( report many if multihit, report "nohit" if not found)
  my ($bdg, $pos) = @_;
  my @hits;
  for my $i(0..$#$bdg){my($chr,$s,$e,$v)=split/[\t ]+/,$bdg->[$i]; if($s <= $pos && $pos <= $e){ push @hits,$i} }
  if(scalar @hits == 0){ push @hits,"nohit"}
  return \@hits;
}

