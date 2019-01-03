sub digitNumber(){
  my $n = shift;
  die "emtpy of <0: $n" if($n < 0 || $n eq "");
  my @box = split//,"$n";
  my $cnt = 0;
  foreach (@box){ if(defined $_ && $_ ne "."){$cnt++}else{last}   }
  return $cnt;
}
