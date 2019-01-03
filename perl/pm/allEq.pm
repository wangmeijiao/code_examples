sub allEq(){
  #check all points equal
  my @data = @_;
  die "allEq data empty" if(scalar @data == 0);
  if(scalar @data == 1){return 1}
  my $temp;
  my $flag = 1;#all equal by default
  foreach(@data){if(!defined $temp){$temp = $_}else{ if($temp != $_){$flag = 0}  }}
  return $flag;
}

