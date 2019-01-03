sub isValid(){
  #check and make sure that hash key is not space, is defined and not empty
  my $v=shift;
  my $flag=1;
  if(!defined $v || $v eq "" || $v=~/^\s+$/){$flag = 0}
  return $flag;
}

