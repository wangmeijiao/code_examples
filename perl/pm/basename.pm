

sub basename(){
  my ($path,$cut)=@_;
  die "your path has no file" if($path =~/\/$/);
  my @path=split/\/+/,$path;
  #if($path=~/^\//){shift @path} #path[0] is empty
  my $base=pop @path;
  my $dir=join("/",@path);
  my @names=split/\./,$base;
  for(my $i=1;$i<=$cut;$i++){pop @names}
  my $prefix=join(".",@names);
  return ($dir,$base,$prefix);
}

