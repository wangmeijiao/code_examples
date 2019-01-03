


sub readInfo(){ #genome info file with spec order
  my %info;
  my @order;
  my $file = shift @_;
  open INF, $file or die"$!";
  $/="#";
  <INF>;
  while(<INF>){
   chomp;
   my @box=split/\n+/,$_;
   my $species=shift @box;
   push @order,$species;
   foreach(@box){
     my ($type,$file)=split/[\t ]+/,$_;
     next if($type =~/^"/);
     if(!exists $info{$species}->{$type}){ $info{$species}->{$type}=$file }else{die "dup $type in species $species\n"}
   }
  }
  close INF;
  $/="\n";
  return (\%info,\@order);
}




sub readInfo(){ #genome info file
  my %info;
  #my @order;
  my $file = shift @_;
  open INF, $file or die"$!";
  $/="#";
  <INF>;
  while(<INF>){
   chomp;
   my @box=split/\n+/,$_;
   my $species=shift @box;
   #push @order,$species;
   foreach(@box){
     my ($type,$file)=split/[\t ]+/,$_;
     next if($type =~/^"/);
     if(!exists $info{$species}->{$type}){ $info{$species}->{$type}=$file }else{die "dup $type in species $species\n"}
   }
  }
  close INF;
  $/="\n";
  return \%info;
}






