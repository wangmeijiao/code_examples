sub gtfRead(){ #gene lines only
  my $file = shift;
  my %attr;
  open GTF, $file or die "$!";
  while(<GTF>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
    my @box=split/\t/,$_;
    next if($box[2] ne "gene");
    #my @attr=split/;/,$box[8];
    $box[8] =~ /gene_id "([^;]+)";/;
    my $id=$1;
    my @attr=split/;/,$box[8];
    foreach my $ctrl(@attr){
      $ctrl=~s/^ +//;
      $ctrl=~s/"//g;
      my ($k,$v)=split/ +/,$ctrl;
      if(!exists $attr{$id}->{$k}){$attr{$id}->{$k} = $v}else{ die "dup key $k in attr $id"}
    }
  }
  close GTF;
  return \%attr;
}
