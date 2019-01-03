sub mfaRead(){
  my $file = shift;
  $/=">";
  my %align;
  my $num_align;
  my $len;
  my @align_order;
  open MFA, $file or die "$1";
  while(<MFA>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp = split/[\t ]+/,$id;
    my $id_real = $temp[0];
    push @align_order,$id_real;
    my $align=join"",@box;
    if(! defined $len){$len = length $align}else{ if($len != length $align){die "align len diff at $id"}}
    if(!exists $align{$id_real}){$align{$id_real}=$align; $num_align++}else{die "dup alignment $id_real\n"}
  }
  close MFA;
  $/="\n";
  return (\%align,\@align_order,$num_align,$len);
}

