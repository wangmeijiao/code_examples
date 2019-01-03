
  sub faRead_gene(){
  #store fa of different versions
  my $file=shift;
  open IN,$file or die "$!";
  $/=">";
  my %seqs;
  while (<IN>){
     chomp;
      next if($_ eq "" || $_ =~/^\s+$/);
     my @box=split/\n+/,$_;
     my $head=shift @box;
     my @temp=split/[\t ]+/ ,$head;
     my $nameString=shift @temp;
     my ($name,$ver)=split/\./,$nameString;
     my $seq=join("",@box);
     $seq=~s/\*$//;
     $seqs{$name}->{$nameString}=$seq;     
  }
  close IN;
  $/="\n";
  return \%seqs;
}



sub faRead_chr(){
  my $file = shift;
  my %fa;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp=split/[\t ]+/ ,$id;
    $id=shift @temp;
    my $seq = join "", @box;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return \%fa;
}

sub faRead_chr_ext(){
  my $file = shift;
  my %fa;
  my @order;
  my $len;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp=split/[\t ]+/ ,$id;
    $id=shift @temp;
    push @order, $id;
    my $seq = join "", @box;
    if(!defined $len){$len = length $seq}else{if($len != length $seq){die "$id length unequal"}}
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return (\%fa,\@order,$len);
}


sub faFormat(){
  my ($seq,$num)=@_;
  my @seq=split//,$seq;
  my ($string,$cnt);
  foreach(@seq){
   $cnt++;
   if($cnt%$num==0){$string.=$_."\n"}else{$string.=$_}
  }
  if($cnt%$num!=0){return $string."\n"}else{return $string}
}



sub RC(){
  my $seq=shift;
  $seq=reverse $seq;
  $seq=~tr/atcgATCG/tagcTAGC/;
  return $seq;
}


sub faRC_coord_gap(){
  #faRC with coordinates and gap
  #start end related to this strand
  my ($seq,$chr_len,$start,$end) = @_;
  $seq = reverse $seq;
  $seq =~tr/atcgn\-ATCGN/tagcn\-TAGCN/;

  my $start_transform = $chr_len - $end;
  my $end_transform = $chr_len - $start;
  die "start <=0 || start > $chr_len" if($start <= 0 || $start > $chr_len);
  die "end <=0 || end > $chr_len" if($end <= 0 || $end > $chr_len);
  die "start == end " if($start == $end);

  return ($seq,$start_transform,$end_transform);

}




sub faExtract(){
    my ($index,$coord)=@_;
    my ($chr,$s,$e,$strand)=split/:|-/,$coord;
    die "not exists $chr in your fa file\n"if(!exists $index->{$chr});
    my $seq_extract=($strand eq "+")?(substr($index->{$chr},$s-1,$e-$s+1)):(&RC(substr($index->{$chr},$s-1,$e-$s+1)));
    #return uc($seq_extract);
    return $seq_extract;
}















