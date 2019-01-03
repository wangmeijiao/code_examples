
sub readTag(){
  #readin tags and store in mem as a density table
  my $tagf=shift;
  my @tags;
  #my $cnt;
  open TAG, $tagf or die "$!";
  while(<TAG>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my($chr,$s,$e)=split/[\t ]+/,$_;
    push @tags,"$chr\t$s\t$e";
    #$cnt++;
    #if($cnt % 1000 ==0){print STDERR "#"}
  }
  #print STDERR "\n";
  close TAG;
  return \@tags;
}


sub tagOverlap(){
  my ($tag,$region)=@_;
  my @hits;
  my ($chr_region,$s_region,$e_region)=split/:|-/,$region;
  die "s >= e : $s_region >= $e_region" if($s_region >= $e_region);
  foreach(@$tag){
    my ($chr,$s,$e)=split/[\t ]+/,$_;
    next if($chr_region ne $chr);
    die "s >= e : $s >= $e" if($s >= $e);
    if($e <= $s_region || $e_region <= $s){}else{push @hits,$_}
  }
  return \@hits;
}


sub tagCount(){
  #make sure that tags are nonOverlap before count
  my ($tag,$region)=@_;
  my ($chr_region,$s_region,$e_region)=split/-|:/,$region;
  die "s >= e: $s_region >= $e_region" if($s_region >= $e_region);
  my $cnt=0;
  my $sum_len=0;
  foreach(@{$tag}){
    my ($chr,$s,$e)=split/\t/,$_;
    die "s >= e: $s >= $e" if($s >= $e);
    die "chr diffs: $chr_region ne $chr" if($chr_region ne $chr);
    if($e <= $s_region || $s >= $e_region){}else{
       $cnt++;
       if($s < $s_region){$s = $s_region}
       if($e > $e_region){$e = $e_region}
       $sum_len+=($e-$s+1);
    }
  }
  my $occupy_perc=sprintf "%.3f",$sum_len/($e_region-$s_region+1);
  my $densi_TE = sprintf "%.3f",$cnt/($e_region-$s_region+1);
  return ($cnt,$densi_TE,$sum_len,$occupy_perc);
}



