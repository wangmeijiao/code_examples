

sub statBed_densi_genome(){
  #don't check feature overlap, since this happens too often
  #only stats one level of feature, gene and LTR, DNATE etc lost inside detail structure
  #stat densi for each segs, gene use exonBed to exclude intron
  #use tabix.so low level c function to speed up

  #need tabix_open first : $k92_densi = tabix_open($info->{$spec}->{"k92"}.".bgz");

  my ($bed,$densi,$range) = @_;

  my ($min,$mid,$max)=split/-/,$range;
 
  #my $genome_len=0;
  #foreach my $chr(keys %$len){$genome_len += $len->{$chr}}

  #calculate accumulated densi sum for all feature genome-wide
  my $cnt_feat;
  my $total_len = 0;
  my $total_densi = 0;
  my $total_ave_densi = 0;

  foreach my $id(keys %{$bed}){
    my ($chr,$start, $end,$id_in, undef, $strand) = split/\t/,$bed->{$id};
    die "start >= end: $start >= $end at $id" if($start > $end);
    next if($start == $end);
    $cnt_feat++;
    print STDERR "#" if($cnt_feat % 1000 == 0);
    #if($end > $e_r){ $end = $e_r}
    #if ($start < $s_r){ $start = $s_r}
    $total_len+=($end-$start+1);

    #extract data by tabix and process at the same time
    my $t = tabix_open($densi);
    my $iter = tabix_query($t, $chr, $start, $end);
    if(!defined $iter){ print EMPTY "can't find data at $chr, $start, $end\n" ;next}
    #my $cnt_hit;
    #my @hits;
    my $sum_densi_feat = 0;
    while (my $line = tabix_read($t, $iter)){
      chomp $line;
      next if($line eq "" );
      #$cnt_hit++;
      my ($c,$s,$e,$v) = split/[\t ]+/,$line;
      die"s>=e :$s >= $e or not defined v\n" if($s >= $e || !defined $v);
      if ($s < $start){ $s = $start}
      if($e > $end){ $e = $end}

      $v-=$min;
      if($v<0){$v=0}
      if($v>$max){$v=$max}

      $sum_densi_feat += $v*($e-$s+1);
      #push @hits,$line;

    }
    $total_densi += $sum_densi_feat;

    tabix_iter_free($iter);
  }#foreach id end
  print STDERR "\n";

  die "region_len($total_len) = 0 " if($total_len == 0);
  $total_ave_densi = sprintf ("%.3f",$total_densi/$total_len);
  return ($total_densi, $total_len,$total_ave_densi);

}#sub end



