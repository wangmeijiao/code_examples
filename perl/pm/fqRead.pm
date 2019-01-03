
sub fqRead(){
  my ($file,$len) = @_;
  $/="\@SRR";
  open FQ, $file or die "$!";
  open OUT24, "> out.24nt.fastq" or die "$!";
  open OUT_big, "> out.big24nt.fastq" or die "$!";
  open OUT_less, "> out.less24nt.fastq" or die "$!";
  my $cnt;
  while(<FQ>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box = split/\n+/,$_;
    die "fastq group num != 4 at $_" if(scalar @box != 4);
    $cnt++;
    my $id = $box[0];
    $id=~s/\/1$//;
    my $seq = $box[1];
    my $qual = $box[3];
    my $len_seq = length $seq;
    my $len_qual = length $qual;
    die "length of seq($len_seq) != length of qual($len_qual) at $_" if($len_seq != $len_qual);


    if($len_seq == $len){
      print OUT24 "\@SRR$id\n$seq\n+\n$qual\n"
    }elsif($len_seq < $len){
       print OUT_less "\@SRR$id\n$seq\n+\n$qual\n"
     }else{print OUT_big "\@SRR$id\n$seq\n+\n$qual\n"}
    if($cnt % 1000000 == 0){print STDERR "#"}

  }
  close FQ;
  close OUT24;
  close OUT_big;
  close OUT_less;
  $/="\n";
  print STDERR "\n";
  return 1;

}



