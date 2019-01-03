


sub readBdg(){
   my $file=shift;
   my @bdg;
   open BDG, $file or die "$! $file";
   while(<BDG>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
     my($chr,$s,$e,$v)=split/[\t ]+/,$_;
     die "$s >= $e at $_" if($s >= $e);
     push @bdg,$_;
   }
   close BDG;

   my $bdg_pack=&bdgPack(\@bdg);
   return $bdg_pack;
}

sub readBdg_advance(){
  my $file = shift;
  my ($min,$max,$mean);
  my $len;
  my ($left,$right);
  my @bdg;
  my $sum;
  my $cnt;
  open BDG, $file or die "$!";
  while(<BDG>){
     chomp;
     next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
     push @bdg,$_;
     my ($chr,$s,$e,$v)=split/[\t ]+/,$_;
     die "$s >= $e at $_" if($s >= $e);
     $len+=($e-$s);
     if(!defined $min && !defined $max){($min,$max)=($v,$v)}else{
       if($min > $v){ $min = $v}
       if($max < $v){ $max = $v}
     }
      if(!defined $left && !defined $right){($left,$right)=($s,$e)}else{
       if($left > $s){ $left = $s}
       if($right < $e){ $right = $e}
     }
     $cnt++;
     $sum+=$v;
  }
  close BDG;
  $mean=$sum/$cnt;
  return (\@bdg,$min,$max,$mean,$len,"$left-$right");
}



