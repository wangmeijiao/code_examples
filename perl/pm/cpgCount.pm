use strict;
use warnings;



my $test = "CTCTTACGGTCGTTCGTTTAGCT";
print &cpgCount($test),"\n";



sub cpgCount(){
  my $seq=shift @_;
  $seq = uc($seq);
  die "seq empty" if($seq eq "");
  my $cnt=0;
  #my @matches = ($seq=~/CG/g);
  while($seq=~/CG/g){$cnt++}
  return $cnt;
}




