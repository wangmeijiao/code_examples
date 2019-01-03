use strict;
use warnings;


print join(" ",&isOverlap_detail("chr1|2|10","chr1|3|8")),"\n";


##sub

sub isOverlap_detail(){
  #input two segs,output overlap flag and overlap bp (if true)
  my ($seg1,$seg2) = @_;
  my ($chr1,$s1,$e1) = split/\|/,$seg1;
  die "$s1 > $e1 at $seg1, $seg2" if($s1 > $e1);
  my $len1 = $e1-$s1+1;
  my ($chr2,$s2,$e2) = split/\|/,$seg2;
  die "$s2 > $e2 at $seg1, $seg2" if($s2 > $e2);
  die "chr diffs when compare $seg1, $seg2" if($chr1 ne $chr2);
  my $len2 = $e2-$s2+1;
  my $isOverlap = 0;
  my $len_overlap = 0;
  my ($perc1,$perc2) = (0,0);
  if($e1 < $s2 || $s1 > $e2){return ($isOverlap,$len_overlap,$perc1,$perc2)}else{
     $isOverlap = 1;
     my @temp = ($s1,$e1,$s2,$e2);
     my ($min,$max) = &minMax(\@temp);      
     my @star;
     my ($cnt_1,$cnt_2) = (0,0);
     die "coordinate err at $s1-$min, $e1-$min, $s2-$min, $e2-$min" if($s1-$min < 0 || $e1-$min < 0 || $s2-$min < 0 || $e2-$min < 0);
     for my $i(($s1-$min)..($e1-$min)){$star[$i]++}
     for my $i(($s2-$min)..($e2-$min)){$star[$i]++}
     for my $i(0..$#star){
        if(!defined $star[$i]){}elsif($star[$i] == 1){$cnt_1++}elsif($star[$i] == 2){ $cnt_2++}elsif($star[$i] >= 2){die "unknow situation at $star[$i]"}
     }
     $len_overlap = $cnt_2;
     $perc1= sprintf ("%.3f",$len_overlap/$len1);
     $perc2= sprintf ("%.3f",$len_overlap/$len2);
     return($isOverlap,$len_overlap,$perc1,$perc2);
  }

}

sub minMax(){
  my $index=shift;
  my ($min,$max)=($index->[0],$index->[0]);
  foreach (@{$index}){
    if($min>$_){$min=$_}
    if($max<$_){$max=$_}
  }
  #print ("$min,$max\n");
  return ($min,$max);
}#end of minMax sub


