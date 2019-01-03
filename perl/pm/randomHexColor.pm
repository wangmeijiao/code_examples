
sub randomHexCol(){
  my $num=shift;
  my @colors;
  for (my $i = 0; $i < $num; $i++) {
    my ($rand,$x);
    my @hex;
    for ($x = 0; $x < 3; $x++) {
      $rand = rand(255);
      $hex[$x] = sprintf ("%x", $rand);
      if ($rand < 9) {
        $hex[$x] = "0" . $hex[$x];
      }
      if ($rand > 9 && $rand < 16) {
        $hex[$x] = "0" . $hex[$x];
      }
    }
    $colors[$i] = "\#" . $hex[0] . $hex[1] . $hex[2];
  }
  #print "@colors\n";
  return \@colors;
}


