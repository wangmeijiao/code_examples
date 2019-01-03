use strict;
use warnings;


  my $group = $ARGV[0];

  my $cnt=0;
  $/ = "\n\n";
  while(<stdin>){
    chomp;
    $cnt++;
    if($cnt == $group+1){ print "Block $cnt\n$_\n";last}
  }
  $/="\n";




