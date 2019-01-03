use strict;
use warnings;


my $test = "ndshufu23jjh4n382387fjhf86fd8d7s5f9g56s6ds7dddf7dfffdddf67sss7s7dddddfdsdddfdffffdddddddss788sddd78ssd7d7dddd";
print "$test\n";

my $cnt=0;
while($test =~/d+/g){
  $cnt++;
  print "$cnt match at $-[0],$+[0] ";
  print "hit is:",substr($test,$-[0],$+[0]-$-[0]),"\n";


}

