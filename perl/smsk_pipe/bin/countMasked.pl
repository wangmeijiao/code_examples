use warnings;
use strict;


$/=">";
while(<stdin>){
   chomp; 
   next if($_ eq "");
   my @box=split/\n/,$_;
   my $name=shift @box;
   my $seq=join"",@box;
   my $len=length $seq;
   my $cnt;
   while($seq=~/(a|t|c|g|n|N)/g){
     $cnt++;

   }
   #my $tail=substr($seq,$len-380);
   #print "$name-$len\n$seq\n$tail\n";
   my $perc=100*$cnt/$len;
   print "total masked nts:$cnt,total nts:$len,percent:$perc\n";

}
$/="\n";
