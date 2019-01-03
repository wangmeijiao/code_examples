use strict;
use warnings;
use Data::Dumper;

my @test = ("abcd3","abcd31","abcd34","abcd3627");


my $basename = &getCommonStr(\@test,"left");
print "$basename\n";

sub getCommonStr(){
   my ($name,$direct)=@_;
   die "empty str" if(scalar @$name == 0);
   my $base;
   my @box;
   foreach my $str(@{$name}){push @box, [split//,$str]   }
   #get the minimal length 
   my $min_len = length $name->[0];
   foreach my $str(@$name){ if(length $str < $min_len ){ $min_len = length $str}}
   die "min_len == 0 @$name" if($min_len == 0);
   #print Dumper \@box;
   if($direct eq "left"){
     #my $c = $box[0][0];
     my $flag = 0;
     for(my $i=0;$i<$min_len;$i++){
       my $c = $box[0][$i];
       for my $j(0..$#box){if( $box[$j][$i] eq $c){}else{$flag =1}}
       if($flag == 1 ){last}
       $base.=$c;
     }
   }elsif($direct eq "right"){$base="NA"}else{$base="NA"}
   if(!defined $base){$base = 'NA'}
   return $base;

}





