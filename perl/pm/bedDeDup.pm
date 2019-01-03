use strict;
use warnings;

my %segs;
my @order;
my $cnt;
while (<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
  $cnt++;
  if($cnt % 1000000 == 0){print STDERR "#"}
  my @box=split/[\t ]+/,$_;
  if(!exists $segs{"$box[0]:$box[1]-$box[2]"}){
     $segs{"$box[0]:$box[1]-$box[2]"}=join("\t",@box);
     push @order,"$box[0]:$box[1]-$box[2]"
  }else{
         my @temp=split/\t/,$segs{"$box[0]:$box[1]-$box[2]"};
         for my $i(3..$#box){next if($i == 5);$temp[$i].="|".$box[$i]}
         $segs{"$box[0]:$box[1]-$box[2]"}=join("\t",@temp);
       }

}

print STDERR "\nread done\n";

$cnt=0;
foreach my $id(@order){
   print "$segs{$id}\n";
   $cnt++;
   if($cnt % 1000000 == 0){print STDERR "#"}
}


print STDERR "\nwrite done\n";
