use strict;
use warnings;
#1, strand transform
#2, coordinate len-start, len-end

my $len=$ARGV[0];
die "len undefined" if(!defined $len);
die "len <= 0" if($len <=0);
while(<stdin>){
 chomp;
 next if($_ eq "" || $_=~/^#/ || $_ =~ /^\s+$/);
 my @box=split/[\t ]+/,$_;
 if($box[6] eq "+"){$box[6]="-"}elsif($box[6] eq "-"){$box[6]="+"}else{die"unknow strand type $box[6]"}
 $box[3]=$len-$box[3];
 $box[4]=$len-$box[4];
 die "start <=0 || start > $len" if($box[3] <= 0 || $box[3] > $len);
 die "end <=0 || end > $len" if($box[4] <= 0 || $box[4] > $len);
 die "start == end " if($box[3] == $box[4]);
 ($box[3] < $box[4])?():(($box[3],$box[4])=($box[4], $box[3]));
 print join("\t",@box),"\n";
}




