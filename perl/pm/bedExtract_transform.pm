use strict;
use warnings;

my $file=$ARGV[0];
my $region=$ARGV[1];
my $prefix=$ARGV[2];

#my $bed=&readBed_extract($file,$region,$prefix);
my $bed=&readBed_region_transform($file,$region,$prefix);
#my $bed=&bedOverlap_transform($file,$region,$prefix);

#above are the same or almost the same(1 and 2 are the same, 3 has one more linethan 1, 2)




foreach(@{$bed}){print "$_\n"}


###sub##

sub readBed_extract(){
   #read in and extract out target region
   my ($file,$region,$prefix)=@_;
   my ($chr,$rstart,$rend)=split/:|-/,$region;
   my @regions;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my @temp=split/\t/,$_;
     die "start >= end\n" if($temp[1]>=$temp[2]);
     next if($temp[0] ne $chr || $temp[2]<=$rstart || $temp[1] >=$rend);
     if($temp[1] < $rstart){$temp[1]=$rstart}
     if($temp[2]>$rend){$temp[2]=$rend}
     $temp[1]=$temp[1]-$rstart;
     $temp[2]=$temp[2]-$rstart;
     $temp[0]=$prefix;
     push @regions,join("\t",@temp);
   }
   close BED;
   return \@regions;

}


sub bedOverlap_transform(){
  my ($file,$region,$prefix)=@_;
  my @hits;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  open BED, $file or die "$!";
  while(<BED>){
   chomp;
   next if($_ eq ""|$_=~/^#/ || $_=~/^\s+$/);
   my @box=split/[\t ]+/,$_;
   my ($chr_in,$start_in,$end_in)=@box;
   next if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{
     #transform and cut ends
     $box[1]-=$start;
     if($box[1] < 0 ){$box[1] = 0}
     $box[2]-=$start;
     if($box[2] <= 0 ){die "wrong at $start_in and $end_in after transform: transformed end <= 0"}
     $box[0]=$prefix;
     push @hits,join("\t",@box);
   }
  }
  close BED;
  return \@hits;
}


sub readBed_region_transform(){
  my ($file,$region,$prefix)=@_;
  my ($chr_region,$s_region,$e_region)=split/:|-/,$region;
  my @bedg;
  open BEDG, $file or die "$!";
  while(<BEDG>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my ($chr,$s,$e,$v)=split/[\t ]+/,$_;
    die "s > e $s,$e" if($s > $e);
    next if($chr ne $chr_region);
    if($e <= $s_region || $e_region <= $s){}else{
       if($s < $s_region){$s = $s_region}
       if($e > $e_region){$e = $e_region}
       $s-=$s_region;
       $e-=$s_region;
       push @bedg,"$prefix\t$s\t$e\t$v";
    }
  }
  close BEDG;
  return \@bedg;

}


