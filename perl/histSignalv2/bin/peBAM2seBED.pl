#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

#this script recieve up stream datalines and
#join or link read pairs to single end bed
#bam -> bed -> joined,linked and sorted bed
#example:  bamToBed -i hits_k92_tigr7_beststrataM1v2.bam |perl peBAM2seBED.pl |sort -k 1,1 -k2,2n > hits_k92_tigr7_beststrataM1v2.se.bed

my $insert=$ARGV[0];
$insert||=600;
my %bed;

while (<STDIN>){
   chomp;
   my @box=split/\t/,$_;
   my $read=$box[3];
   $read=~s/\/[1-2]//;
   #print "$read\t$start\t$end\t$chr\t$score\t$strand\n";
   if(exists $bed{$read}){push @{$bed{$read}},$_ }else{my @temp; $bed{$read}=\@temp; push @{$bed{$read}}, $_;}

}
#print Dumper \%bed;

foreach (keys %bed){
    if (@{$bed{$_}}==2){ #well paired, then link or join them
           my @box0=split /\t/, ${$bed{$_}}[0] ;
           my $start0=$box0[1];         
           my $end0=$box0[2];         
           my @box1=split /\t/, ${$bed{$_}}[1] ;         
           my $start1=$box1[1];         
           my $end1=$box1[2];         
           my($min,$max)=&minMax($start0,$end0,$start1,$end1);
           if($box0[0] eq $box1[0] && $box0[4]==$box1[4] && $box0[5] ne $box1[5] && $max-$min+1<=$insert){print "$box0[0]\t$min\t$max\t$_\t$box0[4]\t+\n";}
               else{print STDERR "failed to join for @{$bed{$_}}\n";} #make sure these two reads are paired or stop
      }elsif(@{$bed{$_}}==1){ #single read with the other failed, remain unchanged
              print @{$bed{$_}},"\n";
        }elsif(@{$bed{$_}}>2||@{$bed{$_}}<1){ print STDERR "pair reads >2 or < 1 at @{$bed{$_}}\n"; } #wrong if this happens



}

sub minMax(){
  my @temp=@_;
  my ($min,$max);
  $min=$temp[0];
  $max=$temp[0];
  foreach (@temp){
    if($_<$min){$min=$_};
    if($_>$max){$max=$_};
  }
  return($min,$max);
}
