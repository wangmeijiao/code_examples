#!/bin/perl
use strict;
use warnings;
use Getopt::Long; 
use Data::Dumper;

###***this file contain general code for files with ## annotated headlines, remove of these headlines before readin is not necessary*****
#usage:  cat final.chainnet.sorted.axt |perl axtFilter.pl > final.chainnet.sorted.filtered.axt 2> filter.stats 
#description:
 # filter axt and do basic statistics
 # axt format are orgnized special for pairwise alignment;
 # for multiplex alignment use maf;
#data example:
 # ##matrix=axtChain 16 91,-114,-31,-123,-114,100,-125,-31,-31,-125,100,-114,-123,-31,-114,91
 # ##gapPenalties=axtChain O=400 E=30
 # ##aligner=blastz.v7 H=2000  Y=3400  L=6000  K=2200  Q=HoxD55.q 
 # ##matrix=blastz 16 91,-90,-25,-100,-90,100,-100,-25,-25,-100,100,-90,-100,-25,-90,91
 # ##gapPenalties=blastz O=400 E=30
 # ##blastzParms=O=400,E=30,K=2200,L=6000,M=0
 # 12 OsjapChr04 45094 45152 chr04 10120927 10120984 + 667
 # TAAcactggtggagaaaccctttgtagtcccggtttgtaaccccc--ctttagtcccggtt
 # TAAGACCGCCGCCGGCGCCCTTTATA---CTGGTTTAGGCCCGCCAGCGCTAATCAAGGTC

 my ($len_refchr,$len_quechr,$showHead,$showStats);
 GetOptions( "refLen=i"=>\$len_refchr,
             "queLen=i"=>\$len_quechr,
             "showHead"=>\$showHead,
             "showStat" => \$showStats
           );

#globe variates
my %blocks; 
my @headlines;
my $cnt; 

##read axt file by block and do several jobs (include head lines)
$/="\n\n";
while(<stdin>){
   chomp;
   #print "$_\n";
   my @box=split/\n/,$_;
   my $line;
   while(($line=shift @box)=~/^#/){  # shift -> =~ -> while(0/1)
       push @headlines,$line."\n";
   }#check headlines, brillient code!
   my $align=$line;
   my ($blockID,$refID,$refStart,$refEnd,$queID,$queStart,$queEnd,$strand,$score)=split/ / ,$align;
   my $ref=shift @box;
   my $query=shift @box;
   my ($align_len,$identity,$gap_perc,$N_perc)=&starCnt($ref,$query);
  #do filtering: align_len >=200,identity%(not include -, N/n)>=0.9, (lowcase+N)%<=0.4, average distance from kx-y=0 <=$dist
   if($align_len>=200 && $identity>=0.9 && &isNear($refStart,$refEnd,$queStart,$queEnd,$len_refchr,$len_quechr,5000000)){ 
     #output to filtered axtfile
     print "$align\n$ref\n$query\n\n";
     #store filtered blocks into hash
     if(exists $blocks{$blockID}){die "duplicated blocks\n"}else{ my @temp;push @temp,$blockID,$refID,$refStart,$refEnd,$queID,$queStart,$queEnd,$strand,$score,$ref,$query,$align_len,$identity,$gap_perc,$N_perc;$blocks{$blockID}=\@temp;}
   }
  #inner clock
   $cnt++;
   #if($cnt==10){last};
}#while end here
$/="\n";
#print STDERR Dumper \%blocks;

#headlines
if($showHead){print STDERR join ("",@headlines)}
#do final statistics:num of blocks, longest block, shortest block, covarage_ref% covarage_que%
if($showStats){
  my @keylist=sort{$a<=>$b} keys %blocks;
  my $first_key=shift @keylist;
  my ($num_blocks,$len_min,$len_max,$coverage_ref,$coverage_que)=(0,${$blocks{$first_key}}[11],${$blocks{$first_key}}[11],0,0);
  my @ref_segs;
  my @que_segs;
  foreach (sort{$a<=>$b} keys %blocks){
    $num_blocks++; 
    if(${$blocks{$_}}[11]<$len_min){$len_min=${$blocks{$_}}[11]}
    if(${$blocks{$_}}[11]>$len_max){$len_max=${$blocks{$_}}[11]}
    push @ref_segs,${$blocks{$_}}[2]."-".${$blocks{$_}}[3];
    push @que_segs,${$blocks{$_}}[5]."-".${$blocks{$_}}[6];
  }#foreach end here
  $coverage_ref=&coverageSegs(\@ref_segs,$len_refchr);
  $coverage_que=&coverageSegs(\@que_segs,$len_quechr);
  print STDERR "Filtered summary:\n------------\ntotal blocks: $num_blocks\nalignLen_min/max: $len_min/$len_max\ncoverage_ref: $coverage_ref %\ncoverage_que: $coverage_que %\n";
}

##############subs ######################
sub starCnt(){ 
  my ($ref,$query)=@_;
  my $align_len=length $ref;
  my @ref=split//,$ref;
  my @query=split//,$query;
  my ($match,$low,$N,$gap)=(0,0,0,0); #print"$match,$low,$N,$gap\n";
  if($#ref == $#query){
   for(0..$#ref){
       if($ref[$_] eq "-"|| $query[$_] eq "-" ){
               $gap++;
            }elsif($ref[$_] eq "N"|| $query[$_] eq "N"|| $ref[$_] eq "n"|| $query[$_] eq "n" ){
                 $N++;
               }elsif($ref[$_] eq $query[$_] || $ref[$_] eq lc($query[$_]) || $ref[$_] eq uc($query[$_]) ){
                   $match++;
                   #if($ref[$_] eq uc($ref[$_]))
                 }
   }#for end here
   return($align_len,$match/$align_len,$gap/$align_len,$N/$align_len);
  }else{die "err, alignment diffs"}
}#sub end here

sub isNear(){
  #average of distance from points ($refStart,$queStart), ($refEnd,$queEnd) to line ($len_query/$len_ref)x-y=0 <= 1k
  #use dist= |kx-y|/sqrt(k^2+(-1)^2)
  my ($refStart,$refEnd,$queStart,$queEnd,$len_refchr,$len_quechr,$dist)=@_;
  my $k=$len_quechr/$len_refchr;
  my $distStart=abs($k*$refStart-$queStart)/sqrt($k*$k+1);
  my $distEnd=abs($k*$refEnd-$queEnd)/sqrt($k*$k+1);
  my $distAve=($distStart+$distEnd)/2;
  if($distAve<=$dist){return 1}else{return 0}
}#sub end here

sub coverageSegs(){
  my $segs=shift;
  my $len=shift;
  my @count;
  my $cnt;
  foreach(@{$segs}){
   my ($start,$end)=split/-/,$_;
   for($start..$end){$count[$_]++}
  }
  for(1..$len){if(defined $count[$_]){$cnt++}}
  if(defined $cnt){return int (100*$cnt/$len)}else{die "count problem!\n"}
}#sub end here
