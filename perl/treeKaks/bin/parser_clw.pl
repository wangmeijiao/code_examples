#!/usr/bin/perl
use strict;
use warnings;

# usage: cat in.alnFa| perl thisfile.pl space_torrence filter_block_size
#        for example cat chr01_Os01g57270.1.txt.alnFa | perl parser_clw.pl 3 50 > chr01_Os01g57270.1.txt.out.parsed
# input file format: seq1	--------------ATGGCT-GGT-----------
#                    seq2	GGCTTATATTATCGATGCCTACGTTAGGATCATAT
# How does it work?
#
#    1), read in raw pairwise alignment and cut suspending ends (gap)
#    2), find conserved blocks
#    3), link small blocks to big one, filter small blocks then
#    4), caculate overall identity,coverage and block identity,coverage
#    5), summary
#

  my %align;
  $/=">";
  my $cnt;
  while(<stdin>){
   chomp;
   next if($_ eq "");
   $cnt++;
   my @box=split /\n/,$_;
   my $id=shift @box;
   my $seq=join "",@box;
   my @fragments=split//,$seq;   
   if($cnt==1){$align{"query"}=\@fragments}elsif($cnt==2){$align{"subject"}=\@fragments}else{die "wrong aligment file $id\n$seq\n"}
  }
  $/="\n";

  my @stars;
  if($#{$align{"query"}}==$#{$align{"subject"}}){
     for(0..$#{$align{"query"}}){
      if(${$align{"query"}}[$_] eq ${$align{"subject"}}[$_]){
          $stars[$_]="*"
         }elsif(${$align{"query"}}[$_] eq "-" || ${$align{"subject"}}[$_] eq "-"){
            $stars[$_]="-"
           }elsif(${$align{"query"}}[$_]=~/C|G/i && ${$align{"subject"}}[$_]=~/A|T/i || ${$align{"query"}}[$_]=~/A|T/i && ${$align{"subject"}}[$_]=~/C|G/i){
             $stars[$_]="v" #for transversion base pairs
             }elsif(${$align{"query"}}[$_]=~/C|G/i && ${$align{"subject"}}[$_]=~/C|G/i || ${$align{"query"}}[$_]=~/A|T/i && ${$align{"subject"}}[$_]=~/A|T/i){
              $stars[$_]="f" #for transformation base pairs
               }else{die("err when parse alignment\n");}
     }
   }else{die "err"}

=pod
  foreach (values %align){
    print @{$_};
    print "\n";
  }
  print @stars;
=cut

&StarCount($align{"query"},$align{"subject"},\@stars,$ARGV[0],$ARGV[1]);
                                                                   


########sub##############

sub StarCount(){ 

    my $star=$_[2]; #@star address
    my $qindex=$_[0];
    my $sindex=$_[1];
    my $link_space=$_[3];#space tolerance 
    my $block_size=$_[4];
    my $query=join("",@{$qindex});
    my $subject=join("",@{$sindex});
    my $star_full=join("",@{$star}); #print "\n\n$star_full\n";
    my $len_full=length($star_full);
    my $star_cut_ends=$star_full;$star_cut_ends=~s/^-+|-+$//g;#print "$star_cut_ends\n";
    my $len_cut=length($star_cut_ends); #print "\n$len_cut\n";
    my $cnt=0;
    my $i;
    my ($start,$end);
    my @segments;
    my @segments_linked;
    my @segments_linked_del;

    #find continual segments
    for $i(0..$#{$star}){
       if($$star[$i] eq "*"){
        $cnt++;
        if($i==$#{$star}){
          $start=$i-$cnt+1;
          $end=$i;
          push @segments,$start."-".$end;
          $cnt=0;
        }
       }else{ 
          $start=$i-$cnt;
          $end=$i-1;
          if($start<=$end){push @segments,$start."-".$end;}
          $cnt=0;
        }
    } 
    #print "@segments\n";
    #&pileup(0,$len,\@segments);

    #link small segments to long ones with gap<=$link_space and block size>$block_size
    for ($i=1;$i<=$#segments;$i++){
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        if(($tmp2[0]-$tmp1[1]-1)<=$link_space){
        $start=$tmp1[0];
        $end=$tmp2[1];
        $segments[$i]=$start."-".$end;
        }else{push @segments_linked, $segments[$i-1];}
        if($i==$#segments){push @segments_linked,$segments[$i];} 
    } 
    if($#segments==0){push @segments_linked,$segments[0];}
    #print "@segments_linked\n";
    #&pileup(0,$len,\@segments_linked);

    foreach (@segments_linked){
         my @tmp=split/-/,$_;
         if(abs($tmp[0]-$tmp[1]+1)>=$block_size){push @segments_linked_del, $_}
    }
    #print "@segments_linked_del\n";
    #&pileup(0,$len,\@segments_linked_del);


    #The overall %identity of this alignment should be caculated as:
    # %identity_total=matches_total/aligment_length_total_after_cut*100
     print "\t\t\t\t\t********begin clustalw2 parsing*********\n\n";
     print"\nOverall:\n.....................\n\n";
     my $coverage_cut=$len_cut/$len_full*100;
     print"alignment_length_full:$len_full (100%)\nalignment_length_remain_after_cut_ends:$len_cut($coverage_cut%)\n";
     my @temp;
     my $matches_total=(@temp)=($star_full=~/\*/g);
     my $identity_total=$matches_total/$len_cut*100;
     print "%identity_total(after_cut_ends):$identity_total%\n";

 
    #in an alignment block, alignment_length=matches+mismatches(namely,G<->C or A<->T or G/C<->A/T)+gaps_total(namely, gaps_q+gaps_s), N neither is mismatch nor gap
    #here we define %identity=matches/alignment_length*100 %coverage=aligment_length_of_this_block/aligment_length_total_after_cut*100
     print"\nconserved blocks:\n....................\n\n";
     my $total_conserved_block_length;
     foreach(@segments_linked_del){
       my @tmp=split/-/,$_;
       my $fragment=substr($star_full,$tmp[0],$tmp[1]-$tmp[0]+1);
       my $fragment_q=substr($query,$tmp[0],$tmp[1]-$tmp[0]+1);
       my $fragment_s=substr($subject,$tmp[0],$tmp[1]-$tmp[0]+1);       
       my $alignment_length=length $fragment;
       $total_conserved_block_length+=$alignment_length;
       my @temp;
       my $matches=(@temp)=($fragment=~/\*/g);#print @temp;
       my $identity=$matches/$alignment_length*100;
       print"\n>alignment_block $tmp[0]-$tmp[1] identity:%$identity alignment_length:$alignment_length\n$fragment_q\n$fragment_s\n$fragment\n";
     }      
     my $percent_of_conserved_parts=$total_conserved_block_length/$len_cut*100;
     print"\n\npercentage of total >$block_size conserved blocks (after remove ends):$percent_of_conserved_parts\n ";
     print"\n\t\t\t\t\t*********end of clustalw2 parsing************\n\n";
}#end of sub


sub pileup(){  # draw coordinate from $min->$max and in this region pile up segments in array @add WITHOUT sorting
 
  my ($min,$max,$add)=@_;
  if($max>9000){die"query too long to pile up"}
 
  #draw_coordinate
  for($min..$max){my $unit=$_%10;print"$unit"}
  print"\n";
  print " "x(10-$min%10);
  for($min..$max){
    if($_%10==0){print"$_";print " "x(10-&digit_num($_)) }
  }
  print"\n";
  print "-"x($max-$min+1),"\n";
 
  #anchor fragments 
  foreach(@{$add}){
    my ($start,$end)=split /-/,$_;
    for($min..$start-1){print" "}
    for($start..$end){print"x"}
    for($end+1..$max){print" "}
    print"\n";
  }
 
}#end of sub

sub digit_num(){
 
  my $number=$_[0];
  my $bits=1;
  while(1){
   if($number/10>=1){$bits++;$number/=10;}else{last}
  }
  return $bits;
}#end of sub


