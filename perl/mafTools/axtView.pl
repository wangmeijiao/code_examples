#!/bin/perl
use strict;
use warnings;
use Getopt::Long;
use SVG;

#usage:  cat final.chainnet.sorted.axt |grep -v "^#" |perl axtView.pl len_refchr len_quechr |display
#                                                                      mjwang, Dec. 11
#draw axt with xyplot and collinear plot
# axt format are orgnized special for pairwise alignment;
# for multiplex alignment use maf;
#data example:
# 12 OsjapChr04 45094 45152 chr04 10120927 10120984 + 667
# TAAcactggtggagaaaccctttgtagtcccggtttgtaaccccc--ctttagtcccggtt
# TAAGACCGCCGCCGGCGCCCTTTATA---CTGGTTTAGGCCCGCCAGCGCTAATCAAGGTC

 my ($len_refchr,$len_quechr);
 GetOptions( "refLen=i"=>\$len_refchr,
             "queLen=i"=>\$len_quechr,
           );

#globe variates
my $ratio=$len_refchr/$len_quechr;# width/height
my $height=2000;  
my $width=$height*$ratio;
my $len_ref=$width;
my $factor=$len_ref/$len_refchr; #zoom factor
my $len_que=$len_quechr*$factor;
my $ref_y=200;
my $que_y=400;

#prepare drawing board
 my $svg=SVG->new(height=>$height,width=>$width);
 $svg->rect(x=>100,y=>$ref_y-5,height=>5,width=>$len_ref); #refchr
 $svg->rect(x=>100,y=>$que_y,height=>5,width=>$len_que); #quechr
 $svg->line(x1=>$width*0.75,y1=>$height*0.75,x2=>$width*0.75+1000000*$factor,y2=>$height*0.75,stroke=>"red"); #1Mb ruler
 $svg->text(x=>$width*0.75+500*$factor,y=>$height*0.75-5,"-cdata"=>"1Mb");
 #$svg->line(x1=>0,y1=>2,x2=>$width-2,y2=>2,stroke=>"black",); #ref_coord
 #$svg->text(x=>$width-100,y=>20,"-cdata"=>"Refer:$len_refchr");
 #$svg->line(x1=>2,y1=>0,x2=>2,y2=>$height-2,stroke=>"black"); #que_coord
 #$svg->text(x=>20,y=>$height-20,"-cdata"=>"Query:$len_quechr","writing-mode"=>"tb");#writing-mode doesn't work in display but works well in chrome
 #fitted line
 #$svg->line(x1=>0,y1=>0,x2=>$len_refchr*$factor,y2=>$len_quechr*$factor,stroke=>"green","stroke-dasharray"=>"1,1");

$/="\n\n";
while(<stdin>){
   chomp;
   #print "$_\n";
   my @box=split/\n/,$_;
   my $align=shift @box;
   my ($blockID,$refID,$refStart,$refEnd,$queID,$queStart,$queEnd,$strand,$score)=split/ / ,$align;
   #my $ref=shift @box;
   #my $query=shift @box;
 
   #draw collinear 
    $svg->line(x1=>$refStart*$factor+100,y1=>$ref_y,x2=>$queStart*$factor+100,y2=>$que_y,stroke=>"grey","stroke-width"=>0.1);     $svg->line(x1=>$refEnd*$factor+100,y1=>$ref_y,x2=>$queEnd*$factor+100,y2=>$que_y,stroke=>"grey","stroke-width"=>0.1);   
   #draw xyplot   
    #$svg->line(x1=>$refStart*$factor,y1=>$queStart*$factor,x2=>$refEnd*$factor,y2=>$queEnd*$factor,stroke=>"black");   
   
}#while end here
$/="\n";
print $svg->xmlify();

