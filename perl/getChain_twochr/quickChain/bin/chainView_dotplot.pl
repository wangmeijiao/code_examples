#!/bin/perl
use strict;
use warnings;
use Getopt::Long;
use SVG;
#usage:  cat in.chain |grep -v "^#" |perl chainView.pl -refLen xxx -queLen xxx|display
 my ($len_refchr,$len_quechr);
 GetOptions( "refLen=i"=>\$len_refchr,
             "queLen=i"=>\$len_quechr,
           );
#globe variates
 my ($width,$height)=(600,600*$len_quechr/$len_refchr);
 my $edge=10;
 my $ppx=$width/$len_refchr;  #x pix/bp;

#prepare drawing board
 my $svg=SVG->new(height=>$height+2*$edge,width=>$width+2*$edge);
 #fitted line
 #$svg->line(x1=>$edge,y1=>$edge,x2=>$edge+$width,y2=>$edge+$height,stroke=>"green","stroke-dasharray"=>"1,1");

$/="\n\n";
my $cnt=0;
my $i=0;
my @colors=@{&colors()};print STDERR @colors;
while(<stdin>){
  chomp;
  $cnt++;
  #next if($cnt==1);
  #print STDERR "$_\n";
  my @box=split/\n/,$_;
  ##draw chain frame
  my $align=shift @box;
  my ($keyword,$score,$refChr,$refLen,$refStrand,$refStart,$refEnd,$queChr,$queLen,$queStrand,$queStart,$queEnd,$blockId)=split/ / ,$align;
  #print STDERR "$keyword,$score,$refChr,$refLen,$refStrand,$refStart,$refEnd,$queChr,$queLen,$queStrand,$queStart,$queEnd\n";
  die"err at $_\n" if($keyword ne "chain" || $refLen!=$len_refchr || $queLen!=$len_quechr);
  #$svg->line(x1=>$edge+$refStart*$ppx,y1=>$edge+$queStart*$ppx,x2=>$edge+$refEnd*$ppx,y2=>$edge+$queEnd*$ppx,stroke=>"grey","stroke-width",2);   
  ##draw details of a certain chain frame
  my($block_refStart,$block_queStart)=($refStart,$queStart);
  foreach(@box){
    my ($blockSize,$dt,$dq)=split/\t/,$_;
    next if(!defined $dt);
    #print STDERR "$blockSize,$dt,$dq\n";
    my($block_refEnd,$block_queEnd)=($block_refStart+$blockSize,$block_queStart+$blockSize);
    $svg->line(x1=>$edge+$block_refStart*$ppx,y1=>$edge+$block_queStart*$ppx,x2=>$edge+$block_refEnd*$ppx,y2=>$edge+$block_queEnd*$ppx,stroke=>$colors[$i],"stroke-width",4);   
    ($block_refStart,$block_queStart)=($block_refEnd+$dt,$block_queEnd+$dq);
    if($i==$#colors){$i=0}else{$i++};
  }
  $svg->line("x1"=>$edge+$refStart*$ppx,"y1"=>$edge+50+$queStart*$ppx,"x2"=>$edge+$refEnd*$ppx,"y2"=>$edge+50+$queEnd*$ppx,"stroke"=>"red");
}#while end here
$/="\n";
print $svg->xmlify();


#####sub###########
sub colors(){

 my @color=("#1B9E77", "#61864B", "#A66F20", "#CE6014", "#A96755" ,"#846D97", "#8D61AA","#B7469B", "#E12C8C", "#BE5067", "#8E7E40", "#6CA61C", "#9BA812", "#CBA907","#DBA206", "#C48F10", "#AC7B1A" ,"#957130", "#7D6B4B" ,"#666666");
 #my @color=("#B7469B", "#7D6B4B" );
 return (\@color);
  
}
