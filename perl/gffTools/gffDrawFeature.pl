#!/usr/bin/env perl
use strict;
use warnings;
use SVG;
use Getopt::Long;
#read gfflist file and draw feature box alone chromosome
#useage: cat gfffile.gff.list |perl thisfile -chrLen xxxx |display
#        cat gfffile.gff.list |perl thisfile -chrLen xxxx --region 230000-240000|display

my ($chrLen,$region);
GetOptions(
         "chrLen=s"=>\$chrLen,
         "region=s"=>\$region,
);

my %genes;
my ($chrStart,$chrEnd);
if(defined $region){($chrStart,$chrEnd)=split/-/,$region; die if($chrStart>=$chrEnd)}else{$chrStart=1;$chrEnd=$chrLen}
$chrLen=$chrEnd-$chrStart+1;
my $width=12000;
my $height=300;
my $factorWH=sprintf ("%.2f",$width/$height);
my $space=50;
my $barThick=10;
my $chrLenDraw=$width-2*$space;
my $factorZ=$chrLenDraw/$chrLen; 
my $chrPos=$height*0.7;
my $featurePosPlus=$height*0.6;
my $featurePosMinu=$height*0.65;

my $svg=SVG->new(width=>$width,height=>$height);
#draw chrom
$svg->rect(x=>$space,y=>$chrPos-$barThick/2,height=>$barThick,width=>$chrLenDraw,fill=>"grey","stroke"=>"none");
#draw metering system
$svg->line(x1=>$space,y1=>$height*0.73,x2=>$space+$chrLenDraw,y2=>$height*0.73,stroke=>"grey","stroke-width"=>0.2);
my $cnt=-1;
for(my $i=$space;$i<=$space+$chrLenDraw;$i+=1000000*$factorZ){
$svg->line(x1=>$i,y1=>$height*0.73,x2=>$i,y2=>$height*0.73+5,stroke=>"grey","stroke-width"=>0.3); 
$cnt++;
if($cnt%5==0){$svg->text(x=>$i-2,y=>$height*0.73+12,"-cdata"=>$cnt,"font-size"=>7);}}
#draw scale ruler
$svg->line(x1=>$space,y1=>$height*0.8,x2=>$space+1000000*$factorZ,y2=>$height*0.8,stroke=>"black");
$svg->text(x=>$space+1000000*$factorZ/3,y=>$height*0.82,"-cdata"=>"1Mb","font-size"=>7);

#1,read in gfflist and prepare coord data. draw feature along chromosome
$/=">";
while(<stdin>){
 chomp;
 next if ($_ eq "");
 my @box=split /\n/, $_;
 #draw gene region
 my $gene=shift @box;
# my ($mark,$string)=split / /,$gene;
# if($mark=~/^(gene)$/i){
#   my ($geneID,$chr,$strand,$region,$desc)=split/:/,$string;
#   my ($start,$end)=split/-/,$region;
#   if($end<$chrStart || $start>$chrEnd){next}
   #if($strand eq "+"){$svg->rect(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosPlus-$barThick/2, width=>($end-$start)*$factorZ,height=>$barThick*0.5,fill=>"grey","fill-optical"=>0.5,stroke=>"none")}
   #elsif($strand eq "-"){$svg->rect(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosMinu-$barThick/2, width=>($end-$start)*$factorZ,height=>$barThick*0.5,fill=>"grey","fill-optical"=>0.5,stroke=>"none")}   
#  } 
 #draw mRNA structure 
 foreach(@box){
  my ($mark,$string)=split / /,$_;
  if($mark=~/^(mRNA)$/i){ 
   my ($mrnaID,$strand,$region,$exon,$CDS)=split/:/,$string;
   my ($start,$end)=split/-/,$region;
   my @exon=split/\|/,$exon;
   my @CDS=split/\|/,$CDS;
   if($end<$chrStart || $start>$chrEnd){next}
   if($strand eq "+"){$svg->rect(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosPlus-$barThick/2, width=>($end-$start)*$factorZ,height=>$barThick*0.5,fill=>"grey","fill-optical"=>0.5,stroke=>"none");$svg->text(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosPlus-$barThick*0.75,"-cdata"=>$mrnaID,);}
   elsif($strand eq "-"){$svg->rect(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosMinu-$barThick/2, width=>($end-$start)*$factorZ,height=>$barThick*0.5,fill=>"grey","fill-optical"=>0.5,stroke=>"none");$svg->text(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosMinu-$barThick*0.75,"-cdata"=>$mrnaID,);}
   #draw details
   foreach(@CDS){
     my ($start,$end)=split/-/,$_;
     if($strand eq "+"){$svg->rect(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosPlus-$barThick*0.75, width=>($end-$start)*$factorZ,height=>$barThick,fill=>"green","fill-optical"=>0.5,stroke=>"none");}
   elsif($strand eq "-"){$svg->rect(x=>($start-$chrStart)*$factorZ+$space,y=>$featurePosMinu-$barThick*0.75, width=>($end-$start)*$factorZ,height=>$barThick,fill=>"red","fill-optical"=>0.5,stroke=>"none")}
    }
  }#if mRNA end here
   $featurePosPlus-=$barThick*2;
   $featurePosMinu-=$barThick*2; 
 }#foreach @box end here
   $featurePosPlus=$height*0.6;
   $featurePosMinu=$height*0.65;
}#while end here
print $svg->xmlify();

