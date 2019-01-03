#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#give me the raw input maf file as it was(better projected and ordered, but it is OK if not,single chromosome file), I will do neccesary things:
#      remove head lines;
#      remove single line blocks;
#      remove length<=minlen blocks(becareful to use this, since you may loose information, e.g small but nearby conserved              blocks);
#      output pairs of MSPcrunch files;
#      output maf_block position files of all species

##maf format:
#a score=21724.0
#s japo    6343 162 + 19604000 CGGAAGCAGGGGCGGAATTTTT---TTTTTGCTATACATCATCTATTTTCACCATGTATGCACCCCCCTCGATACAAAC
#s glab 6429557 161 - 18621000 CAGGAACCTGGTTAAAAATTTT---TTTTAGCTATACCTCATCTATTTTTACCATATATGCACCCTCCTCAATATAAAC
#s brac 3828196 164 - 16420000 CAGGAGTGGAGGTACATATATTGAATTTTAGCTAGAGTTTATTAATTTTCACCATGTACGTA-TCTCCTTAATCTAAAT

##MSPcruch format:use *space* to seprate columns; score must be integer
#   score, %identity, qstart, qend, qid, sstart, send, sid
##  828    95.0   2837  3841  brac  7211 8276 japo

##blast m8 format (not used)
#  $queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore



##gtf format (almost the same with gff)
# ##gff-version 2 
# ##source-version <source> <version text> 
# ##date <date> 
# <seqname>, <source>, <feature>, <start>, <end>, <score>, <strand>, <frame>, [attributes] [comments];
#  10      protein_coding  exon    22495   22586   .       -       .        gene_id "ORGLA10G0000200"; transcript_id "ORGLA10G0000200.1"; exon_number "8"; seqedit "false";
 
my $help=<<HELP;
*Usage:perl bin/maf2act.pl --species "japo glab brac" --minlen 100 --maf all.projected.sorted.maf --mspcrunch-dir msp --gtf-dir gtf
HELP
my %opts;
GetOptions(\%opts,"species=s","minlen=i","maf=s","mspcrunch-dir=s","gtf-dir=s","help!");
if($opts{help}||keys %opts<5){print $help;exit;}
 

#1,prepare variable number of files to write into. Is this a dynamic program topic?? Unnecessary here. Hash tables can play a good role of containers
&checkDir($opts{'mspcrunch-dir'});
&checkDir($opts{'gtf-dir'});
my (%gtf,%msp);
my @species=split/ /,$opts{'species'};
my $gtf_file_num=$#species+1;
my $nspecies=$gtf_file_num;
my $msp_pair_num=&fac($nspecies)/(&fac($nspecies-2)*&fac(2)); 

 #initialize hash_box for data strorage 
  foreach (@species){
    my @tmp;
    if(exists($gtf{$_})){print "err, duplicate species\n";exit}else{$gtf{$_}=\@tmp;}
  }

  my @species_map=@species;
  foreach my $ctrl_out(@species){  
    shift @species_map;
    foreach my $ctrl_in(@species_map){
       my $msp_pair;
       $msp_pair=$ctrl_out.'-'.$ctrl_in;
       my @tmp;
       if(exists($msp{$msp_pair})){print "err, duplicate species\n";exit}else{$msp{$msp_pair}=\@tmp;}
     }#inner foreach end here 
  }#outter foreach end here


  #print "\nspecies:@species,$nspecies; gtf_file_num:$gtf_file_num; msp_pair_num:$msp_pair_num\n";
  #print Dumper \%gtf;
  #print "\n\n";
  #print Dumper \%msp;
  


#2, read maf blocks and caculate coordinates. Push data lines
open MAF, $opts{maf} or die "$!";
$/="\n\n";
my $cnt;
while(<MAF>){
 
  chomp;
  my $flag=0;
  my @box=split /\n/,$_;
  $cnt++;
  #last if ($cnt==10);
  my $block_score;
  my %block;#specises name stand for the key 
  foreach my $line(@box){
    #exit if ($line=~/##eof/);
    last if ($line=~/##eof/); #modified on 12, May 2013
    next if ($line=~/#/);
    if($line=~/a score/){
        my @tmp=split/=/,$line;$block_score=int($tmp[1]);
      }elsif($line=~/s /){
            my @tmp=split/ +/,$line;
            my $name=$tmp[1];
            my $start=$tmp[2];
            my $len=$tmp[3];if ($len<$opts{"minlen"}){$flag=1};
            my $strand=$tmp[4];
            my $chr_len=$tmp[5];
            my $seq_ali=$tmp[6];
            my $seq_ori=$seq_ali; $seq_ori=~s/-//g;
            my $end=$start+$len;
            my @dataline=($cnt,$name,$start,$end,$strand,$block_score,$seq_ali,$seq_ori); 
            if(exists $block{$name}){
                die"duplicated species in maf block $cnt\n";
              }else{
                   $block{$name}=\@dataline;
                 }
       }else{die "\nerr\n";} 
  }#foreach end here  
   
   #do some filering
   if ($flag==1){print STDERR "MESSAGE:Block $cnt has been discarded since short block length, where first line reads \"$box[1]\"\n";next};
   if (keys %block <=1){print STDERR "MESSAGE:Block $cnt has been discarded since single line,where first line reads \"$box[1]\"\n";next};
 
  #beging to write gff and msp lines from %block
   #print "\nBlock $cnt in maf file:\n";
   #print Dumper \%block;
 
   my @gtflist=keys %block;
   my @msplist=&permutation(keys %block); #here need more test....
  
   foreach (@gtflist){
       if (exists $block{$_}){
            my $gtfline=join("\t",(${$block{$_}}[1],"tba","exon",${$block{$_}}[2],${$block{$_}}[3],${$block{$_}}[5],${$block{$_}}[4],'.','block_id'.' '.$cnt.';'));
            push @{$gtf{$_}},$gtfline;           
         }else{die"gtf file name and %block key don't match";}
    }   

   foreach (@msplist){
       my @tmp=split/-/,$_;
       my $query_name=$tmp[0];
       my $subject_name=$tmp[1];
       my $identity=&idenCalc(${$block{$query_name}}[6],${$block{$subject_name}}[6]);
       my $mspline=join(" ",(${$block{$query_name}}[5],$identity,${$block{$query_name}}[2],${$block{$query_name}}[3],$query_name,${$block{$subject_name}}[2],${$block{$subject_name}}[3],$subject_name))."\n";
       if(exists $msp{$_}){
          push @{$msp{$_}},$mspline;
         }else{
              my ($l,$r)=split/-/,$_;
              my $new=$r."-".$l;
              if(exists $msp{$new}){ 
                 my $mspline_rv=join(" ",(${$block{$query_name}}[5],$identity,${$block{$subject_name}}[2],${$block{$subject_name}}[3],$subject_name,${$block{$query_name}}[2],${$block{$query_name}}[3],$query_name))."\n";
                 push @{$msp{$new}},$mspline_rv;
               }else{die"err";}
          }
   }#foreach end here 

}#while end here
$/="\n";
close MAF;

  #print"\n\nfinal gtf and msp file.\n";
  #print Dumper \%gtf;
  #print "\n\n";
  #print Dumper \%msp;

#3, write output files into specified dirs, use hash keys as file names
 
 &writeFile(\%gtf,\%msp,$opts{"gtf-dir"},$opts{"mspcrunch-dir"});
 
 
##########

sub checkDir(){
  #if exsits, ignor; if not, create it.
  my $name=$_[0];
  if (opendir(DIR,$name)){}else{mkdir($name,0755)}
  close DIR;
 }

sub fac(){
  my $num=1;
  my $para=shift;
  if($para==0){
     $num=0;
    }else{
        for(1..$para){$num*=$_;}
      }
  return $num;
}

sub idenCalc(){
  my($query,$subject)=@_;
  #print "query:subject $query:$subject\n";
  $query=~tr/agctn/AGCTN/;
  $subject=~tr/agctn/AGCTN/;
  my @query=split//,$query;
  my @subject=split//,$subject;
  my $cnt;
  for(0..$#query){
    if ($query[$_] eq $subject[$_]){$cnt++;}
  }
  my $iden=$cnt/(length $query)*100;
  return $iden;
 
}

sub permutation(){
   my @elements=@_;
   my @elements_map=@elements;
   my @list;
   foreach my $ctrl_out(@elements){
     shift @elements_map;
     foreach my $ctrl_in(@elements_map){
       my $pair;
       $pair=$ctrl_out.'-'.$ctrl_in;
       push @list,$pair;
     }#inner foreach end here 
   }#outter foreach end here
   return @list;
}

sub writeFile(){
 my ($gtf,$msp,$gtfdir,$mspdir)=@_;
 foreach (keys %{$gtf}){
   my $filename=$gtfdir."/".$_.".gtf";
   open GTF, ">$filename" or die "$!";
   print GTF "##gff-version 2\n";
   print GTF join("\n",@{$$gtf{$_}});
   close GTF;
 }
 foreach(keys %{msp}){
   my $filename=$mspdir."/".$_.".mspcrunch";
   open MSP, ">$filename" or die "$!";
   print MSP join("\n",@{$$msp{$_}});
   close MSP;
 }
}
