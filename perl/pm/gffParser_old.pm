#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

my $address;
$address=&gffParse($ARGV[0]);

######################################
sub gffParse(){
###Read in standard gff 9 columes file in a line by line way, rather than read all in and classify, and add into the gene table immediately; After finish to read, do sort thing ,reconstrute gene detailsi and do some statistics
#gff in, %genes  out
#unsorted will be ok(children feature may come first); 
#only for single chromosome gff annotation
#simplified structure string:1, >gene  ID:Chr:Strand:Region:description
#                            2, mRNA ID:Region exon1:exon2:...:exonN CDS1:CDS2:CDS3:...CDSn 5UTR:3UTR    
#                            3, mRNA ID:Region exon1:exon2:...:exonN CDS1:CDS2:CDS3:...CDSn 5UTR:3UTR    
#                            ...
#                            n, mRNA ID:Region exon1:exon2:...:exonN CDS1:CDS2:CDS3:...CDSn 5UTR:3UTR    
#allowed feature key words: gene->mRNA-> exon/CDS,UTR. Case nonsensitive
#allowed attrib key words: ID=;Note=;Parent=;
  my $cnt=0;
  my %genes;
  open GFF, (shift @_) or die "$!";
  ##1,readin and store gff file by line
  while(<GFF>){
    chomp;
    next if ($_=~/^#|^$/);
    $cnt++;
    my ($chr,$source,$feature,$start,$end,$score,$strand,$frameShift,$attrib)=split/\t/,$_;
    #make sure that start < end, or switch them, but continue
    if($start>=$end){my $temp=$start;$start=$end;$end=$temp;print STDERR "check coords at line $cnt\n";}
    my ($id,$id1,$id2,$desc);
    #get gene name accordingly (this is the most variable part)
    if($feature=~/^(gene)$/i){ 
         ($id1,$desc)=&attribSplit($attrib,"ID","Note"); 
        }elsif($feature=~/^(mRNA)$/i){ 
             ($id2,$id1)=&attribSplit($attrib,"ID","Parent"); 
            }elsif($feature=~/^(exon|CDS|UTR)$/i){ 
                   ($id2)=&attribSplit($attrib,"Parent"); 
               }else{  
                     print STDERR "Unrecognized feature at $attrib\n"; exit;
                    }
    if(defined $id1){$id=$id1}elsif(defined $id2){$id2=~/(.+)\.\d+/;$id=$1}
    #read a single line finished, start to initializa or fill in gene table;
    if(!exists $genes{$id}){ 
        my (@mRNA,@exons,@cds,@UTR);
        my %structure=(   
               "id"=>"", 
               "chr"=>"", #share by gene mRNA exon cds utr
               "strand"=>"", #share the same things
               "gene"=>"",  #region
               "mRNA"=>\@mRNA, #AS models region $id2:start-end
               "exons"=>\@exons, #exon regions
               "CDS"=>\@cds,    #cds regions
               "UTR"=>\@UTR,    #utr regions
               "function"=>"",  
               "linkto"=>"",      #reserved variety to store gene to gene relationship
               "string"=>""      #recode gff gene blocks to a formatted string
         );  
        #start to fill structure table for the first time
        if($feature=~/^(gene)$/i) {
                              $structure{"id"}=$id ;
                              $structure{"chr"}=$chr;
                              $structure{"strand"}=$strand;
                              $structure{"gene"}="$start-$end";
                              $structure{"function"}=$desc;
           }elsif($feature=~/^(mRNA)$/i) {
                                $structure{"id"}=$id;
                                $structure{"chr"}=$chr;
                                $structure{"strand"}=$strand;
                                push @{$structure{"mRNA"}}, "$id2:$start-$end";
              }elsif($feature=~/^(exon)$/i){
                                  $structure{"id"}=$id;
                                  $structure{"chr"}=$chr;
                                  $structure{"strand"}=$strand;
                                  push @{$structure{"exons"}}, "$id2:$start-$end";
                  }elsif($feature=~/^(CDS)$/i){ 
                                  $structure{"id"}=$id;
                                  $structure{"chr"}=$chr;
                                  $structure{"strand"}=$strand;
                                  push @{$structure{"CDS"}}, "$id2:$start-$end";
                       }elsif($feature=~/^(UTR)$/i){
                                   $structure{"id"}=$id;
                                   $structure{"chr"}=$chr;
                                   $structure{"strand"}=$strand;
                                   push @{$structure{"UTR"}}, "$id2:$start-$end";
                          }#fill end
         #bind to $id
         $genes{$id}=\%structure;
    #if exsits, fill in (check and fill)
    }else{ 
        if($feature=~/^(gene)$/i) {
                  if(${$genes{$id}}{"id"} ne $id || ${$genes{$id}}{"chr"} ne $chr || ${$genes{$id}}{"strand"} ne $strand){
                       die "err in $id\n"
                      }else{
                             ${$genes{$id}}{"gene"}="$start-$end";
                             ${$genes{$id}}{"function"}=$desc;
                            }
           }elsif($feature=~/^(mRNA)$/i) {
                      if(${$genes{$id}}{"id"} ne $id || ${$genes{$id}}{"chr"} ne $chr || ${$genes{$id}}{"strand"} ne $strand){
                          die "err in $id\n"
                         }else{
                                push @{${$genes{$id}}{"mRNA"}}, "$id2:$start-$end";
                               }
              }elsif($feature=~/^(exon)$/i){ 
                      if(${$genes{$id}}{"id"} ne $id || ${$genes{$id}}{"chr"} ne $chr || ${$genes{$id}}{"strand"} ne $strand){
                          die "err in $id\n"
                         }else{
                                push @{${$genes{$id}}{"exons"}}, "$id2:$start-$end";
                               }
                  }elsif($feature=~/^(CDS)$/i){ 
                      if(${$genes{$id}}{"id"} ne $id || ${$genes{$id}}{"chr"} ne $chr || ${$genes{$id}}{"strand"} ne $strand){
                          die "err in $id\n"
                          }else{
                                  push @{${$genes{$id}}{"CDS"}}, "$id2:$start-$end";
                                } 
                    }elsif($feature=~/^(UTR)$/i){
                      if(${$genes{$id}}{"id"} ne $id || ${$genes{$id}}{"chr"} ne $chr || ${$genes{$id}}{"strand"} ne $strand){
                           die "err in $id\n"
                           }else{
                                 push @{${$genes{$id}}{"UTR"}}, "$id2:$start-$end";
                                } 
                      }#inner if end here 
         }#fill structure table end here
  }#while ends here, next entry line
  close GFF; 
#  print Dumper \%genes;

  ##2, All lines are filled and stored: caculate "string" 
  foreach my $key(sort keys %genes){
     my($geneString,@mRNAString);
     $geneString="gene"." ".$key.":".${$genes{$key}}{"chr"}.":".${$genes{$key}}{"strand"}.":".${$genes{$key}}{"gene"}.":".${$genes{$key}}{"function"};
     $genes{$key}->{"mRNA"}=&mrnaSort($genes{$key}->{"mRNA"});
     foreach my $ctrl(@{$genes{$key}->{"mRNA"}}){
       my ($id,$region)=split/:/,$ctrl; 
       my @exons=&CDSSort($id,$genes{$key}->{"exons"});
       my $exons=join":",@exons;
       my @CDS=&CDSSort($id,$genes{$key}->{"CDS"});
       my $CDS=join":",@CDS;
       push @mRNAString, "mRNA $id:$region $exons $CDS";
     }#output mRNA string end
     my $mRNAString=join"\n",@mRNAString;
     $genes{$key}->{"string"}=">$geneString\n$mRNAString\n";
     print ">$geneString\n$mRNAString\n";
 }#foreach end here
  #print Dumper \%genes; 

 return \%genes;
}#end of gffparser()

###sorts array with elements like LOC_Os04g47240.1:28049241-28057106
sub mrnaSort(){ # sort by mRNA ID
   my $index=shift;
   #print "before:\n","@{$index}\n";
   my  @list=sort{
         my ($id1,$n1,$region1,$id2,,$n2,$region2);
        ($id1,$n1,$region1)=split/\.|:/,$a;
        ($id2,$n2,$region2)=split/\.|:/,$b;
        if($n1>$n2 && $id1 eq $id2){1}elsif($n1<$n2 && $id1 eq $id2){-1}elsif($n1==$n2 && $id1 eq $id2){0}
       } @{$index};
   #print "after:\n","@list\n";
   return \@list;
}#end of sub mrnaSort

###sorts array with elements like 'LOC_Os04g47240.1:28054414-28054958'
sub CDSSort(){ # ascending sort, without strand consideration
   my $id=shift;
   my $index=shift;
   #print "before:\n","@{$index}\n";
   my @list;
   foreach (@{$index}){
     my ($id1,$region)=split/:/,$_;
     if($id eq $id1){push @list,$region}
   }
   @list=sort{
         my ($s1,$e1,$s2,$e2);
         ($s1,$e1)=split/-/,$a;
         ($s2,$e2)=split/-/,$b;
         if($s1>$s2 && $s1<$e1 && $s2<$e2){1}elsif($s1<$s2 && $s1<$e1 && $s2<$e2){-1}elsif($s1==$s2 && $s1<$e1 && $s2<$e2){0}
        } @list;
   #print "after:\n","@list\n";
   return @list;
}#end of sub CDSSort

#get value(s), according to key(s)
sub attribSplit(){
  my %attrib;
  my $string=shift @_;
  my @fkey=@_;
  my @report;
  my @temp=split/;/,$string;
  foreach(@temp){
     my ($fkey,$fval)=split/=/,$_;
     if(!exists $attrib{$fkey}){ $fval=~s/\s+$//;$attrib{$fkey}=$fval;}else{print"repeat attrib names $fkey\n";}
  }
  foreach(@fkey){
   if(exists $attrib{$_}){push @report, $attrib{$_}}else{die "can't find attrib key in colume 9 at $string\n"}
  }
  return @report; 
}#end of sub attribSplit

