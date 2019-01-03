#!usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
#usage: cat in.rmsk.out |perl thisfile.pl < -filter 0.5(coverage of target repeat hit)> < -join > < -stat > < -bed > <-prefix rigr6> <-outdir outdir>  2> lost.err 
#note: coverage filtered groups output nothing (assuming false discovered by repeatMasker)
#      groups failed in qID, strand or target_hit consistancy checking , report to stderr(lost.gff/bed, without segment linking or coverage filtering, use with caution)
#      groups members have different tleft length, choose the first one to calculate coverage_total%
#      output gff by default

########method#########
#1,readin rmsk.out and group by id (%te)
#2,link fragments in the same group to get the summary line
#3,filter low coverage_total or identity_ave group
#4,classify by TE class (%class)
#5,join nearby( <100bp) groups in the same class and do statistics by class
#6,output gff(nest gff) or bed(simple bed or detail bed) by class


########TE class namenclature####
#repeatMasker repbase classification levels: family->subfamily->repeatID (targetID#family/subfamily/subsubfamily)
#msu rice repeatbase code format: superclass->class->subclass->repeatID (Transposable Elements->Retrotransposons->Ty3-gypsy->repeatID, example :ORSg-TERT002-00703)


#######version#########
#Nov.19 2015 modified method to parse LTR structure and update to version 2.2
#Nov.10 2014 by mjwang


my($filter,$join,$stat,$prefix,$bed,$outdir);
GetOptions("filter:f",\$filter,"join",\$join,"stat",\$stat,"bed",\$bed,"prefix:s",\$prefix,"outdir:s",\$outdir);
my %te;
<stdin>;<stdin>;<stdin>; #skip head 3 lines, for RepeatMasker version open-4.0.3
#1,readin rmsk.out, preprocess and group by id (%te)
while(<stdin>){
  chomp;
  next if($_ eq "" ||$_=~/^#/);
  $_=~s/^ +| +$//;
  my @box=split/ +/,$_;
  print STDERR "column number err at $_\n" if(scalar @box !=16 && scalar @box != 15);
  ##rmsk_out format:
  #score div del ins   que    qs    qe    (left)  strand target_in_lib   class/family  (left) te   ts    id  mark
  # 887  6.6 1.5 0.0 H1_japo 18490 18626 (2781375)   C       Gaijin    DNA/PIF-Harbinger (6)  141   3    22   
  # 883  6.1 1.1 0.0 H1_japo 18290 18326 (2781443)   +       Gaijin    DNA/PIF-Harbinger  1   362 (234)  21   
  # 883  5.0 1.3 0.0 H1_japo 18320 18420 (2781425)   +       Gaijin    DNA/PIF-Harbinger  358  376 (120)  21  *
  #
  #  0    1   2   3     4      5     6       7       8          9             10          11   12  13 14 (15)
  # mark * indicates that there is a higher-scoring match whose domain partly (<80%) includes the domain of this match 
  # id will be the same if match multi-fragments
  my($score,$div,$del,$ins,$qID,$qs,$qe,$qleft,$strand,$tID,$class,$ts,$te,$tleft,$id,$mark)=@box;
  ($ts,$te,$tleft,$strand)=($strand eq "C")?($tleft,$te,$ts,"-"):($ts,$te,$tleft,$strand);
  $qleft=~s/\(|\)//g;
  $tleft=~s/\)|\(//g;
  if(!defined $mark){$mark="NULL"}
  die "check start>end at line @box\n" if($qs>=$qe || $ts>=$te);
  push @{$te{$id}},join(":",($score,$div,$del,$ins,$qID,$qs,$qe,$qleft,$strand,$tID,$class,$ts,$te,$tleft,$id,$mark)); 
}

#2,link fragments in the same group to get the summary line
#3,filter low coverage_total or identity_total groups
#4,classify by TE class (%class: family->id(group)->lines)
my %class;
my %lost;
foreach my $id(sort{$a<=>$b} keys %te){
   my @group=@{$te{$id}};
   #if multi-fragments hits
   if(scalar @group >1){ 
     #check: all fragments must have the same strand, class or qID, output to stderr if not
     my $flag_break=0;
     for(0..$#group-1){
        my @box1=split/:/,$group[$_];
        my @box2=split/:/,$group[$_+1];
        #if($box1[4] ne $box2[4] ||  $box1[10] ne $box2[10]){ #make sure that chr and class/family eq 
        if($box1[4] ne $box2[4] || $box1[8] ne $box2[8]  || $box1[10] ne $box2[10]){ #make sure that chr and class/family eq 
           print STDERR "chr, strand or class/family not the same within a group at group $id(output to lost.gff/bed): \n  @box1\n  @box2\n"; 
           $flag_break=1;
        }    
     }
     if($flag_break){$lost{$id}=$te{$id};next} #next id 
     #sort
     my @group_sorted=sort{ #sort by query start
                         my @box1=split/:/,$a;  
                         my @box2=split/:/,$b;  
                         if($box1[5]>$box2[5]){return 1}elsif($box1[5]<$box2[5]){return -1}elsif($box1[5]==$box2[5]){return 0}
                      } @group;
  
     #find outter region of query & calculate total coverage for target, and filter low coverage groups; finally get summary line
     my @segments_query;
     my @segments_target;
     my $coverage_total;
     my ($identity_ave,$identity_total,$cnt);
     my @firstLine=split/:/,$group_sorted[0]; #get first line to represent the whole group
     foreach(@group_sorted){
       my($score,$div,$del,$ins,$qID,$qs,$qe,$qleft,$strand,$tID,$class,$ts,$te,$tleft,$id,$mark)=split/:/,$_;       
       push @segments_query,"$qs-$qe";
       push @segments_target,"$ts-$te:$tleft";
       $identity_total+=100-$div-$del-$ins;
       $cnt++;
     }
     #$identity_ave=sprintf "%.3f",$identity_total/$cnt;
     my($leftMostQuery,$rightMostQuery)=&minMax(\@segments_query);
     #$coverage_total=&getCoverage(\@segments_target);
     #filter by coverage of target hit repeat
     #if(defined $filter){next if($coverage_total< $filter*100)};
     my $name=$firstLine[9];
     $name=~s/-LTR$|-I$|-I-int$|_LTR$|_I$|_I-int$//i; #rm "Copia-104_SB-I" tail
     my $attribs="ID=".$id."#".$firstLine[10].";"."Note="."$name#$firstLine[10]";
     #my $attribs="ID=".$id."#".$firstLine[10].";"."Note="."$name#$firstLine[10]|coverage_total:$coverage_total%|identity_ave:$identity_ave%";
     my $summary_line="$firstLine[4]\trpmk2gffv2\tTransposon\t$leftMostQuery\t$rightMostQuery\t.\t$firstLine[8]\t.\t$attribs"; 
     push @group_sorted,$summary_line;
     #store into a class hash
     my ($family,$subfamily,$subsubfamily)=split/\//,$firstLine[10];
     if(!defined $family || $family eq""){die "\t$id:can't identify te class at @firstLine\n"}elsif(!defined $subfamily || $subfamily eq ""){$subfamily="unamed_subfamily";$subsubfamily="unamed_subsub"}elsif(!defined $subsubfamily || $subsubfamily eq ""){$subsubfamily="unamed_subsub"}
      $class{$family}->{$subfamily}->{$subsubfamily}->{$id}=\@group_sorted;
  #if single line group       
  }elsif(scalar @group ==1){
        my($score,$div,$del,$ins,$qID,$qs,$qe,$qleft,$strand,$tID,$class,$ts,$te,$tleft,$id,$mark)=split/:/,$group[0];
        my $identity=100-$div-$del-$ins;
        my $total_len=$te+$tleft;
        my $coverage=sprintf "%.3f",100*($te-$ts+1)/$total_len;
        #if(defined $filter){next if($coverage< $filter*100)};
        my $name=$tID;
        #$name=~s/-LTR$|-I$|-I-int$|_LTR$|_I$|_I-int$//i; #rm "Copia-104_SB-I" tail
        my $attribs="ID=".$id."#".$class.";"."Note="."$name#$class|coverage:$coverage%($ts-${te},${total_len})|identity:$identity%";
        my $summary_line="$qID\trpmk2gffv2\tTransposon\t$qs\t$qe\t.\t$strand\t.\t$attribs"; 
        #store into a class hash
        my ($family,$subfamily,$subsubfamily)=split/\//,$class;
        if(!defined $family || $family eq""){die "\t$id:can't identify te class at $group[0]\n"}elsif(!defined $subfamily || $subfamily eq ""){$subfamily="unamed_subfamily";$subsubfamily="unamed_subsub"}elsif(!defined $subsubfamily || $subsubfamily eq ""){$subsubfamily="unamed_subsub"}
        $class{$family}->{$subfamily}->{$subsubfamily}->{$id}=$summary_line;
     }else{die "empty lines at id $id\n"}
}#foreach groups end

#print Dumper \%class;
#exit;

#5,join nearby( <100bp) groups in the same class (target name will be different),add up target names, do some statistics by class 
if($join){




}

if($stat){



}

#6,output gff(nest gff) or bed(simple bed or detail bed) according to class restored hash
  if(defined $prefix){}else{$prefix="out"}
  if(-d $outdir){}else{mkdir $outdir}
  foreach my$family(sort keys %class){
    my $append;
    if(defined $bed){$append=".bed"}else{$append=".gff"}
    my $outfile=$prefix."_".$family.$append;
    #print "$outfile\n";
    open OUT,">$outdir/$outfile" or die "$!";
    if(defined $bed){}else{print OUT "##gff version 3\n"}
    foreach my $subfamily(sort keys %{$class{$family}}){
      print OUT "##$subfamily\n";
      foreach my $subsub(sort keys %{$class{$family}->{$subfamily}}){
        foreach my $id(sort {$a<=>$b} keys %{ $class{$family}->{$subfamily}->{$subsub} }){
          my $index=$class{$family}->{$subfamily}->{$subsub}->{$id};
          if(ref ($index) eq "ARRAY"){# if multilines
             #summay line
             if(defined $bed){
                my @temp=split/\t/,$index->[-1];
                my (undef,$target,undef,$detail)=split/;|=/,$temp[8];
                #$target=~s/#.+//;
                print OUT "$temp[0]\t$temp[3]\t$temp[4]\tte$target|$detail\t.\t$temp[6]","\n";
               }else{
                      print OUT $index->[-1],"\n";
                      #detail lines (only in gff)
                      foreach(0..$#{$index}-1){
                         my($score,$div,$del,$ins,$qID,$qs,$qe,$qleft,$strand,$tID,$class,$ts,$te,$tleft,$id,$mark)=split/:/,$index->[$_];         
                         my ($ident,$coverage);
                         $ident=100-$div-$del-$ins;
                         die "ts >= te at $index->[$_]" if($ts >= $te);
                         my $total_len=$te+$tleft;
                         $coverage=sprintf("%.3f",100*($te-$ts+1)/$total_len);
                         print OUT "$qID\trmsk2gffv2\trepeat_fragment\t$qs\t$qe\t.\t$strand\t.\tParent=$id#$class;Note=$tID|coverage:$coverage%($ts-${te},${total_len})|identity:$ident%\n";
                      }
                    }
          }else{ #if sigle line group
                if(defined $bed){
                  my @temp=split/\t/,$index;
                  my (undef,$target,undef,$detail)=split/;|=/,$temp[8];
                  #$target=~s/#.+//;
                  print OUT "$temp[0]\t$temp[3]\t$temp[4]\tte$target|$detail\t.\t$temp[6]","\n";
                }else{
                       #print $index,"\n";
                       print OUT $index,"\n"
                     }               
               }
       }#id foreach
      }#subsub foreach
    }#subfamily foreach
    close OUT;
  }#family foreach

#7, output lost 
  my $lost_file;
  if (defined $bed){$lost_file=$prefix."_lost.bed"}else{$lost_file=$prefix."_lost.gff"}
  open LOST, ">$lost_file" or die "$!";
  foreach my $id(sort {$a<=>$b} keys %lost){
      my @group=@{$lost{$id}};     
      foreach(@group){
        my($score,$div,$del,$ins,$qID,$qs,$qe,$qleft,$strand,$tID,$class,$ts,$te,$tleft,$id,$mark)=split/:/,$_;
        if(defined $bed){
            print LOST "$qID\t$qs\t$qe\tte$id|$tID#$class\t.\t$strand","\n";            
           }else{
                  print LOST "$qID\trmsk2gffv2\trepeat_fragment\t$qs\t$qe\t.\t$strand\t.\tID=$id#$tID;Note=$tID#$class\n";
                }
      }
  }
  close LOST;


##############subs################
sub minMax(){
  my $index=shift;
  my ($min,$max)=split/-/,${$index}[0];
  foreach(@{$index}){
    my ($start,$end)=split/-/,$_;
    if($start<$min){$min=$start}
    if($end>$max){$max=$end}
  }
  return ($min,$max);
}

#deprecated
sub getCoverage(){
    my @segs=@{shift @_};
    my @segs_true;
    my @len;
    foreach(@segs){
      my ($segs,$left)=split/:/,$_;
      push @segs_true,$segs;
      my ($s,$e)=split/-/,$segs;
      push @len,$e+$left;
    }
    my $len_total=$len[0];
    foreach(@len){
       if($len_total != $_){ print STDERR "length of target not equal at @segs\n";}
    }
    my $len_nonOverlap;
    my @segs_linked=@{&segLink(\@segs_true)}; 
    foreach(@segs_linked){
        my($s,$e)=split/-/,$_;
        $len_nonOverlap+=$e-$s+1;
    }
    return sprintf("%.3f",100*$len_nonOverlap/$len_total); 
}


sub segLink(){ #the same with deOverlap()
   #standard function to sort&link(deOverlap) segments, which belong to a single sequence
   #use 3seg method if (seg_num >=3) or if (seg_num ==2 ) use 2seg method . do nothing if (seg_num <2)
   #input format @segs=("s1-e1","s2-e2","..-..")
   my @segs=@{shift @_};
   if (scalar @segs<2){
     return \@segs
    }elsif(scalar @segs ==2){
       my @segs_sorted=@{&segSort(\@segs)};
       my @segs_sorted_linked=&relate2segs($segs_sorted[0],$segs_sorted[1]);
       return \@segs_sorted_linked;
      }elsif(scalar @segs >2 ){
         my @segs_sorted=@{&segSort(\@segs)};
         #print "sorted: @segs_sorted\n";
         my @result=shift @segs_sorted;
         for(0..$#segs_sorted){
             my @relate;
             if($_ != $#segs_sorted){
               @relate=&relate2segs($segs_sorted[$_],$segs_sorted[$_+1]);
               @relate=&relate2segs($result[-1],$relate[0]);
               pop @result;
               push @result,@relate;
             }elsif($_ == $#segs_sorted){
                   @relate=&relate2segs($result[-1],$segs_sorted[$_]);
                   pop @result;
                   push @result,@relate;
                 }
         }
         return \@result;
       }
}

sub segSort(){
     my @segs=@{shift @_};
     my @segs_sorted=sort{
       my ($s1,$e1)=split/-/,$a;
       my ($s2,$e2)=split/-/,$b;
       if($s1<$s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
     } @segs;
     return \@segs_sorted;
}

sub relate2segs(){
    #input paras ("s1-e1","s2-e2") must be sorted by start
    my ($s1,$e1)=split/-/,(shift @_);
    my ($s2,$e2)=split/-/,(shift @_);
    die "start1>start2 at $s1-$e1|$s2-$e2\n" if ($s1>$s2);
    #three types of relationship between two segments
    #1,separated
    if($s2>$e1){
       return ("$s1-$e1","$s2-$e2");
     #2,overlap
     }elsif($s2<=$e1 && $e1<=$e2){
        return ("$s1-$e2")
      #3,included
      }elsif($e2<$e1 && $e2<=$e1){
         return ("$s1-$e1");
        }else{print STDERR "unknown relationship at $s1-$e1|$s2-$e2, return original\n"; return ("$s1-$e1","$s2-$e2")}
}


sub hspLink(){
   #standard function to sort&link(deOverlap) segments pairs (HSP, high scoring segment pairs)
   #sort & link query and target segments at the same time
   #input format @segs=("s1-e1|s1'-e1'","s2-e2|s2'-e2'","..-..|..-..")

}


exit 0;

