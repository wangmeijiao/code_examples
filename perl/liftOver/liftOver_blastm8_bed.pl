use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#            processed_link_blastm8
#  speciesA  ----------------->  speciesB
#                 

#try to learn from ncbi remap "http://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/alignments"


#link file types for liftover:  blast m8
#usage: perl liftOver_blastm8_bed.pl -tab H1_japo-VS-H1_zs97.blast -inbed japo.mRNA_CDS.embl.gff.bed -gapSize 50000 -iterStep 100 > japo.mRNA_CDS.embl.liftOver.H1_zs97.bed 2> liftOver.bed.err 

#this script is based on liftOver_blastm8_gff.pl
#version 1.0 finished at 2016, 9.18 
#iteLift need test and improvement


## correct minus strand mapping error, Oct. 28, 2016


my($tab,$inbed,$gap,$step);
GetOptions("tab=s",\$tab,"inbed=s",\$inbed,"gapSize=i",\$gap,"iterStep=i",\$step);
$gap||=50000;
$step||=10;

my ($tab_idx,$pair)=&buildLinkTab_blastm8($tab);
my ($chr_q,$chr_t) = split/\|/,$pair;
#print Dumper $tab_idx,$pair; exit; 

open BED, $inbed or die "$!";
my $cnt;
my $cnt_ok;
while(<BED>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
  my ($chr,$start,$end,$id,$score,$strand)=split/[\t ]+/,$_;
  die"empty value at $_\n" if(!defined $chr || !defined $start || !defined $end || !defined $id || !defined $score || !defined $strand );
  next if($chr eq "chrUn" || $chr eq "chrSy" || $chr eq "Pt" || $chr eq "Mt"); #get rid of other chrs
  print STDERR "###begin to lift region $chr:$start-$end..###\n";
  if(!exists $tab_idx->{$chr}){die "can't find $chr in link tab"}
  my $region_lifted = &liftOver("$start-$end",$tab_idx->{$chr},$gap); 
  if($region_lifted !~/^#/ ){
    my ($s,$e) = split/-/,$region_lifted;
    if($s > $e){($s,$e)=($e,$s);$strand=&RC_strand($strand)}
    #($s > $e)?(($s,$e)=($e,$s)):(($s,$e)=($s,$e));
    print "$chr_t\t$s\t$e\t$id\t$score\t$strand\n";
    print STDERR "$chr:$start-$end liftOver sucessfully to $region_lifted\n";
    $cnt_ok++
  }else{print "#$chr\t$start\t$end\t$id\t$score\t$strand\n";print STDERR "$chr:$start-$end liftOver fail: $region_lifted\n"}
  ##iteLift: try to cut step bp inward and try liftOver iterately
  #my $region_lifted = &iteLift("$chr:$start-$end",$tab_idx,$step,$gap);
  #if($region_lifted ne "nohit"){print "$chr:$start-$end -> $region_lifted\n";$cnt_ok++}else{print STDERR "liftOver fail\n"}
  print STDERR "###Finished to lift region $chr:$start-$end###\n\n";
  $cnt++;

}
close BED;
print STDERR "all done, total $cnt regions and $cnt_ok success.";





##sub###


sub buildLinkTab_blastm8(){
#reconstruct blast gapped alignment details  (customed blast+ outfmt 6)
## qID-> query -> subject,  one2one blast only
#with filter alignLen percIdent, and score and evalue? do it in another script
my $file = shift;
my %tab;
my $pair;
open M8, $file or die "$!";
while(<M8>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
  my ($qID,$sID,$percIdent,$alignLen,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$score,$qseq,$sseq) = split/[\t ]+/,$_;
  if(!defined $pair){$pair = "$qID|$sID"}else{if($pair ne "$qID|$sID"){die "blast one2one only:$qID|$sID"}}
  die"qstart >= qend at $_" if($qstart >= $qend);  #blast qstart must < qend
  die"sstart == send at $_" if($sstart == $send);   #blast sstart must not equal send

  die "qseq != alignlen || sseq != alignlen at $_" if(length $qseq != $alignLen || length $sseq != $alignLen);

  my $flip=0;
  if ($sstart > $send){$flip=1}
  #if ($sstart > $send){($sstart,$send)=($send,$sstart);$flip=1}

  my ($qgap,$qgap_open,$qgap_n) = &gapScan($qseq,$qstart);  #blast sstart > send, and sseq already reversed, gapscan need not flip 
  #my ($qgap,$qgap_open,$qgap_n) = &gapScan($qseq,$qstart,0);
  $qgap = join"|",@$qgap;

  my ($sgap,$sgap_open,$sgap_n) = &gapScan($sseq,$sstart);
  #my ($sgap,$sgap_open,$sgap_n) = &gapScan($sseq,$sstart,$flip);
  $sgap = join"|",@$sgap;

  if($qgap eq ""){$qgap = "NA"}
  if($sgap eq ""){$sgap = "NA"}



  $tab{$qID}->{"$qstart-$qend"}="$sID:$sstart-$send:$flip $percIdent:$alignLen:$mismatch:$gapopen:$evalue:$score $qgap_open,$sgap_open:$qgap_n,$sgap_n $qgap $sgap"  
}
close M8;
return (\%tab,$pair);
}


sub gapScan(){
  my ($seq,$start)=@_;
  my $i=0;
  my $cnt=0;
  my @index;
  while($seq=~/-+/g){
    $i++;
    my $up;
    $up=$start+(length $`)-1-$cnt;
    #if($flip == 0){$up=$start+(length $`)-1-$cnt}elsif($flip == 1){$up=$start-((length $`)-1-$cnt)}else{die"unknow flip:$flip"}
    my $hit=length $&;
    $cnt+=$hit;
    #my $down=$up+$hit-1;
    #if($up == $down){$down+=1}
    push @index, "$up,$hit";
  }
  return \@index,$i,$cnt;
}

sub iteLift(){
  my ($region,$tab_idx,$step,$gap)= @_;
  my @box = split/:|-/,$region; #chr:start-end
  die "region err: start >= end:$region" if($box[1] >= $box[2]);
  print STDERR "liftOver for region $box[0]:$box[1]-$box[2]\n";
  if(!exists $tab_idx->{$box[0]}){die "can't find $box[0] in link tab"}
  my $ite=0; # $ite=0
  my $chr = $box[0];
  my $start_coord=$box[1];
  my $end_coord=$box[2];

  my $result;
  while(1){   #add #nohit, flip result start end
    my $region_lifted = &liftOver("$start_coord-$end_coord",$tab_idx->{$chr},$gap);
    if($region_lifted=~/^##/){
       #print join("\t",@box),"\ttarget\t$region_lifted\n";
       #cut region end until result satisified me
       if($region_lifted=~/^##start/){
         if($start_coord - $box[1] == 0 && $end_coord - $box[2] != 0){
                $end_coord-=$step;
               }elsif($start_coord - $box[1] != 0 && $end_coord - $box[2] == 0){
                  $start_coord+=$step;
                }elsif($start_coord - $box[1] == 0 && $end_coord - $box[2] == 0){
                   $start_coord+=$step;
                  }else{die "err: $region_lifted,  $start_coord - $end_coord | $box[1] - $box[2]"}
         #$start_coord+=$step;
       }elsif($region_lifted=~/^##end/){
           if($start_coord - $box[1] == 0 && $end_coord - $box[2] != 0){
                $end_coord-=$step;
               }elsif($start_coord - $box[1] != 0 && $end_coord - $box[2] == 0){
                  $start_coord+=$step;
                }elsif($start_coord - $box[1] == 0 && $end_coord - $box[2] == 0){
                    $end_coord-=$step;
                  }else{die "err: $region_lifted,  $start_coord - $end_coord | $box[1] - $box[2]"}
           #$end_coord-=$step;
          }elsif($region_lifted=~/^##gap too large/){
              if($start_coord - $box[1] == 0 && $end_coord - $box[2] != 0){ 
                $end_coord-=$step;
               }elsif($start_coord - $box[1] != 0 && $end_coord - $box[2] == 0){
                  $start_coord+=$step;
                }elsif($start_coord - $box[1] == 0 && $end_coord - $box[2] == 0){
                     die"$region_lifted, $start_coord - $box[1] == 0 && $end_coord - $box[2] == 0";
                   }else{die "err: $region_lifted,  $start_coord - $end_coord | $box[1] - $box[2]"}
            }elsif($region_lifted=~/^##one end has multiple/){
                 if($start_coord - $box[1] == 0 && $end_coord - $box[2] != 0){ 
                    $end_coord-=$step
                   }elsif($start_coord - $box[1] != 0 && $end_coord - $box[2] == 0){
                      $start_coord+=$step
                    }elsif($start_coord - $box[1] == 0 && $end_coord - $box[2] == 0){
                        die"$region_lifted, $start_coord - $box[1] == 0 && $end_coord - $box[2] == 0";
                      }else{die "err: $region_lifted,  $start_coord - $end_coord | $box[1] - $box[2]"}
               }else{die "unknow circumstence: $region_lifted"}
       if($start_coord >= $end_coord || $start_coord <= 0 || $end_coord <= 0){print STDERR"iter $ite start >= end or < 0 ($start_coord >= $end_coord || $start_coord <= 0 || $end_coord <= 0), stop\n";$result="nohit";last};
    }else{
           #my ($s,$e)=split/-/,$region_lifted;
           #print "$target\t$s\t$e\n";
           $result = "$region_lifted";
           #print "$box[0]\t$start_coord\t$end_coord\ttarget\t$s\t$e\n";
           print STDERR "iter $ite runs(+/- $step bp per run, 0 runs means no need to try iterate) and found the mapped coord\n";
           last;
         }
     $ite++;     
     if($ite !=0 ){print STDERR "searching  $ite runs(+/- $step bp per run) done: $start_coord\t$end_coord\n\n"};
     #if($ite !=0 && $ite % 10 ==0){print STDERR "$ite runs(+/- $step bp per run) passed: $start_coord\t$end_coord\n"};
     #if($ite >3){last}
  }##while loop: try hard to map both ends by iterately -step
  return $result;
}#sub end




sub liftOver(){
 ##map region_query to region_target by searching in link_table
 ## qlen_real+qgap = slen_real+sgap
 ## use range qlen_real+qgap to search sgap (note: this is a approximate method); check and valid the result
 my ($region,$tab,$gap)=@_;
 my ($start,$end)=split/-/,$region;
 my %hits_start;
 my %hits_end;
 #1, map start coordinate and end coordinate seperately
   foreach my$ctrl(keys %{$tab}){
      my ($qstart,$qend)=split/-/,$ctrl;
      my ($target_coord,$blast_stat,$gap_stat,$qgap,$sgap)=split/ / ,$tab->{$ctrl};
      my ($chr_t,$tstart,$tend,$flip) = split/:|-/,$target_coord;
 
      #region start gap search
      if($start >= $qstart && $start <= $qend){  #qstart always < qend
         #get query gaps
         print STDERR "liftOver:start($start):query\n";
         my $qgap_count = &gapQuery($start,$qgap,0);
         my $alignLen_map = abs($start-$qstart)+1+$qgap_count; print STDERR "\nblast link is :$qstart,$qend => $chr_t,$tstart,$tend; alignLen_map is $alignLen_map\n\n";
         #get target gaps: use start_map_esti to get gap and calculate start_map_real
         my $start_map_esti;
         my $sgap_count_esti;
         my $start_map_real;
         if($flip == 0){ 
           $start_map_esti = $tstart+$alignLen_map;
            print STDERR "liftOver:start($start):target:flip$flip\n";
           $sgap_count_esti = &gapQuery($start_map_esti,$sgap,$flip);
           $start_map_real = $tstart+ ($alignLen_map - $sgap_count_esti)-1;
         }elsif($flip == 1){
              $start_map_esti = $tstart-$alignLen_map; # $tstart > $tend in blast link pair
               print STDERR "liftOver:start($start):taget:flip$flip\n";
              $sgap_count_esti = &gapQuery($start_map_esti,$sgap,$flip);
              $start_map_real = $tstart - ($alignLen_map - $sgap_count_esti)-1; # $tstart > $tend in blast link pair
         }else{die "unknow flip $flip"   }
         #check sgap again to verify the estimate is right
         print STDERR "liftOver:start($start):target_real:flip$flip\n";
         my $sgap_count_real = &gapQuery($start_map_real,$sgap,$flip);
         if($sgap_count_esti == $sgap_count_real){
           $hits_start{"$qstart-$qend"}=$start_map_real;
           print STDERR "good luck, target estimate gap equal real: estimate is $sgap_count_esti and real is $sgap_count_real\n"
         }else{
             $hits_start{"$qstart-$qend"}=$start_map_real;
             print STDERR "estimate inappropriate: estimate gap is $sgap_count_esti while real is $sgap_count_real\n"
         }
      }

      #region end gap search
      if($end >= $qstart && $end <= $qend){  #qstart always < qend
         #get query gaps
         print STDERR "\n\nliftOver:end($end):query\n";
         my $qgap_count = &gapQuery($end,$qgap,0);
         my $alignLen_map = abs($end-$qstart)+1+$qgap_count;print STDERR "\nblast link is :$qstart,$qend => $chr_t,$tstart,$tend; alignLen_map is $alignLen_map\n\n"; 
         #get target gaps: use end_map_esti to get gap and calculate end_map_real
         my $end_map_esti;
         my $sgap_count_esti;
         my $end_map_real;
         if($flip == 0){ 
           $end_map_esti = $tstart+$alignLen_map;
            print STDERR "liftOver:end($end):target:flip$flip\n";
           $sgap_count_esti = &gapQuery($end_map_esti,$sgap,$flip);
           $end_map_real = $tstart+ ($alignLen_map - $sgap_count_esti)-1;
         }elsif($flip == 1){
              $end_map_esti = $tstart-$alignLen_map; # $tstart > $tend in blast link pair
               print STDERR "liftOver:end($end):taget:flip$flip\n";
              $sgap_count_esti = &gapQuery($end_map_esti,$sgap,$flip);
              $end_map_real = $tstart - ($alignLen_map - $sgap_count_esti)-1; # $tstart > $tend in blast link pair
         }else{die "unknow flip $flip"   }
         #check sgap again to verify the estimate is right
         print STDERR "liftOver:end($end):target_real:flip$flip\n";
         my $sgap_count_real = &gapQuery($end_map_real,$sgap,$flip);
         if($sgap_count_esti == $sgap_count_real){
           $hits_end{"$qstart-$qend"}=$end_map_real;
           print STDERR "good luck, target estimate gap equal real: estimate is $sgap_count_esti and real is $sgap_count_real\n"
         }else{
             $hits_end{"$qstart-$qend"}=$end_map_real;
             print STDERR "estimate inappropriate: estimate gap is $sgap_count_esti while real is $sgap_count_real\n"
         }
      }

   }#foreach end
   #print Dumper \%hits_start,\%hits_end;
#=pod
 #2, link mapped start and end
  my @id_hits_start=keys %hits_start;
  my $num_hits_start=scalar @id_hits_start;
  my @id_hits_end=keys %hits_end;
  my $num_hits_end=scalar @id_hits_end;

  my $region_lifted;
  if($num_hits_start == 0 || $num_hits_end == 0){ #1 one end mapping failed
    if($num_hits_start == 0){$region_lifted="##start side has no mapping result"}elsif($num_hits_end == 0){$region_lifted="##end side has no mapping result"}
  }elsif($num_hits_start > 1 || $num_hits_end > 1){ #2 one end multimapped
       $region_lifted="##one end has multiple mapping results"
     }elsif($num_hits_start == 1 && $num_hits_end == 1){ #3 both ends uniqal mapped
          if($id_hits_start[0] eq $id_hits_end[0]){   # start and end mapped onto the same segment
              #die "mapped coordinate start >= end" if($hits_start{$id_hits_start[0]} >= $hits_end{$id_hits_end[0]}  ); 
              $region_lifted="$hits_start{$id_hits_start[0]}-$hits_end{$id_hits_end[0]}";
            }else{    #start and end mapped onto two different segments
                   my $space_q=$end-$start+1;
                   my $space_t=abs($hits_start{$id_hits_start[0]} - $hits_end{$id_hits_end[0]}+1);
                   if(abs($space_q-$space_t) <= $gap){
                     $region_lifted="$hits_start{$id_hits_start[0]}-$hits_end{$id_hits_end[0]}";
                   }else{$region_lifted="##gap too large: $space_q, $space_t"}  
                 }
        }else{die "unknow circument"}
  
  return $region_lifted;
#=cut

}

sub gapQuery(){
   my ($pos,$gap_db,$flip) = @_;
   my @gaps = split/\|/, $gap_db;
   my $gap_open= scalar @gaps;
   my $gap_count=0; #default gap = 0
   if( $gap_open == 0){die "empty gap record:$gap_db"}elsif( $gap_open == 1 && $gap_db eq "NA"){print STDERR "record is NA(0),pos is $pos\n";;return 0}else{
     foreach my $rec(@gaps){
       my ($coord,$n) = split/,/,$rec; 
       print STDERR "record is $rec,pos is $pos\n"; 
       if($flip == 0){ 
          if($coord <= $pos){$gap_count+=$n}        
       }elsif($flip == 1){
           if($coord <= $pos){$gap_count+=$n} #the same?
         }else{die "unknow flip: $flip"}   

     }
   }
   print STDERR "hit gap $gap_count\n";
   return $gap_count;
}




sub RC_strand(){
  my $str = shift;
  if($str eq "+"){return "-"}elsif($str eq "-"){return "+"}else{die "err strand when strand_RC: $str"}



}





