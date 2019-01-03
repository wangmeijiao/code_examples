use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#               link_mummer
#  speciesA  ----------------->  speciesB
#                 

#link file types for liftover:  mummer (can't use blast m8)
#usage: perl liftOver_mummer.pl -tab /home/mjwang/pwdexx/oryza_epiCompara_sixSpecies/stats_knobs/seq_percen/mummer_hmsk/H1_niva.fa.hmsk-H1_japo.fa.hmsk.mums -inbed H1_japo.GIs.bed -gapSize 50000   



my($tab,$inbed,$gap,$step);
GetOptions("tab=s",\$tab,"inbed=s",\$inbed,"gapSize=i",\$gap,"iterStep=i",\$step);
$gap||=50000;
$step||=100;

#1, construct link_table
my $tab_idx=&buildLinkTab($tab);
#print Dumper $tab_idx;

open BED, $inbed or die "$!";
while(<BED>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
  my @box=split/[\t ]+/,$_;
  print STDERR "\nliftOver for region $box[0]:$box[1]-$box[2]\n";
  if(!exists $tab_idx->{$box[0]}){die "can't find $box[0] in link tab"}
  my $ite=0;
  my $start_coord=$box[1];
  my $end_coord=$box[2];
  while(1){
    my $region_lifted = &liftOver("$start_coord-$end_coord",$tab_idx->{$box[0]});
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
       if($start_coord >= $end_coord || $start_coord <= 0 || $end_coord <= 0){print STDERR"iter $ite start >= end or < 0 ($start_coord >= $end_coord || $start_coord <= 0 || $end_coord <= 0), stop\n";last};
    }else{
           my ($s,$e)=split/-/,$region_lifted;
           print "$box[0]\t$start_coord\t$end_coord\ttarget\t$s\t$e\n";
           print STDERR "$ite runs(+/- $step bp per run) to find the mapped coord\n";
           last;
         }
     $ite++;     
     if($ite !=0 && $ite % 1000 ==0){print STDERR "$ite runs(+/- $step bp per run) passed: $start_coord\t$end_coord\n"};
     #if($ite >3){last}
  }##try hard to map both ends by iterately -1

}
close BED;



##sub###


sub buildLinkTab(){
#reconstruct ungapped alignment details from mummer result 
##query -> target
my $file = shift;
my %tab;
open MUM, $file or die "$!";
$/=">";
while(<MUM>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
  my @box=split/\n+/,$_; 
  my $head=shift @box;
  $head=~s/ +//g;
  foreach my $line(@box){
    my (undef,$ref_s,$que_s,$match)=split/[ \t]+/,$line;    
    my ($ref_e,$que_e)=($ref_s+$match,$que_s+$match);     
    $tab{$head}->{"$que_s-$que_e"}="$ref_s-$ref_e";
  }
  #last;
}
  $/="\n";
  return \%tab;
}

sub liftOver(){
 ##map region_query to region_target by searching in link_table
 my ($region,$tab)=@_;
 my ($start,$end)=split/-/,$region;
 my %hits_start;
 my %hits_end;
 #1, map start coordinate and end coordinate seperately
   foreach my$ctrl(keys %{$tab}){
      my ($qstart,$qend)=split/-/,$ctrl;
      my ($tstart,$tend)=split/-/,$tab->{$ctrl};
      if($start >= $qstart && $start <= $qend){ 
         my $start_map = $tstart+($start-$qstart);
         $hits_start{"$qstart-$qend"}="$start_map";
         #print "start coord hit:$start_map\n";
      }      
      if($end >= $qstart && $end <= $qend){ 
         my $end_map = $tstart+($end-$qstart);
         $hits_end{"$qstart-$qend"}="$end_map"
      }      
   }
   #print Dumper \%hits_start,\%hits_end;
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
              die "mapped coordinate start >= end" if($hits_start{$id_hits_start[0]} >= $hits_end{$id_hits_end[0]}  ); 
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
}

