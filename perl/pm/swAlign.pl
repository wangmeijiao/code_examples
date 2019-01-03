
use strict;
use warnings;
use Data::Dumper;


my $seq1 = "cggtcgttacggtcgcct";
my $seq2 = "agtttagtcgtacgcggaagttatgtttg";

my ($alignout,$alignlen,$identity,$perc_ident,$range1,$range2) = &align($seq1,$seq2);

print "$seq1\n$seq2\n\n";
print join("\n",($alignout,$alignlen,$identity,$perc_ident,$range1,$range2)),"\n";



## do local alignment using Smith-Waterman Algorithm
#http://etutorials.org/Misc/blast/Part+II+Theory/Chapter+3.+Sequence+Alignment/3.2+Local+Alignment+Smith-Waterman/
#http://etutorials.org/Misc/blast/Part+II+Theory/Chapter+3.+Sequence+Alignment/3.1+Global+Alignment+Needleman-Wunsch/

sub align
{
my ($seq1,$seq2)=@_;
$seq1=~tr/atgc/ATGC/;
$seq2=~tr/atgc/ATGC/;
my $match     =1;
my $mismatch  =-1;
my $gap       =-1;


## iniialization: 
my @matrix;
$matrix[0][0]{score}   =0;
$matrix[0][0]{pointer} ="none";
for (my $j =1 ;$j <= length($seq1) ;$j++){
    #$matrix[0][$j]  =$j;
    $matrix[0][$j]{score}  =0;
    $matrix[0][$j]{pointer}="none"; 
}
for (my $i =1 ;$i <= length($seq2) ;$i++){
    #$matrix[$i][0]  = $i;
    $matrix[$i][0]{score}  =0;
    $matrix[$i][0]{pointer}="none"; 
}


## fill score matrix
my $max_i  =0;
my $max_j  =0;
my $max_score =0;

for (my $i=1;$i <= length($seq2); $i++){
    for (my $j=1;$j <= length($seq1); $j++){
        my ($diagonal_score,$left_score,$up_score);
        
        ## calculate match score
        my $letter1 =substr($seq1,$j-1,1);
        my $letter2 =substr($seq2,$i-1,1);
        if ($letter1 eq $letter2){
           $diagonal_score =$matrix[$i-1][$j-1]{score}+$match;
        }else{
           $diagonal_score =$matrix[$i-1][$j-1]{score}+$mismatch;
        }

        ## calulate gap scores

        $up_score  =$matrix[$i-1][$j]{score}+$gap;
        $left_score=$matrix[$i][$j-1]{score}+$gap;
        if ($diagonal_score <= 0 and $up_score <=0 and $left_score <= 0){
              $matrix[$i][$j]{score}   =0;
              $matrix[$i][$j]{pointer} ="none";
              next; ## no need to choose best score in the next step, go to next iteration.
        }

        ## choose best score
        
        if ($diagonal_score >= $up_score){
             if ($diagonal_score >=$left_score){
                 $matrix[$i][$j]{score}  = $diagonal_score;
                 $matrix[$i][$j]{pointer}= "diagonal";
             }else{
                 $matrix[$i][$j]{score}  = $left_score;
                 $matrix[$i][$j]{pointer}= "left";
             }
        }else{
             if ($up_score >= $left_score){
                 $matrix[$i][$j]{score}  = $up_score;
                 $matrix[$i][$j]{pointer}= "up";
             }else{
                 $matrix[$i][$j]{score}  = $left_score;
                 $matrix[$i][$j]{pointer}= "left";
             }
        }
        
        ## set maximum score

        if ($matrix[$i][$j]{score} > $max_score){
             $max_i    = $i;
             $max_j    = $j;
             $max_score= $matrix[$i][$j]{score};
        }
    }
}

# trace-back

my $align1="";
my $align2="";
my $j =$max_j;
my $i =$max_i;
my $totalscore=0;
while (1){
      last if $matrix[$i][$j]{pointer} eq "none";
      $totalscore+=$matrix[$i][$j]{score};
      if ($matrix[$i][$j]{pointer} eq "diagonal"){
          $align1 .=substr($seq1, $j-1,1);
          $align2 .=substr($seq2, $i-1,1);
          $i--;
          $j--;
      }elsif($matrix[$i][$j]{pointer} eq "left"){
          $align1 .=substr($seq1, $j-1,1);
          $align2 .="-";
          $j--;
      }elsif($matrix[$i][$j]{pointer} eq "up"){
          $align1 .="-";
          $align2 .=substr($seq2, $i-1,1);
          $i--;
      }
}

$align1 = reverse $align1;
$align2 = reverse $align2;

# post calculation

my $alignout="$align1\n$align2";
#print "max_score\tmax_i\tmax_j\n";
#print "$max_score\t$max_i\t$max_j\n";
#print "Score: $totalscore\n";
#print "$align1\n"; 
#print "$align2\n";
my $seq2start=$i+1;
my $seq1start=$j+1;
my $seq2end  =$max_i;
my $seq1end  =$max_j;
#print "offset query:  $offset1\t$max_j\n";
#print "offset target: $offset2\t$max_i\n";
my $alignlen=length $align1;
my $identity;

for (my $i=0;$i < $alignlen; $i++){
     my $word1=substr($align1,$i,1);
	 my $word2=substr($align2,$i,1);
	 if ($word1 eq $word2) {
	      $identity++;
	 }
}
$perc_ident=sprintf("%.3f",100*$identity/$alignlen);

return ($alignout,$alignlen,$identity,"$perc_ident%","$seq1start-$seq1end","$seq2start-$seq2end");
}

