use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

##input: mfa only; >=2 seqs ; alignment length must equal

##output:
 #calculate similarity density data along matched sequence
 #calculate #seg_num, #length and  coverage_overall% and identity_overall% , %identity_match from similarity density data
 #find the best continuous blocks (similar to mablock does )
 #output the filtered alignment:
 #                               cutting the low similarity ends
 #                               filter gaps and cutting the ends
 #deduce the concensous sequence finally


#usage: cat in.mfa |perl thisfile.pl <-quiet> <-mablock>  <-cons> <-cut>  <-degap>   > in.mfa.stats


my ($quiet,$mablock,$cut,$degap,$cons);
GetOptions("quiet!",\$quiet,"mablock!",\$mablock,"cut!",\$cut,"degap!",\$degap,"cons!",\$cons);

$/=">";
my %align;
my $num_align;
while(<stdin>){
  chomp;
  next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
  my @box=split/\n+/,$_; 
  my $id=shift @box;
  my $align=join"",@box;
  if(!exists $align{$id}){$align{$id}=$align; $num_align++}else{die "dup alignment $id\n"}
}
$/="\n";

if(!$quiet){
 foreach my $id(sort keys %align){
  printf("%-20s  %s\n",$id,$align{$id});
 }
}


my $seq_index=&getRealSeq(\%align);
#print Dumper $seq_index;

my ($flag,$len)=&checkAlignLenAll(\%align);
if($flag){
  #my $stars=&getStars(\%align,$len);
  my ($flag_gaps, $densi, $cons)=&getSimilarDensi(\%align,$len,"score");
  #print Dumper $densi;
  #my $stars=&densi2star($densi,50);
  my $stars=&densi2star($densi,0.75*$num_align);
  if(!$quiet){printf "%-20s  %s\n"," ",join("",@{$stars})}
  if(!$quiet){printf "%-20s  %s\n"," ",join("",@{$densi})}
  if(!$quiet){printf "%-20s  %s\n"," ",join("",@{$cons})}
  if(!$quiet){printf "%-20s  %s\n"," ",join("",@{$flag_gaps})}
  if(!$quiet){print "#all equal#\n"}
  my($coverage,$identity,$segs,$segs_linked)=&statBlocks($stars); 
  if(!$quiet){print "@{$segs}\n@{$segs_linked}\n"}
  if(!$quiet){print "total aligned length: $len\ncoverage: $coverage\nidentity: $identity\n"}
  if($quiet){print "$coverage\t$identity\n"}
  if($mablock && !$quiet){
    foreach my $id(sort keys %align){
      print STDERR ">$id\n";
      foreach my $seg(@{$segs}){
         my ($s,$e)=split/-/,$seg;
         next if($e-$s+1 <= 5);# filter for all members, because @segs is the same
         my $extract = substr($align{$id},$s,$e-$s+1);         
         $extract =~s/-//g;
         print STDERR "$extract\n";
      } 
    }
  }

}else{die" length not equal\n"}




###sub### 

sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=length($index->{$id})}else{if($len != length($index->{$id}) ){$flag=0}}
  }
  return ($flag,$len);
}


sub getRealSeq(){
   my $index=shift;
   my %seq=%{$index};
   foreach my $id(keys %seq){
      my $seq=$seq{$id};
      $seq=~s/-//g;
      $seq{$id}=$seq;
   }   
   return \%seq;
}

sub getSimilarDensi(){
  my ($index,$len,$type)=@_; 
  my @densi; #use score or perc;
  my @cons;  
  my @flag_gaps; #len == aligned_len ;
  my %align;
  my @ids;
  foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align{$id}=\@box;
      push @ids,$id;
  }
  for my $i (0.. $len-1){
     my @elements;
     foreach my $id(@ids){
       push @elements, $align{$id}->[$i];   
     }
     my($flag, $score, $perc, $cons)=&checkIdentity_densi(\@elements);
     if($type eq "score"){push @densi, $score}elsif($type eq "perc"){push @densi, $perc}else{die "don't understand densi type $type"}
     push @cons, $cons;
     push @flag_gaps,$flag;
  }
  return(\@flag_gaps, \@densi, \@cons) ;
}

sub densi2star(){
  my ($densi,$cutoff)=@_; 
  my @stars;
  my $len=scalar @{$densi};
  for(0..$len-1){
    die "undefed value in @{$densi}" if(!defined $densi->[$_]);
    if($densi->[$_] >= $cutoff){push @stars, "*"}else{push @stars, " "}
  } 
  return \@stars;
}


sub getStars(){
    my $index=shift;
    my $len=shift;
    my @stars;
    my %align_array;
    my @ids;
    foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align_array{$id}=\@box;
      push @ids,$id;
    } 
    #print Dumper \%align_array;
    for my $i (0..$len-1){
       my @elements;
       foreach my $id(@ids){
         push @elements, $align_array{$id}->[$i];
       }
       my $flag=&checkIdentity(\@elements);
       if($flag == 1){push @stars,"*"}elsif($flag == 0){push @stars," "}elsif($flag == -1){push @stars,"-"}else{die "undefined flag $flag\n"}
    }
    return \@stars;
}

sub statBlocks(){
    #*: the same
    #-: gap
    #(space): diff
    my $star=shift;
    my $len=scalar @{$star};
    my ($cnt_star,$cnt_gap,$cnt_diff)=(0,0,0);
    my ($coverage,$coverage_total,$identity)=(0,0,0);
    my $gap=4; #for link nearby similar blocks
    #find continual segments
    my @segments;
    my ($start,$end,$cnt)=(0,0,0);
    for my $i(0..$len-1){
      if($star->[$i] eq "*"){$cnt_star++}elsif($star->[$i] eq "-"){$cnt_gap++}elsif($star->[$i] eq " "){$cnt_diff++}     
      if($star->[$i] eq "*"){
         $cnt++;
         if($i==($len-1)){ #last element
           $start=$i-$cnt+1;
           $end=$i;
           push @segments,$start."-".$end;
           $cnt=0;
         }
      }else{
            $start=$i-$cnt;
            $end=$i-1;
            if($start<=$end){push @segments,$start."-".$end}
            $cnt=0;
           }
    }#for end
    my @segs=@segments; #@segments changed in next steps
    #print "@segs\n";

    #link small segments to long ones with gap<=4
    my @segments_linked;
    for (my $i=1;$i<=$#segments;$i++){
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        if(($tmp2[0]-$tmp1[1]-1)<= $gap){
        $start=$tmp1[0];
        $end=$tmp2[1];
        $segments[$i]=$start."-".$end;
        }else{push @segments_linked, $segments[$i-1];}
        if($i==$#segments){push @segments_linked,$segments[$i];}
    }
    if($#segments==0){push @segments_linked,$segments[0];}
    #print "@segments_linked\n";

    #start to calculate percentage etc
    $coverage=sprintf "%.3f",($cnt_star+$cnt_diff)/$len; 
    $identity=sprintf "%.3f",$cnt_star/($cnt_star+$cnt_diff);
    return ($coverage,$identity,\@segs,\@segments_linked);
}

sub checkIdentity(){
    #return -1(has gap), 1(all identical), 0 (not all identical)
    my @elements=@{shift @_};
    my $flag=1;
    my $foo;
    foreach (@elements){
      if($_ eq "-"){$flag=-1;return $flag}
      if(!defined $foo){$foo = $_}else{if($foo ne $_){$flag = 0}}
    }
    return $flag;
}


sub checkIdentity_densi(){
   #count -> calculate perc%, similarity score (0-num_seq) and output a cons base 
   my @elements=@{shift @_};
   my $num=scalar @elements;
   my ($score,$flag,$cons,$perc);
   $flag=0;
   my %cnt; #for gap, nucleotide base and aa and others
   foreach my $ele (@elements){
      if($ele=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in $ele"}
      $cnt{$ele}++;
      if($ele eq "-"){$flag=1}
   }
   #print Dumper \%cnt;
   my @sorted=sort{ $cnt{$b} <=> $cnt{$a} }keys %cnt;
   #print "@sorted\n";
   if($sorted[0] ne "-"){$score = $cnt{$sorted[0]}; $cons = $sorted[0]}else{$score = 0; $cons = $sorted[0]}
   #if($sorted[0] ne "-"){$score = $cnt{$sorted[0]}; $cons = $sorted[0]}else{$score = $cnt{$sorted[1]}; $cons = $sorted[1]}
   $perc = sprintf "%.3f",100*$score/$num;
   return ($flag, $score, $perc, $cons);
}

=pod
G - Glycine (Gly)
P - Proline (Pro)
A - Alanine (Ala)
V - Valine (Val)
L - Leucine (Leu)
I - Isoleucine (Ile)
M - Methionine (Met)
C - Cysteine (Cys)
F - Phenylalanine (Phe)
Y - Tyrosine (Tyr)
W - Tryptophan (Trp)
H - Histidine (His)
K - Lysine (Lys)
R - Arginine (Arg)
Q - Glutamine (Gln)
N - Asparagine (Asn)
E - Glutamic Acid (Glu)
D - Aspartic Acid (Asp)
S - Serine (Ser)
T - Threonine (Thr)
=cut

