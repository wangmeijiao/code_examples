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
 #                               degap ?
 #deduce the concensous sequence finally


#usage: cat in.mfa |perl thisfile.pl <-quiet>  [<-mablock mablock.fa>  <-cons cons.fa> <-cut cutend.fa>  <-degap degap.fa>]   > in.mfa.stats
   ##quiet mode output only regions and perc

my ($quiet,$nohead,$cut,$cons,$mablock,$degap);
GetOptions("quiet!",\$quiet,"nohead!",\$nohead,"mablock:s",\$mablock,"cut:s",\$cut,"degap:s",\$degap,"cons:s",\$cons);



$/=">";
my %align;
my $num_align;
my @align_order;
while(<stdin>){
  chomp;
  next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
  my @box=split/\n+/,$_; 
  my $id=shift @box;
  push @align_order,$id;
  my $align=join"",@box;
  if(!exists $align{$id}){$align{$id}=$align; $num_align++}else{die "dup alignment $id\n"}
}
$/="\n";


my $seq_index=&getRealSeq(\%align);
#print Dumper $seq_index;

my ($flag,$len)=&checkAlignLenAll(\%align);
die "length diffs" if($flag == 0);

#my $stars=&getStars(\%align,$len);
my ($flag_gaps, $densi_score, $cons_seq)=&getSimilarDensi(\%align,$len,"score");
#print Dumper $densi;
#my $stars=&densi2star($densi,50);
my $stars=&densi2star($densi_score,0.75*$num_align);
my($coverage,$identity,$segs,$segs_linked)=&statBlocks($stars); 
my ($left_end,$right_end)=&endCut($stars,$flag_gaps);

my $info_idx=&getCoreInfo(\%align,"$left_end-$right_end");



###start to report
if(!$quiet){ ## detail mode
  foreach my $id(sort keys %align){printf("%-20s  %s\n",$id,$align{$id})}
  printf "%-20s  %s\n"," ",join("",@{$stars});
  printf "%-20s  %s\n"," ",join("",@{$densi_score});
  printf "%-20s  %s\n"," ",join("",@{$cons_seq});
  printf "%-20s  %s\n"," ",join("",@{$flag_gaps});
  print "#all equal#\n" if($flag == 1);
  print "@{$segs}\n@{$segs_linked}\n";
  print "ends after cut: $left_end-$right_end\n";
  #print "total aligned length: $len\ncoverage: $coverage\nidentity: $identity\n";

  if($cut){
  
  
  
  }

  if($cons){


  }

  if($mablock){
    open MAB, ">$mablock" or die "$!";
    foreach my $id(sort keys %align){
      print MAB ">$id\n";
      foreach my $seg(@{$segs}){
         my ($s,$e)=split/-/,$seg;
         next if($e-$s+1 <= 5);# filter for all members, because @segs is the same
         my $extract = substr($align{$id},$s,$e-$s+1);         
         $extract =~s/-//g;
         print MAB "$extract\n";
      } 
    }
    close MAB;
  }


  if($degap){


  }
 
}else{ #quiet mode, output summary tab format:

##header 

     if(!$nohead){print "orthID\talign_len\treal_len\treal_GC\trange_core\tlen_core\tcoverage_core\tidentity_core\tidentity_core_real\tGC_core_real\tperc_notgap_core\tleft_eat_up\tright_eat_up\tperc_isN_core\n"}
     print join("|",@align_order),"\t";

##real seqs

    #align_len
    print "$len\t";
    #real_len ref/que
    foreach my $i(0..$#align_order){ 
     my $id=$align_order[$i];
     if($i != $#align_order){print length($seq_index->{$id}),"|"}else{print length($seq_index->{$id}),"\t"}
    }    
    #real_GC ref/que
    foreach my $i(0..$#align_order){    
     my $id=$align_order[$i];
     if($i != $#align_order){print &gcPerct($seq_index->{$id}),"|"}else{print &gcPerct($seq_index->{$id}),"\t"}
    }

##after cut ends

    #range_core, len_core
    print "$left_end-$right_end\t";
    print $right_end-$left_end+1,"\t";
    #coverage_core%
    print sprintf"%.3f\t",($right_end-$left_end+1)/$len;
    
    print "$info_idx->{'identity_core'}\t$info_idx->{'identity_core_real'}\t";

    foreach my $i(0..$#align_order){
     my $id=$align_order[$i];
     if($i != $#align_order){print $info_idx->{$id}->{'GC_core'},"|"}else{print $info_idx->{$id}->{'GC_core'},"\t"}
    }

    foreach my $i(0..$#align_order){
     my $id=$align_order[$i];
     if($i != $#align_order){print $info_idx->{$id}->{'perc_notgap_core'},"|"}else{print $info_idx->{$id}->{'perc_notgap_core'},"\t"}
    }

    foreach my $i(0..$#align_order){
     my $id=$align_order[$i];
     if($i != $#align_order){print $info_idx->{$id}->{'left_eat_up'},"|"}else{print $info_idx->{$id}->{'left_eat_up'},"\t"}
    }
  
    foreach my $i(0..$#align_order){
     my $id=$align_order[$i];
     if($i != $#align_order){print $info_idx->{$id}->{'right_eat_up'},"|"}else{print $info_idx->{$id}->{'right_eat_up'},"\t"}
    }
   
    foreach my $i(0..$#align_order){
     my $id=$align_order[$i];
     if($i != $#align_order){print $info_idx->{$id}->{'perc_isN_core'},"|"}else{print $info_idx->{$id}->{'perc_isN_core'},"\n"}
    }


    #print Dumper $info_idx;

 
   }




###sub### 

sub mfaRead(){
  my $file = shift;
  $/=">";
  my %align;
  my $num_align;
  my @align_order;
  open MFA, $file or die "$1";
  while(<MFA>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp = split/[\t ]+/,$id;
    my $id_real = $temp[0];
    push @align_order,$id_real;
    my $align=join"",@box;
    if(!exists $align{$id_real}){$align{$id_real}=$align; $num_align++}else{die "dup alignment $id_real\n"}
  }
  close MFA;
  $/="\n";
  return (\%align,\@align_order,$num_align);
}



sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=length($index->{$id})}else{if($len != length($index->{$id}) ){$flag=0}}
  }
  return ($flag,$len);
}

sub checkAlignLenAll_array(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=scalar(@{$index->{$id}})}else{if($len != scalar(@{$index->{$id}})){$flag=0}}
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
        push @elements, uc($align{$id}->[$i]);   #modified @ 2018.9.14, to deal with softmask seqs
       #push @elements, $align{$id}->[$i];   
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
      if($ele=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in $ele"} #nucleotide included
      $cnt{$ele}++;
      if($ele eq "-"){$flag=1}
   }
   #print Dumper \%cnt;
   my @sorted=sort{ $cnt{$b} <=> $cnt{$a} }keys %cnt;
   #print "@sorted\n";
   if($sorted[0] ne "-"){$score = $cnt{$sorted[0]}; $cons = $sorted[0]}else{$score = 0; $cons = $sorted[0]}
   #if($sorted[0] ne "-"){$score = $cnt{$sorted[0]}; $cons = $sorted[0]}else{$score = $cnt{$sorted[1]}; $cons = $sorted[1]}
   $perc = sprintf "%.3f",100*$score/$num;
   if($perc <= 50){$score=0; $cons=" "}



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

sub endCut(){
   #cut hang unaligned ends from both direction. interval regions remain untouched
   my ($star,$flag)=@_; 
   die "star, flag length not all equal" if(scalar @{$star} != scalar @{$flag} );
   my $len=scalar @{$flag};

   ##eat from left, try jump 10bp to eat
   my $left_cut=1;
   my $right_cut=$len;
   my $outter_loop_check=0;
   for (my $i=0;$i<=$len-1;$i++){
     last if ($outter_loop_check == 1);
     if($flag->[$i] == 1){next}elsif($flag->[$i] == 0){
        my $cnt_flag=0;
        my $cnt_star=0;
        for (my $j=$i;$j<=$len-1;$j++){
          if($flag->[$j] == 0){
             $cnt_flag++;
             if($star->[$j] eq "*"){$cnt_star++}
           }else{
             #if($cnt_flag <= 10 ){
             if($cnt_flag <= 10 || $cnt_star/$cnt_flag < 0.5){
               $i=$j;
               last;
             }else{
                $left_cut=$i;
                $outter_loop_check=1;
                last;
              }
           }
        }##inner loop
     }else{die "unknow flag value $flag->[$i]"}     
   }#outter loop   

   ##eat from right, try jump 10bp to eat
   $outter_loop_check=0;
   for (my $i=$len-1;$i>=0;$i--){
     #print STDERR "$i\n";
     last if ($outter_loop_check == 1);
     if($flag->[$i] == 1){next}elsif($flag->[$i] == 0){
        my $cnt_flag=0;
        my $cnt_star=0;
        for (my $j=$i;$j>=0;$j--){
          if($flag->[$j] == 0){
             $cnt_flag++;
             if($star->[$j] eq "*"){$cnt_star++}
           }else{
             #if($cnt_flag <= 10 ){
             if($cnt_flag <= 10 || $cnt_star/$cnt_flag < 0.5){
               $i=$j;
               last;
             }else{
                $right_cut=$i;
                $outter_loop_check=1;
                last;
              }
           }
        }##inner loop
     }else{die "unknow flag value $flag->[$i]"}
   }#outter loop   


   return ($left_cut,$right_cut);

}


sub gcPerct(){
   my $seq = shift;
   my $len = length $seq;
   die "empty seq: $seq" if($len == 0);
   my @seqs=split//,$seq;
   die "illegal char found $seq" if($len != scalar @seqs);
   my $cnt=0;
   foreach(@seqs){if($_ eq "c" || $_ eq "C" || $_ eq "g" || $_ eq "G"){$cnt++}}
   return sprintf ("%.3f",$cnt/$len);
}

sub getCoreInfo(){
   my ($align,$range)=@_;
   my ($s,$e)=split/-/,$range;
   die "s >= e in range of core" if($s >= $e);
   my %info;
   my %core;
   foreach my $id(keys %{$align}){
      my @box=split//,$align->{$id};
      #deal with core alignment
      my @core=@box[$s-1..$e-1];
      die "core align of $id is empty" if(scalar @core == 0);
      $core{$id}=\@core;
      my $core_cnt_gap=0;
      my $core_cnt_isN=0;
      my $core_cnt_total=0;
      my $core_cnt_real=0;
      my $core_cnt_CG=0;
      foreach(@core){
        if($_=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in @core at $_"} #nucleotide included
        $core_cnt_total++;
        if($_ ne "-"){$core_cnt_real++}
        if($_ eq "-"){$core_cnt_gap++}
        if($_ eq "N" || $_ eq "n"){$core_cnt_isN++}
        if($_ eq "c" || $_ eq "C" || $_ eq "g" || $_ eq "G"){$core_cnt_CG++}
      }
      die "core_cnt_total == 0" if ($core_cnt_total == 0);
      my $perc_notgap_core=sprintf "%.3f",1-$core_cnt_gap/$core_cnt_total;
      $info{$id}->{"perc_notgap_core"}=$perc_notgap_core;
      my $perc_isN_core=sprintf"%.3f",$core_cnt_isN/$core_cnt_real;
      $info{$id}->{"perc_isN_core"}=$perc_isN_core;
      my $CG_core = sprintf "%.3f",$core_cnt_CG/$core_cnt_real;
      $info{$id}->{"GC_core"}=$CG_core;
      #deal with left eat up
      my @left_eat_up=@box[0..$s-1];
      my $cnt_left_eat_up=0;
      foreach(@left_eat_up){
        if($_=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in @left_eat_up at $_"} #nucleotide included
        if($_ ne "-"){$cnt_left_eat_up++}        
      }
      $info{$id}->{"left_eat_up"}=$cnt_left_eat_up;
      #deal with right eat up
      my @right_eat_up=@box[$e+1..$#box];
      my $cnt_right_eat_up=0;
      foreach(@right_eat_up){
        if($_=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in @right_eat_up at $_"} #nucleotide included
        if($_ ne "-"){$cnt_right_eat_up++}
      }
      $info{$id}->{"right_eat_up"}=$cnt_right_eat_up;
   }#foreach id

   #identity of core alignment
   my ($flag,$len)=&checkAlignLenAll_array(\%core);
   die "core align length diffs at ",keys %core if($flag == 0);
   die "length of core alignment is zero" if($len <= 0);
   my $cnt_core_similar=0;
   my $cnt_core_notgap=0;
   for my $i(0..$len-1){
     my $slice;
     my $flag_similar=1;
     my $flag_notgap=1;
     foreach my $id(keys %core){
       if($core{$id}->[$i] eq "-"){$flag_notgap = 0; $flag_similar = 0}
       if(!defined $slice){
         $slice = $core{$id}->[$i];
       }else{
              if($slice ne $core{$id}->[$i] ){$flag_similar = 0}
            }
     }      
     if($flag_similar == 1){$cnt_core_similar++}
     if($flag_notgap == 1){$cnt_core_notgap++}
   }
   die "cnt_core_notgap = 0" if($cnt_core_notgap <= 0);
   my $identity_core_real=$cnt_core_similar/$cnt_core_notgap;      
   my $identity_core=$cnt_core_similar/$len;      
   $info{"identity_core"}=sprintf "%.3f",$identity_core;
   $info{"identity_core_real"}=sprintf "%.3f",$identity_core_real;

   return \%info;

}


sub isNuc(){
   my $seq=shift;
   if($seq=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in $seq"} #nucleotide and pep alphabet
   my @seqs=split//,$seq;
   my $isNuc=0;
   my $total=0;
   foreach(@seqs){
     $total++;
     if($_ eq "A" ||$_ eq "a" || $_ eq "T" || $_ eq "t" || $_ eq "C" || $_ eq "c" || $_ eq "G" || $_ eq "g" || $_ eq "N" || $_ eq "n"){$isNuc++}
   }
   if($isNuc/$total >= 0.9){return 1}else{return 0} #tolerate for some consolidated characters
}



