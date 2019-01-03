use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor;

#annotate for all alignment groups and for all species 
#method: find the high quality block -> find exon/cds,k43_peak,CGI overlap region -> report for each species

#report good blocks and transform to act and igv bed12 format


my ($cds_addcoord,$mrna_addcoord);
GetOptions("cds=s",\$cds_addcoord,"mrna=s",\$mrna_addcoord);


#my $mrna_index=&readFa_list(&readList($mrna_addcoord));
#print STDERR "read mrna_addcoord done\n";
my $cds_index=&readFa_list(&readList($cds_addcoord));
print STDERR "read cds_addcoord done\n";

#print Dumper $mrna_index,$cds_index; exit;

my @order = ("japo","glab","punc","brac","lper","sorg");

$/="\n\n";
my $cnt=0;

while(<stdin>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/>/,$_;
    die "empty line at $_" if(scalar @box == 0);
    my $num_align;
    my @align_order;
    my %align_header;
    my %align;
    my %align_anno;
    my %marks;
    foreach my $str(@box){
       next if($str eq "" || $str eq "\n");
       my @temp = split/\n+/,$str;
       die "empty line at |$str|" if(scalar @temp == 0);
       my $id= shift @temp;
       my $spec;
       if($id=~/^OSJAP/){$spec = "japo"}elsif($id=~/^OGLAB/){$spec = "glab"}elsif($id=~/^OPUNC/){$spec = "punc"}elsif($id=~/^OBRAC/){$spec = "brac"}elsif($id=~/^LPERR/){$spec = "lper"}elsif($id=~/^Sobic/){$spec = "sorg"}else{die "unknow spec from id $id"}
       push @align_order,$spec;
       my $seq = join"",@temp;
       if(!exists $align{$spec}){$align{$spec}=$seq;$num_align++}else{die "dup spec $spec at $_"}

      ###1, align annotation
      my @align = split//,$seq;
      my $i_real=0;
      my $idx=0;
      for(my $i=0;$i<=$#align;$i++){
        if($align[$i] ne "-"){$i_real++}else{next}
        #my $real_bp = $i_real + $cds_index->{$spec}->{$id}->{"start"};
        if(exists $cds_index->{$spec}->{$id}->{"starts"}->{$i_real}){$idx++;$align[$i] = $idx ;push @{$marks{$spec}},$i}
        #if(exists $cds_index->{$spec}->{$id}->{"ends"}->{$i_real}){$align[$i] = "|" }
        #if($real_bp  == $cds_index->{$spec}->{$id}->{"start"}){$align[$i] = "!" }
      }
    
      my $align_anno = join"",@align;
      if(!exists $align_anno{$spec}){$align_anno{$spec} = $align_anno}else{die "dup spec $spec at $_"}
   
      if(!exists $align_header{$spec}){$align_header{$spec} = "$cds_index->{$spec}->{$id}->{chr}|$cds_index->{$spec}->{$id}->{start}|$cds_index->{$spec}->{$id}->{end}|$id|$cds_index->{$spec}->{$id}->{strand}|$cds_index->{$spec}->{$id}->{regions}" }else{die "dup spec $spec at $_"}
    }#foreach str


  ###2, qualify the alignment block
  my $seq_index=&getRealSeq(\%align);
  #print Dumper $seq_index;
  
  my ($flag,$len)=&checkAlignLenAll(\%align);
  die "length diffs" if($flag == 0);
  
  #my $stars=&getStars(\%align,$len);
  my ($flag_gaps, $densi_score, $cons_seq)=&getSimilarDensi(\%align,$len,"score");
  #print Dumper $densi;
  #my $stars=&densi2star($densi,50);
  my $stars=&densi2star($densi_score,0.75*$num_align);
  my ($coverage,$identity,$segs,$segs_linked)=&statBlocks($stars);
  my ($left_end,$right_end)=&endCut($stars,$flag_gaps);
  
  my $info_idx=&getCoreInfo(\%align,"$left_end-$right_end");
  
  ###3, start to report
  my $width||=6;
  foreach my $id(@order){
     next if(!exists $align_anno{$id});
     printf("%-${width}s  ",$id);
     my @string = split//,$align_anno{$id};
     foreach my $s(@string){
         if($s eq "!"){print color 'bold green';print $s;print color 'reset';next}
         if($s eq "|"){print color 'bold red';print $s;print color 'reset';next}
         if($s =~/\d/){print color 'bold red';print $s;print color 'reset';next}
         if($s eq "<"){print color 'bold red';print $s;print color 'reset';next}
         if($s eq ">"){print color 'bold red';print $s;print color 'reset';next}
         if($s eq "("){print color 'bold blue';print $s;print color 'reset';next}
         if($s eq ")"){print color 'bold blue';print $s;print color 'reset';next}
         if($s eq "{"){print color 'bold magenta';print $s;print color 'reset';next}
         if($s eq "}"){print color 'bold magenta';print $s;print color 'reset';next}
         print $s;
     }
     print "\n";
   }
  
  #foreach my $id(@order){printf("%-${width}s  %s\n",$id,$align_anno{$id})}  
  #foreach my $id(sort keys %align){printf("%-${width}s  %s\n",$id,$align{$id})}  
  printf "%-${width}s  %s\n"," ",join("",@{$stars});
  printf "%-${width}s  %s\n"," ",join("",@{$densi_score});
  printf "%-${width}s  %s\n"," ",join("",@{$cons_seq});
  printf "%-${width}s  %s\n"," ",join("",@{$flag_gaps});
  print "#all equal# $len\n" if($flag == 1);
  #print join"\n",@align_header,"\n";
  foreach my $id(@order){next if(!exists $align_header{$id});print "$align_header{$id}\n"}
  foreach my $id(@order){next if(!exists $marks{$id});print "$id marks in alignment: @{$marks{$id}}\n"}
  print "@{$segs}\n@{$segs_linked}\n";
  print "ends after cut: $left_end-$right_end\n";
  my $extract = &alignExtract(\%align,"$left_end-$right_end",0);
  foreach my $id(@order){next if(!exists $extract->{$id});print "$id extract(degap): $extract->{$id}\n"}
  #print "total aligned length: $len\ncoverage: $coverage\nidentity: $identity\n";
  my $ends_real = &alignTransform2real(\%align,\%align_header,"$left_end-$right_end");
  foreach my $id(@order){next if(!exists $ends_real->{$id});print "$id real start-end: $ends_real->{$id}\n"}

  print "\n\n";
  $cnt++;
  #if($cnt == 1){last}

}#while end

$/="\n";



##sub

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

sub readList(){
  my $file=shift;
  my %list;
  open LIST,$file or die"$!";
  while(<LIST>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my($spec,$file)=split/[\t ]+/,$_;
    if(!exists $list{$spec}){$list{$spec}=$file}else{die "dup species $spec\n"}
 }
  close LIST;
  return \%list;
}

sub readFa_list(){
    my $list=shift;
    my %fa;
    $/=">";
    foreach my $spec(sort keys %{$list}){
      open FA, $list->{$spec} or die "$!";
      while(<FA>){
        chomp;
        next if($_ eq "" || $_=~/^#/);
        my @box=split/\n+/,$_;
        my $headline = shift @box;
        my ($chr,$start,$end,$regions,$id,$strand,$type)  = split/\|+/,$headline;
        my @segs = split/:/,$regions; #strand eq +
        my @segs_sorted = sort{
                                my ($s1,$e1) = split/-/,$a; 
                                die "s1 > e1 in segs sort" if($s1 > $e1);
                                my ($s2,$e2) = split/-/,$b; 
                                die "s2 > e2 in segs sort" if($s2 > $e2);
                                if($e1 <= $s2 || $e2<=$s1){}else{die "overlap segs $a,$b at $headline"}
                                if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                               }@segs;
        my %seg_starts;
        my %seg_ends;

        my $len=0;
        my $len_sum=0;
        foreach my $seg(@segs_sorted){
          my ($s,$e) = split/-/,$seg;
          $len=($e-$s+1);
          $seg_starts{$len_sum+1}=1;
          $seg_ends{$len_sum+$len+1}=1;
          $len_sum+=$len;
        }
        
        $id=~s/\.\d+$//;
        if($id=~/^Sobic/){}else{$id=uc($id)}
        #my $seq=join"",@box;
        if(!exists $fa{$spec}->{$id}){
          #$fa{$spec}->{$id}->{'anno'}=$headline;
          #$fa{$spec}->{$id}->{'seq'}=$seq;
          $fa{$spec}->{$id}->{'starts'}=\%seg_starts;
          $fa{$spec}->{$id}->{'start'}=$start;
          $fa{$spec}->{$id}->{'ends'}=\%seg_ends;
          $fa{$spec}->{$id}->{'end'}=$end;
          $fa{$spec}->{$id}->{'chr'}=$chr;
          $fa{$spec}->{$id}->{'strand'}=$strand;
          $fa{$spec}->{$id}->{'regions'}=$regions;
        }else{die "dup $id in file $list->{$spec}\n"}
      }
      close FA;
    }
    $/="\n";
    return \%fa;
}

  
sub alignTransform2real(){
  #transform aligment coordinates to real sequence coordinates
  #considering CDS coordinates jumpping 
  my ($align,$align_header,$region) = @_;
  my ($start,$end) = split/-/,$region; 
  my %region_real;
  foreach my $spec(keys %{$align}){
     my @seq_align = split//,$align->{$spec};
     my $cnt=0;
     my $cnt_real=0;
     my $start_real = $start;
     my $end_real = $end;
     foreach my$base(@seq_align){
       if($base ne "-"){$cnt_real++;$cnt++}else{$cnt++}
       if($cnt == $start){$start_real = $cnt_real}
       if($cnt == $end){$end_real = $cnt_real}
     }
     
     my ($chr,$s,$e,$id,$strand,$segs) = split/\|/,$align_header->{$spec};
     my $coord_tab = &realTransform2Align($segs,$strand);
     foreach my $ctrl(keys %{$coord_tab}){
        my ($s,$e) = split/-/,$ctrl;
        my ($s_bp_real,$e_bp_real) = split/-/,$coord_tab->{$ctrl};
        if($start_real + 1 >= $s && $start_real <= $e){my $d = $start_real - $s +1;die "d < 0" if($d < 0);$start_real = $s_bp_real + $d  }
        if($end_real >= $s && $end_real <= $e){my $d = $end_real - $s +1;die "d < 0" if($d < 0);$end_real = $s_bp_real + $d  }
     }


     $region_real{$spec} = "$start_real-$end_real";
  }
  return \%region_real;
}


sub realTransform2Align(){
   #for CDS joined seqence aligment
   # if strand is +, check for overlap and sorted, if strand is -, do not sort or reverse
   my ($regions,$strand) = @_;
   my @segs = split/:/,$regions;
   my @segs_sorted = sort{
                           my ($s1,$e1) = split/-/,$a;
                           die "s1 > e1 in segs sort" if($s1 > $e1);
                           my ($s2,$e2) = split/-/,$b;
                           die "s2 > e2 in segs sort" if($s2 > $e2);
                           if($e1 <= $s2 || $e2<=$s1){}else{die "overlap segs $a,$b "}
                           if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                          }@segs;
   my %segs_transform;
   my $len=0;
   my $len_sum=0;
   foreach my $seg(@segs){
     my ($s,$e) = split/-/,$seg;
     $len=($e-$s+1);
     my $start_align = $len_sum+1;
     my $end_align = $len_sum+$len;
     $segs_transform{"$start_align-$end_align"} = "$s-$e";
     $len_sum+=$len;
   }
   return \%segs_transform;
}

sub alignExtract(){
  #aln must length checked
  my ($aln,$region,$degap) = @_;
  my ($s,$e) = split/-/,$region;
  my $len = $e-$s+1;
  die "extract aln length <= 0 at $region" if ($len <=0);
  my %extract;
  foreach my $spec(keys %{$aln}){
    my $seq = substr($aln->{$spec},$s,$len);
    if($degap == 1){$seq=~s/-//g;$extract{$spec}=$seq}elsif($degap == 0){$extract{$spec}=$seq}else{die "don't know degap or not"}
  }
  return \%extract;
}

sub refDegap(){
  #degap for given species, note coordinates of the others will be lost
  my ($align,$ref) = @_;
  foreach my $id(keys %{$align}){
    
  }

  die "ref $ref not exist in mfa file" if(!exists $align->{$ref});
    
  print STDERR "ref is $ref, start to degap\n";
  my %align_real;
  my $cnt;
  for my $i(0..$len-1){
    my $ref_nuc = substr($align->{$ref},$i,1);
    if($i % 1000000 == 0){print STDERR "#"}
    if($ref_nuc ne "-"){
      $cnt++;
      foreach my $id(@$order){
        if($id eq $ref){
          $align_real{$id}.=$ref_nuc;
        }else{
          my $nuc_extract = substr($align->{$id},$i,1);
          $align_real{$id}.=$nuc_extract;
        }
        if($cnt % 60 == 0){$align_real{$id}.="\n"}
      }
  
    }
  
  }
  
  
  print  STDERR "\nstart to output\n";
  foreach my $id(@$order){
    print ">$id\n$align_real{$id}\n";
  
  }
  

}#refdegap end

sub mablock(){
  my $align = shift;
  foreach my $id(sort keys %align){
      print ">$id\n";
      foreach my $seg(@{$segs}){
         my ($s,$e)=split/-/,$seg;
         next if($e-$s+1 <= 5);# filter for all members, because @segs is the same
         my $extract = substr($align{$id},$s,$e-$s+1);
         $extract =~s/-//g;
         print MAB "$extract\n";
      }
    }

}



=pod
sub alnStat(){
     my ($align_order,$seq_index) = @_;
     my %info_idx;
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

}#stat end

=cut


