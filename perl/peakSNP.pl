
use strict;
use warnings;
use Data::Dumper;

#check segments pair number match
#align identity
#detect SNP (use method in mfa_call_snp_indel_multi.pl)
#CpG type SNP


my ($mfa,$map) = &readMFA_multi($ARGV[0]);

my $orthpeak = &readOrthPeak_simple($ARGV[1]);


#print Dumper $mfa,$map,$orthpeak;exit;

my $cnt;
foreach my $id(sort {$a<=>$b} keys %$orthpeak){
  $cnt++;
 
  my $flag = $orthpeak->{$id}->{'flag'};

  my $id_japo = $orthpeak->{$id}->{'japo'}->{'geneId'}; 
  my $id_brac = $orthpeak->{$id}->{'brac'}->{'geneId'}; 

  my $segs_japo = $orthpeak->{$id}->{'japo'}->{'core_segs'}; 
  my $segs_brac = $orthpeak->{$id}->{'brac'}->{'core_segs'}; 

  my $strand_japo = $orthpeak->{$id}->{'japo'}->{'geneStrand'}; 
  my $strand_brac = $orthpeak->{$id}->{'brac'}->{'geneStrand'}; 

  my $orthregion_japo = $orthpeak->{$id}->{'japo'}->{'orth'}; 
  my $orthregion_brac = $orthpeak->{$id}->{'brac'}->{'orth'}; 

  #check overall aln length
  my $len_aln_japo = length ($mfa->{$id_japo}->{'aln'});
  my $len_aln_brac = length ($mfa->{$id_brac}->{'aln'});
  die "aln length diffs at orthpeak tab" if($len_aln_japo != $len_aln_brac);
  print STDERR "\n$id  with $id_japo-$id_brac aln_len:$len_aln_japo|$len_aln_brac $orthregion_japo|$orthregion_brac\n$mfa->{$id_japo}->{'aln'}\n$mfa->{$id_brac}->{'aln'}\n";

  #filter flag
  if($flag eq "1|1|1" || $flag eq "1|1|NF" || $flag eq "NF|1|1" || $flag eq "1|1|0" || $flag eq "1|0|1" || $flag eq "1|0|NF" || $flag eq "NF|0|1" ){}else{print STDERR "flag($flag) filtered\n";next}  

  #filter only gap alns
  my $real_aln_japo = $mfa->{$id_japo}->{'aln'}; 
  $real_aln_japo =~s/-+//g;
  my $real_aln_brac = $mfa->{$id_brac}->{'aln'}; 
  $real_aln_brac =~s/-+//g;
  if(length $real_aln_japo  == 0 ){print STDERR "\nGap_only: aln_japo only gap at $id_japo-$id_brac\n";next}
  if(length $real_aln_brac  == 0 ){print STDERR "\nGap_only: aln_brac only gap at $id_japo-$id_brac\n";next}

  #filter aln with too low identity;
  my $cnt_ident_all =0;
  my @aln1_all = split//,$mfa->{$id_japo}->{'aln'};  
  my @aln2_all = split//,$mfa->{$id_brac}->{'aln'};  
  for my $i(0..$#aln1_all){if($aln1_all[$i] eq $aln2_all[$i]){$cnt_ident_all++}}
  if($cnt_ident_all/$len_aln_japo < 0.2){print STDERR "Filter for aln_all too low identity\n";next}

  #check strand
  my $flag_rev_japo = 0;
  my $flag_rev_brac = 0;
  if($strand_japo eq $strand_brac){}else{
    if($strand_japo eq "+" && $strand_brac eq "-"){$flag_rev_japo = 0; $flag_rev_brac = 1}
    if($strand_japo eq "-" && $strand_brac eq "+"){$flag_rev_japo = 1; $flag_rev_brac = 0}
  }

  if(scalar @$segs_japo != 0 && scalar @$segs_brac != 0){
    if(scalar @$segs_japo !=  scalar @$segs_brac){die "segments number not equal at $id_japo"} 
    if(!exists $mfa->{$id_japo}){die "not found japo id $id_japo"}
    if(!exists $mfa->{$id_brac}){die "not found brac id $id_brac"}

    my $orth_1 = "$mfa->{$id_japo}->{'orth'}->[0]:$mfa->{$id_japo}->{'orth'}->[1]-$mfa->{$id_japo}->{'orth'}->[2]:$mfa->{$id_japo}->{'orth'}->[3]";
    my $orth_2 = "$mfa->{$id_brac}->{'orth'}->[0]:$mfa->{$id_brac}->{'orth'}->[1]-$mfa->{$id_brac}->{'orth'}->[2]:$mfa->{$id_brac}->{'orth'}->[3]";
    if($orth_1 ne $orthregion_japo){die "mfa region diff:$orth_1 ne $orthregion_japo"}#else{print STDERR "$orth_1 eq $orthregion_japo\n"}
    if($orth_2 ne $orthregion_brac){die "mfa region diff:$orth_2 ne $orthregion_brac"}#else{print STDERR "$orth_2 eq $orthregion_brac\n"}

    print ">$id_japo-$id_brac\n";
    for my $i(0..$#$segs_japo){
       my $aln_japo = &alnExtract($mfa->{$id_japo}->{'aln'},$segs_japo->[$i],$orthregion_japo);
       my $aln_brac = &alnExtract($mfa->{$id_brac}->{'aln'},$segs_brac->[$i],$orthregion_brac);
       #if(length $aln_japo != length $aln_brac){die "aln len diffs at $id_japo, $id_brac $segs_japo->[$i] $segs_brac->[$i] of \n$aln_japo\n$aln_brac\n"}
       if(length $aln_japo != length $aln_brac){
         print STDERR "Warning!!! aln len diffs at $id_japo, $id_brac $segs_japo->[$i] $segs_brac->[$i]: \n$aln_japo\n$aln_brac\n";
         ($aln_japo,$aln_brac,$segs_japo->[$i],$segs_brac->[$i]) = &endCut($aln_japo,$aln_brac,$segs_japo->[$i],$segs_brac->[$i]);print STDERR "endCut as: \n$segs_japo->[$i]: $aln_japo\n$segs_brac->[$i]: $aln_brac\n"
         #($aln_japo,$aln_brac) = &endRepair($aln_japo,$aln_brac);print STDERR "\nendRepair as: \n$aln_japo\n$aln_brac\n"
       }
       if (length $aln_japo == 0|| length $aln_brac == 0){print "$segs_japo->[$i]|$segs_brac->[$i]  NA NA NA NA\n\n\n";next}
       #my %aln;
       #$aln{'japo'} = $aln_japo;
       #$aln{'brac'} = $aln_brac;
       #my @order = ("japo","brac");
       #my ($flag,$len)=&checkAlignLenAll(\%aln);
       #if($flag != 1){die "aln len diffs at $id_japo, $id_brac $segs_japo->[$i] $segs_brac->[$i] of \n$aln_japo\n$aln_brac\n"}
       #my ($densi_score, $cons_seq, $snp, $flagN, $flagGap)=&getSimilarDensi(\%aln,\@order,$len,"score");
       #my $stars=&densi2star($densi_score,0.75*2);
       my ($ident,$snp,$CpGsnp,$perc_CpGsnp) = &alnStat($aln_japo,$aln_brac);
       #real_seq_length;
       my $real_japo = $aln_japo;
       $real_japo=~s/-+//g;
       my $len_real_japo = length $real_japo;
       my $real_brac = $aln_brac;
       $real_brac=~s/-+//g;
       my $len_real_brac = length $real_brac;
       #nCpG
       my $ncpg_japo = &cpgCount($real_japo);
       my $ncpg_brac = &cpgCount($real_brac);

       #k43 data by bp  
       my $extract_region_japo = "$mfa->{$id_japo}->{'orth'}->[0]".":$segs_japo->[$i]";
       my $k43_japo = &extractDensi("/home/mjwang/pwdexx/epi-oryza-data-sixSpecies/japo-brac_rept1/japo/k43-920-2/k43_tigr6_bowtiev0m1/hits_k43_tigr6_bowtiev0m1.density.gz",$extract_region_japo,$flag_rev_japo);
       my $k43_japo_str = join",",@$k43_japo;
       my $extract_region_brac = "$mfa->{$id_brac}->{'orth'}->[0]".":$segs_brac->[$i]";
       my $k43_brac = &extractDensi("/home/mjwang/pwdexx/RTline-seedling/ffseedling/k43-927-8/k43_bracRS2_bowtiev0m1/bed2density/hits_k43_bracRS2_bowtiev0m1.density.gz",$extract_region_brac,$flag_rev_brac);
       my $k43_brac_str = join",",@$k43_brac;



       #summary line
       print "$segs_japo->[$i]|$segs_brac->[$i]  $ident $snp  $CpGsnp $perc_CpGsnp\n$aln_japo $len_real_japo $ncpg_japo $k43_japo_str\n$aln_brac $len_real_brac $ncpg_brac $k43_brac_str\n";

    }
    print "\n";

  }else{
    print STDERR "segs empty:japo,@$segs_japo;or brac:@$segs_brac\n"
  }


}

print STDERR "all done\n";


###sub

sub readMFA_multi(){
  #read in file like out.pecan.3
  my $file = shift;
  open MFA, $file or die "$!";
  $/ = "\n\n";
  my %align;  
  my %map;
  while(<MFA>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
    my @temp = split/>/,$_;
    shift @temp;
    die "not two seqs at $_" if(scalar @temp != 2);
    my @pair;
    foreach my $ctrl(@temp){
      my @box=split/\n+/,$ctrl;
      my $head=shift @box;
      my @tmp = split/[\t ]+/,$head;
      my ($chr,$range,$id,$strand) = split/\|/,$tmp[2];
      my $align=join"",@box;
      die "dup id at $_"if(exists $align{$id});
      push @pair,$id;
      $align{$id}->{'aln'} = $align;
      $align{$id}->{'alnLen'} = $tmp[1];
      $align{$id}->{'gene'} = $tmp[2];
      $align{$id}->{'peak'} = $tmp[3];
      #link multi-orth-regions
      if($tmp[4]=~/\|/){
        print STDERR "link multi-orth-regions when readin mfa at $tmp[2] with orth $tmp[4]\n";
        my @orth = split/\|/,$tmp[4];
        my $result_link = &segLink_all(\@orth);print STDERR "@$result_link\n";
        if(scalar @$result_link == 1){
           my ($chr,$range,$strand) = split/:/,$result_link->[0];
           my ($start,$end) = split/-/,$range;
           my @orth = ("$chr","$start","$end","$strand");
           $align{$id}->{'orth'} = \@orth;
        }else{die "link failed at $tmp[2] with orth $tmp[4]"}
      }else{
             my ($chr,$range,$strand) = split/:/,$tmp[4];
             my ($start,$end) = split/-/,$range;
             my @orth = ("$chr","$start","$end","$strand");
             $align{$id}->{'orth'} = \@orth;
           }
    }
    die "not two genes at $_" if(scalar @pair !=2);
    $map{$pair[0]} = $pair[1];
    $map{$pair[1]} = $pair[0];

  }#while end 
  return \%align,\%map;

}#sub end

sub segLink_all(){
  #link and delete  original overlap segs
  my $idx = shift;
  my $gap ||= 5;
  my ($chr,$strand);
  die "segs empty when link orth" if(scalar @$idx == 0);
  die "segs == 1 when link orth" if(scalar @$idx == 1);
  my @segs_sort = sort { #sort and check
                         my ($chr1,$range1,$strand1) = split/:/,$a;
                         my ($s1,$e1) = split/-/,$range1;
                         my ($chr2,$range2,$strand2) = split/:/,$b;
                         my ($s2,$e2) = split/-/,$range2;
                         die "chr diffs when link orth $a and $b" if($chr1 ne $chr2);
                         die "strand diffs when link orth $a and $b" if($strand1 ne $strand2);
                         die "s1 > e1 when link orth $a and $b" if($s1 > $e1);
                         die "s2 > e2 when link orth $a and $b" if($s2 > $e2);
                         if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                         #$chr = $chr1;print STDERR "see??? $chr\n"; #not valid in sort {}
                         #$strand = $strand1;
                       } @$idx;
  for(my $i = 0;$i<= $#segs_sort-1;$i++){
     my ($chr1,$range1,$strand1) = split/:/,$segs_sort[$i];
     my ($s1,$e1) = split/-/,$range1;
     my ($chr2,$range2,$strand2) = split/:/,$segs_sort[$i+1];
     my ($s2,$e2) = split/-/,$range2;
     if(!defined $chr && !defined $strand){$chr = $chr1; $strand= $strand1}else{if($chr1 eq $chr && $strand1 eq $strand){}else{die "chr or strand diffs when link orth $segs_sort[$i] and $segs_sort[$i+1]"}    }
     if($s2 - $e1 <= $gap){$segs_sort[$i] = "NA";$segs_sort[$i+1] = "$chr:$s1-$e2:$strand"}else{die "broken seg when link orth at @segs_sort"}

  }

  my @seg_link;
  foreach my $ctrl(@segs_sort){if($ctrl ne "NA"){push @seg_link,$ctrl}}
  return \@seg_link;

}



sub readOrthPeak_simple(){
   my $file = shift;
   $/=">";
   my %tab;
   my $cnt;
   open TAB, $file or die "$!";
   while(<TAB>){
      chomp;
      next if($_ eq "" || $_=~/^\s+$/);
      $cnt++;
      my @box = split/\n+/,$_;
      my $flag = shift @box;
      #if($flag eq "1|1|1" || $flag eq "1|1|NF" || $flag eq "NF|1|1" || $flag eq "1|1|0" || $flag eq "1|0|1" || $flag eq "1|0|NF" || $flag eq "NF|0|1" ){
        foreach my $line (@box){
          my ($gene,$peak,$orth,$orth_core,$segs_core) = split/[\t ]+/,$line;
          if($gene=~/^#/){$gene=~s/^#//}
          my ($chr_gene,$range_gene,$id,$strand_gene) = split/\|/,$gene;#no ChrUN-Contig
          my ($s_gene,$e_gene) = split/-/,$range_gene;
          my $spec;
          if($chr_gene =~/^chr/){$spec = "japo"}elsif($chr_gene=~/^Chr/){$spec = "brac"}else{die "unknow chr $chr_gene"}
          #if($orth_core eq "NA"){print STDERR "orth_core eq NA at $gene\n";next}
          #my ($s_orth_core,$e_orth_core) = split/-/,$orth_core;
          #die "$s_orth_core >= $e_orth_core at $line" if($s_orth_core >= $e_orth_core);
          #my $len_orth_core = $e_orth_core - $s_orth_core +1;
          my @segs_core = split/\|/,$segs_core;
          my @segs_core_ori = @segs_core;
          #die "segs_core empty at $line" if(scalar @segs_core == 0);
          #modify segments range to orth_core
          if($orth_core ne "NA" && scalar @segs_core != 0){
            my ($min,$max) = split/-/,$orth_core;
            my ($left1,$right1) = split/-/,$segs_core[0];
            my ($left2,$right2) = split/-/,$segs_core[-1];
            if($min != $left1){$segs_core[0] = "$min-$right1"}
            if($max != $right2){$segs_core[-1] = "$left2-$max"}
          }

          $tab{$cnt}->{'flag'} = $flag;
          $tab{$cnt}->{$spec}->{'core_region'} = "$chr_gene:$orth_core";
          $tab{$cnt}->{$spec}->{'core_segs_ori'} = \@segs_core_ori;
          $tab{$cnt}->{$spec}->{'core_segs'} = \@segs_core;
          $tab{$cnt}->{$spec}->{'gene'} = $gene;
          $tab{$cnt}->{$spec}->{'peak'} = $peak;
          $tab{$cnt}->{$spec}->{'orth'} = $orth;
          $tab{$cnt}->{$spec}->{'geneId'} = $id;
          $tab{$cnt}->{$spec}->{'geneStrand'} = $strand_gene;
       }#foreach box end
     #}else{print STDERR "flag($flag) filtered at $box[0]\n";next}
     
   }#while end
   close TAB;
   $/="\n";

   return \%tab;
}




sub alnExtract(){
##extract aln with real coord in a given aln block
 my ($aln,$region,$range) = @_;
 my @aln = split//,$aln;
 my ($start,$end) = split/-/,$region;
 my ($chr_r,$start_r,$end_r,$strand_r) = split/:|-/,$range;
 if($start < $start_r){print STDERR "seg start < region start at $region and $range, cut as range  start\n" ;$start = $start_r}
 if($end > $end_r){print STDERR "seg end > region end at $region and $range, cut as range  end\n ";$end = $end_r}
 #if($start < $start_r){print STDERR "seg start < region start at $region and $range, use seg  start\n" ;$start_r = $start }
 #if($end > $end_r){print STDERR "seg end > region end at $region and $range, use seg  end\n ";$end_r = $end}
 #real coordinates to aln
 my $cnt_real = 0;
 #my $cnt =0;
 my ($start_aln,$end_aln);# = (0,$#$aln);


 foreach my $j(0..$#aln){
 #   if($aln[$j] ne "-"){$cnt_real++}
   if($start_r+$cnt_real == $start){$start_aln = $j}
   if($start_r+$cnt_real == $end ){$end_aln = $j}

=pod
   if($j == 0){ 
      if($start_r == $start){$start_aln = 0}
      if($start_r+$cnt_real == $end ){$end_aln = $j}
   }elsif($j == $#aln){
     if($start_r+$cnt_real == $start){$start_aln = $j}
     if($end_r == $end ){$end_aln = $j-1}
    }else{
       if($start_r+$cnt_real == $start){$start_aln = $j}
       if($start_r+$cnt_real == $end ){$end_aln = $j}
     }
=cut   

#   if($j == $#aln){
#      if($end_r == $end ){$end_aln = $j}
#   }else{
#     if($start_r+$cnt_real == $end ){$end_aln = $j}
#   }


   if($aln[$j] ne "-"){$cnt_real++}

 }

 die "detect start_aln fail at $aln\n$region\n$range\n" if(!defined $start_aln);
 die "detect end_aln fail at $aln\n$region\n$range\n" if(!defined $end_aln);


 my $aln_extract = substr($aln,$start_aln,$end_aln-$start_aln+1); #must +1 and start == end is ok
 print STDERR "#info seg:$region,seg_transformed:$start_aln-$end_aln,aln_extract:$aln_extract\n";
 return $aln_extract;

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

sub endRepair(){
  #fill gap for shorter one of 'not equal' pair alignments, try from both ends
  #do not do this if you have other choices
  #must warning out if you have done this step
  #autoguess left-end or right-end with better identity
  my ($aln1,$aln2) = @_;
  #$aln1 = uc($aln1);
  #$aln2 = uc($aln2);
  my @aln1 = split//,$aln1;
  my @aln2 = split//,$aln2;
  my $len1 = length $aln1;
  my $len2 = length $aln2;
  die "empty aln" if($len1 == 0 || $len2 == 0);
  die "equal aln" if($len1 == $len2);
  my $len_max;
  ($len1 > $len2)?($len_max = $len1):($len_max = $len2);
  my $len_min;
  ($len1 < $len2)?($len_min = $len1):($len_min = $len2);
  my $len_diff = abs($len1-$len2);
  #try left end align
  my $cnt_ident_left = 0;
  for my $i(0..$len_min-1){if(uc($aln1[$i]) eq uc($aln2[$i])){$cnt_ident_left++}}
  #try right end align
  my $cnt_ident_right = 0;
  for my $i($len_diff-1..$len_max-1){
    if($len1 > $len2){
      if(uc($aln1[$i]) eq uc($aln2[$i-$len_diff+1])){$cnt_ident_right++}
    }else{
       if(uc($aln1[$i-$len_diff+1]) eq uc($aln2[$i])){$cnt_ident_right++}
     }
  }

  if($cnt_ident_left >= $cnt_ident_right){ #repair for right end
    if($len1 > $len2){for (1..$len_diff){push @aln2,"-"} }elsif($len1 < $len2){for (1..$len_diff){push @aln1,"-"}}
  }elsif($cnt_ident_left < $cnt_ident_right){ #repair for left end
      if($len1 > $len2){for (1..$len_diff){unshift @aln2,"-"} }elsif($len1 < $len2){for (1..$len_diff){unshift @aln1,"-"}}
    }

  my $aln1_new = join"",@aln1;
  my $aln2_new = join"",@aln2;

  die "length diffs after end repair:$cnt_ident_left,$cnt_ident_right\naln1:$aln1\naln2:$aln2\naln1_new:$aln1_new\naln2_new:$aln2_new\n" if(length $aln1_new != length $aln2_new);

 return ($aln1_new,$aln2_new);

}


sub endCut(){
  #cut overhang end for longer one of 'not equal' pair alignments, try from both ends
  #do not do this if you have other choices
  #must warning out if you have done this step
  #autoguess left-end or right-end with better identity
  my ($aln1,$aln2,$range1,$range2) = @_;
  #$aln1 = uc($aln1);
  #$aln2 = uc($aln2);
  my @aln1 = split//,$aln1;
  my @aln2 = split//,$aln2;
  my $len1 = length $aln1;
  my $len2 = length $aln2;
  die "empty aln" if($len1 == 0 || $len2 == 0);
  die "equal aln" if($len1 == $len2);

  my ($s1,$e1) = split/-/,$range1;
  my ($s2,$e2) = split/-/,$range2;

  my $len_max;
  ($len1 > $len2)?($len_max = $len1):($len_max = $len2);
  my $len_min;
  ($len1 < $len2)?($len_min = $len1):($len_min = $len2);
  my $len_diff = abs($len1-$len2);
  #try left end align
  my $cnt_ident_left = 0;
  for my $i(0..$len_min-1){if(uc($aln1[$i]) eq uc($aln2[$i])){$cnt_ident_left++}}
  #try right end align
  my $cnt_ident_right = 0;
  for my $i($len_diff-1..$len_max-1){
    if($len1 > $len2){
      if(uc($aln1[$i]) eq uc($aln2[$i-$len_diff+1])){$cnt_ident_right++}
    }else{
       if(uc($aln1[$i-$len_diff+1]) eq uc($aln2[$i])){$cnt_ident_right++}
     }
  }
 
  my $cnt_real;
  my $base;
  #my ($s1_new,$e1_new);
  #my ($s2_new,$e2_new);
  if($cnt_ident_left >= $cnt_ident_right){ #cut longer one for right end
    if($len1 > $len2){
      for (1..$len_diff){ $base = pop @aln1;if($base ne "-"){$cnt_real++}} 
      $e1-=$cnt_real;
    }elsif($len1 < $len2){
       for (1..$len_diff){ $base = pop @aln1;if($base ne "-"){$cnt_real++}}
       $e2-=$cnt_real;
     }
  }elsif($cnt_ident_left < $cnt_ident_right){ #cut longer one for left end
      if($len1 > $len2){
         for (1..$len_diff){$base = shift @aln1;if($base ne "-"){$cnt_real++}} 
         $s1+=$cnt_real;
       }elsif($len1 < $len2){
          for (1..$len_diff){$base = shift @aln2;if($base ne "-"){$cnt_real++}}
          $s2+=$cnt_real;
        }
    }

  my $aln1_new = join"",@aln1;
  my $aln2_new = join"",@aln2;
  my $range1_new = "$s1-$e1";
  my $range2_new = "$s2-$e2";

  die "length diffs after end repair:$cnt_ident_left,$cnt_ident_right\naln1:$aln1\naln2:$aln2\naln1_new:$aln1_new\naln2_new:$aln2_new\n" if(length $aln1_new != length $aln2_new);

 return ($aln1_new,$aln2_new,$range1_new,$range2_new);

}


sub getSimilarDensi(){
  my ($index,$order,$len,$type)=@_;
  my @densi; #use score or perc;
  my @cons;
  my @Tsnp;
  my @flag_gaps; #len == aligned_len ;
  my @flag_N; #len == aligned_len ;
  my %align;
  my @ids;
  foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align{$id}=\@box;
      push @ids,$id;
  }
  for my $i (0.. $len-1){
     my @elements;
     foreach my $id(@$order){
       push @elements, $align{$id}->[$i];
     }
     my($score,$cons, $Tsnp,$flag_N, $flag_gap)=&checkIdentity_pair(\@elements);
     #my($flag, $score, $perc, $cons)=&checkIdentity_densi(\@elements);
     #if($type eq "score"){push @densi, $score}elsif($type eq "perc"){push @densi, $perc}else{die "don't understand densi type $type"}
     push @densi, $score;
     push @cons, $cons;
     push @Tsnp,$Tsnp;
     push @flag_gaps,$flag_gap;
     push @flag_N, $flag_N;
  }
  return(\@densi, \@cons, \@Tsnp, \@flag_N, \@flag_gaps ) ;
}


sub checkIdentity_pair(){
    #only compare two species
    #upcase and check all then compare    
    #use bracRS2 as ref (has very few N, only ATGC-)
    #return 2(all identical), 0 (not identical)
    #return flag_gap flag_N cons_base score
    my @elements=@{shift @_};
    die "more than 2 or empty at @elements" if(scalar @elements != 2);
    my ($ref,$que) = @elements;
    $ref = uc($ref); #first element as ref
    $que = uc($que);
    if($ref=~/[^-ATGCN]/i){die "illegal character found in ref $ref"} #nucleotide only
    if($que=~/[^-ATGCN]/i){die "illegal character found in que $que"} 

    die "both gap -" if($ref eq "-" && $que eq "-");
    die "both N " if($ref eq "N" && $que eq "N");

    my ($score,$cons,$Tsnp,$flag_N,$flag_gap);
    # score: 0(not match),2(match)
    # Tsnp : S(A<->G or C<->T transition,Ts),V(other,transversion,Tv), relate to ref
    # CpGsnp : mismatch in a CpG dinucleotide (Yes/No/Na)
    # flag_N: 0(none) , 1 (ref), 2(que)
    # flag_gap: 0(none) , 1 (ref), 2(que), 3(both N or -)
    $Tsnp = " ";
    if($ref=~/[ATGC]/){ 
       if($ref eq $que){ $score = 2; $cons=$ref; $flag_N = 0; $flag_gap = 0 } #identical and not N or -
       elsif($ref ne $que && $que eq "-"){  # not identical and is  -, #InDel, N similar to -
         $score = 0; $cons=" "; $flag_N = 0; $flag_gap = 2             
       }elsif( $ref ne $que && $que eq "N"){ # not identical and is  N
          $score = 0; $cons=" "; $flag_N = 2; $flag_gap = 2
        }elsif($ref ne $que && $que=~/[ATGC]/){ # not identical and is ATCG, SNP
           $score = 0; $cons=" "; 
           if($ref eq "A" && $que eq "G" || $ref eq "G" && $que eq "A" || $ref eq "C" && $que eq "T" || $ref eq "T" && $que eq "C"){$Tsnp = "S"}else{$Tsnp = "V"};            
           $flag_N = 0; $flag_gap = 0 
         }else{die "unknow situation at ref:$ref and que:$que"}

    }elsif($ref eq "-"){ 
       if($ref eq $que){ die "both gap -"}
       elsif($ref ne $que && $que=~/[ATGC]/){  # not identical and is  ATCG
         $score = 0; $cons=" "; $flag_N = 0; $flag_gap = 1
       }elsif( $ref ne $que && $que eq "N"){ # not identical and is  N
          $score = 0; $cons=" "; $flag_N = 2; $flag_gap = 3
         }else{die "unknow situation at ref:$ref and que:$que"}

     }elsif($ref eq "N"){ #very few
       if($ref eq $que){ die "both gap N"}
       elsif($ref ne $que && $que=~/[ATGC]/){  # not identical and is  ATCG
         $score = 0; $cons=" "; $flag_N = 1; $flag_gap = 1
       }elsif( $ref ne $que && $que eq "-"){ # not identical and is  -
          $score = 0; $cons=" "; $flag_N = 1; $flag_gap = 3
         }else{die "unknow situation at ref:$ref and que:$que"}

      }else{die "ref err: $ref:$que"} #should not happen

    return ($score,$cons,$Tsnp,$flag_N,$flag_gap);
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
 
sub alnStat(){
  #simple stats for two alns with few gaps
  #detect snp type: Tv,Ts and CpG,nonCpG,Unknown
  my ($aln1,$aln2) = @_;
  my @aln1 = split//,$aln1;
  my @aln2 = split//,$aln2;
  my $n1 = scalar @aln1;
  my $n2 = scalar @aln2;
  die "n1 != n2, $n1 != $n2" if($n1 != $n2);

  my @Ts;
  my @Tv;
  my @CpGsnp;
  my $cnt_ident = 0;
  for my $i(0..$#aln1){
    my ($ref,$que) = ($aln1[$i],$aln2[$i]);
    $ref = uc($ref); #first element as ref
    $que = uc($que);
    if($ref=~/[^-ATGCN]/i){die "illegal character found in ref $ref"} #nucleotide only
    if($que=~/[^-ATGCN]/i){die "illegal character found in que $que"}
    #die "both gap -" if($ref eq "-" && $que eq "-");
    die "both N " if($ref eq "N" && $que eq "N");

    my $flag_cpg=0;
    if($i!=$#aln1){
      my ($ref_next,$que_next) = ($aln1[$i+1],$aln2[$i+1]);
      $ref_next = uc($ref_next); #first element as ref
      $que_next = uc($que_next);
      if($ref_next=~/[^-ATGCN]/i){die "illegal character found in ref $ref_next"} #nucleotide only
      if($que_next=~/[^-ATGCN]/i){die "illegal character found in que $que_next"}
      #die "both gap -" if($ref_next eq "-" && $que_next eq "-");
      die "both N " if($ref_next eq "N" && $que_next eq "N");
      if($ref.$ref_next eq "CG" || $que.$que_next eq "CG"){$flag_cpg = 1}


      if($i!=0){
        my ($ref_pre,$que_pre) = ($aln1[$i-1],$aln2[$i-1]);
        $ref_pre = uc($ref_pre); #first element as ref
        $que_pre = uc($que_pre);
        if($ref_pre=~/[^-ATGCN]/i){die "illegal character found in ref $ref_pre"} #nucleotide only
        if($que_pre=~/[^-ATGCN]/i){die "illegal character found in que $que_pre"}
        #die "both gap -" if($ref_pre eq "-" && $que_pre eq "-");
        die "both N " if($ref_pre eq "N" && $que_pre eq "N");
        if($ref_pre.$ref eq "CG" || $que_pre.$que eq "CG"){$flag_cpg = 1}
      }

    }
    
    if($ref ne $que && $que=~/[ATGC]/ && $ref=~/[ATGC]/){ # not identical and is ATCG, SNP       
      if($ref eq "A" && $que eq "G" || $ref eq "G" && $que eq "A" || $ref eq "C" && $que eq "T" || $ref eq "T" && $que eq "C"){      push @Ts, "$i:S";
       if($flag_cpg == 1){push @CpGsnp,"$i"}
      }else{
        push @Tv, "$i:V";
        if($flag_cpg == 1){push @CpGsnp,"$i"}
      }
    }elsif($ref eq $que && $que=~/[ATGC]/ && $ref=~/[ATGC]/){$cnt_ident++}

  }#for end
  my $n_CpGsnp = scalar @CpGsnp;
  my $n_Ts = scalar @Ts;
  my $n_Tv = scalar @Tv;
  my $n_snp = $n_Ts + $n_Tv;
  my $perc_cpgsnp;
  if($n_snp != 0){$perc_cpgsnp = sprintf("%.3f",$n_CpGsnp/$n_snp)}else{$perc_cpgsnp = 0}
  return (sprintf("%.3f",$cnt_ident/$n1),$n_snp,$n_CpGsnp,$perc_cpgsnp);
  #return ($cnt_ident,\@Ts,\@Tv,\@CpGsnp);
  #return (sprintf("%.3f",$cnt_ident/$n1),\@Ts,\@Tv,\@CpGsnp);

}#sub end



sub cpgCount(){
   my $seq = shift;
   $seq = uc($seq);
   die "seq empty" if($seq eq "");
   my $cnt=0;
   while($seq=~/CG/g){$cnt++}
   return $cnt;
}

sub extractDensi(){
   my ($file,$region,$flag_rev)=@_;
   die "file $file not exists\n" if(! -f $file);
   my ($chr,$s,$e)=split/:|-/,$region;
   my $len = $e-$s+1;
   die "len <= 0 when extract densi" if($len <=0);
   my @data_region;
   if($file=~/bedgraph\.bgz$/ || $file=~/density\.gz$/){ @data_region=`tabix $file $region `}elsif($file=~/bw$/){
     @data_region=`bigWigToBedGraph -chrom=$chr -start=$s -end=$e $file stdout `
   }
  
   if(scalar @data_region == 0){die "too few data at region $region in file $file\n"}
   my @results;
   foreach my $ctrl(@data_region){
      chomp $ctrl;
      my ($chr,$s,$e,$v) = split/[\t ]+/,$ctrl;
      push @results,$v
   }

   #if($flag_rev == 1){@results = reverse @results}
   return (\@results);

}


