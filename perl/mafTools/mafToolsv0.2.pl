use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#output to stdout or output file without headlines
#stderr give necessery message

#Phast: MAFs used in phastCons need to have the first (reference) sequence on the '+' strand,be non-overlapping and sorted by their chromStart coordinate.
#mafstrander: when the block is flipped,  all start coordinates are transformed, and all sequence fields are reverse-complemented. https://github.com/dentearl/mafTools/tree/master/mafStrander


#result maf:  1,must has ref in each block, ref is the first line in block
#             2,species order agree with -order (empty species allowed)
#             3,ref start coordinate sorted and nonOverlap_checked
#             4,ref strand must be + (flip whole block when not)
             

#process procedure:  (stat) -> filter -> sort,overlapcheck -> strander -> refDegap -> write with order,ref_checked_again 
#                    each step standalone                 
#v0.2 Dec. 1 2016



my ($maf_file,$ref,$filter,$order,$sort,$overlap,$stat,$colinear, $refDegap, $strander, $prefix);
GetOptions("maf=s",\$maf_file,"ref=s",\$ref,"filter!",\$filter,"order=s",\$order,"sort!",\$sort,"overlap!",\$overlap,"stat!",\$stat,"colinear!",\$colinear,"refDegap!",\$refDegap, "strander!",\$strander,"prefix=s",\$prefix);
$prefix ||= "out-";
die "maf_file is a must" if(!defined $maf_file || $maf_file eq ""); # maf_file is a must
die "ref is a must" if(!defined $ref || $ref eq ""); # ref is a must


# 1.readin maf
my $maf = &mafRead($maf_file); #exit;
my @IDs = sort {$a <=> $b} keys %$maf;
my $groupID = \@IDs;
#print Dumper $maf;exit;


# 2,simple statistics #maf hash changed after this step

if($stat){
  print STDERR "stats..(this step cost much time)\n";
  my ($blockn,$min_degree,$max_degree,$min_aliLen,$max_aliLen,$mean_aliLen,$maf_out,$groupID_out) = &mafStat($maf,$groupID);
  print "blockn,min_degree,max_degree,min_aliLen,max_aliLen,mean_aliLen:\n$blockn,$min_degree,$max_degree,$min_aliLen,$max_aliLen,$mean_aliLen\n";
}
#exit;#print Dumper $maf;exit;


# 3,filter   
  my $groupID_filter;
  if ($filter) { 
    
     print STDERR "filter.. "; 
     ($groupID,$groupID_filter) = &mafFilter( $maf, $groupID, $ref, 3 , 2000, 1000, 0.1, 0.1, 0.8);     
                           #$maf, $groupID,$ref, $degree, $aliLen, $score, $perc_gapAny, $perc_nAny, $perc_ident_big60
     my $n = scalar @$groupID;
     print STDERR "$n blocks remained after filter\n";
   }


# 4,sort by start position of given species, check overlap(remove shorter, low degree ones)

  if($sort){ 
    print STDERR "refsort..\n"; 
    $groupID = &mafSort($maf,$groupID, $ref, $overlap)
  }


# 5,refDegap 
  if ($refDegap){    
     print STDERR "refDegap..\n"; 
     $groupID = &maf_refDegap($maf, $groupID, $ref);
     #print Dumper $maf; exit;
 
  }


# 6, ref strander
  my $groupID_strander;
  if ($strander){    
     print STDERR "refstrander..\n"; 
     ($groupID, $groupID_strander) = &mafStrander($maf, $groupID, $ref);
     #print Dumper $maf,$groupID_strander; exit;
     print STDERR "these maf groups have been faRC_coord_gap:@$groupID_strander\n";
  }


# 7, write with order or not
#=pod
   if ($order) {
     print STDERR "write with order ($order)..\n";
     my @specOrder = split/-/,$order;
     die "order not include ref or not at the begining"if ($specOrder[0] ne $ref);
     if($refDegap){
       print STDERR "write refDegaped alignment..\n";
       &mafWrite_subset_order_refdegap($maf,$groupID,\@specOrder,$prefix,$ref);
     }else{
       &mafWrite_subset_order($maf,$groupID,\@specOrder,$prefix);     
     }
   } else {
            print STDERR "write maf in original..\n";
            if($refDegap){
                print STDERR "write refDegaped alignment..\n";
                &mafWrite_subset_order_refdegap($maf,$groupID,"asis",$prefix,$ref)
            }else{
              &mafWrite_subset_order($maf,$groupID,"asis",$prefix);
            }
          }

# 8, write filtered groups  into $prefix.dump.maf  file
  if($filter){
     print STDERR "output filtered maf groups\n";
     &mafWrite($maf,$groupID_filter,"out.filtered.dump");
  }


#=cut





##sub

sub mafRead(){
  #read and check maf file
  #s japoChr01      1554 120 + 43268879 CTATCTAGGCATCCATCCGATATTTGGAGTATGGAGGAGAAAAACAGTGCTCCAGCAGAGTCTCCATCACATGCTTCATTTTTGG
  #s OpuncChr01 30196520 120 + 46096743 CTATCCCACCCTTCATATGAGAAATAGAGTATGTAAGCAAAAAAAGAGACTCCAGCAGACACTCCAAAATATCCTCCAAAAATAG
  #s LperrChr01  1183910 106 + 32922458 CCATCCCATACTCCATCCTATATTTGGTATATATGGAAGGAAAAATGGGCTCCAGTA----------TATATACCCATAAACTAG
  my $file = shift;
  my %maf;
  my $cnt=0;
  open MAF, $file or die "$!";
  $/ = "\n\n";
  while(<MAF>){
    chomp;
    last if ($_=~/##eof/); #modified on 28, Dec 2016
    my @box=split /\n/,$_;
    my $block_score;
    my %block;
    my @order; 
    my $degree = 0;
    foreach my $line(@box){
      #exit if ($line=~/##eof/);
      #last if ($line=~/##eof/); #modified on 12, May 2013
      next if ($line=~/#/);
      if($line=~/^a score/){
          my @tmp=split/=/,$line;$block_score=int($tmp[1]);
          if(!exists $block{"score"}){$block{"score"} = $block_score}else{die "dup score line: $line"}
        }elsif($line=~/^s /){
              my @tmp=split/ +/,$line;
              my ($name,$start,$len,$strand,$chr_len,$seq_ali) = ($tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6]); 
              push @order, $name;
              $degree++;
              my $seq_ori=$seq_ali; $seq_ori=~s/-//g;
              my $end=$start+$len;
              my @dataline=($name,$start,$end,$len,$strand,$chr_len,$seq_ali,$seq_ori);
              die "length $len diff with seq_rmgap at block $cnt" if($len != length $seq_ori);
              my $alignLen = length $seq_ali;
              if(exists $block{$name}){die"duplicated species in maf block $cnt\n"}else{$block{$name}=\@dataline}
              if(exists $block{"alignLen"}){die "align length diff at block $cnt" if($block{"alignLen"} != $alignLen)}else{$block{"alignLen"} = $alignLen}
         }else{print STDERR "\nomit line at group $cnt\n: unknow maf line type:$line\n"}
    }#foreach end here 
    if($degree == 0){die "degree eq 0 at group $cnt: $_"}else{$block{"degree"} = $degree}
    $block{"order"} = \@order;   
    $maf{$cnt} = \%block;
    $cnt++;
    #last if($cnt == 10000);
    print STDERR "#" if($cnt % 5000 == 0);
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks read in\n";
  return \%maf;
}


sub mafStat(){
  # warning: maf hash add info
  my ($maf,$groupID) = @_;
  my $cnt=0;
  my ($degree_min,$degree_max);
  my ($min_aliLen,$max_aliLen);
  my $sum_aliLen = 0;
  foreach my $id(@$groupID){
    $cnt++;
    my $degree = $maf->{$id}->{"degree"};
    if( !defined $degree_min){ $degree_min = $degree }else{if($degree_min > $degree){ $degree_min = $degree}}
    if( !defined $degree_max){$degree_max = $degree}else{if($degree_max < $degree){ $degree_max = $degree}}
    my $aliLen = $maf->{$id}->{"alignLen"};
    if( !defined $min_aliLen){$min_aliLen = $aliLen}else{if( $min_aliLen > $aliLen ){$min_aliLen = $aliLen}}
    if( !defined $max_aliLen){$max_aliLen = $aliLen}else{if($max_aliLen < $aliLen){$max_aliLen = $aliLen}}
    $sum_aliLen += $aliLen; 
    #getPerc
    my %align;
    foreach my $spec(@{$maf->{$id}->{"order"}}){ $align{$spec} = $maf->{$id}->{$spec}->[6]}
    my ($perc_gapAny,$perc_nAny,$perc_ident_big60) = &getPerc(\%align);
    $maf->{$id}->{"perc_gapAny"} = $perc_gapAny;
    $maf->{$id}->{"perc_nAny"} = $perc_nAny;
    $maf->{$id}->{"perc_ident_big60"} = $perc_ident_big60;
  }
  return ($cnt,$degree_min,$degree_max,$min_aliLen,$max_aliLen,sprintf("%.1f",$sum_aliLen/$cnt));
}


sub getPerc(){
    #check align length
    #check illegal character
    #global gap_any(-), N_any and ident(>60%) percent
    my $align = shift;
    my ($flag , $len) = &checkAlignLenAll($align);
    if($flag != 1){print STDERR "align length diff when refDegap\n";print Dumper $align;exit}
    my %align_split;
    my ($cnt_gapAny,$cnt_nAny,$cnt_ident_big60)=(0,0,0);
    foreach my $spec(keys %$align){
       if($align->{$spec}=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in $align->{$spec}"} #nucleotide and pep alphabet
       if(! &isNuc($align->{$spec}) ){die "not in nucleotide name space(peptide?):$align->{$spec}"}
       my @box = split//,$align->{$spec};
       $align_split{$spec} = \@box;
    }
    
    
    for my $i(0..$len-1){
      my ($flag_gap,$flag_N,$cnt_simiBase)=(0,0,0);
      my $base;
      my $cnt;
      foreach my $spec(keys %align_split){
          $cnt++;
          if($align_split{$spec}->[$i] eq "-"){ $flag_gap=1} 
          if( $align_split{$spec}->[$i] eq "n" || $align_split{$spec}->[$i] eq "N"){ $flag_N=1}
          if(!defined $base){$base = $align_split{$spec}->[$i]; $cnt_simiBase++ }else{if($base eq $align_split{$spec}->[$i] || $base eq uc($align_split{$spec}->[$i])){ $cnt_simiBase++  }}
      }
      if($flag_gap == 1){$cnt_gapAny++}
      if($flag_N == 1){$cnt_nAny++}
      if($cnt_simiBase/$cnt >=0.6 ){$cnt_ident_big60++}
    }

    return(sprintf("%.1f",$cnt_gapAny/$len),sprintf("%.1f",$cnt_nAny/$len),sprintf("%.1f",$cnt_ident_big60/$len));

}

sub isNuc(){
   #not pep 
   my $seq=shift;
   if($seq=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in $seq"} #nucleotide and pep alphabet
   my @seqs=split//,$seq;
   my $isNuc=0;
   my $total=0;
   foreach(@seqs){
     $total++;
     if($_ eq "A" ||$_ eq "a" || $_ eq "T" || $_ eq "t" || $_ eq "C" || $_ eq "c" || $_ eq "G" || $_ eq "g" || $_ eq "N" || $_ eq "n" || $_ eq "-"){$isNuc++}
   }
   if($isNuc/$total >= 0.9){return 1}else{return 0} #tolerate for some consolidated characters
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



sub mafSort(){
  #sort by given species name start position, output block id;
  #if given species not exists, append to tail 
  my ($maf, $groupID, $ref) = @_;
  my @orderID;
  my @hasRef;
  my @tail;
  foreach my $id(@$groupID){if(!exists $maf->{$id}->{$ref}){push @tail,$id}else{push @hasRef,$id}}

  @orderID = sort {
                    my $start1 = $maf->{$a}->{$ref}->[1]; 
                    my $start2 = $maf->{$b}->{$ref}->[1]; 
                    if($start1 > $start2){ return 1}elsif($start1 < $start2){return -1}else{return 0}
                   } @hasRef;
  my @resultID = (@orderID,@tail);
  return \@resultID;
 
}


sub mafFilter(){
   #filter by hasref, degree, aliLen, score and gap_perc without considering the group order
   my ($maf, $groupID,$ref, $degree, $aliLen, $score, $perc_gapAny, $perc_nAny, $perc_ident_big60) = @_;
   my @groupId_good;
   my @groupId_filter;
   foreach my $id (@$groupID){
      my $flag = 1;
      if($maf->{$id}->{"degree"} >= $degree){}else{$flag = 0}
      if($maf->{$id}->{"alignLen"} >= $aliLen){}else{$flag = 0}
      if($maf->{$id}->{"score"} >= $score){}else{$flag = 0}
      if(!exists $maf->{$id}->{$ref}){$flag = 0}
      if(exists $maf->{$id}->{"perc_gapAny"}){if($maf->{$id}->{"perc_gapAny"} >= $perc_gapAny ){$flag = 0}}
      if(exists $maf->{$id}->{"perc_nAny"}){if($maf->{$id}->{"perc_nAny"}  >= $perc_nAny ){$flag = 0}}
      if(exists $maf->{$id}->{"perc_ident_big60"}){if($maf->{$id}->{"perc_ident_big60"} < $perc_ident_big60 ){ $flag = 0}}  
      if($flag == 1){push @groupId_good, $id}else{push @groupId_filter, $id}
   }
   return \@groupId_good,\@groupId_filter;
}



sub mafStrander(){
  my ($maf, $groupID, $ref) = @_;
  #mafstrander: when the block is flipped,  all start coordinates are transformed, and all sequence fields are reverse-complemented. https://github.com/dentearl/mafTools/tree/master/mafStrander
  my @groupID_strander;
  foreach my $id (@$groupID){
    die "not exists $ref" if( !exists $maf->{$id}->{$ref});
    if($maf->{$id}->{$ref}->[4] ne "+"  && $maf->{$id}->{$ref}->[4] eq "-"){
      push @groupID_strander,$id;
      foreach my $spec(@{$maf->{$id}->{'order'}}){
        #print STDERR "spec: $spec\n";
        my ($seq_RC,$start_transform,$end_transform,$strand_transform) = &faRC_coord_gap($maf->{$id}->{$spec}->[6],$maf->{$id}->{$spec}->[5],$maf->{$id}->{$spec}->[1],$maf->{$id}->{$spec}->[2],$maf->{$id}->{$spec}->[4]);      
        #$name,$start,$end,$len,$strand,$chr_len,$seq_ali,$seq_ori
        $maf->{$id}->{$spec}->[6] = $seq_RC;
        $maf->{$id}->{$spec}->[1] = $start_transform;
        $maf->{$id}->{$spec}->[2] = $end_transform;
        $maf->{$id}->{$spec}->[4] = $strand_transform;

      }
    }elsif($maf->{$id}->{$ref}->[4] ne "-"  && $maf->{$id}->{$ref}->[4] eq "+"){}else{die "unknow strand $maf->{$id}->{$ref}->[4] at $id of ref $ref"} 


  }
  return ($groupID,\@groupID_strander);
}



sub faRC_coord_gap(){
  #faRC with coordinates and gap
  #start end related to this strand
  my ($seq,$chr_len,$start,$end,$strand) = @_;
  #print STDERR "$seq,$chr_len,$start,$end,$strand\n";
  $seq = reverse $seq;
  $seq =~tr/atcgn\-ATCGN/tagcn\-TAGCN/;

  my $start_transform = $chr_len - $end;
  my $end_transform = $chr_len - $start;
  die "start <=0 || start > $chr_len" if($start <= 0 || $start > $chr_len);
  die "end <=0 || end > $chr_len" if($end <= 0 || $end > $chr_len);
  die "start == end " if($start == $end);

  if($strand eq "+"){$strand = "-"}elsif($strand eq "-"){$strand = "+"}else{die "unknow strand:$strand"};

  return ($seq,$start_transform,$end_transform,$strand);

}


sub maf_refDegap(){
    #warning: maf hash add information
    my ($maf, $groupID, $ref) = @_;
    my @groupID_good;
    foreach my $id(@$groupID){
       if( !exists $maf->{$id}->{$ref}){next}else{push @groupID_good, $id};
       my %align;
       foreach my $spec (@{$maf->{$id}->{"order"}}){ $align{$spec} = $maf->{$id}->{$spec}->[6] }       
       my ($align_refDegap,$align_deNum) = &refDegap(\%align,$ref);
       $maf->{$id}->{"align_refDegap"} = $align_refDegap;
       $maf->{$id}->{"align_deNum"} = $align_deNum;
    }

    return ($maf,\@groupID_good);
}



sub refDegap(){
  #malign in, refDegap_malign out
  my ($align,$ref) = @_; #align hash
  my ($flag,$len) = &checkAlignLenAll($align);
  if($flag != 1){print STDERR "align length diff when refDegap\n";print Dumper $align;exit}
  my %align_split;
  my %align_refDegap;
  my %align_deNum;
  foreach my $spec(keys %$align){
     my @box = split//,$align->{$spec};
     $align_split{$spec} = \@box;
     $align_deNum{$spec}=0;
  }

  for my $i(0..$len-1){
     if($align_split{$ref}->[$i] ne "-"){
        foreach my $spec(keys %align_split){ $align_refDegap{$spec}.=$align_split{$spec}->[$i] }
     }else{
           foreach my $spec(keys %align_split){
              if($spec eq $ref){ $align_deNum{$ref}++;next};
              if($align_split{$spec}->[$i] ne "-"){$align_deNum{$spec}++}

           }
          }
  }


  return (\%align_refDegap,\%align_deNum);

}


sub mafWrite(){ 
   #write with given order, as is if 'order' not supplied
   my ($maf,$groupID,$prefix) = @_;
   open OUT, ">${prefix}.maf" or die "$!";
   print OUT "##maf version 12\n";
   my $cnt;
   foreach my $id(@$groupID){ 
      die "empty id $id when write maf" if($id eq "" || !defined $id);
      $cnt++;
      print OUT "a score=$maf->{$id}->{'score'}\n";
      foreach my $spec(@{$maf->{$id}->{'order'}}){
         if(exists $maf->{$id}->{$spec}){
           printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3], $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{$spec}->[6]);
         }else{die "not exist $spec in group $id"}
      } 
      print OUT "\n";
   }   
   print OUT "##eof maf\n";
   close OUT;
   print STDERR "finally output $cnt maf blocks into $prefix.maf\n";
   return 1;
}

sub mafWrite_subset_order(){
   #write with given order, as is if 'order' not supplied
   my ($maf,$groupID, $ref,$specOrder,$prefix) = @_;
   open OUT, ">$prefix.maf" or die "$!";
   print OUT "##maf version 12\n";
   my $cnt;
   foreach my $id(@$groupID){ #output given groupIDs  only
      die "empty id $id when write maf" if($id eq "" || !defined $id);
      if(!exists $maf->{$id}->{$ref}){ print STDERR "group$id has no ref:$ref, skip\n";next}
      $cnt++;
      print OUT "a score=$maf->{$id}->{'score'}\n";
      if( $specOrder ne "asis" && scalar @$specOrder != 0){
        foreach my $spec(@$specOrder){ #use given species order
           if(exists $maf->{$id}->{$spec}){
             printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3], $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{$spec}->[6]);
           }
        }
        print OUT "\n";
      }elsif($specOrder eq "asis"){  #use within group order 
         foreach my $spec(@{$maf->{$id}->{'order'}}){
           if(exists $maf->{$id}->{$spec}){
             printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3], $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{$spec}->[6]);
           }else{die "not exist $spec in group $id"}
         } 
         print OUT "\n";
       }
   }   
   print OUT "##eof maf\n";
   close OUT;
   print STDERR "finally output $cnt maf blocks into $prefix.maf\n";
   return 1;
}



sub mafWrite_subset_order_refdegap(){
   #write with given order, as is if 'order' not supplied
   my ($maf,$groupId, $specOrder,$prefix,$ref) = @_;
   #print"ref is $ref\n";
   open OUT, ">$prefix.maf" or die "$!";
   print OUT "##maf version 12\n";
   my $cnt;
   foreach my $id(@$groupId){ #output given groupIDs  only
      die "empty id $id when write maf" if($id eq "" || !defined $id);
      if(!exists $maf->{$id}->{$ref}){ print STDERR "group$id has no ref:$ref, skip\n";next}
      $cnt++;
      print OUT "a score=$maf->{$id}->{'score'}\n";
      if( $specOrder ne "asis" && scalar @$specOrder != 0){
        foreach my $spec(@$specOrder){ #use given species order
           if(exists $maf->{$id}->{$spec}){
             if(exists $maf->{$id}->{'align_refDegap'}->{$spec} && exists $maf->{$id}->{'align_deNum'}->{$spec}){
               if($spec eq $ref){
                 printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3], $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{'align_refDegap'}->{$spec})
               }else{
                  printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3]-$maf->{$id}->{'align_deNum'}->{$spec}, $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{'align_refDegap'}->{$spec})
                }
             }else{print STDERR "not exist align_refDegap or align_deNum at $id in spec $spec\n"}#;next}
           }else{print STDERR "not exist $spec in group $id\n"}#;next}
        }
        print OUT "\n";
      }elsif($specOrder eq "asis"){  #use within group order 
         foreach my $spec(@{$maf->{$id}->{'order'}}){
           if(exists $maf->{$id}->{$spec}){
             if(exists $maf->{$id}->{'align_refDegap'}->{$spec} && exists $maf->{$id}->{'align_deNum'}->{$spec}){
               if($spec eq $ref){
                 printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3], $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{'align_refDegap'}->{$spec})
               }else{
                  printf OUT ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $maf->{$id}->{$spec}->[1], $maf->{$id}->{$spec}->[3]-$maf->{$id}->{'align_deNum'}->{$spec}, $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $maf->{$id}->{'align_refDegap'}->{$spec})
                }
             }else{die "not exist align_refDegap or align_deNum at $id in spec $spec"}
           }else{die "not exist $spec in group $id"}
         }
         print OUT "\n";
       }
   }
   print OUT "##eof maf\n";
   close OUT;
   print STDERR "finally output $cnt maf blocks into $prefix.maf\n";
   return 1;
}




