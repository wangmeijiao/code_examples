use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#output to stdout or output file without headlines
#stderr give necessery message

my ($maf_file,$ref,$filter,$order,$sort,$stat,$colinear, $refDegap, $prefix);
GetOptions("maf=s",\$maf_file,"ref=s",\$ref,"filter!",\$filter,"order=s",\$order,"sort!",\$sort,"stat!",\$stat,"colinear!",\$colinear,"refDegap!",\$refDegap, "prefix=s",\$prefix);
$prefix ||= "out-";


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



# 3,sort by start position of given species
if($sort){ 
   print STDERR "sort..\n"; 
   $groupID = &mafSort($maf,$ref)}


# 4,filter 
#
  if ($filter) { 
    
     print STDERR "filter.. "; 
     $groupID = &mafFilter( $maf, $groupID, 3 , 2000, 1000, 0.1, 0.1, 0.8);     
                           #$maf, $groupID,$degree, $aliLen, $score, $perc_gapAny, $perc_nAny, $perc_ident_big60
     my $n = scalar @$groupID;
     print STDERR "$n blocks remained after filter\n";
   }

# 5,refDegap 
  my $maf_refDegap;
  if ($refDegap){    
     print STDERR "refDegap..\n"; 
     ($maf_refDegap,$groupID) = &maf_refDegap($maf, $groupID, $ref);
     #print Dumper $maf; exit;
 

  }


# 6,write with order or not
#=pod
   if ($order) {
     print STDERR "write with order ($order)..\n";
     my @specOrder = split/-/,$order;
     if($refDegap){
       print STDERR "write refDegaped alignment..\n";
       &mafWrite_subset_order_refdegap($maf_refDegap,$groupID,\@specOrder,$prefix,$ref);
     }else{
       &mafWrite_subset_order($maf,$groupID,\@specOrder,$prefix);     
     }
   } else {
            print STDERR "write maf in original..\n";
            if($refDegap){
                print STDERR "write refDegaped alignment as is..\n";
                &mafWrite_subset_order_refdegap($maf_refDegap,$groupID,"asis",$prefix,$ref)
            }else{
              &mafWrite($maf,$groupID,$prefix)
            }
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
   #filter by degree, aliLen, score and gap_perc without considering the group order
   my ($maf, $groupID,$degree, $aliLen, $score, $perc_gapAny, $perc_nAny, $perc_ident_big60) = @_;
   my @groupId_good;
   foreach my $id (@$groupID){
      my $flag = 1;
      if($maf->{$id}->{"degree"} >= $degree){}else{$flag = 0}
      if($maf->{$id}->{"alignLen"} >= $aliLen){}else{$flag = 0}
      if($maf->{$id}->{"score"} >= $score){}else{$flag = 0}
      if(exists $maf->{$id}->{"perc_gapAny"}){if($maf->{$id}->{"perc_gapAny"} >= $perc_gapAny ){$flag = 0}}
      if(exists $maf->{$id}->{"perc_nAny"}){if($maf->{$id}->{"perc_nAny"}  >= $perc_nAny ){$flag = 0}}
      if(exists $maf->{$id}->{"perc_ident_big60"}){if($maf->{$id}->{"perc_ident_big60"} < $perc_ident_big60 ){ $flag = 0}}  
      if($flag == 1){push @groupId_good, $id}
   }
   return \@groupId_good;
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
   my ($maf,$groupID, $specOrder,$prefix) = @_;
   open OUT, ">$prefix.maf" or die "$!";
   print OUT "##maf version 12\n";
   my $cnt;
   foreach my $id(@$groupID){ #output given groupIDs  only
      die "empty id $id when write maf" if($id eq "" || !defined $id);
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
   print"ref is $ref\n";
   open OUT, ">$prefix.maf" or die "$!";
   print OUT "##maf version 12\n";
   my $cnt;
   foreach my $id(@$groupId){ #output given groupIDs  only
      die "empty id $id when write maf" if($id eq "" || !defined $id);
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




