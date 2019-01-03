use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#readin mafExtract maf blocks, output concatenanted mfa align file, need to has ref chr, fill missing spec with *
#no need to fill gaps between blocks because there are extract regions in the mfa headline

my($maf_file,$width,$ref,$fill);
GetOptions("maf=s",\$maf_file,"width:i",\$width,"ref=s",\$ref);
$width||=50;

my ($ref_chr,$region,$ref_strand) =split/:/,$ref;
my ($ref_start,$ref_end) = split/-/,$region;
my $ref_base = $ref_chr;
$ref_base=~s/\d{1,2}$//;

my ($maf,$maf_cnt)=&mafRead_refchr($maf_file,$ref_base);
print STDERR "read maf done with refbase $ref_base\n";
#print Dumper $maf_cnt,$maf;exit; 


#my ($maf_concated,$max_order,$aliLen,$concat_regions) = &mafConcat_withcut($maf->{$ref_chr},$ref,$ref_chr);
my ($maf_concated,$max_order,$aliLen,$concat_regions) = &mafConcat_withcut($maf->{$ref_chr},$ref,"no");
#print Dumper $maf_concated,$max_order,$aliLen,$concat_regions;

#output mfa
my $cnt;
foreach my $spec(@$max_order){
   $cnt++;
   my $faStr_format = &faFormat($maf_concated->{$spec},$width);
   #if($cnt == 1){
     print ">$spec $aliLen ".join("|",@{$concat_regions->{$spec}})."\n$faStr_format"
   #}else{print ">$spec $aliLen\n$faStr_format"} 

}

#output SS
#my $ss = &align2SS($maf_concated,$max_order,$aliLen);
#foreach my $ctrl(0..$#$ss){ print $ctrl,"\t",$ss->[$ctrl],"\n"}



#sub

sub mafRead_refchr(){
  #read and check maf file
  #s japoChr01      1554 120 + 43268879 CTATCTAGGCATCCATCCGATATTTGGAGTATGGAGGAGAAAAACAGTGCTCCAGCAGAGTCTCCATCACATGCTTCATTTTTGG
  #s OpuncChr01 30196520 120 + 46096743 CTATCCCACCCTTCATATGAGAAATAGAGTATGTAAGCAAAAAAAGAGACTCCAGCAGACACTCCAAAATATCCTCCAAAAATAG
  #s LperrChr01  1183910 106 + 32922458 CCATCCCATACTCCATCCTATATTTGGTATATATGGAAGGAAAAATGGGCTCCAGTA----------TATATACCCATAAACTAG
  # store as chr, index block num for each chrs separatively
  # must has ref species
  # return alignLen checked maf and chr_maf_count
  # maf hash has keys: order, ref_region, aliLen, score, degree, and species align lines
  my ($file,$ref_chr_base) = @_;
  my %maf;
  my %cnt; # count for each chrs
  my $cnt_block;
  open MAF, $file or die "$!";
  $/ = "\n\n";
  while(<MAF>){
    chomp;
    last if ($_=~/##eof/); #modified on 28, Dec 2016
    my @box=split /\n/,$_;
    my @box_valid;
    my $has_aline = 0;
    my $has_sline = 0;
    foreach my $ctrl(@box){if($ctrl=~/^#/){}else{push @box_valid,$ctrl}; if($ctrl=~/^a score/){$has_aline = 1};if($ctrl=~/^s /){$has_sline = 1} }
    if($has_aline == 0 || $has_sline == 0){print STDERR "block not has_aline or has_sline at :@box\n ";next};
    if(scalar @box_valid <= 2){print STDERR "block has too few lines at @box\n"; next}
    my $block_score;
    my $ref_chr;
    my %block;
    my @order; 
    my $degree = 0;
    $cnt_block++;
    foreach my $line(@box){
      #exit if ($line=~/##eof/);
      #last if ($line=~/##eof/); #modified on 12, May 2013
      next if ($line=~/#/);
      if($line=~/^a score/){
          my @tmp = split/[\t ]+/,$line;
          my (undef,$block_score) = split/=/,$tmp[1];$block_score=int($block_score);
          if(!exists $block{"score"}){$block{"score"} = $block_score}else{die "dup score line: $line"}
        }elsif($line=~/^s /){
              my @tmp=split/ +/,$line;
              my ($name,$start,$len,$strand,$chr_len,$seq_ali) = ($tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6]); 
              if($name=~/^$ref_chr_base/){$ref_chr = $name}
              push @order, $name;
              $degree++;
              my $seq_ori=$seq_ali; $seq_ori=~s/-//g;
              my $end=$start+$len;
              my @dataline=($name,$start,$end,$len,$strand,$chr_len,$seq_ali,$seq_ori);
              die "length $len diff with seq_rmgap at block $cnt_block" if($len != length $seq_ori);
              my $alignLen = length $seq_ali;
              if(exists $block{$name}){die"duplicated species in maf block $cnt_block\n"}else{$block{$name}=\@dataline}
              if(exists $block{"alignLen"}){die "align length diff at block $cnt_block" if($block{"alignLen"} != $alignLen)}else{$block{"alignLen"} = $alignLen}
              if(exists $block{'align'}->{$name}){die"duplicated species in maf block $cnt_block\n"}else{$block{'align'}->{$name}=$seq_ali}
              if(exists $block{'all_region'}->{$name}){die"duplicated species in maf block $cnt_block\n"}else{$block{'all_region'}->{$name}="$name:$start-$end:$strand"}
         }else{print STDERR "\nomit line at group $cnt_block\n: unknow maf line type:$line\n"}
    }#foreach end here 
    if($degree == 0){die "degree eq 0 at group $cnt_block: $_"}else{$block{"degree"} = $degree}
    $block{"order"} = \@order;   
    if(defined $ref_chr ){
       $cnt{$ref_chr}++;
       my ($ref_start,$ref_end,$ref_strand)= ($block{$ref_chr}->[1],$block{$ref_chr}->[2]-1,$block{$ref_chr}->[4]);
       $block{"ref_region"} = "$ref_chr:$ref_start-$ref_end:$ref_strand";
       $maf{$ref_chr}->{$cnt{$ref_chr}} = \%block; 
    }else{print STDERR "ref chr not found, skip at block $cnt_block\n"}
    #last if($cnt_block >= 100);
    print STDERR "#" if($cnt_block % 10000 == 0);
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt_block maf blocks read in\n";
  return \%maf,\%cnt;
}


sub mafConcat_withcut(){
  #input mafs are strandered to +
  #sort->checkOverlap->concatenant are must, refDegap and fill_region_missing_with_* are optional
  #if empty, then warning
  my ($idx, $ref_cutregion,$flag_refdegap) = @_;
  my @ids = sort{$a<=>$b} keys %{$idx};  
  my $n_idx = scalar @ids;
  if($n_idx == 0){die "empty align in mafConcat"}
  #find the max spec_order of all maf blocks, sort by ref coords and check overlap
  my $max_nspec; 
  my $max_order;
  my $isOverlap = 0;
  my @ids_sorted = sort{
                          my ($chr1,$range1,$strand1) = split/:/,$idx->{$a}->{'ref_region'};                          
                          my ($s1,$e1) = split/-/,$range1;
                          die "s1($s1) > e1($e1) in segs sort" if($s1 > $e1);
                          my $n1 = scalar @{$idx->{$a}->{'order'}};
                          my $order1 = $idx->{$a}->{'order'};
                          my ($chr2,$range2,$strand2) = split/:/,$idx->{$b}->{'ref_region'};
                          my ($s2,$e2) = split/-/,$range2;
                          die "s2($s2) > e2($e2) in segs sort" if($s2 > $e2);
                          my $n2 = scalar @{$idx->{$b}->{'order'}};
                          my $order2 = $idx->{$b}->{'order'};
                          if($e1 <= $s2 || $e2<=$s1){}else{print STDERR "overlap segs found at ids $a,$b\n";$isOverlap = 1}
                          my $bigger = ($n1 >= $n2)?($n1):($n2);
                          my $bigger_order = ($n1 >= $n2)?($order1):($order2);
                          if(!defined $max_nspec){ $max_nspec = $bigger; $max_order = $bigger_order}else{ if($max_nspec<$bigger){$max_nspec = $bigger; $max_order = $bigger_order} }
                          if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                        } @ids ;
  #print "@ids @ids_sorted\n";
  #concatenant and report then append alignments, calculate missing regions and fill, refdegap if asked
  my %concat;
  my %concat_regions;
  if(scalar @ids_sorted == 1){ $max_order = $idx->{'1'}->{'order'}   } #for only one block
  if($isOverlap == 0){
    foreach my $id(@ids_sorted){
        #my $aliLen = &getAlignLen($idx->{$id}->{'align'});        
        my ($align_cut,$region_cut,$order_cut,$aliLen_cut) = &alignExtract_withcut($idx->{$id},$ref_cutregion);        
        #print Dumper $isOverlap,$align_cut,$region_cut,$order_cut,$aliLen_cut;
        #push @concat_regions,$region_cut;
        foreach my $spec(@$max_order){
           #print Dumper $max_order,$spec;
           if(exists $align_cut->{$spec}){
             my $aln = $align_cut->{$spec};
             $concat{$spec}.=$align_cut->{$spec}
           }else{$concat{$spec}.="*" x $aliLen_cut  }
           if(exists $region_cut->{$spec}){
              push @{$concat_regions{$spec}},$region_cut->{$spec}
           }else{push @{$concat_regions{$spec}},"fillgap_${aliLen_cut}"   }

        }
    }
  }else{die "overlap found" }

  #print Dumper %concat,$max_order,$aliLen,\%concat_regions;
  if($flag_refdegap eq "no" || $flag_refdegap eq "" ){my $aliLen = &getAlignLen(\%concat); return \%concat,$max_order,$aliLen,\%concat_regions}else{ my $concat_refdegap =  &aliRefDegap(\%concat,$flag_refdegap); my $aliLen = &getAlignLen($concat_refdegap); return $concat_refdegap, $max_order, $aliLen,\%concat_regions}

}





sub alignExtract_withcut(){
  #input ref real_coord, transform to ali_coord_all and extract
  my ($idx,$region) = @_;
  my ($chr,$range,$strand) = split/:/,$region; #ref region need to be cut
  my ($start,$end) = split/-/,$range;
  die "alignExract start($start) >= end($end)" if($start >= $end);
  my ($chr_ref,$range_ref,$strand_ref) = split/:/,$idx->{'ref_region'};  
  my ($start_ref,$end_ref) = split/-/,$range_ref;
  die "ref start($start_ref) >= end($end_ref)" if($start_ref >= $end_ref);
  die "chr($chr) ne chr_ref($chr_ref)" if($chr ne $chr_ref);

  my %align;

  if($start_ref >= $start && $end_ref <= $end){ #maf block included in cutregion, keep as is
    return $idx->{'align'},$idx->{'all_region'},$idx->{'order'},$idx->{'alignLen'}
    #return $idx->{'align'},$idx->{'ref_region'},$idx->{'order'},$idx->{'alignLen'}
  }elsif($start_ref < $start || $end_ref > $end){ #with overhangs
    my $region_cut;
    if($start_ref < $start && $end_ref <= $end){$region_cut = "$chr:$start-$end_ref:$strand" } #left overhangs
    if($start_ref >= $start && $end_ref > $end){$region_cut = "$chr:$start_ref-$end:$strand"} #right overhangs
    if($start_ref < $start && $end_ref > $end){$region_cut = "$chr:$start-$end:$strand"} #both overhangs
    #1 real_coord to align_coord
    my @align_ref = split//,$idx->{'align'}->{$chr};
    my ($start_align,$end_align);
    my $cnt =0;
    foreach my $i(0..$#align_ref){
      if($align_ref[$i] ne "-"){$cnt++}
      if($start_ref + $cnt == $start+1){$start_align = $i}
      if(!defined $start_align){$start_align=0}
      if($start_ref + $cnt == $end){$end_align = $i}
      if(!defined $end_align){$end_align=$#align_ref}
      
    }
    #2 extract all alignment 
    foreach my $spec(@{$idx->{'order'}}){
      $align{$spec}=substr($idx->{'align'}->{$spec},$start_align,$end_align-$start_align+1);
    }
    my $alilen = &getAlignLen(\%align);
    #3 align_coord to real coord for all species
    my %region_cut;
    foreach my $spec(@{$idx->{'order'}}){
        #print "spec: $spec\n";
        my ($chr_spec,$range_spec,$strand_spec) = split/:/,$idx->{'all_region'}->{$spec};
        my ($start_spec,$end_spec) = split/-/,$range_spec;
        my @align_spec = split//,$idx->{'align'}->{$spec};
        my $cnt_spec=0;
        my ($start_real,$end_real);
        foreach my $i(0..$#align_spec){
           if($align_spec[$i] ne "-" && $align_spec[$i] ne "*"){$cnt_spec++}
           if($i == $start_align){$start_real = $start_spec + $cnt_spec}#;print  "start_align($start_align),start_spec($start_spec),cnt($cnt_spec),start_real($start_real)\n"}
           if($i == $end_align){$end_real = $start_spec + $cnt_spec}#;print "end_align($end_align),end_spec($end_spec),cnt($cnt_spec),end_real($end_real)\n"}           
        }
        $region_cut{$spec}="$spec:$start_real-$end_real:$strand_spec";
    }

    return \%align,\%region_cut,$idx->{'order'},$alilen;

  }else{
         die "impossible: region: $start,$end | in_maf: $start_ref,$end_ref"
        }

}#sub end

sub getAlignLen(){
   my $idx = shift;
   my $len;
   foreach my $spec(keys %{$idx}){
     if(!defined $len){ $len = length $idx->{$spec}}else{
       if($len != length $idx->{$spec}){print Dumper $idx; die "align length not all eq"} #with check
     }
   }
   return $len;
}


sub aliRefDegap(){
  #malign in, refDegap_malign out
  my ($align,$ref) = @_; #align hash
  my $len = &getAlignLen($align);
  die "refId $ref not exists when refDegap" if(!exists $align->{$ref});

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


  return \%align_refDegap;
  #return (\%align_refDegap,\%align_deNum);

}

sub align2SS(){
 my ($align,$order,$len) = @_;

 my %align_split;
 my %ss;
 foreach my $spec(@$order){
   my @box = split//,$align->{$spec};
   $align_split{$spec} = \@box;
 }

 for my $i(0..$len-1){
   my $column;
   foreach my $spec(@$order){ $column.=uc($align_split{$spec}->[$i]) } #upcase
   $ss{$column}->{'cnt'}++;
   push @{$ss{$column}->{'order'}}, $i;
 }

 #sort %ss
 my @id_sorted = sort {
                        my $order1 = $ss{$a}->{'order'}->[0];
                        my $order2 = $ss{$b}->{'order'}->[0];
                        if($order1 > $order2){return 1}elsif($order1 < $order2){return -1}else{return 0}
                      } keys %ss;
 my @ss_sorted;
 foreach my $id(@id_sorted){ push @ss_sorted,"$id\t$ss{$id}->{cnt}"}

 return \@ss_sorted;

}


sub faFormat(){
  my ($seq,$num)=@_;
  my @seq=split//,$seq;
  my ($string,$cnt);
  foreach(@seq){
   $cnt++;
   if($cnt%$num==0){$string.=$_."\n"}else{$string.=$_}
  }
  if($cnt%$num!=0){return $string."\n"}else{return $string}
}




