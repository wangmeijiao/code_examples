use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#perl mafExtract_seg.pl -maf oryza_5way_kent.chr03.filter.degree3.blocklen10.score100.sort.order.strander.maf -width 2000 -seg "japoChr03:122040-123040:+" > test.mafExtract_seg.concat.mfa

my($seg,$maf_file,$width);
GetOptions("seg=s",\$seg,"maf=s",\$maf_file,"width:i",\$width); #seg japoChr03:135895-136900:+
$width||=50;

my $ref;
($ref,undef,undef) = split/:/,$seg;
my $ref_base = $ref;
$ref_base=~s/\d{1,2}$//;

my ($maf,$maf_cnt)=&mafRead_refchr($maf_file,$ref_base);
print STDERR "read maf done\n";
#print Dumper $maf_cnt,$maf; exit;

die "not exists ref chr $ref in your maf file\n" if(!exists $maf->{$ref});
my $result = &mafExtract_seg($maf->{$ref},$seg);
#print Dumper $result;#exit;

if(!%$result){die "result of mafExtract is empty"}else{
  my ($align_concat,$max_order,$align_len,$concat_regions) = &mafConcat($result,"no","no"); 
  #my ($align_concat,$max_order,$align_len,$concat_regions) = &mafConcat($result,$ref,"no"); 
  #print Dumper $align_concat,$max_order,$align_len;
  #output mfa
 
  #print Dumper $align_concat,$max_order,$align_len,$concat_regions; 
  my $cnt;
  foreach my $spec(@$max_order){
     $cnt++;
     my $faStr_format = &faFormat($align_concat->{$spec},$width);
     #if($cnt == 1){
       print ">$spec $align_len ".join("|",@{$concat_regions->{$spec}})."\n$faStr_format"
     #}else{print ">$spec $align_len\n$align_concat->{$spec}\n"}
  }
  
  #output SS
  #my $ss = &align2SS($align_concat,$max_order,$align_len);
  #foreach my $ctrl(0..$#$ss){ print $ctrl,"\t",$ss->[$ctrl],"\n"}

}

###subs###

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
       my ($ref_start,$ref_end,$ref_strand)= ($block{$ref_chr}->[1],$block{$ref_chr}->[2],$block{$ref_chr}->[4]);
       #my ($ref_start,$ref_end,$ref_strand)= ($block{$ref_chr}->[1],$block{$ref_chr}->[2]-1,$block{$ref_chr}->[4]);
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



sub mafExtract_seg(){
    # with ref species coordinates
    # one segment per time
    # don't be duplicated in maf
    # may be lost or shorter than segment
    # pay attention to strand
    # return align with trim, add missing with * and refdegap or not
    my ($index,$coord)=@_;
    my ($chr,$region,$strand)=split/:/,$coord; 
    my ($s,$e) = split/-/,$region;
    die "s($s) >= e($e) when query maf" if($s >= $e);
    my %results;
    my @hits; #block ids and block ref_regions, in case of multihits.
    foreach my $id(sort {$a<=>$b} keys %{$index}){ #loop over all ids of a chromosome maf
       my ($ref_chr,$ref_range,$ref_strand) = split/:/,$index->{$id}->{'ref_region'};
       my ($ref_start,$ref_end) = split/-/,$ref_range; 
       die "$chr ne $ref_chr at $index->{$id}->{'ref_region'}" if($chr ne $ref_chr);
       die "ref_start($ref_start ) >= ref_end ($ref_end)" if($ref_start > $ref_end);
       #die "ref_start($ref_start ) >= ref_end ($ref_end)" if($ref_start >= $ref_end);

       #match feature_strand and maf_strand
       if($strand eq "+" && $ref_strand eq "+"){ 
         if($s >= $ref_start && $e <= $ref_end){ # feature is included in the ref block, extract alignment from the ref block
           #print STDERR "ref_block($ref_start-$ref_end) includes feature($s-$e) :use $ref_chr:$s-$e:$ref_strand\n";
           my ($ali_extract,$order,$region_cut,$alilen) = &alignExtract($index->{$id},"$ref_chr:$s-$e:$ref_strand");
           $results{$id}->{'ref_region'}="$ref_chr:$s-$e:$strand";
           $results{$id}->{'ref_region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
           $results{$id}->{'align'}=$ali_extract;
           $results{$id}->{'order'}=$order;
           $results{$id}->{'region_cut'}=$region_cut;
           $results{$id}->{'alilen'}=$alilen;
         }elsif($ref_start >= $s && $ref_end <= $e){ # feature covers the ref block, report the whole ref block
             #print STDERR "ref_block($ref_start-$ref_end) covered by feature($s-$e) :use $ref_chr:$ref_start-$ref_end:$ref_strand\n";
             my ($ali_extract,$order,$region_cut,$alilen) = &alignExtract($index->{$id},"$ref_chr:$ref_start-$ref_end:$ref_strand");
             $results{$id}->{'ref_region'}="$ref_chr:$ref_start-$ref_end:$strand";
             $results{$id}->{'ref_region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
             $results{$id}->{'align'}=$ali_extract;
             $results{$id}->{'order'}=$order;
             $results{$id}->{'region_cut'}=$region_cut;
             $results{$id}->{'alilen'}=$alilen;
          }elsif($s >= $ref_start && $s  <= $ref_end && $e > $ref_end ){#overlap with feature left part, report the right part
             #print STDERR "ref_block($ref_start-$ref_end) overlap with feature($s-$e) left part:use $ref_chr:$s-$ref_end:$ref_strand\n";
             my ($ali_extract,$order,$region_cut,$alilen) = &alignExtract($index->{$id},"$ref_chr:$s-$ref_end:$ref_strand");
             $results{$id}->{'ref_region'}="$ref_chr:$s-$ref_end:$strand";
             $results{$id}->{'ref_region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
             $results{$id}->{'align'}=$ali_extract;
             $results{$id}->{'order'}=$order;
             $results{$id}->{'region_cut'}=$region_cut;
             $results{$id}->{'alilen'}=$alilen;
            }elsif($e >= $ref_start && $e  <= $ref_end && $s < $ref_start){#overlap with feature right part, report the left part
               #print STDERR "ref_block($ref_start-$ref_end) overlap with feature($s-$e) right part:use $ref_chr:$ref_start-$e:$ref_strand\n";
               my ($ali_extract,$order,$region_cut,$alilen) = &alignExtract($index->{$id},"$ref_chr:$ref_start-$e:$ref_strand");
               $results{$id}->{'ref_region'}="$ref_chr:$ref_start-$e:$strand";
               $results{$id}->{'ref_region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
               $results{$id}->{'align'}=$ali_extract;
               $results{$id}->{'order'}=$order;
               $results{$id}->{'region_cut'}=$region_cut;
               $results{$id}->{'alilen'}=$alilen;
             }

      }elsif($strand eq "+" && $ref_strand eq "-"){
          die "strand diffs:$strand, $ref_strand at $ref_chr,$ref_range,$ref_strand";

         }elsif($strand eq "-" && $ref_strand eq "+"){
                 die "strand diffs:$strand, $ref_strand at $ref_chr,$ref_range,$ref_strand";

              }elsif($strand eq "-" && $ref_strand eq "-"){
                    die "strand are -:$strand, $ref_strand at $ref_chr,$ref_range,$ref_strand";

               }
    }#foreach id end

    return \%results;

}



sub alignExtract(){
  my ($idx,$region) = @_;
  my ($ref,$range,$strand) = split/:/,$region;
  my ($start,$end) = split/-/,$range;
  #print STDERR "extract align for $region ";
  die "alignExract start($start) > end($end)" if($start > $end);
  my %align;
  #real_coord to align_coord
  my @align_ref = split//,$idx->{$ref}->[6];
  my ($start_ref,$end_ref) = ($idx->{$ref}->[1],$idx->{$ref}->[2]);
  #print STDERR " with ref_block $start_ref,$end_ref ";
  my ($start_align,$end_align); 
  my $cnt =0;
  foreach my $i(0..$#align_ref){
    if($align_ref[$i] ne "-"){$cnt++}
    if($start_ref + $cnt == $start){$start_align = $i}#;print STDERR "start hit: cnt $cnt "} 
    if($start_ref + $cnt == $end){$end_align = $i}#; print STDERR "end hit: cnt $cnt "}
  }
  if($start == $start_ref){$start_align = 0}
  if($end == $end_ref){$end_align = $#align_ref}
  if($start == $end){$start_align = 0; $end_align = 0}
  #extract all alignment 
  #print STDERR " with align coord: $start_align,$end_align\n";
  foreach my $spec(@{$idx->{'order'}}){
    $align{$spec}=substr($idx->{$spec}->[6],$start_align,$end_align-$start_align+1);
  }
  my $alilen = &getAlignLen(\%align);
  #align_coord to real coord for all species
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

  return \%align,$idx->{'order'},\%region_cut,$alilen;
}

sub mafConcat(){
  #sort->checkOverlap->concatenant are must, refDegap and fill_missing_with_* are optional
  #if empty, then warning
  my ($idx, $flag_refdegap, $flag_fill_missing) = @_;
  my @ids = sort{$a<=>$b} keys %{$idx};  
  my $n_idx = scalar @ids;
  if($n_idx == 0){print STDERR "empty align in mafConcat"}
  #find the max spec_order of all maf blocks, sort by ref coords and check overlap
   my $max_nspec; 
   my $max_order;
   my $isOverlap = 0;
   my @ids_sorted = sort{
                          my ($chr1,$s1,$e1,$strand1) = split/:|-/,$idx->{$a}->{'ref_region'};                          
                          die "s1 > e1 in segs sort" if($s1 > $e1);
                          my $n1 = scalar @{$idx->{$a}->{'order'}};
                          my $order1 = $idx->{$a}->{'order'};
                          my ($chr2,$s2,$e2,$strand2) = split/:|-/,$idx->{$b}->{'ref_region'};
                          die "s2 > e2 in segs sort" if($s2 > $e2);
                          my $n2 = scalar @{$idx->{$b}->{'order'}};
                          my $order2 = $idx->{$b}->{'order'};
                          if($e1 <= $s2 || $e2<=$s1){}else{print STDERR "overlap segs found at ids $a,$b\n";$isOverlap = 1}
                          my $bigger = ($n1 >= $n2)?($n1):($n2);
                          my $bigger_order = ($n1 >= $n2)?($order1):($order2);
                          if(!defined $max_nspec){ $max_nspec = $bigger; $max_order = $bigger_order}else{ if($max_nspec<$bigger){$max_nspec = $bigger; $max_order = $bigger_order} }
                          if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
                        } @ids ;
  #concatenant and report then append alignments, calculate missing regions and fill, refdegap if asked
  my %concat;
  my %concat_regions;
  if(scalar @ids_sorted == 1){ $max_order = $idx->{$ids_sorted[0]}->{'order'}   } #for only one block
  if($isOverlap == 0){
    foreach my $id(@ids_sorted){
        my $aliLen = $idx->{$id}->{'alilen'};
        foreach my $spec(@$max_order){
           if(exists $idx->{$id}->{'align'}->{$spec}){
             my $aln = $idx->{$id}->{'align'}->{$spec};
             $concat{$spec}.=$idx->{$id}->{'align'}->{$spec}
           }else{$concat{$spec}.="*" x $aliLen  }
           if(exists $idx->{$id}->{'region_cut'}->{$spec}){
              push @{$concat_regions{$spec}},$idx->{$id}->{'region_cut'}->{$spec}
           }else{push @{$concat_regions{$spec}},"fillgap_${aliLen}" }
        }
    }
  }else{die "overlap found when concat maf"  }

  if($flag_refdegap eq "no" || $flag_refdegap eq "no" ){my $aliLen = &getAlignLen(\%concat); return \%concat,$max_order,$aliLen,\%concat_regions;}else{ my $concat_refdegap =  &aliRefDegap(\%concat,$flag_refdegap); my $aliLen = &getAlignLen($concat_refdegap); return $concat_refdegap, $max_order, $aliLen,\%concat_regions}

}


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


