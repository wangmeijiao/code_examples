use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#readin mafExtract maf blocks, output concatenanted mfa align file, need to has ref chr, fill missing spec with *
#no need to fill gaps between blocks because there are extract regions in the mfa headline

my($maf_file,$width,$ref,$fill);
GetOptions("maf=s",\$maf_file,"width:i",\$width,"ref=s",\$ref);
$width||=50;


my $ref_base = $ref;
$ref_base=~s/\d{1,2}$//;

my ($maf,$maf_cnt)=&mafRead_refchr($maf_file,$ref_base);
print STDERR "read maf done with refbase $ref_base\n";
#print Dumper $maf_cnt,$maf; exit;



my ($maf_concated,$max_order,$aliLen,$concat_regions) = &mafConcat($maf->{$ref},$ref);
#print Dumper $maf_concated,$max_order,$aliLen;

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


sub mafConcat(){
  #sort->checkOverlap->concatenant are must, refDegap and fill_missing_with_* are optional
  #if empty, then warning
  my ($idx, $flag_refdegap) = @_;
  my @ids = sort{$a<=>$b} keys %{$idx};  
  my $n_idx = scalar @ids;
  if($n_idx == 0){print STDERR "empty align in mafConcat"}
  #if($n_idx == 1){ die "only one maf block found" } one block also applicable?
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
  if($isOverlap == 0){
    foreach my $id(@ids_sorted){
        my $aliLen = &getAlignLen($idx->{$id}->{'align'});
        #push @concat_regions,$idx->{$id}->{'ref_region'};
        foreach my $spec(@$max_order){
           if(exists $idx->{$id}->{'align'}->{$spec}){
             my $aln = $idx->{$id}->{'align'}->{$spec};
             $concat{$spec}.=$idx->{$id}->{'align'}->{$spec}
           }else{$concat{$spec}.="*" x $aliLen  }
           if(exists $idx->{$id}->{'all_region'}->{$spec}){
              push @{$concat_regions{$spec}},$idx->{$id}->{'all_region'}->{$spec}
           }else{push @{$concat_regions{$spec}},"fillgap_$aliLen"   }
        }
    }
  }else{ die "overlap found" }

  if($flag_refdegap eq "no" || $flag_refdegap eq "" ){my $aliLen = &getAlignLen(\%concat); return \%concat,$max_order,$aliLen,\%concat_regions}else{ my $concat_refdegap =  &aliRefDegap(\%concat,$flag_refdegap); my $aliLen = &getAlignLen($concat_refdegap); return $concat_refdegap, $max_order, $aliLen,\%concat_regions}

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




