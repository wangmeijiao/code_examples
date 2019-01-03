use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

# read in gff and maf,  then  extract cds codons and choose 4d sites (code3), output as concatenated mfa file

#4d codons:
#A       GCA     GCC     GCG     GCT     
#G       GGA     GGC     GGG     GGT     
#P       CCA     CCC     CCG     CCT     
#T       ACA     ACC     ACG     ACT     
#V       GTA     GTC     GTG     GTT     
#
##4d-2d
#L       CTA     CTC     CTG     CTT ,    TTA     TTG     
#R       AGA     AGG     CGA     CGC ,    CGG     CGT     
#S       AGC     AGT     TCA     TCC ,    TCG     TCT     
#


my($type,$gff_file,$maf_file,$width);
GetOptions("type=s",\$type,"gff=s",\$gff_file,"maf=s",\$maf_file,"width:i",\$width);
$width||=50;
#type = gene, exon, cds, codon, 4dsite

my ($gff_idx,$gene_idx,$mrna_idx) = &readGff($gff_file);
print STDERR "read gff done\n";
#print Dumper $gff_idx,$gene_idx,$mrna_idx; exit;



my ($maf,$maf_cnt)=&mafRead_refchr($maf_file,"japoChr");
print STDERR "read maf done\n";
#print Dumper $maf_cnt,$maf; exit;


#test for mafExtract_seg;
die "not exists chr japoChr03 in your maf file\n" if(!exists $maf->{"japoChr03"});
my $result = &mafExtract_seg($maf->{"japoChr03"},"japoChr03:8500-9000:+");
#print Dumper $result;

#test for mafConcat;
if(!%$result){die "result of mafExtract is empty"}else{
  my ($align_concat,$max_order,$align_len,$concat_regions) = &mafConcat($result,"japoChr03","no"); #sort->checkOverlap->concatenant, refDegap and fill_missing_with_* are optional
  #print Dumper $align_concat,$max_order,$align_len;
  #output mfa
  #foreach my $spec(@$max_order){print ">$spec $align_len\n$align_concat->{$spec}\n"}
  my $cnt;
  foreach my $spec(@$max_order){
     $cnt++;
     if($cnt == 1){
       print ">$spec $align_len ".join("|",@$concat_regions)."\n$align_concat->{$spec}\n"
     }else{print ">$spec $align_len\n$align_concat->{$spec}\n"}
  }



  
  #output SS
  my $ss = &align2SS($align_concat,$max_order,$align_len);
  foreach my $ctrl(0..$#$ss){ print $ctrl,"\t",$ss->[$ctrl],"\n"}

}


=pod

#start to extract and output
#$outdir||=".";
#if(-d $outdir){}else{mkdir $outdir}
#$prefix||="out";
print STDERR "start to extract from maf: $maf_file according to gff $gff_file :\n";

foreach my$chr(sort keys %{$gff_idx}){
 print STDERR "$chr ..";
 foreach my $geneId(sort keys %{$gff->{$chr}}){

    #1,output gene region
    my ($chr,$start,$end,$ID,$strand)=split/\t/,$gene->{$geneId};
    if($type eq "gene"){
     #my $seq = &faFormat(&faExtract($fa_index,"$chr:$start-$end:$strand"),$width);
     #print GENE ">$chr|$start|$end|$ID|$strand|gene\n$seq";

    }
    my $max_len_mrna;
    foreach my $mrna(keys %{$gff->{$chr}->{$gene}}){
      my @box=split/\t/,$mrna{$mrna};
      if($box[6] eq "longest_AS"){$max_len_mrna=$mrna}
    }
    #mrna (the longest one, which all of the following caculation are based on)
    my $index=$gff{$chr}->{$gene}->{$max_len_mrna};
    if(!defined $index || $index eq ""){print STDERR "the longest mrna of $gene is lost:$max_len_mrna\n"; next}
    my ($chr_mrna,$start_mrna,$end_mrna,$id_mrna,$strand_mrna,$pID,$mrna_note)=split/\t/,$mrna{$max_len_mrna};
    #$id_mrna=~s/\.v2\.1$//;

    #2,output mrna region longest
    #if($mrna){
    # my $seq = &faFormat(&faExtract($fa_index,"$chr_mrna:$start_mrna-$end_mrna:$strand_mrna"),50);
    # print MRNA ">$alias_mrna\n$seq";
    #}

    #3,output mrna transcritps (exons with utrs)
    if(exists $index->{"exon"}){
      my @exon=@{&segSort($index->{"exon"})};
      my @exon_seq;
      my @exon_coords;
      foreach(@exon){
        my($exon_start,$exon_end)=split/-/,$_;
        #push @exon_coords,"$exon_start-$exon_end";
        my $exon_seq=&faExtract($fa_index,"$chr_mrna:$exon_start-$exon_end:$strand_mrna");
        if($strand_mrna eq "+"){push @exon_seq, $exon_seq;push @exon_coords,"$exon_start-$exon_end"}elsif($strand_mrna eq "-"){unshift @exon_seq, $exon_seq;unshift @exon_coords,"$exon_start-$exon_end"}
        #($strand_mrna eq "+")?(push @exon_seq, $exon_seq):(unshift @exon_seq, $exon_seq);
      }
      if($mrna_flag){
         my $seq = &faFormat(join("",@exon_seq),$width);
         my $coords = join":",@exon_coords;
         print MRNA ">$chr_mrna|$start_mrna|$end_mrna|$coords|$id_mrna|$strand_mrna|mRNA\n$seq";
      }
    }else{print STDERR "no exon records for $id_mrna, mrna transcripts output nothing\n"}


    #4,output CDS and/or pep
    if(exists $index->{"CDS"}){
       my @CDS=@{&segSort($index->{"CDS"})};
       #my $i;
       my @CDS_seq;
       my @CDS_coords;
       my ($cds_left,$cds_right);
       #my @PEP_seq;
       foreach(@CDS){
         #$i++;
         my($cds_start,$cds_end,$frameshift)=split/-/,$_;
         #if(!defined $cds_left){$cds_left = $cds_start}else{if($cds_left > $cds_start){$cds_left = $cds_start}}
         #if(!defined $cds_right){$cds_right = $cds_end}else{if($cds_right < $cds_end){$cds_right = $cds_end}}
         my $cds_seq=&faExtract($fa_index,"$chr_mrna:$cds_start-$cds_end:$strand_mrna");
         ($strand_mrna eq "+")?(push @CDS_seq, $cds_seq):(unshift @CDS_seq, $cds_seq);
         push @CDS_coords,"$cds_start-$cds_end"; #frameshift is not so important as we may thought
       }

       if($type eq "cds"){
         my $seq=&faFormat(join("",@CDS_seq),$width);
         my $coords = join":",@CDS_coords;
         #my $flag_strand;
         #if($strand_mrna eq "+"){$flag_strand = 1}else{$flag_strand = -1}
         print CDS ">$chr_mrna|$cds_left|$cds_right|$coords|$id_mrna|$strand_mrna|CDS\n$seq";
       }
    }else{die "no CDS at $ID\n"}

 }#foreach gene
  print STDERR " done\n";
}#foreach chr end

=cut

###subs###

sub readGff(){
  my $gff_f = shift;
  my %gff;#main structure table:  chr->gene->mRNAs->(elements:CDS,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)
  my %gene;#gene info
  my %mrna;#mrna info
  
  #1,read&fill without order consideration
  open GFF, $gff_f or die "$!";
  while(<GFF>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
    die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
    next if($chr eq "chrUn" || $chr eq "chrSy"); #get rid of other chrs
    next if($ftype eq "chromosome" || $ftype eq "Chromosome" || $ftype eq "contig");
    die "start > end :$start > $end at $_" if($start > $end);
    #$chr=~s/^Chr/chr/;
    #($start,$end)=($start<$end)?($start,$end):($end,$start);
    my ($ID,$pID,$note,$anno);
    if($ftype eq "gene"){
      $attr=~/ID=([^;]+)/;
      $ID=$1;
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand";
        #$gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$note\t$alias";
      }else{die "dup gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+);Note=([^;]+)/;
        ($ID,$pID,$note)=($1,$2,$3);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID\t$note";
          #if(exists $gene{$pID}){ $gene{$pID}.="\t$note"}else{die "mrna $ID come before gene $pID at $_"}
        }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"|| $ftype eq 'tss' || $ftype eq 'tts'){
           $attr=~/Parent=([^;]+);Note=([^;]+)/;
           ($pID,$note)=($1,$2);
           my $geneID;
           if(exists $mrna{$pID}){
            my @box=split/\t/,$mrna{$pID};
            if(defined $box[5]){$geneID=$box[5]}else{die"pID empty at $pID\n"}
           }else{die "elements comes first, I can't assign it to gene\n"}
           if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
              my @temp;
              if($ftype eq "CDS"){push @temp,"$start-$end-$fshift"}else{push @temp,"$start-$end"};
              $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
           }else{ if($ftype eq "CDS"){push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$fshift"}else{push  @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end"} }
         }else{die "unknow feature type $ftype\n"}
  }#while end
  
  close GFF;
  #print STDERR "read gff done\n";  
  return (\%gff,\%gene,\%mrna);
}


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
         }else{print STDERR "\nomit line at group $cnt_block\n: unknow maf line type:$line\n"}
    }#foreach end here 
    if($degree == 0){die "degree eq 0 at group $cnt_block: $_"}else{$block{"degree"} = $degree}
    $block{"order"} = \@order;   
    if(defined $ref_chr ){
       $cnt{$ref_chr}++;
       my ($ref_start,$ref_end,$ref_strand)= ($block{$ref_chr}->[1],$block{$ref_chr}->[2],$block{$ref_chr}->[4]);
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


sub faRead(){
  my $file = shift;
  my %fa;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp = split/[\t ]+/,$id;
    $id = $temp[0];
    my $seq = join "", @box;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return \%fa;
}



sub segSort(){
 my @segs=@{shift @_};
 my @sorted=sort{
      my ($s1,$e1,$s2,$e2);#$s1<=$e1 && $s2<=$e2 and deOverlaped beforehand by default
      ($s1,$e1)=split/-/,$a;
      ($s2,$e2)=split/-/,$b;
      if($e2<=$s1){1}elsif($e1<=$s2){-1}else{print STDERR "unknown circumstance at $s1,$e1:$s2,$e2\n"}
    } @segs;
# print "@CDS\n@sorted\n";
 return \@sorted;
}


sub mafExtract_seg(){
    # with ref species coordinates
    # one segment per time
    # don't be duplicated in maf
    # may be lost or shorter than segment
    # pay attention to strand
    # return align with trim, add missing with * and refdegap
    my ($index,$coord)=@_;
    my ($chr,$s,$e,$strand)=split/:|-/,$coord; 
    die "s($s) >= e($e) when query maf" if($s >= $e);
    my %results;
    my @hits; #block ids and block ref_regions, in case of multihits.
    foreach my $id(sort {$a<=>$b} keys %{$index}){ #loop over all ids of a chromosome maf
       my ($ref_chr,$ref_start,$ref_end,$ref_strand) = split/:|-/,$index->{$id}->{'ref_region'};
       die "$chr ne $ref_chr at $index->{$id}->{'ref_region'}" if($chr ne $ref_chr);
       die "ref_start($ref_start ) >= ref_end ($ref_end)" if($ref_start >= $ref_end);

       #match feature_strand and maf_strand
       if($strand eq "+" && $ref_strand eq "+"){ 
         if($s >= $ref_start && $e <= $ref_end){ # feature is included in the ref block, extract alignment from the ref block
           my ($ali_extract,$order) = &alignExtract($index->{$id},"$ref_chr:$s-$e:$ref_strand");
           $results{$id}->{'region'}="$ref_chr:$s-$e:$strand";
           $results{$id}->{'region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
           $results{$id}->{'align'}=$ali_extract;
           $results{$id}->{'order'}=$order;
         }elsif($ref_start >= $s && $ref_end <= $e){ # feature covers the ref block, report the whole ref block
             my ($ali_extract,$order) = &alignExtract($index->{$id},"$ref_chr:$ref_start-$ref_end:$ref_strand");
             $results{$id}->{'region'}="$ref_chr:$ref_start-$ref_end:$strand";
             $results{$id}->{'region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
             $results{$id}->{'align'}=$ali_extract;
             $results{$id}->{'order'}=$order;
          }elsif($s >= $ref_start && $s  <= $ref_end && $e > $ref_end ){#overlap with feature left part, report the right part
             my ($ali_extract,$order) = &alignExtract($index->{$id},"$ref_chr:$s-$ref_end:$ref_strand");
             $results{$id}->{'region'}="$ref_chr:$s-$ref_end:$strand";
             $results{$id}->{'region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
             $results{$id}->{'align'}=$ali_extract;
             $results{$id}->{'order'}=$order;
            }elsif($e >= $ref_start && $e  <= $ref_end && $s < $ref_start){#overlap with feature right part, report the left part
               my ($ali_extract,$order) = &alignExtract($index->{$id},"$ref_chr:$ref_start-$e:$ref_strand");
               $results{$id}->{'region'}="$ref_chr:$ref_start-$e:$strand";
               $results{$id}->{'region_ori'}="$ref_chr:$ref_start-$ref_end:$ref_strand";
               $results{$id}->{'align'}=$ali_extract;
               $results{$id}->{'order'}=$order;
             }

      }elsif($strand eq "+" && $ref_strand eq "-"){


         }elsif($strand eq "-" && $ref_strand eq "+"){


              }elsif($strand eq "-" && $ref_strand eq "-"){


               }
    }#foreach id end

    return \%results;

}



sub alignExtract(){
  my ($idx,$region) = @_;
  my ($ref,$start,$end,$strand) = split/:|-/,$region;
  die "alignExract start($start) >= end($end)" if($start >= $end);
  my %align;
  #real_coord to align_coord
  my @align_ref = split//,$idx->{$ref}->[6];
  my ($start_ref,$end_ref) = ($idx->{$ref}->[1],$idx->{$ref}->[2]);
  my ($start_align,$end_align); 
  my $cnt =0;
  foreach my $i(0..$#align_ref){
    if($align_ref[$i] ne "-"){$cnt++}
    if($start_ref + $cnt == $start){$start_align = $i} 
    if($start_ref + $cnt == $end){$end_align = $i} 
  }

  #extract all alignment 
  foreach my $spec(@{$idx->{'order'}}){
    $align{$spec}=substr($idx->{$spec}->[6],$start_align,$end_align-$start_align+1);
  }

  return \%align,$idx->{'order'};
}

sub mafConcat(){
  #sort->checkOverlap->concatenant are must, refDegap and fill_missing_with_* are optional
  #if empty, then warning
  my ($idx, $flag_refdegap, $flag_fill_missing) = @_;
  my @ids = sort{$a<=>$b} keys %{$idx};  
  my $n_idx = scalar @ids;
  if($n_idx == 0){print STDERR "empty align in mafConcat"}
  if($n_idx == 1){ return 0  }
  #find the max spec_order of all maf blocks, sort by ref coords and check overlap
   my $max_nspec; 
   my $max_order;
   my $isOverlap = 0;
   my @ids_sorted = sort{
                          my ($chr1,$s1,$e1,$strand1) = split/:|-/,$idx->{$a}->{'region'};                          
                          die "s1 > e1 in segs sort" if($s1 > $e1);
                          my $n1 = scalar @{$idx->{$a}->{'order'}};
                          my $order1 = $idx->{$a}->{'order'};
                          my ($chr2,$s2,$e2,$strand2) = split/:|-/,$idx->{$b}->{'region'};
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
  my @concat_regions;
  if($isOverlap == 0){
    foreach my $id(@ids_sorted){
        my $aliLen = &getAlignLen($idx->{$id}->{'align'});
        push @concat_regions, $idx->{$id}->{'region'};
        foreach my $spec(@$max_order){
           if(exists $idx->{$id}->{'align'}->{$spec}){
             my $aln = $idx->{$id}->{'align'}->{$spec};
             $concat{$spec}.=$idx->{$id}->{'align'}->{$spec}
           }else{$concat{$spec}.="*" x $aliLen  }
        }
    }
  }else{return 0  }

  if($flag_refdegap eq "no" || $flag_refdegap eq "" ){my $aliLen = &getAlignLen(\%concat); return \%concat,$max_order,$aliLen,\@concat_regions;}else{ my $concat_refdegap =  &aliRefDegap(\%concat,$flag_refdegap); my $aliLen = &getAlignLen($concat_refdegap); return $concat_refdegap, $max_order, $aliLen,\@concat_regions}

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


sub faExtract(){
    my ($index,$coord)=@_;
    my ($chr,$s,$e,$strand)=split/:|-/,$coord;
    die "not exists $chr in your fa file\n"if(!exists $index->{$chr});
    my $seq_extract=($strand eq "+")?(substr($index->{$chr},$s-1,$e-$s+1)):(&RC(substr($index->{$chr},$s-1,$e-$s+1)));
    #return uc($seq_extract);
    return $seq_extract;
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

sub RC(){
  my $seq=shift;
  $seq=reverse $seq;
  $seq=~tr/atcgATCG/tagcTAGC/;
  return $seq;
}

sub translate_CDS_simple(){
#standard code table:http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
#    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#  Starts = ---M---------------M---------------M----------------------------
#  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
   my $orf=shift;
   #my $frameshift=shift; #0,1,2: delete 0,1,2 nucleiotide before translating
   print STDERR "Warning! orf undefined\n" if(!defined $orf);
   #if($frameshift eq "." || $frameshift eq "" || !defined $frameshift){$frameshift=0}
   $orf=uc($orf);
   my $len = length $orf;
   if($len % 3 != 0){my $num = $len % 3; $orf=~s/.{$num}$//;    print STDERR "\nWarning! ##not3x at orf $orf, chomp tails to fit 3x\n"}
   if($len < 3){die "err orf <3 :$orf"}
   #$orf=substr($orf,$frameshift);
   my %code=(
             'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S',
             'TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
             'TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*',
             'TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W',
             'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
             'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
             'CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q',
             'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
             'ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M',
             'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
             'AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
             'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R',
             'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
             'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
             'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E',
             'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G'
            );
   my $aa="";
   my $codon;
   for(my$i=1;$i<=$len-2;$i+=3){
     $codon=substr($orf,$i-1,3);
     if($codon=~/n+/i){$aa.="X";next} #X for gap
     if($i!=(length $orf)-2 && $codon=~/TAA|TAG|TGA/){$aa.="*";print STDERR "\nWarning! ##earlier_stop at orf $orf\n";next}
     if(exists $code{$codon}){$aa.=$code{$codon}}else{$aa.="X";print STDERR "\nWarning! ##not_in_codon_tab $codon at orf $orf\n";next} #X for unknown
     #if(exists $code{$codon}){$aa.=$code{$codon}}else{return "##wrong_codon"}
   }
   return $aa;
}

sub checkPEP(){
   my $pep=shift;
   my @comments;
   if($pep!~/^M/){push @comments,"##start codon not M"}
   if($pep!~/\*$/){push @comments,"##stop codon not exists"}
   if(length $pep <= 10){push @comments,"##short than 50 aa"}
   chop($pep);
   if($pep=~/\*/){push @comments,"##stop codon in pep"}
   if(scalar @comments == 0){ return "OK"}else{return join",",@comments}

}



