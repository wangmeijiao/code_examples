use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

##map ref annotation feature to refDegapped maf
##extract any coordination-based subsets, like msa_view in phast package
## 29. Nov 2016




my ($maf_file,$gff_file,$ref,$outfmt);
GetOptions("maf=s",\$maf_file,"gff=s",\$gff_file,"ref=s",\$ref,"outfmt=s",\$outfmt);



# 1.readin maf(check ref is degapped) and gff
my $maf = &mafRead_plus_withcheck($maf_file,$ref); 
my @IDs = sort {$a <=> $b} keys %$maf;
my ($gff_idx,$gene_idx,$mrna_idx)=&gffRead($gff_file);


#print Dumper $gff_idx,$gene_idx,$mrna_idx,$maf;


if(defined $outfmt && $outfmt eq "maf"){print "##maf version 12\n"}
foreach my $chr(sort keys %{$gff_idx}){
  next if($chr ne $ref);
  foreach my $gene(sort keys %{$gff_idx->{$chr}}){
    foreach my $mrna(sort keys %{$gff_idx->{$chr}->{$gene}}){
        my $idx = $gff_idx->{$chr}->{$gene}->{$mrna};
        die "not exists exon in $chr->$gene->$mrna" if(!exists $idx->{'exon'});
        foreach my $region (@{$idx->{'exon'}}){
           my ($s,$e) = split/-/,$region;
           die "s >= e :$s >= $e at $region of $chr,$gene,$mrna" if($s >= $e);
           foreach my $id(@IDs){
             die "maf not exists ref $ref at $id" if(!exists $maf ->{$id}->{$ref});
             #next if(!exists $maf ->{$id}->{$chr});
             my ($s_maf,$e_maf) = ($maf->{$id}->{$chr}->[1],$maf->{$id}->{$chr}->[2]);
             die "s_maf >= e_maf :$s_maf >= $e_maf at maf id $id" if($s_maf >= $e_maf);
             if($e <= $s_maf || $e_maf <= $s){}else{
                my ($start,$end);
                if($s <= $s_maf){$start = 0}else{$start = $s-$s_maf-1}  
                if($e >= $e_maf){$end = $e_maf-$s_maf}else{$end = $e-$s_maf}  
                die "not found order at $id of maf" if(!exists $maf ->{$id}->{'order'});
                if(defined $outfmt && $outfmt eq "maf"){print "#exon $mrna:$s-$e\na score=$maf->{$id}->{'score'}\n"}elsif(defined $outfmt && $outfmt eq "mfa"){print "#exon $mrna:$s-$e\n"}
                foreach my $spec(@{$maf ->{$id}->{'order'}}){
                   if($spec eq $ref){    
                      my $len = $end-$start;
                      my $extract = substr($maf->{$id}->{$spec}->[6],$start,$len);
                      if(defined $outfmt && $outfmt eq "maf"){
                        printf ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $s_maf+$start, $len, $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $extract);
                      }elsif(defined $outfmt && $outfmt eq "mfa"){
                           my $real_start = $s_maf+$start;
                           my $real_end = $s_maf+$end;
                           print ">$maf->{$id}->{$spec}->[0] $real_start-$real_end:$maf->{$id}->{$spec}->[4]\n$extract\n";
                        }else{die "outfmt not defined"}
                   }else{
                          #print STDERR "spec is not ref:$spec, chr is $maf->{$id}->{$spec}->[0], s_maf is $maf->{$id}->{$spec}->[1]\n";
                          ($s_maf,$e_maf) = ($maf->{$id}->{$spec}->[1],$maf->{$id}->{$spec}->[2]);                         
                          my $len = $end-$start;
                          my $extract = substr($maf->{$id}->{$spec}->[6],$start,$len);
                          my $cnt_real_extract = &countReal($extract);
                          my $extract_left = substr($maf->{$id}->{$spec}->[6],0,$start);                          
                          my $cnt_real_extract_left = &countReal($extract_left);                           
                          if(defined $outfmt && $outfmt eq "maf"){ 
                             printf ("s %-12s %-9s %-6s %s %-9s %s\n",$maf->{$id}->{$spec}->[0], $s_maf+$cnt_real_extract_left, $cnt_real_extract, $maf->{$id}->{$spec}->[4], $maf->{$id}->{$spec}->[5], $extract);
                          }elsif(defined $outfmt && $outfmt eq "mfa"){
                             my $real_start = $s_maf+$cnt_real_extract_left;
                             my $real_end = $s_maf+$cnt_real_extract;
                             print ">$maf->{$id}->{$spec}->[0] $real_start-$real_end:$maf->{$id}->{$spec}->[4]\n$extract\n";
                           }else{die "outfmt not defined"}       

                        }

                }                
                print "\n"

             }             
           }


        }
    }
  }
}


if(defined $outfmt && $outfmt eq "maf"){print "##eof maf\n"}



###sub


sub gffRead(){
  ##standard nest gff3 parsing code v1.0,  @ Dec. 26, 2015
  my $file=shift;
  my %gff;#main structure table:  chr->gene->mRNAs->(elements:CDS,exon,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)
  my %gene;#gene info
  my %mrna;#mrna info
  #1,read&fill without order consideration
  open GFF, $file or die "$!";
  while(<GFF>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
    die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
    next if($chr eq "chrUn" || $chr eq "chrSy" || $chr eq "Pt" || $chr eq "Mt"); #get rid of other chrs
    next if($ftype eq "chromosome" || $ftype eq "Chromosome");
    ($start,$end)=($start<$end)?($start,$end):($end,$start);
    my ($ID,$pID,$name);
    if($ftype eq "gene"){
      $attr=~/ID=([^;]+);Name=([^;]+)/;
      ($ID,$name)=($1,$2);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$name";
      }else{die "dup gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+)/;
        ($ID,$pID)=($1,$2);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID";
        }else{die "dup mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR" || $ftype eq "tss" || $ftype eq "tts"){
           $attr=~/Parent=([^;]+)/;
           $pID=$1;
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
         }else{die "unknow feature type\n"}
  }#while end
  close GFF;
  return (\%gff,\%gene,\%mrna);

}



sub mafRead_plus_withcheck(){
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
    my $has_ref =0;
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
              if($name eq $ref){die "ref $ref align not degapped: $seq_ali at $_\n" if ($seq_ori ne $seq_ali);$has_ref=1}
              my $end=$start+$len;
              my @dataline=($name,$start,$end,$len,$strand,$chr_len,$seq_ali,$seq_ori);
              die "length $len diff with seq_rmgap at block $cnt" if($len != length $seq_ori);
              my $alignLen = length $seq_ali;
              if(exists $block{$name}){die"duplicated species in maf block $cnt\n"}else{$block{$name}=\@dataline}
              if(exists $block{"alignLen"}){die "align length diff at block $cnt" if($block{"alignLen"} != $alignLen)}else{$block{"alignLen"} = $alignLen}
         }else{print STDERR "\nomit line at group $cnt\n: unknow maf line type:$line\n"}
    }#foreach end here 
    if($degree == 0){die "degree eq 0 at group $cnt: $_"}else{$block{"degree"} = $degree}
    die "not has ref:$ref at $_" if($has_ref != 1);
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

sub countReal(){
  my $seq = shift;
  if($seq eq ""){return 0}
  my @box = split//,$seq;
  my $cnt = 0;
  foreach (@box){
      if($_ ne "-"){$cnt++}
  }
  return $cnt;
}




