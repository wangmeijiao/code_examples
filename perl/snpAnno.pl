use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;




#total 3k species: 3046?
#uniq_lzk_snp header: 
#   CHR     POS     Major   Major_freq      Minor   MAF     DEP     Allele_number   Ref_nt  Ref_Num  Alt_nt  Alt_Num
#only statistic for biallele 

#lifted snp.biallele 
#CHR     START END  TYPE(Ref-Alt) . STRAND(always+) Major   Major_freq      Minor   MAF     DEP     Allele_number   Ref_nt  Ref_Num Alt_nt  Alt_Num


my ($snpf, $ancef,$genomef,$chr);
$chr||="chr03";
GetOptions("snp=s",\$snpf,"ance=s",\$ancef,"genome=s",\$genomef,"chr=s",\$chr);



my $genome = &faRead_chr($genomef);
my $ance = &faRead_chr($ancef);

die "$chr not exists either in genome or ances" if(!exists $genome->{$chr} || !exists $ance->{$chr});
die "length of $chr of genome and ancestral not equal" if(length $genome->{$chr} != length $ance->{$chr});

my $snp = &readSNP_uniq_bed($snpf);

print STDERR "readin genome, ancestral and snp file done\n";
#print Dumper $genome,$ance,$snp; exit;

my $cnt;
foreach my $s (sort {$a<=>$b} keys %{$snp->{$chr}}){
   $cnt++;
   my ($name,$major,$majorP,$minor,$minorP,$depth,$totalAlle,$refNuc,$refNum,$varNuc,$varNum) = split/\t/,$snp->{$chr}->{$s};
   #check valid
   die "major is lowcase" if($major=~/[agct]/);
   die "minor is lowcase" if($minor=~/[agct]/);
   die "ref is lowcase" if($refNuc=~/[agct]/);
   die "var is lowcase" if($varNuc=~/[agct]/);  
   die "major eq minor at snp position $s of chr $chr" if($major eq $minor);
   die "refNuc eq varNuc: $refNuc eq $varNuc at snp position $s of chr $chr" if($refNuc eq $varNuc);

   die "varNuc is not valid:$varNuc at snp position $s of chr $chr" if($varNuc ne "A" && $varNuc ne "G" && $varNuc ne "C" && $varNuc ne "T");  #, - etc
   die "refNuc is not valid:$refNuc at snp position $s of chr $chr" if($refNuc ne "A" && $refNuc ne "G" && $refNuc ne "C" && $refNuc ne "T");  #, - etc


   my ($ref_nuc, $var_nuc) = split/-/,$name;
   my $ref_nuc_extract = uc(substr($genome->{$chr},$s-1,1)); #upcase only
   my $nuc_ance = substr($ance->{$chr},$s-1,1); #can be upcase lowcase . N -


   if($ref_nuc ne $refNuc){die "ref_nuc($ref_nuc) ne refNuc($refNuc) at snp position $s"}
   if($var_nuc ne $varNuc){die "ref_nuc($ref_nuc) ne refNuc($refNuc) at snp position $s"}   
   if($ref_nuc_extract ne $refNuc){print STDERR "!!ref_nuc_extract(tigr6 $ref_nuc_extract) ne refNuc($refNuc) at snp position $s\n"; next}

   #print "snp position is $s, ref_nuc is $ref_nuc, varNuc is $varNuc; ances_nuc is $nuc_ance\nmajor($major,$majorP),minor($minor,$minorP)\n";
   

   my ($flag_depth_low,$flag_totalAlle_low) ;
   if($depth <= 10000){$flag_depth_low = "yes"}else{$flag_depth_low = "no"}
   if($totalAlle <= 1500){$flag_totalAlle_low = "yes"}else{$flag_totalAlle_low = "no"}


   my ($mut_type,$DAF,$ances_freq);
   if($nuc_ance eq "-" ){
     #print "ancestral is -(gap), can't determine the derived allele, skip..\n";
     #next;
   }elsif($nuc_ance eq "." ){
     #print "ancestral is .(dot), can't determine the derived allele, skip..\n";
     #next;
    }elsif($nuc_ance eq "N" || $nuc_ance eq "n" ){
       #print "ancestral is N, can't determine the derived allele, skip..\n";
       #next;
     }elsif($nuc_ance=~/[ATCG]/){
         if($nuc_ance eq $refNuc && $nuc_ance ne $varNuc){
           print "good! derived nuc is varNuc: $nuc_ance(ance) eq $refNuc(ref) && $nuc_ance(ance) ne $varNuc(var);";
           $DAF = $varNum/$totalAlle;
           $ances_freq = $refNum/$totalAlle;
           if($refNuc=~/[GC]/ && $varNuc=~/[AT]/){ $mut_type="strong->weak"  
            }elsif($refNuc=~/[AT]/ && $varNuc=~/[GC]/){ $mut_type="weak->strong"}else{ $mut_type="unknow"}   
           print "snp $s: ancestral->derived=$refNuc->$varNuc, DAF=$DAF, mutation_type: $mut_type, low_depth:$flag_depth_low, low_totalAllele:$flag_totalAlle_low\n";     

         }elsif($nuc_ance ne $refNuc && $nuc_ance eq $varNuc){
            print "good! derived nuc is ref: $nuc_ance(ance) ne $refNuc(ref) && $nuc_ance(ance) eq $varNuc(var);";         
            $DAF = $refNum/$totalAlle;
            $ances_freq = $varNum/$totalAlle;
            if($varNuc=~/[GC]/ && $refNuc=~/[AT]/){ $mut_type="strong->weak"  
              }elsif($varNuc=~/[AT]/ && $refNuc=~/[GC]/){ $mut_type="weak->strong"}else{$mut_type="unknow"}
            print "snp $s: ancestral->derived=$varNuc->$refNuc, DAF=$DAF, mutation_type: $mut_type, low_depth:$flag_depth_low, low_totalAllele:$flag_totalAlle_low\n";

          }elsif($nuc_ance ne $refNuc && $nuc_ance ne $varNuc){
              print STDERR "Oops, both refNuc and varNuc ne nuc_ance: $nuc_ance ne $refNuc && $nuc_ance ne $varNuc\n";
              #next;
            }else{die "unknow circurrence at snp position $s of chr $chr"}
             
      }elsif($nuc_ance=~/[atcg]/){
            #print "ancestral is locase, can't determine the derived allele, skip..\n";
            $nuc_ance = uc($nuc_ance);  
            if($nuc_ance eq $refNuc && $nuc_ance ne $varNuc){
            print "good(lowcase)! derived nuc is varNuc: $nuc_ance(ance) eq $refNuc(ref) && $nuc_ance(ance) ne $varNuc(var);";
            $DAF = $varNum/$totalAlle;
            $ances_freq = $refNum/$totalAlle;
            if($refNuc=~/[GC]/ && $varNuc=~/[AT]/){ $mut_type="strong->weak"
             }elsif($refNuc=~/[AT]/ && $varNuc=~/[GC]/){ $mut_type="weak->strong"}else{ $mut_type="unknow"}
            print "snp $s: ancestral->derived=$refNuc->$varNuc, DAF=$DAF, mutation_type: $mut_type, low_depth:$flag_depth_low, low_totalAllele:$flag_totalAlle_low\n";

             }elsif($nuc_ance ne $refNuc && $nuc_ance eq $varNuc){
                print "good(lowcase)! derived nuc is ref:$nuc_ance(ance) ne $refNuc(ref) && $nuc_ance(ance) eq $varNuc(var);";
                $DAF = $refNum/$totalAlle;
                $ances_freq = $varNum/$totalAlle;
                if($varNuc=~/[GC]/ && $refNuc=~/[AT]/){ $mut_type="strong->weak"
                  }elsif($varNuc=~/[AT]/ && $refNuc=~/[GC]/){ $mut_type="weak->strong"}else{$mut_type="unknow"}
                print "snp $s: ancestral->derived=$varNuc->$refNuc, DAF=$DAF, mutation_type: $mut_type, low_depth:$flag_depth_low, low_totalAllele:$flag_totalAlle_low\n";

               }elsif($nuc_ance ne $refNuc && $nuc_ance ne $varNuc){
                  print STDERR "Oops, both refNuc and varNuc ne nuc_ance: $nuc_ance ne $refNuc && $nuc_ance ne $varNuc\n";
                  #next;
                }else{die "unknow circurrence at snp position $s of chr $chr"}


         }
   

   #print STDERR "#" if($cnt % 100000 == 0);
   #last if($cnt >= 5);
}


print STDERR "\ntotal $cnt snp processed\n";






##sub



sub faRead_chr(){
  my $file = shift;
  my %fa;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp=split/[\t ]+/ ,$id;
    $id=shift @temp;
    my $seq = join "", @box;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return \%fa;
}

sub readSNP_uniq(){
   my $file = shift;
   my %snp;
   open IN, $file or die "$!";
   my $cnt=0;
   while(<IN>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my ($chr,$pos,$major,$majorP,$minor,$minorP,$depth,$total,$refNuc,$refNum,$varNuc,$varNum) = split/[\t ]+/,$_;
     $chr=~s/^Chr//;
     if(!exists $snp{$chr}->{$pos}){$snp{$chr}->{$pos} = "$major-$minor"}else{die"dup position at $chr:$pos"}
     $cnt++;
     if($cnt % 200000 == 0){print STDERR "#"}
   }
   close IN;
   print STDERR "\n$file, total:$cnt\n";
   return \%snp;
}


sub readSNP_uniq_bed(){
   my $file = shift;
   my %snp;
   open IN, $file or die "$!";
   my $cnt=0;
   while(<IN>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my ($chr,$s,$e,$name,undef,$strand,$major,$majorP,$minor,$minorP,$depth,$total,$refNuc,$refNum,$varNuc,$varNum) = split/[\t ]+/,$_;
     $chr=~s/^Chr//;
     if(!exists $snp{$chr}->{$s}){$snp{$chr}->{$s} = join("\t",($name,$major,$majorP,$minor,$minorP,$depth,$total,$refNuc,$refNum,$varNuc,$varNum)) }else{die"dup position at $chr:$s"}
     $cnt++;
     #if($cnt % 200000 == 0){print STDERR "#"}
   }
   close IN;
   #print STDERR "\n$file, total:$cnt\n";
   return \%snp;
}


sub readSNP_bim(){
   my $file = shift;
   my %snp;
   open IN, $file or die "$!";
   my $cnt=0;
   while(<IN>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my ($chr,$pos,$minor,$major) = split/[\t ]+/,$_;
     if(!exists $snp{$chr}->{$pos}){$snp{$chr}->{$pos} = "$major-$minor"}else{die"dup position at $chr:$pos"}
     $cnt++;
     if($cnt % 200000 == 0){print STDERR "#"}
   }
   close IN;
   print STDERR "\n$file, total:$cnt\n";
   return \%snp;
}







