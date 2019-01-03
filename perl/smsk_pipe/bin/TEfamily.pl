use strict;
use warnings;
use Data::Dumper;

#read in gff and classify subfamilies of LTRs, DNATEs etc
#the next step of  rpmk2gffv2.2.pl


$/="##";
my %subFamily;
while (<stdin>){
  chomp;
  next if($_ eq "" || $_ =~/^\s+$/ || $_ =~/gff version 3/i);
  my @box=split/\n+/,$_;
  my $familyID=shift @box;
  foreach my $line(@box){
    my @box_in=split/[\t ]+/,$line; 
    if($box_in[2] eq "Transposon"){
      $box_in[8]=~/^ID=([^;]+);Note=([^#]+)#/;
      my $id=$1;
      my $subfamilyID=$2;
      push @{$subFamily{$familyID}->{$subfamilyID}},$line;
    }elsif($box_in[2] eq "repeat_fragment"){
        $box_in[8]=~/^Parent=([^;]+);Note=([^#\|]+)\|/;
        my $id=$1;
        my $subfamilyID=$2;
        $subfamilyID=~s/-LTR$|-I$|-I-int$|_LTR$|_I$|_I-int$//i; #rm "Copia-104_SB-I" tail
        push @{$subFamily{$familyID}->{$subfamilyID}},$line;
      }else{die "unknown feature $box_in[2]"}    

  } 


}

$/="\n";

#print Dumper \%subFamily;
#exit;
print "##gff version 3\n";
foreach my $family(keys %subFamily){
   print "#$family\n";
   foreach my $subfam(sort keys %{$subFamily{$family}}){
     print "##$subfam\n";
     foreach(@{$subFamily{$family}->{$subfam}}){
       print "$_\n";

     }

   }
}


