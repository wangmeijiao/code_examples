use strict;
use warnings;
use Data::Dumper;

#read in gff and classify subfamilies of LTRs, DNATEs etc  
#the next step of  rpmk2gffv2.2.pl
#modified  @2016,  18 Jan. 
#the ID=xxx may not the same in a certain subfamily gff file

$/="##";
my %family;
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
      $subfamilyID=~s/[-_]OS$//i; 
      $subfamilyID=~s/-LTR$|-I$|-I-int$|_LTR$|_I$|_I-int$|_int-int$//i; #rm "Copia-104_SB-I" tail
      push @{$family{$familyID}->{$subfamilyID}->{$id}},$line;
    }elsif($box_in[2] eq "repeat_fragment"){
        $box_in[8]=~/^Parent=([^;]+);Note=([^#\|]+)\|/;
        my $id=$1;
        my $subfamilyID=$2;
      $subfamilyID=~s/[-_]OS$//i; 
        $subfamilyID=~s/-LTR$|-I$|-I-int$|_LTR$|_I$|_I-int$|_int-int$//i; #rm "Copia-104_SB-I" tail
        push @{$family{$familyID}->{$subfamilyID}->{$id}},$line;
      }else{die "unknown feature $box_in[2]"}    

  } 


}

$/="\n";




#print Dumper \%family;
#exit;
#print "family|subfamily|id\t##gff version 3\n";
foreach my $familyID(sort keys %family){
   #print "#$familyID\n";
   if(! -d $familyID){mkdir $familyID}
   foreach my $subfamID(sort keys %{$family{$familyID}}){
     #print "##$subfamID\n";
     open OUT, ">$familyID/$subfamID.gff" or die "$!";
     foreach my $id(sort keys %{$family{$familyID}->{$subfamID}}){
       #print "###$id\n";
       my $idx=$family{$familyID}->{$subfamID}->{$id};
       foreach(@{$idx}){print OUT "$_\n"}

=pod
       if(defined $family && $family ne ""){
          if($subfamID=~/$family/i){
             foreach(@{$idx}){
                 my @box=split/[\t ]+/,$_;
                 if(defined $seg){}              
                               print "$familyID|$subfamID|$id\t$_\n"

             }

          }
       }else{

            foreach(@{$idx}){print "$familyID|$subfamID|$id\t$_\n"}        
        }
=cut
    

    }#foreach id
    close OUT;
  }#foreach subfamID
}




