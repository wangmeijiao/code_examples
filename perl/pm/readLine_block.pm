use strict;
use warnings;
use Data::Dumper;

#&read_mfaStat_view($ARGV[0]);
print Dumper &read_mfaStat_view($ARGV[0]);

###

sub read_mfaStat_view(){
   #segs may be empty or few
   my $file = shift;
   open IN, $file or die "$!"; 
   my $cnt;
   my @box;
   my %segs;
   while(<IN>){
     $cnt++;
     if($cnt % 11 == 0 ){
       push @box,$_;
       #print "cnt$cnt\n",join"\n",@box,"\n\n";    
       chomp @box;
       die "box not 11 lines at $box[0]" if(scalar @box != 11);
       die "not #all equal# at $box[0]" if($box[6] ne "#all equal#");
       my ($gene1,$aln1) = split/[\t ]+/,$box[0]; #empty aln had been removed
       my @temp1 = split/\|/,$gene1;
       my ($gene2,$aln2) = split/[\t ]+/,$box[1]; 
       my @temp2 = split/\|/,$gene2;
       my $id;
       if($temp1[0] =~/^chr/){$id = $temp1[2]}elsif($temp2[0] =~/^chr/){$id=$temp2[2]}else{die "not matched japo id at $gene1,$gene2"}       
       my @segs = split/[\t ]+/,$box[8];
       if($box[8] eq "" || scalar @segs == 0){print STDERR "empty segs at $gene1 and $gene2\n"}
       if(!exists $segs{$id}){$segs{$id}=\@segs}else{die "dup id $id at $gene1, $gene2"}
       @box = ();
     }else{push @box,$_ }

   }

   close IN;
   print STDERR "total cnt $cnt\n";
   return \%segs;
}




