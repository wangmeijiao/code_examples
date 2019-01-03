use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#              link_orthtable
#  speciesA  ----------------->  speciesB
#                 

#link file types for liftover:  orthoTab
#usage: cat in.bed |perl thisfile.pl -tab blastm8.tab -gap xxx  > out.liftOver.bed


my $tab;
GetOptions("tab=s",\$tab);

#1, construct link_table
my ($tab_idx,$pair_idx)=&buildLinkTab($tab);
#print Dumper $tab_idx,$pair_idx;

#2, map anchor_query to anchor_target by searching in link_table
while(<stdin>){
   chomp;
   next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
   #my ($spec,$id)=split/[\t ]+/,$_;
   my $id=$_;
   $id=~s/\s+//g;
   print STDERR "start to lift for $id\n";
   if(exists $tab_idx->{$pair_idx->{$id}}){
       print "$tab_idx->{$pair_idx->{$id}}\n" 
   }else{print STDERR "unmapped $id\n"}
} 





##sub###

sub buildLinkTab(){
  my $file = shift;
  my %tab;
  my %pair;
  open ORTHO, $file or die "$!";
  while(<ORTHO>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
    my @box=split/[\t ]+/,$_;
    if(!exists $tab{$box[4]}){$tab{$box[4]}="$box[0]\t$box[1]\t$box[2]\t$box[3]\t$box[4]\t$box[5]" }else{die "dup $box[4] at $_"}
    if(!exists $tab{$box[10]}){$tab{$box[10]}="$box[6]\t$box[7]\t$box[8]\t$box[9]\t$box[10]\t$box[11]" }else{die "dup $box[10] at $_"}
    $pair{$box[4]}=$box[10];
  }
  return \%tab, \%pair;
}




sub liftOver(){




}



