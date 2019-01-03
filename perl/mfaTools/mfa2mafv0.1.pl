use strict;
use warnings;

#  no need to refDegap
#  no need to supply chr_len etc
#

#perl mfa2maf.pl Ochr05.japo-VS-glab-VS-punc-VS-brac-VS-lper.mfa |gzip > Ochr05.japo-VS-glab-VS-punc-VS-brac-VS-lper.mfa.maf.gz


my ($fa,$order,$len) = &faRead_chr($ARGV[0]);

print "##maf version=1 scoring=mfa2maf\n";

my $step = int $len/500; #500 segments



my $len_all = &getLenReal($fa,$order);


my $cnt;

my %index_real;
for(my $i=0;$i<=$len;$i+=$step){
  $cnt++;
  print "a score=$cnt\n";
  foreach my $id(@$order){
    my $seq_extract = substr($fa->{$id},$i,$step);
    my $seq_real = $seq_extract;
    $seq_real=~s/-+//g;
    my $len_real = length $seq_real;
    if(!exists $index_real{$id}){ printf("s %8s 1 $len_real + $len_all->{$id} $seq_extract\n",$id) }else{printf("s %8s $index_real{$id} $len_real + $len_all->{$id} $seq_extract\n",$id)}
    $index_real{$id}+=$len_real;
  }
  print "\n";

}


##sub

sub faRead_chr(){
  my $file = shift;
  my %fa;
  my @order;
  my $len;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp=split/[\t ]+/ ,$id;
    $id=shift @temp;
    if($id =~/japo/){$id="japo.".$id}elsif($id =~/glab/){$id="glab.".$id}elsif($id =~/punc/){$id="punc.".$id}elsif($id =~/brac/){$id="brac.".$id}elsif($id =~/Lper/){$id="lper.".$id}else{die"match err at $id"}
    push @order, $id;
    my $seq = join "", @box;
    if(!defined $len){$len = length $seq}else{if($len != length $seq){die "$id length unequal"}}
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return (\%fa,\@order,$len);
}



sub getLenReal(){

  my ($fa,$order) = @_;
  my %len_real;
  foreach my $id(@$order){
    my $seq_real = $fa->{$id};
    $seq_real=~s/-+//g;
    $len_real{$id} = length $seq_real;
  }
  return \%len_real;

}

