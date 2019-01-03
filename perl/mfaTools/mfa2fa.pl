use strict;
use warnings;



my ($fa,$order) = &faRead_chr($ARGV[0]);
foreach my $id(@$order){
  print ">$id\n";
  print &faFormat($fa->{$id},50);


}








##sub

sub faRead_chr(){
  my $file = shift;
  my %fa;
  my @order;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp=split/[\t ]+/ ,$id;
    $id=shift @temp;
    push @order,$id;
    my $seq = join "", @box;
    $seq=~s/-+//g;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return (\%fa,\@order);
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






