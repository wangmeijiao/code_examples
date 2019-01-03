use warnings;
use strict;

my ($align,$order,$num,$len) = &mfaRead($ARGV[0]);
print STDERR "readin mfa done, $len columns\n";

my $ref = $ARGV[1];

die "ref $ref not exist in mfa file" if(!exists $align->{$ref});


print STDERR "ref is $ref, start to degap\n";
my %align_real;
my $cnt;
for my $i(0..$len-1){
  my $ref_nuc = substr($align->{$ref},$i,1);
  if($i % 1000000 == 0){print STDERR "#"}
  if($ref_nuc ne "-"){
    $cnt++;
    foreach my $id(@$order){
      if($id eq $ref){
        $align_real{$id}.=$ref_nuc;  
      }else{
        my $nuc_extract = substr($align->{$id},$i,1);
        $align_real{$id}.=$nuc_extract;
      }
      if($cnt % 60 == 0){$align_real{$id}.="\n"}
    }

  }

}


print  STDERR "\nstart to output\n";
foreach my $id(@$order){
  print ">$id\n$align_real{$id}\n";
    
}



###sub


sub mfaRead(){
  my $file = shift;
  $/=">";
  my %align;
  my $num_align;
  my $len;
  my @align_order;
  open MFA, $file or die "$1";
  while(<MFA>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp = split/[\t ]+/,$id;
    my $id_real = $temp[0];
    push @align_order,$id_real;
    my $align=join"",@box;
    if(! defined $len){$len = length $align}else{ if($len != length $align){die "align len diff at $id"}}
    if(!exists $align{$id_real}){$align{$id_real}=$align; $num_align++}else{die "dup alignment $id_real\n"}
  }
  close MFA;
  $/="\n";
  return (\%align,\@align_order,$num_align,$len);
}



