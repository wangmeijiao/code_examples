use strict;
use warnings;


# igv can read the result maf and also accept tba.maf without any modification
# pecan global align and chainet pipeline, seems pecan is better ?

#  no need to refDegap
#  no need to supply chr_len etc
#  only for chromosome level mfa 
#  each maf block has the same length

# perl mfa2mafv0.2.pl chr03.japo-VS-zs97-VS-glab-VS-meri.mfa 100000 > chr03.japo-VS-zs97-VS-glab-VS-meri.mfa.maf

my ($fa,$order,$len) = &faRead_chr($ARGV[0]);
my $step = $ARGV[1];
$step||=10000;

print "##maf version=1 scoring=mfa2maf\n";

#my $step = int $len/$win; #500 segments



my $len_all = &getLenReal($fa,$order);
my @len;
foreach my $id (keys %$len_all){push @len, $len_all->{$id}}
my $maxlen_totallen = &maxStringLen(@len);

my $maxlen_name = &maxStringLen(@$order);


my $cnt;
my %index_real;
for(my $i=0;$i<=$len;$i+=$step){
  $cnt++;
  print "a score=$cnt\n";

  #prepare block
  my %block;
  foreach my $id(@$order){
    my $seq_extract = substr($fa->{$id},$i,$step);
    $block{$id}->{'seq_extract'} = $seq_extract;
    my $seq_real = $seq_extract;
    $seq_real=~s/-+//g;
    my $len_real = length $seq_real;
    $block{$id}->{'len_real'} = $len_real ;
    #if(!exists $index_real{$id}){ 
    #   printf("s %8s 1 $len_real + $len_all->{$id} $seq_extract\n",$id) 
    #}else{
    #      printf("s %8s $index_real{$id} $len_real + $len_all->{$id} $seq_extract\n",$id)
    # }
    if(!exists $index_real{$id}){$index_real{$id} = 1}
    $block{$id}->{'index_real'} = $index_real{$id};
    $index_real{$id}+=$len_real;
  }

  my @index_real;
  foreach my $id (keys %block){push @index_real, $block{$id}->{'index_real'} }
  my $maxlen_indexreal = &maxStringLen(@index_real);

  my @len_real;
  foreach my $id (keys %block){push @len_real, $block{$id}->{'len_real'} }
  my $maxlen_lenreal = &maxStringLen(@len_real);

  #output maf format
  foreach my $id(@$order){
     printf("s %-${maxlen_name}s %${maxlen_indexreal}s %${maxlen_lenreal}s + %${maxlen_totallen}s %s\n",$id,$block{$id}->{'index_real'}-1,$block{$id}->{'len_real'},$len_all->{$id},$block{$id}->{'seq_extract'})

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
    if($id =~/^chr/){$id="japo.".$id}elsif($id =~/^zschr/){$id="zs97.".$id}elsif($id =~/^Oglab/){$id="glab.".$id}elsif($id =~/^Omeri/){$id="meri.".$id}else{die"match err at $id"}
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

sub maxStringLen(){
  my @str = @_;
  die "empty str @str" if(scalar @str == 0);
  my $maxlen = 0;
  foreach (@str){if (!defined $maxlen){$maxlen = length $_ }else{if($maxlen < length $_){ $maxlen = length $_}      }    }
  return $maxlen;
}






