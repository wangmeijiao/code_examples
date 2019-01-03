use strict;
use warnings;
use Getopt::Long;

my $degap;

GetOptions("degap!",\$degap);

my @nameList;
my %align;
$/=">";
while(<stdin>){
   chomp;
   next if($_ eq "");
   my @box=split/\n+/,$_;
   my $id=shift @box;
   my @tmp=split/ +/,$id;
   $id=shift @tmp;
   push @nameList,$id;
   my $seq=join "",@box;
   if(!exists $align{$id}){$align{$id}=$seq}else{die "dup $id\n"}
}
$/="\n";

my $pairs=&pair(\@nameList);

if(&checkAlignLenAll(\%align)){
   foreach(@{$pairs}){
     my ($f,$l)=split/\|/,$_;
     my ($align_f,$align_l);
     if($degap){
       ($align_f,$align_l)=&degap_pair($align{$f},$align{$l});
       print"${f}_VS_${l}\n$align_f\n$align_l\n\n";
     }else{ print"${f}_VS_${l}\n$align{$f}\n$align{$l}\n\n"}
   }
}else{die "align length diffs\n"}



###sub### 

sub pair(){ #choose(\@elements,2)
    my $elements=shift;
    my @choose;
    for(my $i=0;$i<$#{$elements};$i++){
       for (my $j=$i+1;$j<=$#{$elements};$j++){
          push @choose,$elements->[$i]."|".$elements->[$j];
       }
    }
    return \@choose;
}


sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=length($index->{$id})}else{if($len != length($index->{$id}) ){$flag=0}}
  }
  return ($flag,$len);
}


sub choose(){
    my ($element,$n)=@_;
    #unfinished   
}


sub degap_pair(){
   #check len%3==0 before and after
   my ($f,$l)=@_;
   my $f_degap="";
   my $l_degap="";
   my @f=split//,$f;
   my $len_f=scalar @f;
   my @l=split//,$l;
   my $len_l=scalar @l;
   if($len_f != $len_l || $len_f % 3 != 0 || $len_l % 3 != 0){die "either len not equal or not mod 3 before degap: \n$f\n$l"}
   for my $i(0..$len_f-1){
     if($f[$i] eq "-" || $l[$i] eq "-"){next}else{$f_degap.=$f[$i]; $l_degap.=$l[$i]}
   }
   if(length $f_degap != length $l_degap || (length $f_degap) % 3 != 0 || (length $l_degap) % 3 != 0){die "either len not equal or not mod 3 after degap: \n$f_degap\n$l_degap"}
  return ($f_degap,$l_degap);

}


