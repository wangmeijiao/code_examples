use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#usage: cat in.mfa |perl thisfile.pl  > in.mfa.stats
#input: mfa alignment file format only (two or more seqs aligned)
#output: #seg_num, #length, %coverage, %identity_overall, %identity_match, and best continuous blocks( similar to mablock does)

my ($quiet,$mablock);
GetOptions("quiet!",\$quiet,"mablock!",\$mablock);

$/=">";
my %align;
while(<stdin>){
  chomp;
  next if($_ eq "" || $_ =~/^#/);
  my @box=split/\n+/,$_; 
  my $id=shift @box;
  my $align=join"",@box;
  if(!exists $align{$id}){$align{$id}=$align}else{die "dup alignment $id\n"}
}
$/="\n";

if(!$quiet){
 foreach my $id(sort keys %align){
  printf("%-20s  %s\n",$id,$align{$id});
 }
}


my $seq_index=&getRealSeq(\%align);
#print Dumper $seq_index;

my ($flag,$len)=&checkAlignLenAll(\%align);
if($flag){
  my $stars=&getStars(\%align,$len);
  if(!$quiet){printf "%-20s  %s\n"," ",join("",@{$stars})}
  if(!$quiet){print "#all equal#\n"}
  my($coverage,$identity,$segs,$segs_linked)=&statBlocks($stars); 
  if(!$quiet){print "@{$segs}\n@{$segs_linked}\n"}
  if(!$quiet){print "total aligned length: $len\ncoverage: $coverage\nidentity: $identity\n"}
  if($quiet){print "$coverage\t$identity\n"}
  if($mablock && !$quiet){
    foreach my $id(sort keys %align){
      print STDERR ">$id\n";
      foreach my $seg(@{$segs}){
         my ($s,$e)=split/-/,$seg;
         next if($e-$s+1 <= 5);# filter for all members, because @segs is the same
         my $extract = substr($align{$id},$s,$e-$s+1);         
         $extract =~s/-//g;
         print STDERR "$extract\n";
      } 
    }
  }

}else{die" length not equal\n"}




###sub### 

sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=length($index->{$id})}else{if($len != length($index->{$id}) ){$flag=0}}
  }
  return ($flag,$len);
}


sub getRealSeq(){
   my $index=shift;
   my %seq=%{$index};
   foreach my $id(keys %seq){
      my $seq=$seq{$id};
      $seq=~s/-//g;
      $seq{$id}=$seq;
   }   
   return \%seq;
}

sub getStars(){
    my $index=shift;
    my $len=shift;
    my @stars;
    my %align_array;
    my @ids;
    foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align_array{$id}=\@box;
      push @ids,$id;
    } 
    #print Dumper \%align_array;
    for my $i (0..$len-1){
       my @elements;
       foreach my $id(@ids){
         push @elements, $align_array{$id}->[$i];
       }
       my $flag=&checkIdentity(\@elements);
       if($flag == 1){push @stars,"*"}elsif($flag == 0){push @stars," "}elsif($flag == -1){push @stars,"-"}else{die "undefined flag $flag\n"}
    }
    return \@stars;
}

sub statBlocks(){
    #*: the same
    #-: gap
    #(space): diff
    my $star=shift;
    my $len=scalar @{$star};
    my ($cnt_star,$cnt_gap,$cnt_diff)=(0,0,0);
    my ($coverage,$coverage_total,$identity)=(0,0,0);
    my $gap=4; #for link nearby similar blocks
    #find continual segments
    my @segments;
    my ($start,$end,$cnt)=(0,0,0);
    for my $i(0..$len-1){
      if($star->[$i] eq "*"){$cnt_star++}elsif($star->[$i] eq "-"){$cnt_gap++}elsif($star->[$i] eq " "){$cnt_diff++}     
      if($star->[$i] eq "*"){
         $cnt++;
         if($i==($len-1)){ #last element
           $start=$i-$cnt+1;
           $end=$i;
           push @segments,$start."-".$end;
           $cnt=0;
         }
      }else{
            $start=$i-$cnt;
            $end=$i-1;
            if($start<=$end){push @segments,$start."-".$end}
            $cnt=0;
           }
    }#for end
    my @segs=@segments; #@segments changed in next steps
    #print "@segs\n";

    #link small segments to long ones with gap<=4
    my @segments_linked;
    for (my $i=1;$i<=$#segments;$i++){
        my @tmp1=split /-/,$segments[$i-1];
        my @tmp2=split /-/,$segments[$i];
        if(($tmp2[0]-$tmp1[1]-1)<= $gap){
        $start=$tmp1[0];
        $end=$tmp2[1];
        $segments[$i]=$start."-".$end;
        }else{push @segments_linked, $segments[$i-1];}
        if($i==$#segments){push @segments_linked,$segments[$i];}
    }
    if($#segments==0){push @segments_linked,$segments[0];}
    #print "@segments_linked\n";

    #start to calculate percentage etc
    $coverage=sprintf "%.3f",($cnt_star+$cnt_diff)/$len; 
    $identity=sprintf "%.3f",$cnt_star/($cnt_star+$cnt_diff);
    return ($coverage,$identity,\@segs,\@segments_linked);
}

sub checkIdentity(){
    my @elements=@{shift @_};
    my $flag=1;
    my $foo;
    foreach (@elements){
      if($_ eq "-"){$flag=-1;return $flag}
      if(!defined $foo){$foo = $_}else{if($foo ne $_){$flag = 0}}
    }
    return $flag;
}


