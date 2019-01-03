use strict;
use warnings;
use Data::Dumper;

my %tab;
my $cnt;
while(<stdin>){
  chomp;
  my @box = split/[\t ]+/,$_;
  $cnt++;
  #push @tab,$_;
  $tab{$cnt}->{"summary"}=$_;
}

#print Dumper \%tab;
my $tab_sorted = &colSort_multi(\%tab,"0,3","increase");
#my $tab_sorted = &colSort_multi(\%tab,"0,1","decrease");
#my $tab_sorted = &colSort_multi(\%tab,"0,1","increase");
#print Dumper $tab_sorted;

foreach my $id(@$tab_sorted){
  print $tab{$id}->{'summary'},"\n";


}



###sub


sub colSort_multi(){
    #sort columns zero based
    #from left to right
    #detect characteric and numeric automatically
    my ($idx,$cols,$method)=@_;
    my $flag;
    if ($cols eq ''){die "no sort columes"}
    if ($method eq 'increase'){$flag=1}elsif($method eq 'decrease'){$flag=-1}else{die "wrong sort method"}
    my @list;
    my @cols=split/,/,$cols;
    @list= sort {
                  foreach my $ctrl(0..$#cols){
                       #print STDERR "sort $ctrl column\n";
                       my @box1= split/\t/,$idx->{$a}->{"summary"};
                       my $col1 = $box1[$cols[$ctrl]];
                       my @box2= split/\t/,$idx->{$b}->{"summary"};
                       my $col2 = $box2[$cols[$ctrl]];

                       if(&isNumeric($col1) ==1 && &isNumeric($col2) == 1){ #numeric column    
                          #print STDERR "col1,col2 are numeric: $col1,$col2\n";
                          if( $col1 > $col2){return 1*$flag}elsif($col1 < $col2){return -1*$flag}elsif($col1 == $col2){
                             if($ctrl==$#cols){return 0}else{ next}
                          }
                        }elsif(&isNumeric($col1) == 0 && &isNumeric($col2) == 0){ #characteric column
                            #print STDERR "col1,col2 are characteric: $col1,$col2\n";
                            if( $col1 gt $col2){return 1*$flag}elsif($col1 lt $col2){return -1*$flag}elsif($col1 eq $col2){
                               if($ctrl==$#cols){return 0}else{ next}
                            }
                         }else{die "a and b not the same type data at :$a,$b"}
                  }
              }keys %{$idx};
    return \@list;
}



sub isNumeric(){
  my $idx=shift;
  my $flag=1;
  #foreach my $ctrl (@{$idx}){
    #signed    
    $idx=~s/^-//;
    #decimals
     if($idx =~/[^0-9\.]/){$flag = 0}

    #scientific notation

  #}
  return $flag;
}
