
use strict;
use warnings;

my @test=(5,6,7,1,3,19,18,20,21,90,21,21,21,90,91,100,25);
my @segs=@{&getContinuRegions(\@test)};

print "@segs\n";

sub getContinuRegions(){
  my $index=shift;
  my @sorted = sort {$a <=> $b} @{$index};
  print "@sorted\n";
  my $f;
  my @segs;
  foreach(0..$#sorted){
    if(!defined $f){$f=$sorted[$_]}
    if($_ == $#sorted){ #last one
      if($sorted[$_-1] != $sorted[$_]-1){push @segs,"$sorted[$_]"}
    }else{
          if($sorted[$_] == $sorted[$_+1]-1){
             if($_ == $#sorted-1){push @segs,"$f-$sorted[$_+1]"}
           }elsif($sorted[$_] == $sorted[$_+1]){
              #dedup numbers 
             }else{
                    if($f == $sorted[$_]){
                      push @segs,"$sorted[$_]" ;
                      $f=$sorted[$_+1] ;
                    }else{
                          push @segs,"$f-$sorted[$_]" ;
                          $f=$sorted[$_+1] ;
                         }
                  }
         }
  }#foreach end
  return \@segs;
}

