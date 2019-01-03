use strict;
use warnings;
use Data::Dumper;
#use Statistics::Basic; #can't find functions
#use Statistics::Basic::Correlation; #can't find functions
use Statistics::Basic qw(:all); #use this

my @test = (1,2,3);
my $test1 = \@test;
#my $test1 = [1,2,3];
my $test2 = [1,2,3];  #equal to \@test

my $test3 = vector(1,2,3);

my $corr1 = corr($test1,$test2);

my $covariance  = covariance(  [1 .. 3], [1 .. 3] );
my $correlation = correlation( [1 .. 3], [1 .. 3] );

print "$test3,$corr1,$covariance,$correlation\n"


