use strict;
use warnings;
use Tabix;
use Data::Dumper;
use feature qw( say );
#say "hello";

my $file = $ARGV[0];
my ($chrom, $start, $end) = ("H1_japo","100000","200000");


#####high level use as perl pm package

my $t = Tabix->new('-data' => $file); 
my $iter = $t->query($chrom, $start, $end);

while(my $n = $t->read($iter)){print "$n\n"}

my @names = $t->getnames;
print @names,"\n";


#####directly us low level c functions

        my $t = tabix_open($file);
        my $iter = tabix_query($t,$chrom, $start, $end);
        
        while (my $line = tabix_read($t, $iter)){print "$line\n"}
        tabix_iter_free($iter);
        
        my @names = tabix_getnames($t);
        print "@names\n";

