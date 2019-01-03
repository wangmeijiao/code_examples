
use strict;
use warnings;
use GO::Parser;
use Data::Dumper;


#see https://metacpan.org/pod/distribution/go-perl/go-perl.pod
#https://metacpan.org/pod/GO::Model::Term

my $go_file = $ARGV[0];
my $go_id = $ARGV[1];

my $parser=new GO::Parser ({format=>'obo_text',
                            handler=>'obj'});
$parser -> parse($go_file);
print STDERR "parse go file $go_file done\n";


#print Dumper $parser;



##start to query go id

my $graph=$parser->handler->graph;
my $term = $graph->get_term($go_id);

my $id = $term->acc;
my $name = $term->name;
my $name_space = $term->namespace();
#my $code = $term->get_code_from_namespace;
my $def = $term->definition;


printf "Got term: %s\n%s\n%s\n%s\n", $id, $name, $name_space, $def;
#print "$term->acc, $term->name\n" ;


=pod
my $path=$graph->paths_to_top($go_id);
#print Dumper $path;

#my @level;
foreach my $idx(@$path){
  print Dumper $idx;
#  my $len=$_->length+1;
#  push (@level,$len);

}

#print "@level\n";
=cut


