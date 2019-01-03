#!/usr/bin/perl

# parse repeatmasker out file and convert result to embl format

my %repeat;
open IN, "$ARGV[0].txt.out" or die "can not open my infile";
<IN>;
<IN>;
<IN>;
while(<IN>){
  my @unit=split(" ", $_);
  #print "$unit[9]\n";
  $part1="FT   repeat_region   $unit[5]\.\.$unit[6]\n";
  $part2="FT                   /rpt_family=\"$unit[10]\"\n";
  $part3="FT                   /rpt_name=\"$unit[9]\"\n";
  $part="$part1$part2$part3";
  $repeat{$unit[5]}=$part;
  #print "Key $unit[5]\nValue $part\n";
}
close IN;

 
my $genenumber=0;
my @gene;
open EMBL, "$ARGV[0].embl" or die "can not open my infile2";
while (<EMBL>) {
	
	if ($_=~/\sCDS\s/) {
        $genenumber++;
		
		if (length $gene[$genenumber] > 0) {
           $gene[$genenumber].=$_;
		}else{
		   #print "ok";
		   $gene[$genenumber]=$_;
		}
       # print "$gene[$genenumber]";
	}elsif($_=~/^SQ\s/){
		$genenumber++;
		if (length $gene[$genenumber] > 0) {
           $gene[$genenumber].=$_;
		}else{
			#print "yes\n";
		   $gene[$genenumber]=$_;
		}
	}else{
		if (length $gene[$genenumber] > 0) {
           $gene[$genenumber].=$_;
		}else{
			#print "yes\n";
		   $gene[$genenumber]=$_;
		}
        #print "$gene[$genenumber]";
	}
}
close EMBL;
$l=@gene;
print "$l\n";
#print "$gene[0]\n";

$header=shift @gene;
$tail=pop @gene;

#print "$header\t$tail\n";
my %cds;
foreach (@gene) {
    if ($_=~/FT   CDS             join\((\d+)\.\..*\d+\)/ or $_=~/FT   CDS             (\d+)\.\.\d+/) {
        $cds{$1}=$_;
		#print "Key $1\nValue $_\n";
	}
}

%hash=(%repeat,%cds);
my @start=sort {$a <=> $b} keys %hash;
#print @start, "\n";

open ME, ">$ARGV[0].merge" or die "can not open my outfile";
print ME "$header";
foreach  (@start) {
	print ME "$hash{$_}";
}
print ME "$tail";
close ME;