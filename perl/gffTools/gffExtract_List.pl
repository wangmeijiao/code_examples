#!/usr/bin/perl
use strict;
use warnings;
#usage :perl extractList.pl loci.list all.gene >loci.fa

 open LIST, "$ARGV[0]" or die "$!";
 open FA, "$ARGV[1]" or die "$!";

 my @list;
 while (<LIST>){
  chomp;
  push @list, $_;
  
 }
 close LIST;

 $/=">";
 my %seqs;
 while (<FA>){
  chomp;
  next if ($_ eq "");
  my @box=split /\n/,$_;
  my $head=shift @box;
  my @tmp=split / /,$head;  
  my $id=shift @tmp;
  my $seq=join "\n", @box;
  $seqs{$id}=">$head\n$seq\n";
 }
 $/="\n";
 close FA;

 
 foreach(@list){
   if(exists $seqs{$_}){print $seqs{$_};}else{die "not exists!"}
 }



