
  use strict;
  use warnings;



  $/="\n\n";
  while(<stdin>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_ =~/^\s+$/);
    #print "$_\n";
    my @temp=split/>/,$_;
    shift @temp;
    foreach my $ctrl(@temp){
      my @box = split/\n+/,$ctrl;
      my $id=shift @box;  
      my $align=join("",@box);
      printf("%-40s%s\n",$id,$align);
    }
    print "\n";
  }
  $/="\n";



