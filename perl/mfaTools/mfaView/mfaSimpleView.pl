
  use strict;
  use warnings;

  
  $/='>';
  while(<stdin>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_ =~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;  
    my $align=join("",@box);
    printf("%-30s%s\n",$id,$align);
  }
  $/='\n';
