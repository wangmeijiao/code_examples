

sub readBed(){
   my $file=shift;
   my %bed;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my @box=split/[\t ]+/,$_;
     #die "trunct bed file $file at $_\n" if(scalar @box < 6);
     if(!exists $bed{$box[3]}){$bed{$box[3]}=join"\t",@box   }else{die "dup $box[3] at $_\n"}
   }
   close BED;
   return \%bed;
}


sub readBed_array(){
   my $file=shift;
   my @bed;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my @box=split/[\t ]+/,$_;
     #die "trunct bed file $file at $_\n" if(scalar @box < 6);
     push @bed,$_;
   }
   close BED;
   return \@bed;
}



sub readBed_hash(){
   my $file=shift;
   die "file $file not exists\n" if(! -f $file);
   my %bed;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my @box=split/[\t ]+/,$_;
     #die "trunct bed file $file at $_\n" if(scalar @box < 6);
     my $id = $box[3];
     $id=~s/\.\d+$//; $id = uc($id);
     if(!exists $bed{$id}){$bed{$id}=join"\t",@box }else{die "dup id $box[3] found in $_ of file $file"}
   }
   close BED;
   return \%bed;
}



sub readBed_list(){
  my ($info,$list) = @_ ;
  my %gene;
  foreach my $spec(@$list){
     my $file = $info->{$spec}->{'geneBed'};
     $gene{$spec} = &readBed($file);

  }
  return \%gene;
}






