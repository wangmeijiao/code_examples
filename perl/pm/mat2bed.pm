sub mat2bed(){
   #mat table without header and rowid
   my $file = shift;
   my %pairs;
   open MAT, $file or die "$!";
   my ($ncol,$nrow);
   while(<MAT>){
     chomp;
     die "empty line found" if($_ eq "" || $_=~/^\s+$/);
     my @data = split/[\t ]+/,$_;
     if(!defined $ncol){$ncol = scalar @data}else{if($ncol != scalar @data){die "unequal(!=$ncol) data points found at $_"}}
     $nrow++;
     #for my $i(0..$#data){my $real_i = $i + 1; print "$nrow-$real_i","\t",$data[$i],"\n"      }
     for my $i(0..$#data){my $real_i = $i + 1; $pairs{$nrow}->{"$nrow-$real_i"} = $data[$i]      }
     if($nrow % 100 == 0){print STDERR "#"}
   }
   die "ncol($ncol) != nrow($nrow) at last" if($ncol != $nrow);
   print STDERR "\n";
   close MAT;
   return \%pairs;
   #return 1;

}



