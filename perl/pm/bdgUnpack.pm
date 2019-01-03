sub bdgUnpack_plus(){
   #must sort bedgraph first
   #with strict coordinate check
   #cut more and fill less at begin and end
   my ($data,$range)= @_;
   my ($begin,$end)=split/-/,$range;
   my @data_unpacked;
   die "empty dataline" if (scalar @$data == 0 );
   for my $i(0..$#{$data}){
     my($chr1,$s1,$e1,$v1)=split/[\t ]+/,$data->[$i];
     die "s1 >= e1 at $s1, $e1" if($s1 >= $e1);
     if($i == 0){
        if($s1 != $begin){
        #if($s1 != $begin && $begin+1 != $s1){
           if($s1 >= $begin){
             print STDERR "not continuous at begining: $begin-$s1,fill $begin to $s1 with 0\n";
             for my $i($begin..$s1-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t0"}
             for my $i($s1..$e1-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t$v1"}
           }else{
                 print STDERR "extra data points at begining: $s1-$begin-$e1,cut $s1 to $begin\n";
                 if($e1 <= $end){        
                   for my $i($begin..$e1-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t$v1"}
                 }else{
                        for my $i($begin..$end){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t$v1"}
                      }

                }

        }else{for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}}
     }elsif($i != 0 && $i!=$#{$data}){
         my ($chr2,$s2,$e2,$v2)=split/[\t ]+/,$data->[$i+1];
         die "s2 >= e2 at $s2, $e2" if($s2 >= $e2);
         if($e1 <= $s2){}else{
           print STDERR "not sorted at $s1-$e1, $s2-$e2 \n" if($s1 >= $s2 || $e1 >= $e2);
           if($e1 >= $s2){print STDERR "overlap at $s1-$e1, $s2-$e2 \n"};
           exit;
         }
         die "chr diffs at: $chr1, $chr2" if($chr1 ne $chr2);
         if($e1 != $s2 ){
         #if($e1 != $s2 && $e1+1 != $s2){
           for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}
           print STDERR "not continuous at $s1-$e1 and $s2-$e2,fill $e1 to $s2 with 0\n";
           for my $i($e1..$s2-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t0"}
         }else{
                #print STDERR "$e1 == $s2 at  $s1-$e1, $s2-$e2\n ";
                for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}
              }
     #}
#=pod
     }elsif($i == $#{$data}){
         if($e1 != $end){
            if($e1 <= $end){
              for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}
              print STDERR "not continuous at end: $e1-$end,fill $e1 to $end with 0\n";
              for my $i($e1..$end-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t0"}
            }else{
                   print STDERR "extra data at end: $s1-$end-$e1,cut $end to $e1 with 0\n";
                   for my $i($s1..$end-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t$v1"}
                 }

         }else{for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}}
       }else{die "unknow situation"}
#=cut
  }#for end
  return \@data_unpacked;
}



sub bdgUnpack(){
   #must sort bedgraph first
   #with strict coordinate check
   #fill omitted segments with 0, including the begining and end segment
   #be careful when bdg segment is 1, use bdgUnpack_new instead
   my ($data,$range)= @_;
   my ($begin,$end)=split/-/,$range;
   my @data_unpacked;
   for my $i(0..$#{$data}){
     my($chr1,$s1,$e1,$v1)=split/[\t ]+/,$data->[$i];
     die "s1 >= e1 at $s1, $e1" if($s1 >= $e1);
     if($i == 0){
        if($s1 != $begin){
        #if($s1 != $begin && $begin+1 != $s1){
           print STDERR "not continuous at begining: $begin-$s1,fill $begin to $s1 with 0\n";
           for my $i($begin..$s1-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t0"} 
        }
     }
     if($i!=$#{$data}){
         my ($chr2,$s2,$e2,$v2)=split/[\t ]+/,$data->[$i+1];
         die "s2 >= e2 at $s2, $e2" if($s2 >= $e2);
         if($e1 <= $s2){}else{
           print STDERR "not sorted at $s1-$e1, $s2-$e2 \n" if($s1 >= $s2 || $e1 >= $e2);
           if($e1 >= $s2){print STDERR "overlap at $s1-$e1, $s2-$e2 \n"};
           exit;
         }
         die "chr diffs at: $chr1, $chr2" if($chr1 ne $chr2);
         if($e1 != $s2 ){
         #if($e1 != $s2 && $e1+1 != $s2){
           for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}
           print STDERR "not continuous at $s1-$e1 and $s2-$e2,fill $e1 to $s2 with 0\n";
           for my $i($e1..$s2-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t0"}
         }else{
                print STDERR "$e1 == $s2 at  $s1-$e1, $s2-$e2\n ";
                for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}
              }
     }else{
         if($e1 != $end){
            for my $i($s1..$e1-1){my $j=$i+1; push @data_unpacked,"$chr1\t$i\t$j\t$v1"}
            print STDERR "not continuous at end: $e1-$end,fill $e1 to $end with 0\n";
            for my $i($e1..$end-1){my $j=$i+1;push @data_unpacked, "$chr1\t$i\t$j\t0"} 
         }
      }
  }#for end
  return \@data_unpacked;
}




sub bdgUnpack_new(){
   #check chr
   #check sort
   #check overlap
   #ok if not continuous
   #unpack region by directly fill in 
   my ($data,$range)= @_;
   my ($chr,$begin,$end)=split/:|-/,$range;
   die "begin >= end: $begin >= $end" if($begin >= $end);
   my @fill;
   for my $i(0..$end-$begin-1){$fill[$i] = 0} #initialization  
   die "empty dataline" if (scalar @$data == 0 );
   if(scalar @$data == 1){
     my($chr1,$s1,$e1,$v1)=split/[\t ]+/,$data->[0];
     die "chr unequal: $chr,$chr1" if($chr ne $chr1);
     die "s1 >= e1 at $s1, $e1" if($s1 >= $e1);     
     my $cut_s;
     if($s1 <= $begin ){$cut_s = 0}else{$cut_s = $s1-$begin-1}   
     my $cut_e;
     if($e1 >= $end ){$cut_e = $end-$begin-1 }else{$cut_e = $e1-$begin-1}         
     for my $i($cut_s..$cut_e){$fill[$i] = $v1}
   }else{#check and fill
       for my $i(0..$#{$data}-1){
         my($chr1,$s1,$e1,$v1)=split/[\t ]+/,$data->[$i];
         my($chr2,$s2,$e2,$v2)=split/[\t ]+/,$data->[$i+1];
         die "s1 >= e1 at $s1, $e1" if($s1 >= $e1);
         die "s2 >= e2 at $s2, $e2" if($s2 >= $e2);
         die "chr diffs at: $chr1, $chr2 or $chr" if($chr1 ne $chr2 || $chr1 ne $chr);
         if($e1 <= $s2){}else{
           print STDERR "not sorted at $s1-$e1, $s2-$e2 \n" if($s1 >= $s2 || $e1 >= $e2);
           if($e1 >= $s2){print STDERR "overlap at $s1-$e1, $s2-$e2 \n"};
           exit;
         }
         if($i == $#{$data}-1){
             my $cut_s1;
             if($s1 <= $begin ){$cut_s1 = 0;}else{$cut_s1 = $s1-$begin-1}
             my $cut_e1;
             if($e1 >= $end ){$cut_e1 = $end-$begin-1 }else{$cut_e1 = $e1-$begin-1}
             for my $i($cut_s1..$cut_e1){$fill[$i] = $v1}
             my $cut_s2;
             if($s2 <= $begin ){$cut_s2 = 0;}else{$cut_s2 = $s2-$begin-1}
             my $cut_e2;
             if($e2 >= $end ){$cut_e2 = $end-$begin-1 }else{$cut_e2 = $e2-$begin-1}
             for my $i($cut_s2..$cut_e2){$fill[$i] = $v2}
         }else{
                my $cut_s;
                if($s1 <= $begin ){$cut_s = 0;}else{$cut_s = $s1-$begin-1}
                my $cut_e;
                if($e1 >= $end ){$cut_e = $end-$begin-1 }else{$cut_e = $e1-$begin-1}
                for my $i($cut_s..$cut_e){$fill[$i] = $v1}
              }
       }
     }   
   #transform to bedgraph
   my @data_unpack;
   foreach my $i(0..$#fill){push @data_unpack,"$chr\t".($begin+$i)."\t".($begin+$i+1)."\t".$fill[$i]}
   die "data_unpack ne fill:",scalar @data_unpack,"\t",scalar @fill if(scalar @data_unpack != scalar @fill);
   return \@data_unpack;
}



