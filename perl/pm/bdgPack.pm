
sub bdgPack(){
   #bedgraph coords must be in properly ascending order (no need to be continuously, warning when continuity fails)
   my $data=shift @_;
   my @data_packed;
   for (my $i=0;$i<=$#{$data}-1;$i++){
    #if($i%100000 == 0){print STDERR "#"}
    my($chr,$s,$e,$v)=split/[\t ]+/,$data->[$i];
    die "$s >= $e at $data->[$i]" if($s >= $e);
    my $cnt=0;
    for( my $j=$i+1;$j<=$#{$data};$j++){
       my($chr1,$s1,$e1,$v1)=split/[\t ]+/,$data->[$j-1];
       die "$s1 >= $e1 at $data->[$j-1]" if($s1 >= $e1);
       my($chr2,$s2,$e2,$v2)=split/[\t ]+/,$data->[$j];
       die "$s2 >= $e2 at $data->[$j]" if($s2 >= $e2);
       die "$chr1 ne $chr2 at $data->[$j-1] and $data->[$j]" if($chr1 ne $chr2);
       die "$e1 > $s2 at $data->[$j-1] and $data->[$j]" if($e1 > $s2); #must in ascending order
       if($e1 != $s2){print STDERR "warning: continuity fails at $data->[$j-1]  and  $data->[$j]"};
       if($v1 == $v2 && $e1 == $s2){
         $cnt++;
         if($j==$#{$data}){push @data_packed,"$chr\t$s\t$e2\t$v";$i=$#{$data}}
         next;
       }else{
             if($cnt == 0){
               push @data_packed,"$chr\t$s\t$e\t$v";
               if($i==$#{$data}-1){push @data_packed,"$chr2\t$s2\t$e2\t$v2"}
             }else{push @data_packed,"$chr\t$s\t$e1\t$v"}
             $i=$j;
             $i--;
             $cnt=0;
             last;
        }
     }#inner for
    }#outter for
    #print STDERR "\n";
   return \@data_packed;
}




