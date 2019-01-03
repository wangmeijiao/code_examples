


sub segQuery_hash(){  #query a given segment in a table, report all overlap ones
   my ($table,$seg)=@_;
   my %segs_hit;
   my ($s,$e)=split/-/,$seg;
   die "$s >= $e at $seg" if($s >= $e);
   foreach my $id(sort keys %{$table}){
     my (undef,$start,$end)=split/\t/,$table->{$id};
=pod
     if($start >= $s && $end <= $e){
       $segs_hit{$id}=$table->{$id};
      }elsif($start < $s && $end > $s && $end <= $e){
        $segs_hit{$id}=$table->{$id}
       }elsif($end > $e && $start > $s && $start < $e){
          $segs_hit{$id}=$table->{$id}
         }elsif($start <= $s && $end >= $e){$segs_hit{$id}=$table->{$id}}
=cut
     if($e < $start || $end < $s){}else{$segs_hit{$id}=$table->{$id} }  #segment overlap 

   }
   return \%segs_hit;
}



