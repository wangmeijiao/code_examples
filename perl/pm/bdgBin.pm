

sub bdgBin(){
  #bin the bdg (unpack first if packed format, coord must continuous) with bin length win
  #output bdg 
  my ($bdg,$win) = @_;
  $win||=1000; #1kb
  my $cnt;
  my @data; #win lines
  my @bdg_bin;
  foreach my $ctrl(@$bdg){
     $cnt++;
     #print STDERR "#" if($cnt % 200000 == 0);
     if($cnt % $win == 0 || $ctrl eq $bdg->[-1]){  #when reach last but not $cnt%win == 0
      push @data,$ctrl; #the last line of win
      #print STDERR "eof! $_\n" if(eof(STDIN));
      #print  "cnt:$cnt,number of array:",scalar @data,"\n";
      my $sum;
      my $cnt_in;
      for (my $i=0;$i<=$#data;$i++){
         #check that neighboring two lines must be ascendingly linked  and continuous with step 1
         my($chr1,$s1,$e1,$v1)=split/[\t ]+/,$data[$i];
         die "$s1 >= $e1 at $data[$i]" if($s1 >= $e1);
         die "$s1+1 != $e1 at $data[$i]" if($s1+1 != $e1);
         if($i == $#data){
           $sum+=$v1;
           $cnt_in++;
           last;
         }
         my($chr2,$s2,$e2,$v2)=split/[\t ]+/,$data[$i+1];
         die "$s2 >= $e2 at $data[$i+1]" if($s2 >= $e2);
         die "$s2+1 != $e2 at $data[$i+1]" if($s2+1 != $e2);
         die "$chr1 ne $chr2 at $data[$i] and $data[$i+1]" if($chr1 ne $chr2);
         die "$e1 > $s2 at $data[$i] and $data[$i+1]" if($e1 > $s2); #must in ascending order
         if($e1 != $s2){die "\nerr: continuity fails at $data[$i]  and  $data[$i+1]\n"};
         $sum+=$v1;
         $cnt_in++;
      }#outter for
      print STDERR "\nwarning, cnt_in ne win: $cnt_in != $win between $data[0] and $data[-1]\n" if($cnt_in != $win  && !($ctrl eq $bdg->[-1]));
      my $ave = sprintf "%.3f", $sum/$cnt_in;
      my($chr_f,$s_f,$e_f,$v_f)=split/[\t ]+/,$data[0]; 
      my($chr_l,$s_l,$e_l,$v_l)=split/[\t ]+/,$data[-1]; 
      $s_f-=1;
      $e_l-=1;
      #print "$chr_f\t$s_f\t$e_l\t$ave\n";
      push @bdg_bin,"$chr_f\t$s_f\t$e_l\t$ave";
      @data=();
    }else{ push @data, $ctrl}
  }
  
  #print STDERR "\n";
  return \@bdg_bin;
}#sub end




