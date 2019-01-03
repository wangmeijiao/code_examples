sub permutation_list(){ #permutation and combination elements
   my $index=shift;
   my @elements=@{$index};
   my @elements_map=@elements;
   my @list;
   foreach my $ctrl_out(@elements){
     shift @elements_map;
     foreach my $ctrl_in(@elements_map){
       my $pair;
       $pair=$ctrl_out.'-'.$ctrl_in;
       push @list,$pair;
     }#inner foreach end here 
   }#outter foreach end here
   return @list;
}

sub factorial(){ #jie chen calculation
   my $k=shift;
   my $f=1;
   for (1..$k){$f*=$_}
   return $f;
}

sub choose(){ #zu he calculation
#  c(k,n)=p(k,n)/factorial(k)
#         =factorial(n)/(factorial(n-k)*factorial(k))
   my ($k,$n)=@_;
   return &factorial($n)/(&factorial($n-$k)*&factorial($k));
}

sub permutation(){ # pai lie calculation
   my ($k,$n)=@_;
   return &choose($k, $n) * &factorial($k);
}


sub minMax(){
  my $index=shift;
  my ($min,$max)=($index->[0],$index->[0]);
  foreach (@{$index}){
    if($min>$_){$min=$_}
    if($max<$_){$max=$_}
  }
  #print ("$min,$max\n");
  return ($min,$max);
}#end of minMax sub

sub meanSd(){
  #use n-1 methods, to describe the overall dispersivity relative to mean 
  # sqrt(n*total(x^2)-total(x)^2/n/(n-1)) #excel, one step
  # sqrt( total((x-mean)^2)/(n-1) ) #multi steps
  my @list=@{shift @_};
  my $n=scalar(@list);
  my ($total,$total_squa);
  foreach(@list){
   $total+=$_;
   $total_squa+=$_*$_;
  }
  my $mean=$total/$n;
  my $sum;
  foreach(@list){
   $sum+=($_-$mean)*($_-$mean);
  }
  return ($mean,sqrt($sum/($n-1)));#n-1, one step #use this
  #return ($mean,sqrt(($total_squa*$n-$total*$total)/($n*($n-1)))); #multi-step?
}


sub sum(){
   my $data =shift;
   die "empty data" if(scalar @$data == 0);
   my $sum;
   foreach (@$data){
     die"NULL found @$data" if(!defined $_ || $_ eq "");
     $sum+=$_;
   }   
   return $sum;
}


sub mean(){
   my $data =shift;
   die "empty data" if(scalar @$data == 0);
   my $sum;
   my $cnt;
   foreach (@$data){
     die"NULL found @$data" if(!defined $_ || $_ eq "");
     $sum+=$_;
     $cnt++;
   }
   return $sum/$cnt;
}



sub sd(){
   #the same with meanSd()
   my $data =shift;
   if(scalar @$data == 1){return 0}
   my $ave = &mean($data);
   my $sqtotal = 0;
   foreach(@$data) {
      $sqtotal += ($ave-$_) ** 2;
   }
   my $std = ($sqtotal / (@$data-1)) ** 0.5;
   return $std;
}


sub quanTile(){
#matlab defined quantile
#let
#w = (N+1)/4, 2(N+1)/4, 3(N+1)/4
#y = the truncated integer value of w
#z = the fraction component of w that was truncated away
#Q1 = x(y) + z(x(y+1) - x(y))
#Note: when w is an integer, y = w, z = 0, and Q1 = x(y)
  my @list=@{shift @_};
  my $n=scalar(@list);
  @list=sort{$a<=>$b} @list;
  #print STDERR "@list\n";
  #my $min=$list[0];
  #my $max=$list[$#list];
  my ($Q1,$Q2,$Q3,$i,$f);
  ($i,$f)=&trunct(($n+1)*0.25);
  $Q1=$list[$i-1]+$f*($list[$i-1+1]-$list[$i-1]);
  #print"$n,$i,$f,$Q1\n";
  ($i,$f)=&trunct (($n+1)*0.5);
  $Q2=$list[$i-1]+$f*($list[$i-1+1]-$list[$i-1]); #also median
  #print"$n,$i,$f,$Q2\n";
  ($i,$f)=&trunct( ($n+1)*0.75);
  $Q3=$list[$i-1]+$f*($list[$i-1+1]-$list[$i-1]);
  #print"$n,$i,$f,$Q3\n";
  return ($Q1,$Q2,$Q3)
}

sub trunct(){
 my $number=shift @_;#print "$number\n";
 my $integer=int ($number);
 my $fraction=$number-$integer;
 return ($integer,$fraction);
}


sub ceiling(){
  my $num = shift;
  return (int($num) + 1);

}


sub round(){
   my $n=shift;
   my $i=int $n;
   my $d=$n-$i;
   if($d>0.6){return $i+1}else{return $i}
}


sub log2(){
  my $v = shift;
  die "illegal number $v" if($v <=0);
  return log($v)/log(2);
}

sub log10(){
  my $v = shift;
  die "illegal number $v" if($v <=0);
  return log($v)/log(10);
}


sub isNumeric(){
  my $idx=shift;
  my $flag=1;
  foreach my $ctrl (@{$idx}){
    #signed    
    $ctrl=~s/^-//;
    #decimals
     if($ctrl =~/[^0-9\.]/){$flag = 0}

    #scientific notation

  }
  return $flag;
}




sub notNA(){
    my $idx=shift;
    my $flag=1;
    foreach my $ctrl (@{$idx}){if($ctrl eq "NA" || $ctrl eq "NULL" || $ctrl eq "Na" || $ctrl eq "n/a" || $ctrl eq "N/A"){$flag
    return $flag;
}



sub readTable(){
   my $file=shift;
   my $head=shift;
   my @header;
   my @data;
   my $cnt;
   my $ncol;
   my ($min,$max);
   open DATA, $file or die "$!";
   while(<DATA>){
     chomp;
     next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
     $cnt++;
     if($cnt==1){
       if($head eq "T"){
           @header=split/[\t ]+/,$_;
           $ncol=scalar @header;
           next;
       }elsif($head eq "F"){
           my @box=split/[\t ]+/,$_;
           $ncol=scalar @box;
           if(&notNA(\@box) == 0){die "has NA at $_"}
           if(&isNumeric(\@box) == 0){die "not numeric at $_"}
           push @data,$_;
           ($min,$max)=($box[0],$box[0]);
         }else{die "unknow header swith $head"}
     }
     my @box_line=split/[\t ]+/,$_;
     my $ncol_line=scalar @box_line;
     if(!defined $ncol || $ncol_line != $ncol){die "ncol err $ncol != $ncol_line at $_"}
     if(&notNA(\@box_line) == 0){die "has NA at $_"}
     if(&isNumeric(\@box_line) == 0){die "not numeric at $_"}
     push @data,$_;
     if(! defined $min && !defined $max){
       ($min,$max)=($box_line[0],$box_line[0])
     }elsif(defined $min && defined $max){
         foreach my $ctrl(@box_line){
           if($ctrl < $min){$min = $ctrl}
           if($ctrl > $max){$max = $ctrl}
         }
       }else{die "err min max $min , $max"}
   }#while end
   close DATA;
   if($head eq "T"){$cnt-=1}
   return (\@data,\@header,$cnt,$ncol,$min,$max);

}



