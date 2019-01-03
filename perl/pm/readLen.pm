sub readLen(){
    my $file=shift @_;
    my %len;
    open LEN, $file or die "$!";
    while(<LEN>){
      chomp;
      next if ($_ eq "" || $_=~/^#/);
      my ($chr,$len)=split/\t/,$_;
      if(!exists $len{$chr}){$len{$chr}=$len}else{die "dup chr $chr\n"}
    }
    close LEN;
    return \%len;
}


sub getLenAll(){
  my ($info,$archi)=@_;
  my @len;
  foreach my $item(@{$archi}){
   my($species,$chr)=split/:/,$item; #don't use $_
   die "species empty in archi @{$archi}" if($species eq "" || $species eq "NA" || $species eq "NONE");
   my $sizef=$info->{$species}->{"chrSize"};
   if ( -e $sizef ){}else{print "\n$species\t$sizef\n"; die "not exists $sizef\n"}
   #if($chr eq "NONE"){$chr = $refChr}
   #if($chr eq "NONE"){$chr = $refChr}
   my $len=&readLen($sizef)->{$chr} ;  # here $_ will be unpredictable when call another sub in a sub
   push @len, $len;
  }
  return \@len;
}


sub getLenAll_hash(){
  my ($info,$archi)=@_;
  my %len;
  foreach my $species(@{$archi}){
   my $sizef=$info->{$species}->{"chrSize"};
   if ( -e $sizef ){}else{print "\n$species\t$sizef\n"; die "not exists $sizef\n"}
   my $len=&readLen($sizef);
   $len{$species} =  $len;
  }
  return \%len;
}


sub readLen_list(){
    my $list=shift;
    my %len;
    foreach my $spec(keys %{$list}){
      open LEN, $list->{$spec} or die "$!";
      while(<LEN>){
        chomp;
        next if($_ eq "" || $_=~/^#/);
        my ($chr,$length) = split/[\t ]+/,$_;
        if(!exists $len{$spec}->{$chr}){
          $len{$spec}->{$chr}=$length;
        }else{die "dup chr $chr in file $list->{$spec}\n"}
      }

      close LEN;
    }
    return \%len;
}



