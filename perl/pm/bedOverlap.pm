
sub bedOverlap(){
  my ($file,$region)=@_;
  my @hits;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  open BED, $file or die "$!";
  while(<BED>){
   chomp;
   next if($_ eq ""|$_=~/^#/ || $_=~/^\s+$/);
   my @box=split/[\t ]+/,$_;
   my ($chr_in,$start_in,$end_in)=@box;
   next if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{  
     #cut ends
     if($box[1] < $start){ $box[1] = $start}
     if($box[2] > $end){ $box[2] = $end}
     push @hits,join("\t",@box); 
   }
  }
  close BED;
  return \@hits;
}

sub bedOverlap_array(){
  my ($bed,$region)=@_;
  my @hits;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  foreach(@$bed){
   my @box=split/[\t ]+/,$_;
   my ($chr_in,$start_in,$end_in)=@box;
   next if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{
     #cut ends
     if($box[1] < $start){ $box[1] = $start}
     if($box[2] > $end){ $box[2] = $end}
     push @hits,join("\t",@box);
   }
  }
  return \@hits;
}

sub bedOverlap_hash(){
  my ($bed,$region)=@_;
  my @hits;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  foreach(keys %$bed){
   my @box=split/[\t ]+/,$bed->{$_};
   my ($chr_in,$start_in,$end_in)=@box;
   next if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{
     #cut ends
     #if($box[1] < $start){ $box[1] = $start}
     #if($box[2] > $end){ $box[2] = $end}
     push @hits,join("\t",@box);
   }
  }
  return \@hits;
}


sub bedNotOverlap(){
  my ($bed,$region,$id_self)=@_;
  $id_self=~s/\.\d+$//;
  my $flag = 1;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  foreach my $id(keys %$bed){
   my @box=split/[\t ]+/,$bed->{$id};
   my ($chr_in,$start_in,$end_in)=@box;
   die "chr diffs: $chr ne $chr_in at @box" if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{ if($id ne $id_self){$flag = 0;last }   }
  }
  return $flag;
}

sub bedIsOverlap(){
  my ($bed,$region)=@_;
  my $flag = 0;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  foreach my $ctrl(@$bed){
   my @box=split/[\t ]+/,$ctrl;
   my ($chr_in,$start_in,$end_in)=@box;
   die "chr diffs: $chr ne $chr_in at @box" if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{ $flag = 1}
  }
  return $flag;
}



sub bedOverlap_transform(){
  my ($file,$region)=@_;
  my @hits;
  my ($chr,$start,$end)=split/:|-/,$region;
  die "start>=end :$start >= $end" if($start>=$end);
  open BED, $file or die "$!";
  while(<BED>){
   chomp;
   next if($_ eq ""|$_=~/^#/ || $_=~/^\s+$/);
   my @box=split/[\t ]+/,$_;
   my ($chr_in,$start_in,$end_in)=@box;
   next if($chr ne $chr_in);
   die"start>=end :$start_in >= $end_in or less than 3 columns\n" if($start_in >= $end_in || scalar @box <3);
   if($end_in <= $start || $start_in >= $end ){}else{
     #transform and cut ends
     $box[1]-=$start;
     if($box[1] < 0 ){$box[1] = 0}
     $box[2]-=$start;
     if($box[2] <= 0 ){die "wrong at $start_in and $end_in after transform: transformed end <= 0"}
     push @hits,join("\t",@box);
   }
  }
  close BED;
  return \@hits;
}




