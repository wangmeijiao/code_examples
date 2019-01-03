
sub readOrthTab(){
#readin orthTab detailed +peak format
  my $file = shift;
  my @order;
  my %tab;
  $/="//";
  open ORTH, $file or die "$!";
  while(<ORTH>){
    chomp;
    next if($_ eq "" ||  $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    @box=@{&boxFilter(\@box)};
    my $id=shift @box;
    next if($id =~/^"/);
    my ($id_real,$comment)=split/[\t ]+/,$id;
    push @order,$id_real;
    #push @order,$id;
    my %temp;
    foreach my $ctrl(@box){
      next if($ctrl eq "-");
      my @line=split/[\t ]+/,$ctrl;
      my ($spec)=split/\|/,$line[0];
      if(!exists $temp{$spec}){
         $temp{$spec}->{"gene"}=$line[0];
         $temp{$spec}->{"peak"}=$line[1];
         #$temp{$spec}->{"data"}=$line[2];
      }else{die "dup $spec at @box"}
      $temp{"comment"}=$comment;
    }
    if(!exists $tab{$id_real}){$tab{$id_real}=\%temp}else{die "dup $id_real at @box"}
    #if(!exists $tab{$id}){$tab{$id}=\%temp}else{die "dup $id at @box"}
    #print Dumper \@box;

  }
  close ORTH;
  $/="\n";
  #return \%tab;
  return (\@order,\%tab);
}


sub readOrthTab_advanced(){
#readin orthTab detailed +peak format
  my $file = shift;
  my @order;
  my %tab;
  $/="//";
  open ORTH, $file or die "$!";
  while(<ORTH>){
    chomp;
    next if($_ eq "" ||  $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    @box=@{&boxFilter(\@box)};
    my $id=shift @box;
    next if($id =~/^"/);
    my ($id_real,$comment)=split/[\t ]+/,$id;
    push @order,$id_real;
    #push @order,$id;
    my %temp;
    foreach my $ctrl(@box){
      next if($ctrl eq "-");
      my @line=split/\t+/,$ctrl;
      my ($spec)=split/\|/,$line[0];
      if(!exists $temp{$spec}){
         $temp{$spec}->{"gene"}=$line[0];
         $temp{$spec}->{"peak"}=$line[1];
         $temp{$spec}->{"data"}=$line[2];
      }else{die "dup $spec at @box"}
      $temp{"comment"}=$comment;
    }
    if(!exists $tab{$id_real}){$tab{$id_real}=\%temp}else{die "dup $id_real at @box"}
    #if(!exists $tab{$id}){$tab{$id}=\%temp}else{die "dup $id at @box"}
    #print Dumper \@box;

  }
  close ORTH;
  $/="\n";
  #return \%tab;
  return (\@order,\%tab);
}

sub readOrthTab_advanced_multi(){
#readin orthTab detailed +peak format, check line number and guess spec order then report
  my $file = shift;
  my @order;
  my %tab;
  my @spec;
  my $cnt_spec;
  $/="//";
  open ORTH, $file or die "$!";
  while(<ORTH>){
    chomp;
    next if($_ eq "" ||  $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    @box=@{&boxFilter(\@box)};
    my $id=shift @box;
    next if($id =~/^"/);
    my ($id_real,$comment)=split/[\t ]+/,$id;
    push @order,$id_real;
    #push @order,$id;
    if(!defined $cnt_spec){$cnt_spec = scalar @box}else{die "spec number diffs with last record at $id_real:\n @box" if($cnt_spec != scalar @box)}
    my %temp;
    my $i=-1;
    foreach my $ctrl(@box){
      $i++;
      if($ctrl eq "-"){if(!defined $spec[$i]){$spec[$i] = "-"};next};
      my @line=split/[\t ]+/,$ctrl;
      #my @line=split/\t+/,$ctrl; #must be \t
      my ($specID)=split/\|/,$line[0];
      if(!defined $spec[$i]){$spec[$i] = $specID}else{if($spec[$i] eq "-"){$spec[$i] = $specID;next};die "spec order diffs with last record at $id_real:\n @box" if($spec[$i] ne $specID)}
      if(!exists $temp{$specID}){
         $temp{$specID}->{"gene"}=$line[0];
         $temp{$specID}->{"peak"}=$line[1];
         $temp{$specID}->{"data"}=$line[2];
      }else{die "dup $specID at @box"}
      $temp{"comment"}=$comment;
    }
    if(!exists $tab{$id_real}){$tab{$id_real}=\%temp}else{die "dup $id_real at @box"}
    #if(!exists $tab{$id}){$tab{$id}=\%temp}else{die "dup $id at @box"}
    #print Dumper \@box;

  }
  close ORTH;
  $/="\n";
  return (\@order,\@spec,\%tab);
}

sub readOrthTab_region(){
#readin orthTab detailed +peak format
  my $file = shift;
  my @order;
  my %tab;
  $/="//";
  open ORTH, $file or die "$!";
  while(<ORTH>){
    chomp;
    next if($_ eq "" ||  $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    @box=@{&boxFilter(\@box)};
    my $id=shift @box;
    next if($id =~/^"/);
    my ($id_real,$comment)=split/[\t ]+/,$id;
    push @order,$id_real;
    #push @order,$id;
    my %temp;
    foreach my $ctrl(@box){
      next if($ctrl eq "-");
      my @line=split/[\t ]+/,$ctrl;
      my $spec=$line[0];
      if(!exists $temp{$spec}){
         $temp{$spec}->{"chr"}=$line[1];
         $temp{$spec}->{"start"}=$line[2];
         $temp{$spec}->{"end"}=$line[3];
         $temp{$spec}->{"geneID"}=$line[4];
      }else{die "dup $spec at @box"}
    }
    if(!exists $tab{$id_real}){$tab{$id_real}=\%temp}else{die "dup $id_real at @box"}
    #if(!exists $tab{$id}){$tab{$id}=\%temp}else{die "dup $id at @box"}
    #print Dumper \@box;

  }
  close ORTH;
  $/="\n";
  #return \%tab;
  return (\@order,\%tab);
}

sub readOrthTab_region_multi(){
#readin orthregion format, check line number and guess spec order then report
  my $file = shift;
  my @order;
  my %tab;
  my @spec;
  my $cnt_spec;
  $/="//";
  open ORTH, $file or die "$!";
  while(<ORTH>){
    chomp;
    next if($_ eq "" ||  $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    @box=@{&boxFilter(\@box)};
    my $id=shift @box;
    next if($id =~/^"/);
    my ($id_real,$comment)=split/[\t ]+/,$id;
    push @order,$id_real;
    #push @order,$id;
    if(!defined $cnt_spec){$cnt_spec = scalar @box}else{die "spec number diffs with last record at $id_real\n  @box" if($cnt_spec != scalar @box)}
    my %temp;
    my $i=-1;
    foreach my $ctrl(@box){
      $i++;
      if($ctrl eq "-"){if(!defined $spec[$i]){$spec[$i] = "-"};next};
      my @line=split/[\t ]+/,$ctrl;
      if(!defined $spec[$i]){$spec[$i] = $line[0]}else{if($spec[$i] eq "-"){$spec[$i] = $line[0];next};die "spec order diffs with last record at $id_real\n @box" if($spec[$i] ne $line[0])}
      if(!exists $temp{$line[0]}){
         $temp{$line[0]}->{"chr"}=$line[1];
         $temp{$line[0]}->{"start"}=$line[2];
         $temp{$line[0]}->{"end"}=$line[3];
         $temp{$line[0]}->{"geneID"}=$line[4];
      }else{die "dup spec $line[0] at @box"}
    }
    if(!exists $tab{$id_real}){$tab{$id_real}=\%temp}else{die "dup $id_real at @box"}
    #if(!exists $tab{$id}){$tab{$id}=\%temp}else{die "dup $id at @box"}
    #print Dumper \@box;
  }
  close ORTH;
  $/="\n";
  return (\@order,\@spec,\%tab);
}


sub boxFilter(){
  #filter empty or other illegal elements from array, keep the order.
  my $idx=shift;
  my @filtered;
  foreach(@{$idx}){
    if(!defined $_ || $_ eq "" || $_=~/^\s+$/ ||$_=~/^#/){}else{push @filtered,$_}
  }
  return \@filtered;
}


