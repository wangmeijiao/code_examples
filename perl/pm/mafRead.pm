
#mafRead(), mafRead_plus , mafRead_fordraw
#use mafRead_plus;

sub mafRead(){
  my $file = shift;
  my %maf;
  my $cnt=0;
  open MAF, $file or die "$!";
  $/ = "\n\n";
  while(<MAF>){
    chomp;
    last if ($_=~/##eof/); #modified on 28, Dec 2016
    my @box=split /\n/,$_;
    my $block_score;
    my %block;#specises name stand for the key, lost order information 
    foreach my $line(@box){
      #exit if ($line=~/##eof/);
      #last if ($line=~/##eof/); #modified on 12, May 2013
      next if ($line=~/#/);
      if($line=~/^a score/){
          my @tmp=split/=/,$line;$block_score=int($tmp[1]);
          if(!exists $block{"score"}){$block{"score"} = $block_score}else{die "dup score line: $line"} 
        }elsif($line=~/^s /){
              my @tmp=split/ +/,$line;
              my $name=$tmp[1];
              my $start=$tmp[2];
              my $len=$tmp[3];
              my $strand=$tmp[4];
              my $chr_len=$tmp[5];
              my $seq_ali=$tmp[6];
              my $seq_ori=$seq_ali; $seq_ori=~s/-//g;
              my $end=$start+$len;
              my @dataline=($name,$start,$end,$strand,$seq_ali,$seq_ori);
              die "length $len diff with seq_rmgap at block $cnt" if($len != length $seq_ori);
              my $alignLen = length $seq_ali;
              if(exists $block{$name}){die"duplicated species in maf block $cnt\n"}else{$block{$name}=\@dataline}
              if(exists $block{"alignLen"}){die "align length diff at block $cnt" if($block{"alignLen"} != $alignLen)}else{$block{"alignLen"} = $alignLen}
         }else{die "\nunknow maf line:$line\n"}
    }#foreach end here 
    $maf{$cnt} = \%block;
    $cnt++;
    #last if($cnt == 10000);
    print STDERR "#" if($cnt % 100000 ==0);
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks\n";
  return \%maf;
}



sub mafRead_plus(){
  #read and check maf file
  #s japoChr01      1554 120 + 43268879 CTATCTAGGCATCCATCCGATATTTGGAGTATGGAGGAGAAAAACAGTGCTCCAGCAGAGTCTCCATCACATGCTTCATTTTTGG
  #s OpuncChr01 30196520 120 + 46096743 CTATCCCACCCTTCATATGAGAAATAGAGTATGTAAGCAAAAAAAGAGACTCCAGCAGACACTCCAAAATATCCTCCAAAAATAG
  #s LperrChr01  1183910 106 + 32922458 CCATCCCATACTCCATCCTATATTTGGTATATATGGAAGGAAAAATGGGCTCCAGTA----------TATATACCCATAAACTAG
  my $file = shift;
  my %maf;
  my $cnt=0;
  open MAF, $file or die "$!";
  $/ = "\n\n";
  while(<MAF>){
    chomp;
    last if ($_=~/##eof/); #modified on 28, Dec 2016
    my @box=split /\n/,$_;
    my $block_score;
    my %block;
    my @order; 
    my $degree = 0;
    foreach my $line(@box){
      #exit if ($line=~/##eof/);
      #last if ($line=~/##eof/); #modified on 12, May 2013
      next if ($line=~/#/);
      if($line=~/^a score/){
          my @tmp=split/=/,$line;$block_score=int($tmp[1]);
          if(!exists $block{"score"}){$block{"score"} = $block_score}else{die "dup score line: $line"}
        }elsif($line=~/^s /){
              my @tmp=split/ +/,$line;
              my ($name,$start,$len,$strand,$chr_len,$seq_ali) = ($tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6]); 
              push @order, $name;
              $degree++;
              my $seq_ori=$seq_ali; $seq_ori=~s/-//g;
              my $end=$start+$len;
              my @dataline=($name,$start,$end,$len,$strand,$chr_len,$seq_ali,$seq_ori);
              die "length $len diff with seq_rmgap at block $cnt" if($len != length $seq_ori);
              my $alignLen = length $seq_ali;
              if(exists $block{$name}){die"duplicated species in maf block $cnt\n"}else{$block{$name}=\@dataline}
              if(exists $block{"alignLen"}){die "align length diff at block $cnt" if($block{"alignLen"} != $alignLen)}else{$block{"alignLen"} = $alignLen}
         }else{print STDERR "\nomit line at group $cnt\n: unknow maf line type:$line\n"}
    }#foreach end here 
    if($degree == 0){die "degree eq 0 at group $cnt: $_"}else{$block{"degree"} = $degree}
    $block{"order"} = \@order;   
    $maf{$cnt} = \%block;
    $cnt++;
    #last if($cnt == 10000);
    print STDERR "#" if($cnt % 5000 == 0);
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks read in\n";
  return \%maf;
}

sub mafRead_draw(){
  # guess and return species_num, species_order&species_chrs
  my $file = shift;
  my %maf;
  my ($maxLen,%len,@order,%range);
  my $cnt=0;
  open MAF, $file or die "$!";
  $/ = "\n\n";
  while(<MAF>){
    chomp;
    my @box=split /\n/,$_;
    my $block_score;
    my %block;#specises name stand for the key, lost order information 
    my $cnt_spec=-1;
    foreach my $line(@box){
      #exit if ($line=~/##eof/);
      last if ($line=~/##eof/); #modified on 12, May 2013
      next if ($line=~/#/);
      if($line=~/^a score/){
          my @tmp=split/=/,$line;$block_score=int($tmp[1]);
          if(!exists $block{"score"}){$block{"score"} = $block_score}else{die "dup score line: $line"}
        }elsif($line=~/^s /){
              my @tmp=split/ +/,$line;
              my $name=$tmp[1];
              $cnt_spec++;
              if(!defined $order[$cnt_spec] || $order[$cnt_spec] ne $name){$order[$cnt_spec] = $name}  
              my $start=$tmp[2];
              if(!exists $range{$name}->{"begin"}){$range{$name}->{"begin"} = $start}else{if($range{$name}->{"begin"} > $start){$range{$name}->{"begin"} = $start}}
              my $len=$tmp[3];
              my $strand=$tmp[4];
              my $chr_len=$tmp[5];
              if(!defined $maxLen){$maxLen = $chr_len}else{if($maxLen < $chr_len){$maxLen = $chr_len}}
              if (!exists $len{$name}){$len{$name} = $chr_len}else{if($len{$name} != $chr_len){die "len of chr of $name diffs:$len{$name} != $chr_len at block $cnt:$line"}}
              my $seq_ali=$tmp[6];
              my $seq_ori=$seq_ali; $seq_ori=~s/-//g;
              my $end=$start+$len;
              if(!exists $range{$name}->{"end"}){$range{$name}->{"end"} = $end}else{if($range{$name}->{"end"} < $end){$range{$name}->{"end"} = $end}}
              my @dataline=($name,$start,$end,$strand,$seq_ali,$seq_ori);
              die "length $len diff with seq_rmgap at block $cnt" if($len != length $seq_ori);
              my $alignLen = length $seq_ali;
              if(exists $block{$name}){die"duplicated species in maf block $cnt\n"}else{$block{$name}=\@dataline}
              if(exists $block{"alignLen"}){die "align length diff at block $cnt" if($block{"alignLen"} != $alignLen)}else{$block{"alignLen"} = $alignLen}
         }else{die "\nunknow maf line:$line:$_\n"}
    }#foreach end here 
    $maf{$cnt} = \%block;
    $cnt++;
    #last if($cnt == 10000);
    print STDERR "#" if($cnt % 100000 ==0);
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks\n";
  my $n_spec = scalar @order;
  return \%maf,$maxLen,\%len,$n_spec,\@order,\%range;
}


