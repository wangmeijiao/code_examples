use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#readin mfa alignment (>=2 seqs), scan %identity for each bp or by sliding window
#output densi_like data for R plot

  #usage: perl thisfile.pl -mfa  in.mfa  -mode  depth|window #don't use window mode 

  my ($mode,$mfa);  

  GetOptions("mode=s",\$mode,"mfa=s",\$mfa);


  my ($align,$align_order,$num_align) = &mfaRead($mfa);
  my ($flag,$len)=&checkAlignLenAll($align);die "length diffs" if($flag == 0);
  #print Dumper $align; exit;


  my $densi_score;
  if($mode eq "depth"){
    (undef, $densi_score, undef)=&getSimilarDensi($align,$len,"perc");
    #my ($flag_gaps, $densi_score, $cons_seq)=&getSimilarDensi($align,$len,"score");
    #print Dumper $densi_score;
    die "score len err" if(scalar @$densi_score != $len);
  }elsif($mode eq "window"){
      my $window = 5;
      #my $window = int (0.1*$len);
      $densi_score=&getSimilarDensi_win($align,$len,"perc",$window);
      #print Dumper $densi_score; exit;
     }else{die "unknow mode: $mode, use depth|window"}



  #report
  print "##msa $num_align $len\n";
  for(0..$#$densi_score){print sprintf("%.1f",$densi_score->[$_]),"\n"}



###sub



sub mfaRead(){
  my $file = shift;
  $/=">";
  my %align;
  my $num_align;
  my @align_order;
  open MFA, $file or die "$1";
  while(<MFA>){
    chomp;
    next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp = split/[\t ]+/,$id;
    my $id_real = $temp[0];
    push @align_order,$id_real;
    my $align=join"",@box;
    if(!exists $align{$id_real}){$align{$id_real}=$align; $num_align++}else{die "dup alignment $id_real\n"}
  }
  close MFA;
  $/="\n";
  return (\%align,\@align_order,$num_align);
}



sub checkAlignLenAll(){
  my $index=shift;
  my $len;
  my $flag=1;
  foreach my $id(keys %{$index}){
     if(!defined $len){$len=length($index->{$id})}else{if($len != length($index->{$id}) ){$flag=0}}
  }
  return ($flag,$len);
}


sub getSimilarDensi(){
  my ($index,$len,$type)=@_;
  my @densi; #use score or perc;
  my @cons;
  my @flag_gaps; #len == aligned_len ;
  my %align;
  my @ids;
  foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align{$id}=\@box;
      push @ids,$id;
  }
  for my $i (0.. $len-1){
     my @elements;
     foreach my $id(@ids){
       push @elements, $align{$id}->[$i];
     }
     my($flag, $score, $perc, $cons)=&checkIdentity_densi(\@elements);
     if($type eq "score"){push @densi, $score}elsif($type eq "perc"){push @densi, $perc}else{die "don't understand densi type $type"}
     push @cons, $cons;
     push @flag_gaps,$flag;
  }
  return(\@flag_gaps, \@densi, \@cons) ;
}


sub getSimilarDensi_win(){
  my ($index,$len,$type,$win)=@_;
  #my @densi; #use score or perc;
  #my @cons;
  #my @flag_gaps; #len == aligned_len ;
  my @densi_win;
  my %align;
  my @ids;
  foreach my $id(keys %{$index}){
      my @box=split//,$index->{$id};
      $align{$id}=\@box;
      push @ids,$id;
  }
  for my $i ($win.. $len-1-$win){
    my @densi; #use score or perc;
    #my @cons;
    #my @flag_gaps; #len == aligned_len ;
    for my $j($i-$win..$i+$win){
        #my @densi; #use score or perc;
        #my @cons;
        #my @flag_gaps; #len == aligned_len ;
        my @elements;
        foreach my $id(@ids){
          push @elements, $align{$id}->[$j];
        }
        my($flag, $score, $perc, $cons)=&checkIdentity_densi(\@elements);
        if($type eq "score"){push @densi, $score}elsif($type eq "perc"){push @densi, $perc}else{die "don't understand densi type $type"}
        #print "$flag, $score, $perc, $cons\n";
        #push @cons, $cons;
        #push @flag_gaps,$flag;
    }
    push @densi_win, &ave(\@densi);
  }
  return(\@densi_win) ;
}



sub checkIdentity_densi(){
   #count -> calculate perc%, similarity score (0-num_seq) and output a cons base 
   my @elements=@{shift @_};
   my $num=scalar @elements;
   my ($score,$flag,$cons,$perc);
   $flag=0;
   my %cnt; #for gap, nucleotide base and aa and others
   foreach my $ele (@elements){
      if($ele=~/[^-GPAVLIMCFYWHKRQNEDST]/i){die "illegal character found in $ele"} #nucleotide included
      $ele = uc($ele);  #in case of softmasked repeat
      $cnt{$ele}++;
      if($ele eq "-"){$flag=1}
   }
   #print Dumper \%cnt;
   my @sorted=sort{ $cnt{$b} <=> $cnt{$a} }keys %cnt;
   #print "@sorted\n";
   #if($sorted[0] ne "-"){$score = $cnt{$sorted[0]}; $cons = $sorted[0]}else{$score = 0; $cons = $sorted[0]}
   if($sorted[0] ne "-"){$score = $cnt{$sorted[0]}; $cons = $sorted[0]}else{$score = $cnt{$sorted[1]}; $cons = $sorted[1]}
   $perc = sprintf "%.3f",100*$score/$num;
   #if($perc <= 50){$score=0; $cons=" "}



   return ($flag, $score, $perc, $cons);
}


sub ave(){
  my $ind = shift;
  my $cnt=0;
  my $sum=0;
  die "empty: @$ind" if(scalar @$ind == 0);
  foreach (@$ind){ $cnt++;$sum+=$_}
  die "empty cnt: @$ind" if($cnt == 0);
  my $ave = $sum/$cnt;
  return $ave;

}




