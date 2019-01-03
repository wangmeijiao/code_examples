use strict;
use warnings;



&mafRead($ARGV[0],$ARGV[1]);
  #perl maf2mfa.pl tba.maf "japoChr03-zs97Chr03-glabChr03-meriChr03" > tba.maf.mfa







##sub

sub mafRead(){
  my ($file,$order) = @_;
  my @order = split/-/,$order;
  die "order empty: $order" if(scalar @order == 0);
  my %maf;
  my $cnt=0;
  open MAF, $file or die "$!";
  $/ = "\n\n";
  while(<MAF>){
    chomp;
    my @box=split /\n/,$_;
    my $block_score;
    my %block;#specises name stand for the key, lost order information 
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
    print "#maf_group $cnt length:$block{'alignLen'} score:$block{'score'}\n";
    foreach my $spec(@order){if(exists $block{$spec}){print ">$spec $block{$spec}->[1]-$block{$spec}->[2] $block{$spec}->[3]\n",&faFormat($block{$spec}->[4],100)}}
    print "\n\n"; 
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks\n";
  return \%maf;
}


sub faFormat(){
  my ($seq,$num)=@_;
  my @seq=split//,$seq;
  my ($string,$cnt);
  foreach(@seq){
   $cnt++;
   if($cnt%$num==0){$string.=$_."\n"}else{$string.=$_}
  }
  if($cnt%$num!=0){return $string."\n"}else{return $string}
}



