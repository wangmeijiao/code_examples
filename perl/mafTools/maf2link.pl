use strict;
use warnings;
use Data::Dumper;


  my $maf = &mafRead($ARGV[0]);
  print STDERR "readin done\n";
  #print Dumper $maf;exit;

  my @species = ("japoChr01","OglabChr01","OpuncChr01","bracChr01","LperrChr01");
  my $refSpec = "japoChr01";
  open JG, ">japoChr01.OglabChr01.links" or die "$!";  
  open JP, ">japoChr01.OpuncChr01.links" or die "$!";  
  open JB, ">japoChr01.bracChr01.links" or die "$!";  
  open JL, ">japoChr01.LperrChr01.links" or die "$!";  

  foreach my $id(sort {$a <=> $b} keys %$maf){
       if(!exists $maf->{$id}->{$refSpec}){print "ref $refSpec not exists at maf group $id\n"; next}
       if($maf->{$id}->{"score"} <=100 ){print "score $maf->{$id}->{'score'}<=100 at maf group $id\n"; next}
       if($maf->{$id}->{"alignLen"} <=10 ){print "alignLen $maf->{$id}->{'alignLen'}<=10 at maf group $id\n"; next}
       if(exists $maf->{$id}->{"OglabChr01"}){print JG "$maf->{$id}->{$refSpec}->[0] $maf->{$id}->{$refSpec}->[1] $maf->{$id}->{$refSpec}->[2] $maf->{$id}->{'OglabChr01'}->[0] $maf->{$id}->{'OglabChr01'}->[1] $maf->{$id}->{'OglabChr01'}->[2] color=red\n"}    
       if(exists $maf->{$id}->{"OpuncChr01"}){print JP "$maf->{$id}->{$refSpec}->[0] $maf->{$id}->{$refSpec}->[1] $maf->{$id}->{$refSpec}->[2] $maf->{$id}->{'OpuncChr01'}->[0] $maf->{$id}->{'OpuncChr01'}->[1] $maf->{$id}->{'OpuncChr01'}->[2] color=blue\n"}    
       if(exists $maf->{$id}->{"bracChr01"}){print JB "$maf->{$id}->{$refSpec}->[0] $maf->{$id}->{$refSpec}->[1] $maf->{$id}->{$refSpec}->[2] $maf->{$id}->{'bracChr01'}->[0] $maf->{$id}->{'bracChr01'}->[1] $maf->{$id}->{'bracChr01'}->[2] color=orange\n"}    
       if(exists $maf->{$id}->{"LperrChr01"}){print JL "$maf->{$id}->{$refSpec}->[0] $maf->{$id}->{$refSpec}->[1] $maf->{$id}->{$refSpec}->[2] $maf->{$id}->{'LperrChr01'}->[0] $maf->{$id}->{'LperrChr01'}->[1] $maf->{$id}->{'LperrChr01'}->[2] color=black\n"}    



  }

  close JG;
  close JP;
  close JB;
  close JL;



###sub

sub mafRead(){
  my $file = shift;
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
  }
  close MAF;
  $/ = "\n";
  print STDERR "\ntotal $cnt maf blocks\n";
  return \%maf;
}







