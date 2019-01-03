
sub readtable
{
my ($file,$col)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[$col]}=1;
}
close IN;
return \%hash;
}



sub sliding
{
my ($line,$win,$step,$hash,$meth)=@_;
my @unit=split("\t",$line);
my $bin=int ($unit[1]/$step) + 1;
for(my $i=$bin;$i>0;$i--){
   my $start=($i-1)*$step;
   my $end=$start+$win;
   #print "$i\t$unit[1]\t$start\t$end\n";
   if ($unit[1] >= $start and $unit[1] < $end) {
            $hash->{$start}++;
            if ($unit[3] >= $unit[4] and $unit[3] >= 2){  ### if this loci is methylated
               $meth->{$start}++;
            }
   }elsif($unit[1] >= $end or $unit[1] < $start){
            return;
   }
}
}


sub ci
{
my ($num)=@_;
my $loop=0;
my $total=0;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        #my $temp=log($_);
        #print "$_\t$temp\n";
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}


my $number=$loop;
return (0,0,0) if ($number < 2);
my $mean=$total/$number;
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return ($mean,$se,$number);
}


sub mean
{
my ($num)=@_;
my $loop=0;
my $total;
foreach  (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
return 0 if ($loop == 0);
my $mean=$total/$loop;
return $mean;
}

sub log10 {
    my ($n) = shift;
    return log($n)/log(10);
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return \%hash;
}



sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq,$step)=@_;
my $length=length $seq;
my $run=int ($length/$step);
my $newseq;
for(my $i=0;$i<=$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}


sub estimateGAP{
my ($seq)=@_;
my $length=length $seq;
my $n=$seq=~tr/Nn/Nn/;
my $gapfrq=sprintf("%02f",$n/$length);
return $gapfrq;
}

sub estimateTE{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
my $n=$seq=~tr/Nn/Nn/;
my $len=$length-$n;
my $te=$seq=~tr/atcg/atcg/;
if ($len > 0){
  my $tefrq=sprintf("%02f",$te/$len);
  return $tefrq;
}else{
  my $tefrq="Na";
  return $tefrq;
}
}



sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}


sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $index="$seq"."_"."$id";
        $hash{$index}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$index}.="$_\n";
    }

}
close IN;

return \%hash;
}


sub estimateGC{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
my $n=$seq=~tr/Nn/Nn/;
my $len=$length-$n;
my $c=$seq=~tr/Cc/Cc/;
my $g=$seq=~tr/Gg/Gc/;
if ($len > 0){
  my $gc=($g+$c)/$len;
  return $gc;
}else{
  my $gc="Na";
  return $gc;
}
}


sub blastpair
{
####read a blasttable file, return a ref of hash contain pair of homologous gene that have e-value <=1e-15,
####and coverage of > 70% for both protein sizes;
my ($blast)=@_;
my %record;
open IN, "$blast" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
   if ($unit[13] > 1e-15){
       next;
   }
   my $pair=$unit[0]."-".$unit[4];
   if (exists $record{$pair}){
        my $refarray=$record{$pair};
        push (@$refarray,[@unit]);
        $record{$pair}=$refarray;
   }else{  
        my @array;
        push (@array,[@unit]);
        $record{$pair}=\@array;
   }
}
close IN;

my %homologous;
foreach (keys %record){
      #print "$_\n";
      #print Dumper($record{$_}),"\n";
      my $refarray=$record{$_};
      my $num=@$refarray;
      #print "$num\n";
      my $hsp=1;      
      my $queryid =$$refarray[0][0];
      my $hitid   =$$refarray[0][4];
      my $bitscore=$$refarray[0][12];
      my $identity=$$refarray[0][8];
      my $querylen=$$refarray[0][1];
      my $hitlen  =$$refarray[0][5];
      my $longestq=$$refarray[0][3]-$$refarray[0][2]+1;
      my $longesth=$$refarray[0][7]-$$refarray[0][6]+1;
      my $matchq  =$$refarray[0][3]-$$refarray[0][2]+1;
      my $matchh  =$$refarray[0][7]-$$refarray[0][6]+1;
      if ($num > 1){  
          my $hspstartq=$$refarray[0][2];
          my $hspstarth=$$refarray[0][6];
          my $hspendq  =$$refarray[0][3];
          my $hspendh  =$$refarray[0][7];
          for(my $i=1;$i<=$num-1;$i++){
                if ($$refarray[$i][2] < $hspendq or $$refarray[$i][6] < $hspendh){
                      next;
                }else{
                      $hsp++;
                      $hspendq=$$refarray[$i][3];
                      $hspendh=$$refarray[$i][7];
                      $bitscore+=$$refarray[$i][12];
                      $identity+=$$refarray[$i][8];
                      $matchq  +=$$refarray[$i][3]-$$refarray[$i][2]+1;
                      $matchh  +=$$refarray[$i][7]-$$refarray[$i][6]+1;
                } 
          } 
          $longestq=$hspendq-$hspstartq+1;
          $longesth=$hspendh-$hspstarth+1;
          $identity=$identity/$hsp;
      }
      my $qhspcoverage=$matchq/$querylen;
      my $hhspcoverage=$matchh/$hitlen;
      my $qlencoverage=$longestq/$querylen;
      my $hlencoverage=$longesth/$hitlen;
      #if ($qhspcoverage >= 0.7 and $hhspcoverage >= 0.7){
      #if ($qlencoverage >= 0.7 and $hlencoverage >= 0.7){
      if ($qhspcoverage >= 0.3 and $hhspcoverage >= 0.3 and $qlencoverage >= 0.5 and $hlencoverage >= 0.5){
          my $temp="$queryid\t$hitid";
          #print "$temp\n";
          $homologous{$temp}=$identity;
      }
}
return \%homologous;
}


