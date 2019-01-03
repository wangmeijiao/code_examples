use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

my $fa;
GetOptions("fa=s",\$fa);


my %gff;#main structure table:  chr->mRNAs->(elements:CDS,intron,UTR,TSS,TSS_flanking,TTS,TTS_flanking)

#1,read&fill %gff without order consideration
#no gene line, no alternative mrna records
#single mrna with exon and cds
while(<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
  my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/[\t ]+/,$_;
  die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
  if($ftype eq "mRNA"){
    $attr=~/ID=([^;]+);/;
    my $ID=$1;
    my $note="adjust_by_mannual_examination";
    die "undefined Id or Note at $_" if(! defined $ID || !defined $note);
    if(!exists $gff{$chr}->{$ID}){
      $gff{$chr}->{$ID}->{'mrna'}="$chr\t$start\t$end\t$ID\t$strand\t$note"
    }else{die "dup mRNA name $ID\n"}
  }elsif($ftype eq "CDS" || $ftype eq "intron"|| $ftype eq "exon" || $ftype eq "five_prime_UTR" ||  $ftype eq "three_prime_UTR"){
       $attr=~/Parent=([^;]+)/;
       my $pID=$1;
       if(!exists $gff{$chr}->{$pID}->{'mrna'}){die "elements comes first, I can't assign it to gene\n"}
       if($ftype eq "CDS"){push @{$gff{$chr}->{$pID}->{$ftype}},"$start-$end-$fshift"}else{push @{$gff{$chr}->{$pID}->{$ftype}},"$start-$end"}
   }else{die "unknow feature type\n"}
}#while end

#print Dumper \%gff;exit;


#2, check and join exon lines; extract and output exon fa
my $fa_idx=&faRead_chr($fa);

foreach my$chr(sort keys %gff){
  foreach my $mrnaID(sort keys %{$gff{$chr}}){
    my ($chr,$start,$end,$id,$strand,$note)=split/\t/,$gff{$chr}->{$mrnaID}->{"mrna"};
    die "mrnaID $mrnaID ne id $id" if($mrnaID ne $id);
    my @exons;
    if(exists $gff{$chr}->{$mrnaID}->{"exon"}){ foreach(@{$gff{$chr}->{$mrnaID}->{"exon"}}){push @exons, $_} }else{die "not exist exon in mrna $mrnaID"}

    #start to check, join exons and print exon gff, extract fa
    if(&segCheck(\@exons)){
      my $seq_exon;
      foreach(@exons){
        my ($s,$e)=split/-/,$_;
        #print "$chr\tgff_simplify\texon\t$s\t$e\t.\t$strand\t.\tParent=$id.longest.NOLP\n";
        if(exists $fa_idx->{$chr}){
          $seq_exon.=substr($fa_idx->{$chr},$s-1,$e-$s+1);
        }else{die "can't find chr in file $fa at $id"}
      }
      #if($strand eq "-"){$seq_exon=&RC($seq_exon)}
      die "empty seq_exon at $mrnaID" if(length $seq_exon == 0 );
      print ">$mrnaID longestNOLP $chr:$start-$end:$strand ",join(' ',@exons),"\n",&faFormat($seq_exon,60);
    }else{die "unsorted or Overlap at exon: @exons of $mrnaID"}

  }#foreach mrna 
}#foreach chr end




#sub


sub segLink(){ #the same with deOverlap()
   #standard function to sort&link(deOverlap) segments, which belong to a single sequence
   #use 3seg method if (seg_num >=3) or if (seg_num ==2 ) use 2seg method . do nothing if (seg_num <2)
   #input format @segs=("s1-e1","s2-e2","..-..")
   my @segs=@{shift @_};
   if (scalar @segs<2){
     return \@segs
    }elsif(scalar @segs ==2){
       my @segs_sorted=@{&segSort(\@segs)};
       my @segs_sorted_linked=&relate2segs($segs_sorted[0],$segs_sorted[1]);
       return \@segs_sorted_linked;
      }elsif(scalar @segs >2 ){
         my @segs_sorted=@{&segSort(\@segs)};
         #print "sorted: @segs_sorted\n";
         my @result=shift @segs_sorted;
         for(0..$#segs_sorted){
             my @relate;
             if($_ != $#segs_sorted){
               @relate=&relate2segs($segs_sorted[$_],$segs_sorted[$_+1]);
               @relate=&relate2segs($result[-1],$relate[0]);
               pop @result;
               push @result,@relate;
             }elsif($_ == $#segs_sorted){
                   @relate=&relate2segs($result[-1],$segs_sorted[$_]);
                   pop @result;
                   push @result,@relate;
                 }
         }
         return \@result;
       }
}


sub segSort(){
     my @segs=@{shift @_};
     my @segs_sorted=sort{
       my ($s1,$e1)=split/-/,$a;
       my ($s2,$e2)=split/-/,$b;
       if($s1<$s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
     } @segs;
     return \@segs_sorted;
}

sub relate2segs(){
    #input paras ("s1-e1","s2-e2") must be sorted by start
    my ($s1,$e1)=split/-/,(shift @_);
    my ($s2,$e2)=split/-/,(shift @_);
    die "start1>start2 at $s1-$e1|$s2-$e2\n" if ($s1>$s2);
    #three types of relationship between two segments
    #1,separated
    if($s2>$e1){
       return ("$s1-$e1","$s2-$e2");
     #2,overlap
     }elsif($s2<=$e1 && $e1<=$e2){
        return ("$s1-$e2")
      #3,included
      }elsif($e2<$e1 && $e2<=$e1){
         return ("$s1-$e1");
        }else{print STDERR "unknown relationship at $s1-$e1|$s2-$e2, return original\n"; return ("$s1-$e1","$s2-$e2")}
}



sub segCheck(){
   #sorted , nonOverlap or not
   my @segs=@{shift @_};
   my @points;
   my $flag=1;
   foreach (@segs){my ($s,$e)=split/-/,$_;push @points,($s,$e)}
   return("err:points <=1") if(scalar @points <= 1);
   for my $i(0 .. $#points-1){
      if($points[$i] <= $points[$i+1]){}else{$flag = 0;print STDERR "check failed at $points[$i] and $points[$i+1]\n"}
   }
   return $flag;
}



sub faRead_chr(){
  my $file = shift;
  my %fa;
  open FA, $file or die "$!";
  $/=">";
  while(<FA>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my @temp=split/[\t ]+/ ,$id;
    $id=shift @temp;
    my $seq = join "", @box;
    if(!exists $fa{$id}){$fa{$id}=$seq}else{die "dup $id\n"}
  }
  close FA;
  $/="\n";
  return \%fa;
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


sub RC(){
  my $seq=shift;
  $seq=reverse $seq;
  $seq=~tr/atcgATCG/tagcTAGC/;
  return $seq;
}





