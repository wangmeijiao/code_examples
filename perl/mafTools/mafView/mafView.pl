use strict;
use warnings;
use Data::Dumper;
use SVG;


my ($maf,$maxLen,$len,$n_spec,$order,$range)=&mafRead($ARGV[0]);
#print Dumper $maf,$maxLen,$len,$n_spec,$order,$range; exit;

my $space=100;
my $barH=10;
my ($height,$width)=(600,1600);
my $dist_chr=$height/($n_spec-1);  
my ($len_region,$maxLen_region)=&getLenAll($range);
my $ppx=$width/$maxLen_region; 

print STDERR "Guessed and calculated from maf file $ARGV[0]: $n_spec species, maxLen=$maxLen_region, with ppx=$ppx; for @{$order}\n ";
my $svg=SVG->new('height'=>$height+4*$space,'width'=>$width+2*$space);
#1, draw overall architechture
  my $y_next=$space;
  foreach my $spec(@{$order}){
     my $len_chr= $len_region->{$spec};
     $svg->rect("x",$space,"y",$y_next,"height",$barH,"width",$len_chr*$ppx,"rx",8,"ry",8,"fill","white","stroke","none","stroke-width",2,"fill-opacity",0.2);
     $svg->text("x",$space/4,"y",$y_next+5,"height",$barH,"width",$len_chr*$ppx,"-cdata",$spec,);
     $y_next+=$dist_chr;
  } 

#2, draw anchors & links
  foreach my $i( sort {$a<=>$b} keys %{$maf}){
    $y_next=$space;
    my (@x_s,@x_e,@y_s,@y_e);
    #my (@mids_x,@mids_y);
    #my @anchors=split/\t/,$_; 
    foreach my $spec(@$order){
      if(!exists $maf->{$i}->{$spec}){$y_next+=$dist_chr;next}
      my $s=$maf->{$i}->{$spec}->[1]-$range->{$spec}->{'begin'}; 
      my $e=$maf->{$i}->{$spec}->[2]-$range->{$spec}->{'begin'}; 
      my $strand=$maf->{$i}->{$spec}->[3]; 
      #$svg->rect("x",$space+$s*$ppx,"y",$y_next,"height",$barH,"width",($e-$s)*$ppx,"fill","black","stroke","black");
      $svg->rect("x",$space+$s*$ppx,"y",$y_next,"height",$barH,"width",($e-$s)*$ppx,"fill","black","stroke","black","stroke-width",0.1);
      push @x_s,$space+$s*$ppx;
      push @x_s,$space+$s*$ppx;
      push @x_e,$space+$e*$ppx;
      push @x_e,$space+$e*$ppx;
      push @y_s,$y_next;
      push @y_s,$y_next+$barH;
      push @y_e,$y_next;
      push @y_e,$y_next+$barH;

      $y_next+=$dist_chr;
    }

=pod
    my $points=$svg->get_path(
       "x" => \@mids_x,
       "y" => \@mids_y,
       "-type" => "polyline",
       "-closed" => "false",
    );
    $svg->polyline(
        %$points,       
        "style" => {
                     "stroke" => "#777777",
                     "fill" => "none"
                   }
    );
=cut

    my @x = (@x_s,reverse(@x_e));
    my @y = (@y_s,reverse(@y_e));


    my $points=$svg->get_path(
       "x" => \@x,
       "y" => \@y,
       "-type" => "polygon",
       "-closed" => "false",
    );
    $svg->polygon(
        %$points,
        "style" => {
                     "stroke" => "#777777",
                     "fill" => "navy",
                      "fill-opacity"=>0.6,
                      "stroke-width"=>0.1
                   }
    );



  }

  




=pod
#3, draw block shadows
#format :  japo|H1_japo|1230100|1310000	brac|H1_brac|748297|852229	brac|H1_brac|548297|652229
#          japo|H1_japo|1430100|1560000	brac|H1_brac|548997|656229	brac|H1_brac|848297|852229

my @test=("japo|H1_japo|530100|610000\tbrac|H1_brac|348297|552229\tbrac|H1_brac|548297|652229","japo|H1_japo|1430100|1560000\tbrac|H1_brac|548997|656229\tbrac|H1_brac|748297|852229");
#print STDERR "@test\n";
foreach(@test){
  my (@x,@y);
  my $y_next=$space;
  my @species=split/\t/,$_;
  foreach(@species){
    my($species,$chr,$s,$e)=split/\|/,$_;
    unshift @x,$space+$s*$ppx;
    push @x, $space+$e*$ppx;
    unshift @y, $y_next;
    push @y, $y_next;
    $y_next+=$dist_chr;
  }
  my $points=$svg->get_path(
       "x" => \@x,
       "y" => \@y,
       "-type" => "polyline",
       "-closed" => "true",
  );
  $svg->polyline(
        %$points,
        "style" => {
                     "stroke" => "none",
                     "fill" => "navy",
                     "fill-opacity" => 0.2
                   }
  );

}
=cut


print $svg -> xmlify();




###sub###


sub mafRead(){
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



 

sub readLen(){
  my %len;
  my $file=shift @_;
  open LEN, $file or die "$!";
  while(<LEN>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my ($chr,$len)=split/\t/,$_;
    if(!exists $len{$chr}){$len{$chr}=$len}else{die"dup $chr\n"}
  }
  close LEN;
  return \%len;
}

sub getLenAll(){
  my $region=shift;
  my %len;
  my $max = 0;
  foreach my $spec(keys %{$region}){
   my $begin=$region->{$spec}->{"begin"};
   my $end=$region->{$spec}->{"end"};
   die "begin >= end : $begin >= $end at $spec" if($begin >= $end);
   my $len= $end-$begin+1; 
   if($len > $max){$max = $len}
   $len{$spec}=$len;
  }
  return (\%len,$max);
}




sub drawOrthLink(){
  my ($tab,$order,$y_spec,$dist)=@_;

  foreach my $id(sort keys %$tab){
    for(my $i=0;$i<=$#$order-1;$i++){
       next if(!exists $tab->{$id}->{$order->[$i]} || !exists $tab->{$id}->{$order->[$i+1]});
       my ($spec1,$chr1,$start1,$end1,$id1,$strand1)=split/-|\|/,$tab->{$id}->{$order->[$i]};
       my ($spec2,$chr2,$start2,$end2,$id2,$strand2)=split/-|\|/,$tab->{$id}->{$order->[$i+1]};
       die "spec not eq: $spec1, $order->[$i]" if( $spec1 ne $order->[$i]);
       die "spec not eq: $spec2, $order->[$i+1]" if( $spec2 ne $order->[$i+1]);
       my $mid1=0.5*($start1+$end1);
       my $mid2=0.5*($start2+$end2);
       my $y1=$y_spec->{$spec1}+2;
       my $y2=$y1+$dist;
       $svg->line("x1",$space+$start1*$ppx,"y1",$y1,"x2",$space+$end1*$ppx,"y2",$y1,"stroke","black","stroke-width",0.5);
       $svg->line("x1",$space+$start2*$ppx,"y1",$y2,"x2",$space+$end2*$ppx,"y2",$y2,"stroke","black","stroke-width",0.5);
       $svg->line("x1",$space+$mid1*$ppx,"y1",$y1,"x2",$space+$mid2*$ppx,"y2",$y2,"stroke","black","stroke-width",0.5);
       $svg->text("x",$space+$mid1*$ppx,"y",$y1,"-cdata",$id1,style=>{"font-family","Arial","font-size",8,"text-anchor","middle"},"stroke","black");
       $svg->text("x",$space+$mid2*$ppx,"y",$y2,"-cdata",$id2,style=>{"font-family","Arial","font-size",8,"text-anchor","middle"},"stroke","black");
    }
  }

  return 1;
}





