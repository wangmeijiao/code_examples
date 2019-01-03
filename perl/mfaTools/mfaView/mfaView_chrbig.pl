use strict;
use warnings;
use Data::Dumper;
use SVG;



my $space=150;
my $barH=10;
my ($height,$width)=(200,4000);
my $ppx;


my $svg=SVG->new('height'=>$height,'width'=>$width+2*$space);
my $y_next=0.2*$height;


my $alignLen;
$/=">";
while(<stdin>){
  chomp;
  next if($_ eq "" || $_ =~/^\s+$/);
  my @box =split/\n+/,$_;
  my $id = shift @box;
  print STDERR "read and draw chr $id\n";
  my @temp = split/[\t ]+/,$id;  
  my $id_real = shift @temp;
  my $seq = join("",@box);
  my $len = length $seq;
  if(!defined $alignLen){$alignLen = $len}else{if($len != $alignLen){die "align len diffs at $id"}}
  if(!defined $ppx){$ppx=$width/$len}
  #draw chr
  $svg->rect("x",$space,"y",$y_next,"height",$barH,"width",$len*$ppx,"rx",8,"ry",8,"fill","white","stroke","none","stroke-width",2,"fill-opacity",0.2);
  $svg->text("x",$space/4,"y",$y_next+5,"height",$barH,"width",$space,"-cdata",$id_real);
  #annotion and paint chr

  my $cnt=0;
  while($seq =~/[^-nN]+/g){
    #$cnt++;
    #print STDERR "$cnt match at $-[0],$+[0] ";
    #print STDERR "hit is:",substr($seq,$-[0],$+[0]-$-[0]),"\n";
    $svg->rect("x",$space+$-[0]*$ppx,"y",$y_next,"height",$barH,"width",($+[0]-$-[0])*$ppx,"fill","navy","stroke","none","stroke-width",2,"fill-opacity",0.8);
    #last if($cnt  >= 1000);

  }
  print STDERR "$id done\n";
  $y_next+=1.2*$barH;

}
$/="\n";

 &drawBdg_heatmap(&readBdg($ARGV[0]),&greyMono256(),$y_next,"none","phastP.gerp");


print $svg -> xmlify();


##sub

sub readBdg(){
   my $file=shift;
   my @bdg;
   open BDG, $file or die "$! $file";
   while(<BDG>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
     my($chr,$s,$e,$v)=split/[\t ]+/,$_;
     die "$s >= $e at $_" if($s >= $e);
     #next if($v == 0);     
     push @bdg,$_;
   }
   close BDG;

   #my $bdg_pack=&bdgPack(\@bdg);
   return \@bdg;
}

sub drawBdg_heatmap(){
  #fill real bed region with color 
  my ($bdg_idx,$color,$y,$range,$mark)=@_;
  my $ncol=scalar @{$color} - 1;
  my $idx_color=&data2color($bdg_idx,$ncol,$range);
  my $barH=10;  #heatmap bar height
  foreach(@{$idx_color}){
    my ($chr,$s,$e,$idx)=split/\t/,$_;
    #$svg->rect("x",$space+$s*$ppx,"y",$y,"height",$barH,"width",($e-$s+1)*$ppx,"fill",$color->[$idx],"stroke","none");
    $svg->rect("x",$space+$s*$ppx,"y",$y,"height",$barH,"width",($e-$s+1)*$ppx,"fill",$color->[$idx],"stroke",$color->[$idx]);
      #very sharp
  }
  my $txt_col;
  #if ($mark =~/H3k9me2/i){$txt_col = "red"}else{$txt_col = "black"}
  $svg->text("x",$space*0.9,"y",$y+0.6*$barH,"-cdata",$mark,"font-family","ArrialNarrow","font-size",8,"stroke",$txt_col,"fill",$txt_col,"text-anchor","end");
  return $y+$barH+3;
}


sub data2color(){
  #map data from min-max (or customed min-max) to col
  my ($data,$ncol,$range)=@_;
  my @index;
  if($range ne "none"){
    my ($min_cut,$max_cut)=split/-/,$range; #mid_cut not used
    foreach (@{$data}){
       my ($chr,$s,$e,$v)=split/[\t ]+/,$_;
       $v-=$min_cut; #remove background noise 
       if($v < 0){$v = 0}
       if($v > $max_cut){ $v = $max_cut}
       #print ($v,"\t",&round($v*$ncol/$max_cut),"\n"); 
       push @index,join("\t",($chr,$s,$e,&round($v*$ncol/$max_cut))  );
    }
  }else{
      my ($min,$max)=&minMax($data);
      foreach (@{$data}){
         my ($chr,$s,$e,$v)=split/[\t ]+/,$_;
         #print ($v,"\t",&round(($v-$min)*$ncol/($max-$min)),"\n");
         push @index,join("\t",($chr,$s,$e,&round(($v-$min)*$ncol/($max-$min))) );
      }
    }
  return \@index;
}


sub round(){
   my $n=shift;
   my $i=int $n;
   my $d=$n-$i;
   if($d>0.6){return $i+1}else{return $i}
}

sub minMax(){
   my @data=@{shift @_};
   my ($min,$max)=(0,0);
   foreach(@data){
     my ($chr,$s,$e,$v)=split/[\t ]+/,$_;
     if($min> $v){$min= $v}
     if($max< $v){$max= $v}
   }
   return($min,$max);
}




sub redMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb(255,$i,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub redMono9(){
   my @color=("#FFFFFF", "#FFDFDF", "#FFBFBF", "#FF9F9F" ,"#FF7F7F" ,"#FF5F5F", "#FF3F3F" ,"#FF1F1F", "#FF0000");
   return \@color;
}

sub greenMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb($i,255,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub greyMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb($i,$i,$i)";
   }
   #print "@color\n";
   return \@color;
}



