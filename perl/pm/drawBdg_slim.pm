
sub drawBdg_slim(){
   my ($data,$peaks,$region,$y,$mark,$range,$mode)=@_;
   my ($chr_r,$s_r,$e_r) = split/:|-/,$region;
   my $len = $e_r-$s_r+1;
   my ($min,$mid,$max)=split/-/,$range;
   $min||=0;
   $max||=15;
   my $boxH=15;
   my $f=15/$max;
   my $mid_height=$boxH*$mid/$max;
   #my $x_left=$space;
   #my $x_right=$width;
   #my $x=$x_left;
   #my $switch=shift;
   $svg->line("x1",$space,"y1",$y+$boxH,"x2",$space+$len*$ppx,"y2",$y+$boxH,"stroke","grey","stroke-width",0.5);
   $svg->text("x",$space-15,"y",$y+$boxH-2,"width",10,"height",5,"-cdata",$mark,"text-anchor","end","stroke","black");

  #y axis ticks and text
   $svg->line("x1",$space,"y1",$y+$boxH+5,"x2",$space,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
   for(my $i=0;$i<=$boxH;$i+=$boxH/3){
     $svg->line("x1",$space,"y1",$y+$boxH-$i,"x2",$space-2,"y2",$y+$boxH-$i,style=>{"stroke-width",1,"stroke","black"});
     if($i == 0){
       $svg->text("x",$space-5,"y",$y+$boxH-$i-2,"-cdata",int $min,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end"})
     }elsif($i==$boxH){
       $svg->text("x",$space-5,"y",$y+$boxH-$i-2,"-cdata",int $max,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end"})
     }


   }

   #draw data with barplot or filled curve
   if($mode eq "barplot"){
     foreach my$line(@$data){
       #$svg->rect("x",$x,"y",$y+15-$value*$scale{"all"}*$scale{$mark},"width",1,"height",$value*$scale{"all"}*$scale{$mark},"fill",$color{$mark});
       my ($chr_in,$s,$e,$value)=split/\t/,$line;
       next if($chr_r ne $chr_in || $e < $s_r || $s > $e_r);
       $s-=$s_r; if ($s < 0){ $s = 0}
       if($e > $e_r){ $e = $e_r}
       $e-=$s_r; if ($e < 0){ $e = 0}
       $value-=$min;
       if($value<0){$value=0}
       if($value>$max){$value=$max}
       if($value <= $mid){
        #$svg->rect("x",$x,"y",$y+$box_height-$value*$f,"width",1,"height",$value*$f,"fill","grey","stroke","none")
        $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH-$mid_height,"width",($e-$s+1)*$ppx,"height",$value*$f,"fill","grey","stroke","none")
        #$svg->rect("x",$x,"y",$y+$box_height-$mid_height,"width",1,"height",$value*$f,"fill","grey","stroke","none")
       }else{
           #$svg->rect("x",$x,"y",$y+$box_height-$value*$f,"width",1,"height",$value*$f,"fill",$color{$mark},"stroke","none")
           $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH-$value*$f,"width",($e-$s+1)*$ppx,"height",$value*$f-$mid_height,"fill",$color{$mark},"stroke",$color{$mark},"stroke-width",0.2)
           #$svg->rect("x",$x,"y",$y+$box_height-$value*$f,"width",1,"height",$value*$f-$mid_height,"fill",$color{$mark},"stroke","none")
          }
     }
   }elsif($mode eq "curve"){
      



    }else{die "unknow plotting type $mode"}

   #draw peak regions
   if($peaks ne "none"){
     foreach my $id(keys %{$peaks}){
       my($chr_in,$s,$e,$strand)=split/-|:/,$peaks->{$id};
       #if($s<$start && $e >$start){$s=$start}
       #if($s<$end && $e > $end){$e=$end}
       $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH+2,"width",($e-$s+1)*$ppx,"height",4,"fill","black");
     }
   }

   return $y+$boxH+0.5*$boxH; #spaced to another track block
}#sub end


