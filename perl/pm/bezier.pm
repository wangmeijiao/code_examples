
sub drawPesudoLoop(){
   my($regions,$y,$flip,$x_begin)=@_;
   die "less than 2\n" if(scalar @{$regions}<2);
   my $loopH=30;
   my $boxH=2;
   my $string;
   if($flip=~/^T/i){$y=$y+$loopH}
   for(my $i=0;$i<=$#{$regions}-1;$i++){
       my($chr1,$start1,$end1)=split/\|/,$regions->[$i]; die "end < start" if($end1<=$start1);
       my($chr2,$start2,$end2)=split/\|/,$regions->[$i+1]; die "end < start" if($end2<=$start2);
       die "not sorted" if( $end1>$start2);
       die "not the same chr" if($chr1 ne $chr2);
       my $mid1=$x_begin+$start1+($end1-$start1+1)/2;
       my $mid2=$x_begin+$start2+($end2-$start2+1)/2;
       my $midLoop=$mid1+($mid2-$mid1+1)/2;
       my $yLoop1;
       if($flip=~/^T/i){$yLoop1=$y}else{$yLoop1=$y+$boxH};
       my $yLoop2;
       if($flip=~/^T/i){$yLoop2=$y-$loopH}else{$yLoop2=$yLoop1+$loopH};
       $string.="M$mid1 $yLoop1, Q$midLoop $yLoop2, $mid2 $yLoop1 ";# Belzier Curve: "M180 10,Q260 280,340 10 ";
       $svg->rect("x",$x_begin+$start1,"y",$y,"height",$boxH,"width",$end1-$start1+1,"fill","black","stroke","black");
       $svg->rect("x",$x_begin+$start2,"y",$y,"height",$boxH,"width",$end2-$start2+1,"fill","black","stroke","black");
   }
   #print STDERR $string,"\n";
   my $tag = $svg->path(
                     d => $string,
                     #id => 'pline_1',
                     style => {
                               'fill' => 'none',
                               'stroke' => 'grey'
                              }
                     );
  return $y+$loopH;
}

