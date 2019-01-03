sub curvature(){
   # three adjacent points curvature approximate :four times the area of the triangle formed by the three points divided by the product of its three sides
   #url http://cn.mathworks.com/matlabcentral/answers/57194-how-to-find-the-sharp-turn-in-a-2d-line-curve
   #http://www.mathopenref.com/coordtrianglearea.html
   my ($x1,$y1,$x2,$y2,$x3,$y3) = @_;
   if($x1 == $x2 && $x1 == $x3 || $y1 == $y2 && $y1 == $y3){print STDERR "three points ($x1,$y1,$x2,$y2,$x3,$y3) can't line up\n";return (0)}
   #if($x1 == $x2 && $x1 == $x3 || $y1 == $y2 && $y1 == $y3){die "three points ($x1,$y1,$x2,$y2,$x3,$y3) can't line up\n"}
   my $areaTriangle1=0.5*abs($x1*($y2-$y3)+$x2*($y3-$y1)+$x3*($y1-$y2));
   my $areaTriangle2= 0.5*abs( ($x2-$x1)*($y3-$y1)-($x3-$x1)*($y2-$y1) );
   print STDERR "area of triangule not eq:$areaTriangle1,$areaTriangle2\n" if($areaTriangle1 != $areaTriangle2);
   my $sqrt = sqrt( (($x2-$x1)**2+($y2-$y1)**2) * (($x3-$x1)**2+($y3-$y1)**2) * (($x3-$x2)**2+($y3-$y2)**2) ); # sqrt of production of three sides
   if( $sqrt  == 0) {die "three points ($x1,$y1,$x2,$y2,$x3,$y3) sqrt eq 0\n",(($x2-$x1)**2+($y2-$y1)**2)," ",(($x3-$x1)**2+($y3-$y1)**2)," ",(($x3-$x2)**2+($y3-$y2)**2),"\n" }
   my $result = 4*$areaTriangle2 / $sqrt;

   return (sprintf ("%.9f",$result));
   #return (sprintf ("%.6f", $areaTriangle1), sprintf("%.6f", $areaTriangle2));
}


