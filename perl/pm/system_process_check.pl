
sub pcheck
{
    my $cmd = "ps -ef | grep " . $pname . " | grep -v grep | grep -v perl | awk '{print \$2}'";
    my $ret = `$cmd` ;
    chomp $ret;
    if( $ret =~ /\n/ ){
        $ret = ( split /\n/, $ret )[0] ; 
    }
    ($ret eq "") ? print 0 : print $ret ;
}




sub gettime
{
   my($sec,$min,$hour,$day,$mon,$year) = (localtime(time));
   $mon +=1 ;
   $day = ($day < 10 )?"0$day":$day;
   $mon = ($mon < 10 )?"0$mon":$mon;
   $min = ($min < 10 )?"0$min":$min;

   return "$mon$day$hour$min" ; 
}







