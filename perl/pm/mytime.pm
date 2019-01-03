

sub mytime(){
  my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
  $year+=1900;
  $mon+=1;
  #print "time is $year,month is $mon,day is $day\n";
  #print "hour is $hour,minutes is $min,second is $sec\n";
  print "$hour,$min,$sec,$year,$mon,$day\n";
  my $string=sprintf ("%s:%2s:%2s  @%s %s.%s\n",$hour,$min,$sec,$year,$mon,$day);
  return $string;
}



