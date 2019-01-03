my @list=("chr04-chr04","chr05-chr05","chr06-chr06","chr08-chr08","chr09-chr09");
foreach(@list){
  my ($t,$q)=split/-/,$_;
  my $flag=system("bash ../getChain_v3.sh ../../${t}_old.fa.smasked ../../${q}_new.fa.smasked");
  if($flag==0){print"$t,$q suceed\n"}else{print "err $t,$q\n";}
}
