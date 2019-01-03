
sub readEmbl(){
  #note that, coord of mRNA and CDS are unchecked,  gene order is not kept
  my ($file,$keyword)=@_;
  my %gene;
  #my @order;
  #my %id_list;
  #my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^mRNA/ || $_=~/^CDS/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          #start to fill gene hash table
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,undef,undef,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){
              #mush has /ID=xxx or /gene=xxx
              $geneID=$v;
              next;
            }
            if(!exists $feature{$k}){$feature{$k}=$v}else{die "dup key $k at $ctrl"}
          }
          $feature{"strand"}=$strand;
          $feature{"coord"}=\@regions;
          die "can't find ID line at $_" if(!defined $geneID || $geneID eq "");
          if($feature_name eq "mRNA"){
            $gene{$geneID}->{"mRNA"}=\%feature
          }elsif($feature_name eq "CDS"){
              $gene{$geneID}->{"CDS"}=\%feature
            }else{die "unknow feature $feature_name"}

          #inistialize for next chunk
          @block=();
          push @block, $_;
        }
      }elsif($_=~/^\//){push @block, $_}else{die "unknown line $_"}

  }#while end
  close EMBL;
  return \%gene;
  #return (\%gene,\@order);
}



#general code for parsing embl feature structure line
sub parseCoord(){
 my $string=shift;
 if($string=~/^join\((.+)\)$/){
        my @temp=split/\.\./,$1;
        return "+",$temp[0],$temp[$#temp],(split/,/,$1);
    }elsif($string=~/^[0-9]+\.\.[0-9]+$/){
         my @temp=split/\.\./,$string;
         return "+",$temp[0],$temp[$#temp],$string;
       }elsif($string=~/^complement\(join\((.+)\)\)$/){
             my @temp=split/\.\./,$1;
             return "-",$temp[0],$temp[$#temp],split/,/,$1;
          }elsif($string=~/^complement\(([0-9]+\.\.[0-9]+)\)$/){
             my @temp=split/\.\./,$1;
             return "-",$temp[0],$temp[$#temp],$1;
            }else{die  "err when parse coordination $string\n"}
}#sub end here



