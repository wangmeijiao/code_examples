

sub readEmbl_gene(){
  my ($file,$keyword)=@_;
  my %gene;
  my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^mRNA/ || $_=~/^CDS/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          if(eof (EMBL)){push @block, $_}; #add final line
          #start to fill gene hash table
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,$start,$end,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){$geneID=$v;next} #mush has /ID=xxx or /gene=xxx
            if(!exists $feature{$k}){$feature{$k}=$v}else{$feature{$k}.="-".$v}#die "dup key $k at $ctrl"}
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

}


sub readEmbl_LTR(){
  my ($file,$keyword)=@_;
  my %gene;
  my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^LTR/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          if(eof (EMBL)){push @block, $_}; #add final line
          #start to fill gene hash table
          $cnt++;
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,$start,$end,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){$geneID=$v;next} #mush has /ID=xxx or /gene=xxx
            if(!exists $feature{$k}){$feature{$k}=$v}else{ $feature{$k}.="-".$v}#$k.="_".&generate_random_string(4);$feature{$k}=$v}
            #if(!exists $feature{$k}){$feature{$k}=$v}else{die "dup key $k at $ctrl"}
          }
          $feature{"strand"}=$strand;
          $feature{"coord"}=\@regions; 
          die "can't find ID line at $_" if(!defined $geneID || $geneID eq "");
          if(exists $gene{$geneID}){
             #die "dup LTR ID $geneID";
             $geneID.="_".&generate_random_string(4);
             $gene{$geneID}=\%feature;
          }else{
            $gene{$geneID}=\%feature;
          }
          #if($feature_name eq "mRNA"){
          #  $gene{$geneID}->{"mRNA"}=\%feature
          #}elsif($feature_name eq "CDS"){
          #    $gene{$geneID}->{"CDS"}=\%feature
          #  }else{die "unknow feature $feature_name"}
          
          #inistialize for next chunk
          @block=();
          push @block, $_;          
        }
      }elsif($_=~/^\//){push @block, $_}else{die "unknown line $_"}

  }#while end
  close EMBL;
  return \%gene;

}


sub readEmbl_DNATE(){
  my ($file,$keyword)=@_;
  my %gene;
  my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^repeat_region/ || $_=~/^gap/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          if(eof (EMBL)){push @block, $_}; #add final line
          #start to fill gene hash table
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,$start,$end,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){$geneID=$v;next} #mush has /ID=xxx or /gene=xxx
            if(!exists $feature{$k}){$feature{$k}=$v}else{ $feature{$k}.="-".$v}#die "dup key $k at $ctrl"}
          }
          $feature{"strand"}=$strand;
          $feature{"coord"}=\@regions; 
          die "can't find ID line at $_" if(!defined $geneID || $geneID eq "");
          if(exists $gene{$geneID}){die "dup LTR ID $geneID"}else{
           $gene{$geneID}=\%feature;
          }
          #if($feature_name eq "mRNA"){
          #  $gene{$geneID}->{"mRNA"}=\%feature
          #}elsif($feature_name eq "CDS"){
          #    $gene{$geneID}->{"CDS"}=\%feature
          #  }else{die "unknow feature $feature_name"}
          
          #inistialize for next chunk
          @block=();
          push @block, $_;          
        }
      }elsif($_=~/^\//){push @block, $_}else{die "unknown line $_"}

  }#while end
  close EMBL;
  #print Dumper \%gene;
  return \%gene;

}

sub generate_random_string{
    #http://stackoverflow.com/questions/13687643/generate-unique-random-strings
    my $length_of_randomstring = shift; # the length of 
                                        # the random string to generate

    my @chars=('a'..'z','A'..'Z','0'..'9','_');
    my $random_string;
    for(1..$length_of_randomstring){
        # rand @chars will generate a random 
        # number between 0 and scalar @chars
        $random_string.=$chars[rand @chars];
    }

    return $random_string;
}


sub parseCoord(){
#general code for parsing embl feature structure line
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


