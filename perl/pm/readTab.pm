
sub readTab_r(){
   #usage: my ($data_idx,$header_idx,$rowname_idx,$nrow,$ncol,$min,$max)=&readTable($data,$header,$rowname);
   my $file=shift;
   my $head=shift;
   my @header;
   my $rowname = shift;
   my @rowname;
   my @data;
   my $cnt;
   my $ncol;
   my ($min,$max);
   open DATA, $file or die "$!";
   while(<DATA>){
     chomp;
     next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
     $cnt++;
     if($cnt==1){
       if($head eq "T"){
           @header=split/[\t ]+/,$_;
           if($rowname eq "T"){ shift @header}
           $ncol=scalar @header;
           next;
       }elsif($head eq "F"){
           my @box=split/[\t ]+/,$_;
           if(&notNA(\@box) == 0){die "has NA at $_"}
           #if(&isNumeric(\@box) == 0){die "not numeric at $_"}
           if($rowname eq "T"){ my $rowid = shift @box; push @rowname, $rowid}; 
           $ncol=scalar @box;
           push @data,join("\t",@box);
           ($min,$max)=($box[0],$box[0]);
           next;#?
         }else{die "unknow header swith $head"}
     }
     my @box_line=split/[\t ]+/,$_;
     if(&notNA(\@box_line) == 0){die "has NA at $_"}
     #if(&isNumeric(\@box_line) == 0){die "not numeric at $_"} 
     if($rowname eq "T"){ my $rowid = shift @box_line; push @rowname, $rowid};
     my $ncol_line=scalar @box_line;
     if(!defined $ncol || $ncol_line != $ncol){die "ncol err $ncol != $ncol_line at $_"}
     push @data,join("\t",@box_line);
     if(! defined $min && !defined $max){
       ($min,$max)=($box_line[0],$box_line[0])
     }elsif(defined $min && defined $max){
         foreach my $ctrl(@box_line){ 
           if($ctrl < $min){$min = $ctrl}
           if($ctrl > $max){$max = $ctrl}
         }
       }else{die "err min max $min , $max"}    
   }#while end
   close DATA;
   if($head eq "T"){$cnt-=1}
   return (\@data,\@header,\@rowname,$cnt,$ncol,$min,$max);

}

sub readTab_array(){
  my $file = shift;
  my @tab;
  open TAB, $file or die "$!";
  my $cnt;
  my $cnt_valid;
  while(<TAB>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^gene_ID/);
    $cnt++;
    my @box = split/\t/,$_;
    die "abnormal line length at $_ in file $file" if(scalar @box != 13 && scalar @box != 4);
    next if($box[3] eq "NA");
    $cnt_valid++;
    for my $i(0..12){
      push @{$tab[$i]},$box[$i];
    }

  }

  close TAB;
  return ($cnt,$cnt_valid,\@tab);

}

sub readTab_hash(){
  my $file = shift;
  my %tab;
  my @names = ("gene_ID","gene_len","mRNA_ID","intron_num","intron_len_str","intron_len_total","intron_len_perc","te_cnt_str","sum_te_cnt","te_len_str","sum_te_len_total","te_occ_perc_str","te_len_total_perc");
  open TAB, $file or die "$!";
  my $cnt;
  my $cnt_valid;
  while(<TAB>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^gene_ID/);
    $cnt++;
    my @box = split/\t/,$_;
    #my ($gene_ID,$gene_len,$mRNA_ID,$intron_num,$intron_len_str,$intron_len_total,$intron_len_perc,$te_cnt_str,$sum_te_cnt,$te_len_str,$sum_te_len_total,$te_occ_perc_str,$te_len_total_perc) = split/\t/,$_;
    die "abnormal line length at $_ in file $file" if(scalar @box != 13 && scalar @box != 4);
    next if($box[3] eq "NA");
    $cnt_valid++;
    for my $i(0..12){
      push @{$tab{$names[$i]}},$box[$i];
    }

  }

  close TAB;
  return ($cnt,$cnt_valid,\%tab);

}





