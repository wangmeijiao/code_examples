
sub pslRead() {
  # -nohead
  my $file = shift;
  my $cnt;
  my %anchors;
  open PSL, $file or die "$!";
  while(<PSL>){
   chomp;
   next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
   die "use -noHead\n" if($_=~/^psLayout/ || $_=~/^-/);
   my @result = split(/\t/,$_);
   die "psl format err at $_" if(scalar @result != 21);
   foreach my $i(@result){die "undefined element at $_" if(! defined $i)}
   $cnt++;
   my $matches = $result[0];     #Number of bases that matches the query matched to the target
   my $misMatches = $result[1];  # Number of bases that don't match
   my $repMatches = $result[2];  # Number of bases that match but are part of repeats
   my $nCount = $result[3];      # Number of 'N' bases
   my $qNumInsert = $result[4];  # Number of inserts in query
   my $qBaseInsert = $result[5]; # Number of bases inserted in query
   my $tNumInsert = $result[6];  # Number of inserts in target
   my $tBaseInsert = $result[7]; # Number of bases inserted in target
   my $strand = $result[8];      # '+' or '-' for *query* strand.
   my $qName = $result[9];       # Query sequence name
   my $qSize = $result[10];      # Query sequence size
   my $qStart = $result[11];     # Alignment start position in query
   my $qEnd = $result[12];       # Alignment end position in query
   my $tName = $result[13];      # Target sequence name
   #$tName =~ s/\.fa$//;
   my $tSize = $result[14];      # Target sequence size
   my $tStart = $result[15];     # Alignment start position in target
   my $tEnd = $result[16];       # Alignment end position in target
   my $blockCount = $result[17]; # Number of blocks in the alignment (a block contains no gaps)
   my $blockSizes = $result[18]; # Comma-separated list of sizes of each block
   my $qStarts = $result[19];    # Comma-separated list of starting positions of each block in query (always incremental)
   my $tStarts = $result[20];    # Comma-separated list of starting positions of each block in target (always incremental) 

   die "$qStart >= $qEnd || $tStart >= $tEnd at $_" if($qStart >= $qEnd || $tStart >= $tEnd);

   #swap tStart and tEnd
   #if ($strand eq "-"){ my $temp=$tStart;$tStart=$tEnd;$tEnd=$temp;$strand="+"} #coordination transformation minus->plus

   #Calculate score for this match line (based on UCSC blat)
   #score = match - mismatch - Qgapcount - Tgapcount
   my $sizeMul = 1; #This is 1 for dna and 3 for protein
   my $score = ($sizeMul * ($matches + $repMatches)) - ($sizeMul * $misMatches) - $qNumInsert - $tNumInsert;

   #calculate overall perc for this line  
   


   #store full, summary  and detail
   my $full = join("\t",@result); 
   my $summary = join("\t", ($qName,$tName,$qStart,$qEnd,$tStart,$tEnd,$strand,$score));
   my @blockSizes = split/,/,$blockSizes;
   my @qStarts = split/,/,$qStarts;
   foreach my $i(0..$#qStarts-1){if($qStarts[$i] < $qStarts[$i+1]){}else{die "not incremental at $_"}}
   my @tStarts = split/,/,$tStarts;
   foreach my $i(0..$#tStarts-1){if($tStarts[$i] < $tStarts[$i+1]){}else{die "not incremental at $_"}}
   die "block size diffs at $_" if(scalar @blockSizes != $blockCount || scalar @qStarts != $blockCount || scalar @tStarts != $blockCount);
   my @qSegments;
   my @tSegments;
   foreach my $i(0..$#blockSizes){
      my $len = $blockSizes[$i];
      my $q_start = $qStarts[$i];
      my $q_end = $q_start + $len;
      push @qSegments,"$q_start-$q_end";
      my $t_start = $tStarts[$i];
      my $t_end = $t_start + $len;
      push @tSegments,"$t_start-$t_end";
   }
   $anchors{$cnt}->{"full"}=$full;
   $anchors{$cnt}->{"summary"}=$summary;
   $anchors{$cnt}->{"qDetail"}=\@qSegments;
   $anchors{$cnt}->{"tDetail"}=\@tSegments;
   #last if ($cnt == 10);
  }#while end
  close PSL;
  return \%anchors;
}#sub end


