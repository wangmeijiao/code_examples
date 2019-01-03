open HH, ">$outdir/out.orthexon.len.HH.stats" or die "$!";
open HE, ">$outdir/out.orthexon.len.HE.stats" or die "$!";
open EE, ">$outdir/out.orthexon.len.EE.stats" or die "$!";



#...
#...
#...

   my $fh;
   if($type eq "HH"){$fh = \*HH}elsif($type eq "HE"){$fh = \*HE}elsif($type eq "EE"){$fh = \*EE}

   print $fh "\n$gene_japo-$gene_brac\n\n";
   print $fh "#orthGene\texon_num\tintron_num\tgene_len\tintron_len\tintron_perc\texon_len\texon_perc\n";
   foreach my $spec(("japo","brac")){
      my $exon;
      my $len;
      if($spec eq "japo"){$exon=\@segs_japo;$len=$len_gene_japo;print $fh "$gene_japo\t"}else{$exon=\@segs_brac;$len=$len_gene_brac;print $fh "$gene_brac\t"}
      print $fh scalar @$exon,"\t"; #orthexon number
      my $intron=&getIntron($exon);
      print $fh scalar @$intron,"\t"; #intron number
      print $fh "$len\t"; #all orthExon (mrna) length
      my $total_intron=&segSum($intron); #intron length
      print $fh $total_intron,"\t",sprintf("%.3f",$total_intron/$len),"\t";
      my $total_exon=&segSum($exon); #exon length
      print $fh $total_exon,"\t",sprintf("%.3f",$total_exon/$len),"\n";
   }
   print $fh "\n//\n\n";

}
close HH;
close HE;
close EE;

