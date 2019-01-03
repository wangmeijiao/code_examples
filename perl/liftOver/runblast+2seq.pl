opendir DIR, "./" or die "can not open my dir";
foreach my $file (readdir DIR){
    if ($file=~/^(.*)\.hmsk.fa$/){
       print "$file\n";
       push (@seq,$1);
       #push (@name,$1);
       system "cat $file |makeblastdb -dbtype nucl -out $file -title $file";
       #system "formatdb -i $file -p F";
    }
}
close DIR;

for(my $i=0;$i<@seq;$i++){
   for(my $j=0;$j<@seq;$j++){
     system "blastn -task blastn -db $seq[$j].hmsk.fa -out $seq[$i]-VS-$seq[$j].blast -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -num_threads 6 -evalue 1e-5 -query $seq[$i].hmsk.fa";
     #system "blastn -task blastn -db $seq[$j].fa.hmsk -out $seq[$i]-VS-$seq[$j].blast -outfmt 0 -ungapped -num_threads 6 -evalue 1e-5 -query $seq[$i].fa.hmsk";
     #system "blastall -p blastn -U -F F -i $seq[$i].fa.hmsk -d $seq[$j].fa.hmsk -o $seq[$i]-VS-$seq[$j].blast -e 1e-10 -m 8 ";
     print "$i\t$j\n"; 
   }
}








