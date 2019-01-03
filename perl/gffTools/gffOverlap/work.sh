
cat all.simple.allLen_exon.gff |perl gffOverlap.pl chr04:30701378-30902378 |awk -vOFS='\t' '{$4-=30701378;$5-=30701378;if($4<0){$4=1};if($5<0){$5=1};print $0}' > all.simple.allLen_exon.chr04_30701378-30902378.transform.gff


perl gff2embl_gene.pl all.simple.allLen_exon.chr04_30701378-30902378.transform.gff > all.simple.allLen_exon.chr04_30701378-30902378.transform.gff.embl





