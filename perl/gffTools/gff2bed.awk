

 awk -vFS='\t' -vOFS='\t' '($3=="gene"){split($9,a,"=|;| +");print $1,$4,$5,a[length(a)],a[2],$7}' japonica_mRNA_CDS.embl.gff >japonica_mRNA_CDS.embl.gff.tigr6.bed



