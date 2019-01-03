
#must use customed blast m8 (with sequence)



#1. gff to bed
 awk -vOFS='\t' '{split($9,attr,"=|;");print $1,$4,$5,$3"_"attr[2],".",$7   }' japo.mRNA_CDS.embl.gff > japo.mRNA_CDS.embl.gff.bed


#2, filter blast table

cat H1_japo-VS-H1_zs97.blast |perl blastm8filter_quick.pl -len 200 -ident 80 > H1_japo-VS-H1_zs97.blast.filter 2> H1_japo-VS-H1_zs97.blast.lost

perl blastm8filter_slope.pl H1_japo-VS-H1_zs97.blast 500 > H1_japo-VS-H1_zs97.blast.filterslope 2> H1_japo-VS-H1_zs97.blast.lostslope

#3, liftOver_blastm8_bed  liftOver_blastm8_gff

perl liftOver_blastm8_gff.pl -tab H1_japo-VS-H1_zs97.blast -ingff japo.mRNA_CDS.embl.gff -gapSize 50000 -iterStep 100 > japo.mRNA_CDS.embl.liftOver.H1_zs97.gff 2> liftOver.err

perl liftOver_blastm8_bed.pl -tab H1_japo-VS-H1_zs97.blast -inbed japo.mRNA_CDS.embl.gff.bed -gapSize 50000 -iterStep 100 > japo.mRNA_CDS.embl.liftOver.H1_zs97.bed 2> liftOver.bed.err







  perl liftOver_mummer.pl -tab /home/mjwang/pwdexx/oryza_epiCompara_sixSpecies/stats_knobs/seq_percen/mummer_hmsk/H1_niva.fa.hmsk-H1_japo.fa.hmsk.mums -inbed H1_japo.GIs.bed -gapSize 50000 -iterStep 50  > H1_japo.GIs.toniva.hmsk.liftOver.bed 2> H1_japo.GIs.toniva.hmsk.liftOver.err


 perl liftOver_mummer.pl -tab /home/mjwang/pwdexx/oryza_epiCompara_sixSpecies/stats_knobs/seq_percen/mummer/H1_niva.fa-H1_japo.fa.mums -inbed H1_japo.GIs.bed -gapSize 50000 -iterStep 50  > H1_japo.GIs.toniva.liftOver.bed 2> H1_japo.GIs.toniva.liftOver.err


cat test.anchor.bed |perl liftOver_orthGene.pl -tab final.japoVsbrac.simple.orthTab > test.anchor.liftOver.bed 









  #1, solar blast(hmsk) m8 link
  #blast2solar
  perl ~/progs/solar/solar.pl -c -b -m 500 -f m8 H1_japo-VS-H1_zs97.blast > H1_japo-VS-H1_zs97.blast.solar
  #filter stat solar
  cat H1_japo-VS-H1_zs97.blast.solar |perl solarFilterStats.pl 500 > H1_japo-VS-H1_zs97.blast.solar.filter500 2> filter.stats
  #solar2blast
  perl solar2m8.pl H1_japo-VS-H1_zs97.blast.solar.filter500 > H1_japo-VS-H1_zs97.blast.solar.filter500.blast


  #2, liftOver_blastm8_bed
  

  perl liftOver_blastm8_bed.pl -tab H1_japo-VS-H1_zs97.blast.solar.filter500.blast -inbed japo.mRNA_CDS.embl.gff.bed -gapSize 50000 > japo.mRNA_CDS.embl.gff.liftOver.bed  2> liftOver.err


  #3, bed2embl
  perl bed2embl_LTR.pl japo.mRNA_CDS.embl.gff.liftOver.bed > japo.mRNA_CDS.embl.gff.liftOver.bed.embl












