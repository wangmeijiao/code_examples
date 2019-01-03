perl work_smskFa.pl -fa H1_japo.fa -species "rice"
perl work_smskFa.pl -fa H1_glab.fa -species "rice"
perl work_smskFa.pl -fa H1_punc.fa -species "rice"
perl work_smskFa.pl -fa H1_brac.fa -lib "brac"
perl work_smskFa.pl -fa H1_brac500-871264.fa -lib "brac"
perl work_smskFa.pl -fa H1_sorg.fa -species "sorghum bicolor"
perl work_smskFa.pl -fa H1_bd.fa -species "brachypodium"



perl work_smskFa.pl -fa Sbicolor_255_v2.0.fa -species "sorghum bicolor" -keep > log 2>&1 &

bash rpmkout2bedgff_filter.sh Sbicolor_255_v2.0 bed
bash rpmkout2bedgff_filter.sh Sbicolor_255_v2.0 gff





