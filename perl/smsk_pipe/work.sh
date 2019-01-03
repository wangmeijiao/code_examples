fa=$1

cp -r ~/progs/misc-tools/smsk_pipe/bin/ .

cp  ~/progs/misc-tools/smsk_pipe/work_smsk2hmsk.sh .

cp  ~/progs/misc-tools/smsk_pipe/rpmkout2bedgff_filter.sh .

cp  ~/progs/misc-tools/smsk_pipe/work_smskFa.pl .

perl work_smskFa.pl -fa $fa -species "rice" -keep > pipe.log 2>&1 
#!! use .fa file only !  .fasta will fail the pipe
#!! remember fa head id must not be chr04:3400-4500


##########
#perl work_smskFa.pl -fa H1_japo.fa -species "rice"
#perl work_smskFa.pl -fa H1_glab.fa -species "rice"
#perl work_smskFa.pl -fa H1_punc.fa -species "rice"
#perl work_smskFa.pl -fa H1_brac.fa -lib "brac"
#perl work_smskFa.pl -fa H1_brac500-871264.fa -lib "brac"
#perl work_smskFa.pl -fa H1_sorg.fa -species "sorghum bicolor"
#perl work_smskFa.pl -fa H1_bd.fa -species "brachypodium"
#perl work_smskFa.pl -fa H1_zea.fa -species maize -keep




