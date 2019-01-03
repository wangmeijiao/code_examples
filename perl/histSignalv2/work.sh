#bash  work_bracRS2.sh FCC62Y4ACXX_L4_CWHPE0141101005 k363 > pipe_bracRS2.log 2>&1 &

read=$1
mark=$2

#cp bash scripts
#cp -r ~/progs/epiSignalv0.1/histSingal_v0.1/bin/ .
#cp ~/progs/epiSignalv0.1/histSingal_v0.1/work_ChIPseq_fqpe_bowtiev0*.sh .

#bowtiev0M1y
#bash work_ChIPseq_fqpe_bowtiev0M1y.sh $read ${mark}_brachy1.5_bowtiev0M1 /home/mjwang/data/all_genomes/brachy1.5/bowtieidx_all/all_ffv1.5.fa  /home/mjwang/data/all_genomes/brachy1.5/all_ffv1.5.fa.sizes

#bowtiev0m1
bash work_ChIPseq_fqpe_bowtiev0m1.sh $read ${mark}_bracRS2_bowtiev0m1 /home/mjwang/dataex/AGI_data/Genomes/bracRS2/bowtie_idx/ObraRS2.AGAT02.fa  /home/mjwang/dataex/AGI_data/Genomes/bracRS2/ObraRS2.AGAT02.fa.size


