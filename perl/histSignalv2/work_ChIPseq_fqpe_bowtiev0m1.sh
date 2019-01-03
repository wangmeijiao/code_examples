#!/bin/bash
##########################################################
##     mapping pipeline for ChIP-seq fqpe               ##
##     Sep 11, 2013, modified at 6.27, 2018             ## 
##########################################################

echo "Begin at:"
date

#0. get options: input_data_file and prefix
 data=$1 #exclude _1.fq and _2.fq
 prefix=$2
 mkdir $prefix
 bowtieidx=$3
 chrLen=$4

#1. mapping by bowtie
 echo "start bowtie..(uniqual)"
 bowtie --al uniqual.fq --un unmapped.fq --max multi.fq -v 0 -m 1 --sam --threads 14 $bowtieidx -1 ${data}_1.fq -2 ${data}_2.fq hits_${prefix}.sam  2>bowtie.err 

#2. sam2bai
 sam2bai hits_${prefix}.sam
 rm hits_${prefix}.sam
 rm uniqual*.fq 

#3. transform pair mapped reads to single fragments (this process correct small burr points at density peaks, which are caused by double and missed calls. Like issue in FPKM caculation?)
  bamToBed -i hits_${prefix}.bam |perl bin/peBAM2seBED.pl 600|sort -k 1,1 -k2,2n > hits_${prefix}.pe2se.bed

#4, bam2bedgraph_veryquick
echo "start to bam2bedgraph_veryquick.."
mkdir bam2bedgraph
cd bam2bedgraph
ln -s ../hits_${prefix}.pe2se.bed .
cp /home/mjwang/progs/misc-tools/bam2coveragev0.1/bam2bedgraph_veryquick.sh .
cp ~/progs/misc-tools/tabix.sh .
bash  bam2bedgraph_veryquick.sh hits_${prefix}.pe2se.bed $chrLen
bash tabix.sh hits_${prefix}.pe2se.bedgraph
cd ..

##4. bed2density 
# echo "bed2density.. "
# genomeCoverageBed -d -i hits_${prefix}.pe2se.bed -g $chrLen|awk -vOFS='\t' '{print $1,$2,$2+1,$3}' |bgzip >hits_${prefix}.density.gz
# echo "creat tabix index.."
# tabix -s 1 -b 2 -e 3 hits_${prefix}.density.gz

##5.density to bw
# echo "density2bw.."
# bgzip -dc hits_${prefix}.density.gz |awk -vOFS='\t' '{if(!chrom[$1]){print "variableStep chrom="$1; chrom[$1]=1};print $2,$4}' | gzip > hits_${prefix}.wig.gz
# wigToBigWig hits_${prefix}.wig.gz $chrLen hits_${prefix}.bw
# rm hits_${prefix}.wig.gz


#6. clean work dir
 echo "clean.."
 mv  *.bed  *.err *.bw *.bam *.bai *.gz *.tbi unmapped*.fq multi*.fq bam2bedgraph $prefix

#7.summary 
#echo "Summary:"
#cat $prefix/bowtie.err 


echo "end at:"
date


