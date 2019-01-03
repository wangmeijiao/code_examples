

###expression profiling for a given annotation set (and only for this set)
#one condition per time (with reps if possible, use comma to join reps)
#methods(FPKM+density):
#  tophat2->cufflinks
#         ->cuffquant
#         ->bam->density
#              ->bedgraph

echo "Begin at:"
date

#0,get options
data=$1  #<reads1_1[,...,readsN_1]>  <[reads1_2,...readsN_2]> NOT include .fq/.fastq
gff=$2 #reference gff, better has gene lines
prefix=$3
bowtieidx=$4
chrLen=$5
insert=$6
if [ ! -d $prefix  ] ;
  then mkdir $prefix
fi 

#1,tophat  support an exist gff
  ##common use
  #  --segment-length 20
  #  --no-coverage-search
  #  --max-intron-length 10000
  #  -o outdir
  ##for pair-ends 
  #  --no-mixed : when pair-mapped fails, try map one of the pair, this disable the attempt
  #  --no-discordant : if readpair mapped, only report concordant mapped ones( exclude, for example, ones mapped to different chr or with big insertion size)
  #  --mate-inner-dist $insert --mate-std-dev 50 
  ##deal with multi-hit reads
  #  --max-multihits 1 
  ##given gff only
  #  --transcriptome-only  
  #  <--transcriptome-index>  
  #  --transcriptome-max-hits 
  #  --prefilter-multihits

  
echo "start tophat2.."
/home/mjwang/progs/tophat2.0.13/tophat  -p 10  --no-coverage-search  --no-discordant --no-mixed  --mate-inner-dist $insert  --segment-length 20 --max-intron-length 10000 -o $prefix -G $gff $bowtieidx ${data}_1.fq.gz ${data}_2.fq.gz   #can read fq.gz
echo "done.."


#2, cufflinks (-G quantify given gff only); cuffquan for cuffnorm input
cd $prefix
cufflinks -p 10 -G $gff -o cufflink_out accepted_hits.bam

cat cufflink_out/isoforms.fpkm_tracking |perl /home/mjwang/progs/misc-tools/fpkmTools/extractFPKM.pl > cufflink_out/isoforms.fpkm_tracking.bed 

cat cufflink_out/genes.fpkm_tracking |perl /home/mjwang/progs/misc-tools/fpkmTools/extractFPKM.pl -log2 -add 1 > cufflink_out/genes.fpkm_tracking.add1log2.bed
#awk -vOFS='\t' '(NR>1){print $1,$2,$3,$5}' genes.fpkm_tracking.add1log2.bed > genes.fpkm_tracking.add1log2.bedgraph
#awk -vOFS='\t' '{mid=$2+int(($3-$2+1)/2);print $1,mid,mid+1,$4}' genes.fpkm_tracking.add1log2.bedgraph > genes.fpkm_tracking.add1log2.mid.bedgraph

cat cufflink_out/genes.fpkm_tracking |perl /home/mjwang/progs/misc-tools/fpkmTools/extractFPKM.pl > cufflink_out/genes.fpkm_tracking.bed
#awk -vOFS='\t' '(NR>1){print $1,$2,$3,$5}' genes.fpkm_tracking.bed > genes.fpkm_tracking.bedgraph
#awk -vOFS='\t' '{mid=$2+int(($3-$2+1)/2);print $1,mid,mid+1,$4}' genes.fpkm_tracking.bedgraph > genes.fpkm_tracking.mid.bedgraph

#cuffquant -p 12 -o cuffquant_out $gff accepted_hits.bam

#3, bam2bedgraph
echo "sort accepted_hits.bam.."
samtools sort accepted_hits.bam accepted_hits.sort
samtools index accepted_hits.sort.bam accepted_hits.sort.bam.bai
echo "start to bam2bedgraph.."
mkdir bam2bedgraph
cd bam2bedgraph
ln -s ../accepted_hits.sort.bam .
cp ~/progs/misc-tools/bam2coveragev0.1/bam2bedgraph_quick_forRNASEQ.sh .
cp ~/progs/misc-tools/tabix.sh .
bash  bam2bedgraph_quick_forRNASEQ.sh accepted_hits.sort.bam $chrLen 
#bash bam2wigQuick.sh accepted_hits.sort.bam $chrLen 
#bash bam2densi.sh accepted_hits.sort.bam $chrLen 
bash tabix.sh accepted_hits.sort.bedgraph


#4, clean

rm accepted_hits.sort.bedgraph
cd ../

rm accepted_hits.bam
rm unmapped.bam
cd ../

echo "all done, end at:"
date


