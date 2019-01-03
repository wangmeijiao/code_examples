

#if [ ! -f work_gffprofiling_fqpe.sh ] ;then
  #echo "not exsits, cp this file for you"
#  cp /home/mjwang/progs/epiSignalv0.1/mrnaSignalv1.0/work_gffprofiling_fqpe.sh .
#else
  #echo "already exists"
#fi

#if [ ! -d bin ] ; then
#  cp -r /home/mjwang/progs/epiSignalv0.1/mrnaSignalv1.0/bin/ .
#fi


#bash work.sh S10_RRA122016-V S10 > all.log 2>&1 &

data=$1
prefix=$2


bash work_nogff_fqpe.sh $data $prefix-tigr6_noG_g1m1_bowtie-n /home/mjwang/data/all_genomes/rice/tigr6.1/bowtieidx2_all/all_tigr6.fa /home/mjwang/data/all_genomes/rice/tigr6.1/all_tigr6.fa.size 500 > $prefix-tigr6_noG_g1m1_bowtie-n.log 2>&1 #without gff

bash work_gffprofiling_fqpe.sh $data /home/mjwang/data/all_genomes/rice/tigr6.1/all.gff3 ${prefix}_tigr6 /home/mjwang/data/all_genomes/rice/tigr6.1/bowtieidx2_all/all_tigr6.fa /home/mjwang/data/all_genomes/rice/tigr6.1/all_tigr6.fa.size 500 > $prefix-tigr6.log 2>&1  #with gff


