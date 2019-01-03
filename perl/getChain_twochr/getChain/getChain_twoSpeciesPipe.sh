target=$1
query=$2

#pipeline for make chain file of two genomes
#input file: target.genome.smsk.fa query.genome.smsk.fa
#need file: HoxD55, which is a blastz score matrix
#note: change blastz pars in multiBlastzAll2All.sh
#0. init work dirs
echo "start pipe.."
date
echo "check dirs"
if [ ! -d faSplit_T ]; then
  mkdir faSplit_T
fi

if [ ! -d faSplit_Q ]; then
  mkdir faSplit_Q
fi

if [ ! -d lavAxtChain ]; then
  mkdir lavAxtChain
fi

if [ ! -d chainMergeSort ]; then
  mkdir chainMergeSort
fi

if [ ! -d netChain ]; then
  mkdir netChain
fi

if [ ! -f HoxD55.q  ]; then
  echo "HoxD55.q not exists"
  exit
fi

#2. split fa files
echo "split fa..."
faSplit byname $query faSplit_Q/
faSplit byname $target faSplit_T/

#3. run multiBlastz
echo "run multiBlastzAll2All.sh"
bash multiBlastzAll2All.sh > multiblastz.log 2>multiblastz.err

#4. lav->axt->chain
echo "run multiDirCollector.sh"
bash multiDirCollector.sh $target $query

#5. all.chain->net->liftOver
echo "run bash chainHarvest.sh"
bash chainHarvest.sh $target $query

echo
date
echo "all done"
