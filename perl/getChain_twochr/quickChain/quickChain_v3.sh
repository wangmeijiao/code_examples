#get chain file for one chromosome link to the other
target=$1 #old version fa
base_t=`basename $target`
core_t=`basename $target .fa.smasked`
query=$2  #new version fa
base_q=`basename $query`
core_q=`basename $query .fa.smasked`


echo "start...",$target,$query
#1,split query fa(new.fa)
faSplit size $query 3000 ${base_q}.split3K -lift=./${base_q}.lft -oneFile

#2,blat query(new.fa) to targe(old.fa)
blat $target ${base_q}.split3K.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=95 -noHead -minScore=100  ${core_t}.${core_q}.psl

#3,liftup
liftUp -pslQ ${core_t}.${core_q}.liftup.psl ${base_q}.lft warn ${core_t}.${core_q}.psl

#4, chaining and sort chain file
axtChain -linearGap=medium  -faQ -faT -psl ${core_t}.${core_q}.liftup.psl $target $query ${core_t}.${core_q}.chain
#chainMergeSort *.chain | chainSplit . stdin
chainSort ${core_t}.${core_q}.chain ${core_t}.${core_q}.sort.chain


#5, netting
faSize -detailed $target > ${base_t}.size
faSize -detailed $query > ${base_q}.size
chainNet ${core_t}.${core_q}.sort.chain ${base_t}.size ${base_q}.size ${core_t}.${core_q}.net /dev/null

#6, extract liftOver 
netChainSubset ${core_t}.${core_q}.net ${core_t}.${core_q}.chain ${core_t}.${core_q}.liftOver

#7, clean up
echo "done."


