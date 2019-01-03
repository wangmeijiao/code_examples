target=$1
query=$2

chainMergeSort lavAxtChain/*/*.chain >chainMergeSort/all.chain
#rm *.tmp
#chainSplit chainMergeSort chainMergeSort/all.chain

faSize -detailed $target >${target}.size
faSize -detailed $query >${query}.size

chainPreNet chainMergeSort/all.chain ${target}.size ${query}.size netChain/all.prenet.chain

chainNet netChain/all.prenet.chain ${target}.size ${query}.size netChain/all.net /dev/null -minSpace=1 #query net goes to null

netChainSubset netChain/all.net netChain/all.prenet.chain netChain/all.liftOver

chainSort netChain/all.liftOver netChain/all.sorted.liftOver


