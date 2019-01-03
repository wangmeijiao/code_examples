#!/bin/bash
par="H=2000 Y=3400 L=6000 K=2200"
touch startBlastz.stamp
 for i in faSplit_T/* #only use chromosome for target.fa
 do mkdir -p lavAxtChain/`basename $i .fa`
    for j in faSplit_Q/*  #only use chromosome for query.fa
    do #echo $i,$j
     #echo $(ps|grep "blastz"|wc -l)
     while [ $(ps |grep "blastz"|wc -l) -gt 42 ]
      do sleep 1
     done
     echo "start to blastz $j to $i with parameter '$par',  running in background"
     blastz $i $j $par Q=HoxD55.q > lavAxtChain/`basename $i .fa`/`basename $i .fa`-`basename $j .fa`.lav &
    done
 done
touch endBlastz.stamp


