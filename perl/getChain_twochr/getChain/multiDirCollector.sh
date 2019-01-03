target=$1
query=$2
#cnt=0
for i in lavAxtChain/*; do 
  echo "gather lav->axt->chain, running in background.."
  echo $i
  core=`basename $i`
  bash enrollLAV2Chain.sh $i  $target $query >${core}.enroll.log 2>${core}.enroll.err &
  #(($cnt++))
  echo 
  #echo $i "done"
done

#wait until enrollLAV2Chain.sh all completed
while [ $(ps aux |grep enrollLAV2Chain.sh|grep -v grep|wc -l) -ne 0 ]
   do sleep 1
done
echo "wait for ending..."
sleep 5
echo "gathering data all done."
