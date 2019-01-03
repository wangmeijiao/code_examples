dir=$1
target=$2
query=$3

  for  i in $dir/*.lav; do
    echo $i
    lavToAxt $i -tfa $target -fa $query ${i}.axt
    axtChain -linearGap=loose ${i}.axt -faT $target -faQ $query ${i}.axt.chain
  done;
