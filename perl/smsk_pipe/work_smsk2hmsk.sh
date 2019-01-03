  file=$1 
  prefix=`basename $file .smsk.fa`
  cat $file |perl -e 'while(<stdin>){if($_ =~/^>/){print $_; next};$_=~s/[atcgn]/N/g;print $_}' > $prefix.hmsk.fa


