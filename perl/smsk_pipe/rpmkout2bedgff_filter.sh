
prefix=$1  # prefix:  out_{Sbicolor_255_v2.0}-smsk  = {prefix}
type=$2
filter=$3  #not used 

if [ $type = "bed" ]; then
  echo "do $type"
  outdir=beds_$prefix
  cat out_${prefix}-smsk/${prefix}.fa.out |perl bin/rpmk2gffv2.3.3.pl -bed -prefix $prefix -filter $filter -outdir $outdir 2> out2bed.err

  cat $outdir/*.bed > $outdir/${prefix}_Allrepeats.bed
  chmod u+x $outdir/${prefix}_Allrepeats.bed

  mv out2bed.err $outdir
  mv ${prefix}_lost.bed $outdir
  chmod u+x $outdir/${prefix}_lost.bed

elif [ $type = "gff" ]; then
  echo "do $type"
  outdir=gffs_$prefix
  cat out_${prefix}-smsk/${prefix}.fa.out |perl bin/rpmk2gffv2.3.2.pl  -prefix $prefix -filter $filter -outdir $outdir 2> out2gff.err

  cat $outdir/*.gff > $outdir/${prefix}_Allrepeats.gff
  chmod u+x $outdir/${prefix}_Allrepeats.gff

  mv out2gff.err $outdir
  mv ${prefix}_lost.gff $outdir
  chmod u+x $outdir/${prefix}_lost.gff

fi
