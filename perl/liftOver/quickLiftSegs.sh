
  grep "CDS" LOC_Os04g04390.segs| awk -vOFS='\t' '(NF >=9 ){$4-=(1700000+316001);$5-=(1700000+316001);print $0}' LOC_Os04g04390.segs |awk -vORS=',' '{print $4".."$5}'
