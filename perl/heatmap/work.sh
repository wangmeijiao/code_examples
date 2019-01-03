


#  heatmap
#  heatmap.2
#  (heatmap.3)

#  myheatmap_image



#  pheatmap




#  matrixView_heatmap.pl
  
   cat test.matrixView.tab |perl matrixView_heatmap.pl -dimXY 10,20 -zlim 0,10|display


   perl matrixView_heatmap.pl -data s0d_vs_3d.cut.cut.txt -header T -dimXY 100,800 -zlim 0,80 > svg_out/out_s0dVS3d.100X800.filter0-80.svg


