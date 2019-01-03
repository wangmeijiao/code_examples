pheatmap(mat=as.matrix(data),display_numbers= T,cluster_rows = T, cluster_cols = T,clustering_distance_rows = "correlation", clustering_distance_cols ="correlation")

#use annotation to group columns
#annotation is data.frame like this
#          desc
#Murder   type1
#Rape     type1
#Assault  type2
#UrbanPop type3

pheatmap(mat=as.matrix(data),display_numbers= T,cluster_rows = F, cluster_cols = T,annotation=anno1)
pheatmap(mat=as.matrix(data.1),display_numbers= T,cluster_rows = T, cluster_cols = T,cellwidth = 25, cellheight = 25,border_color="black",fontsize=12,clustering_distance_rows="correlation",clustering_distance_cols="correlation")
