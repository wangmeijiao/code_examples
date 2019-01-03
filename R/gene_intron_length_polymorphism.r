


#library(RColorBrewer)
#boxcol <- brewer.pal(11, "Set3")
#boxcol <- boxcol[c(-1,-2)] #like perl shift
#boxcol[5] <- "#5c7a29"


#boxcol <- c( "#BEBADA", "#FB8072" ,"#80B1D3", "#FDB462" ,"#5c7a29", "#FCCDE5" ,"#D9D9D9", "#BC80BD" )
boxcol <- c( "#BEBADA", "#FB8072" ,"#80B1D3", "#FDB462" ,"#b4534b","#5c7a29", "#FCCDE5" ,"#D9D9D9", "#BC80BD")

cex<-2.5

#Gene Length per Species
source("gene_len.mat")
gene.len <- list('ara'=ara,'bd'=bd,'rice' = rice, 'maize' = maize, 'atau'=atau,'fly' = fly, 'zfish'=zfish,'mouse'=mouse,'human' = human)
#dev.new();stripchart(gene.len,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(gene.len,add=T,outline=F);title(main="gene.len")
dev.new();stripchart(gene.len,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(gene.len,add=T,outline=F);title(main="Gene Length per Species");abline(h=100000,lty=2,lwd=1.5,col='grey');text(x=0.9,y=125000,labels='100kb');abline(h=1000000,lty=2,lwd=1.5,col='grey');text(x=0.85,y=1020000,labels='1Mb')
savePlot("gene_len.tiff",type="tiff")

#Total Intron Length per Gene
source("intron_len_total.mat")
intron_len_total <- list('ara'=ara,'bd'=bd,'rice' = rice, 'maize' = maize, 'atau'=atau,'fly' = fly, 'zfish'=zfish,'mouse'=mouse,'human' = human)
#dev.new();stripchart(intron_len_total,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(intron_len_total,add=T,outline=F);title(main="intron_len_total")
dev.new();stripchart(intron_len_total,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(intron_len_total,add=T,outline=F);title(main="Total Intron Length per Gene");abline(h=100000,lty=2,lwd=1.5,col='grey');text(x=0.9,y=125000,labels='100kb');abline(h=1000000,lty=2,lwd=1.5,col='grey');text(x=0.85,y=1020000,labels='1Mb')
savePlot("intron_len_total.tiff",type="tiff")



#Total Exon Length per Gene
ara.exon <- gene.len$ara - intron_len_total$ara
bd.exon <- gene.len$bd - intron_len_total$bd
fly.exon <- gene.len$fly - intron_len_total$fly
human.exon <- gene.len$human - intron_len_total$human
mouse.exon <- gene.len$mouse - intron_len_total$mouse
zfish.exon <- gene.len$zfish - intron_len_total$zfish
maize.exon <- gene.len$maize - intron_len_total$maize
atau.exon <- gene.len$atau - intron_len_total$atau
rice.exon <- gene.len$rice - intron_len_total$rice
exon_len_total <- list('ara'=ara.exon,'bd'=bd.exon,'rice' = rice.exon, 'maize' = maize.exon, 'atau'=atau.exon, 'fly' = fly.exon, 'zfish'=zfish.exon,'mouse'=mouse.exon,'human' = human.exon)
#dev.new();stripchart(exon_len_total,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(exon_len_total,add=T,outline=F);title(main="exon_len_total")
dev.new();stripchart(exon_len_total,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(exon_len_total,add=T,outline=F);title(main="Total Exon Length per Gene");abline(h=100000,lty=2,lwd=1.5,col='grey');text(x=0.9,y=125000,labels='100kb');abline(h=1000000,lty=2,lwd=1.5,col='grey');text(x=0.85,y=1020000,labels='1Mb')
savePlot("exon_len_total.tiff",type="tiff")





#Total Intron Transposon Occupation Length per Gene
source("sum_te_len_total.mat")
sum_te_len_total  <- list('ara'=ara,'bd'=bd,'rice' = rice, 'maize' = maize, 'atau'=atau, 'fly' = fly, 'zfish'=zfish,'mouse'=mouse,'human' = human)
dev.new();stripchart(sum_te_len_total,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(sum_te_len_total,add=T,outline=F);title(main="Total Intron Transposon Occupation Length per Gene");abline(h=100000,lty=2,lwd=1.5,col='grey');text(x=0.9,y=125000,labels='100kb');abline(h=1000000,lty=2,lwd=1.5,col='grey');text(x=0.85,y=1020000,labels='1Mb')
savePlot("sum_te_len_total.tiff",type="tiff")
#mannual add n(gene) and GSIZE under xlab
#note rice is not more intron expansion than ara
#each points represent a gene

#Total Intron Transposon Occupation Percentage per Gene
source("te_len_total_perc.mat")
te_len_total_perc  <- list('ara'=ara,'bd'=bd,'rice' = rice, 'maize' = maize,  'atau'=atau,'fly' = fly, 'zfish'=zfish,'mouse'=mouse,'human' = human)
te_len_total_perc_cut<-list();for(i in names(te_len_total_perc)){  temp = unlist(te_len_total_perc[i]);for(j in 1:length(temp)){if(temp[j] > 1){temp[j]=1}}; te_len_total_perc_cut <-c(te_len_total_perc_cut,list(temp))};names(te_len_total_perc_cut) <- names(te_len_total_perc)
dev.new();stripchart(te_len_total_perc,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(te_len_total_perc,add=T,outline=F);title(main="te_len_total_perc");abline(h=100000,lty=2,lwd=1.5,col='grey');text(x=0.9,y=125000,labels='100kb');abline(h=1000000,lty=2,lwd=1.5,col='grey');text(x=0.85,y=1020000,labels='1Mb')
savePlot("te_len_total_perc.tiff",type="tiff")
dev.new();stripchart(te_len_total_perc_cut,vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(te_len_total_perc_cut,add=T,outline=F);title(main="Total Intron Transposon Occupation Percentage per Gene");abline(h=100000,lty=2,lwd=1.5,col='grey');text(x=0.9,y=125000,labels='100kb');abline(h=1000000,lty=2,lwd=1.5,col='grey');text(x=0.85,y=1020000,labels='1Mb')
savePlot("te_len_total_perc_cut.tiff",type="tiff")

#Transposon occupation percentage in each intron per gene
te_occ_perc_str <- readBigLine(f="te_occ_perc_str.mat",9)
dev.new();stripchart(te_occ_perc_str[c('ara','bd','rice', 'maize', 'atau','fly', 'zfish','mouse','human')],vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(te_occ_perc_str[c('ara','bd','rice', 'maize', 'atau','fly', 'zfish','mouse','human')],add=T,outline=F);title(main="Transposon occupation percentage in each intron per gene")
savePlot("te_occ_perc_str.tiff",type="tiff")

#All Intron Length per Species
intron_len_str <- readBigLine(f="intron_len_str.mat",9)
dev.new();stripchart(intron_len_str[c('ara','bd','rice', 'maize', 'atau', 'fly', 'zfish','mouse','human')],vertical=T, method='jitter',pch='.',col=boxcol, bg = 'bisque',cex=cex);boxplot(intron_len_str[c('ara','bd','rice', 'maize', 'atau','fly', 'zfish','mouse','human')],add=T,outline=F);title(main="All Intron Length per Species")
savePlot("intron_len_str.tiff",type="tiff")



####sub####
 readBigLine <- function(f=NULL,n=NULL){
  #low level readLines() to read very long data by lines, very fast than source
  con <- file(f,'r')
  line <- readLines(con,n=n)
  datlist <- list()
  nameid <- ''
  for (i in 1:n){
    spl <- unlist(strsplit(x=line[i],split=","))
    spl.data <- as.numeric(spl[2:length(spl)])
    spl.name <- spl[1]
    datlist <- c(datlist,list(spl.data))
    nameid <- c(nameid,spl.name)
  }
  close(con)
  #rm(i,con,line,spl,spl.data,spl.name) , no need to clean these local vectors because they will destroied after the function end
  nameid <- nameid[-1]
  names(datlist) <- nameid
  return(datlist)
}


