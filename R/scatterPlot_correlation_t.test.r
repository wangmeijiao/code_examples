

#################readin and filter#######################################

  data <- read.table("orthGene.fpkm.all.1",header=T,row.names=1,comment.char="")

  data.core <- data[,1:4]

  data.core.filter2= data.core[ rowSums(data.core[,1:2] > 2 ) == 2 ,]  #each lib count MUST >2   lost genes like japo=0 brac=210 
  n <- nrow(data.core.filter2) #13459 
  lib1 <-  sum(data.core.filter2$cnt1) #13204060
  lib2 <- sum(data.core.filter2$cnt2) #18166031


########normlization and output#############################

  #normalize by length (or log transform)

  data.core.filter2.normlen <- cbind(data.core.filter2$cnt1/data.core.filter2$len1,data.core.filter2$cnt2/data.core.filter2$len2)
  #data.core.filter2.normlen.log2 <- cbind(log2(data.core.filter2$cnt1/data.core.filter2$len1),log2(data.core.filter2$cnt2/data.core.filter2$len2))


  #norm by lib size ?


  #median norm (between libs,  only change scatter plot xlim/ylim, the point position remains, refer to Perry 2012)  

  med1 <- median(data.core.filter2.normlen[,1])
  med2 <- median(data.core.filter2.normlen[,2])
  medSyn <-  median(rowSums(data.core.filter2.normlen)/2)
  data.core.filter2.normlen.mednorm <- cbind( data.core.filter2.normlen[,1]*(medSyn/med1),data.core.filter2.normlen[,2]*(medSyn/med2) )

  #log transformation
  data.core.filter2.normlen.mednorm.log2 <- log2(data.core.filter2.normlen.mednorm)
  data.core.filter2.normlen.log2 <- log2(data.core.filter2.normlen)


  #add colnames and rownames
 
  colnames(data.core.filter2.normlen.mednorm) <- c("japo","brac")
  row.names(data.core.filter2.normlen.mednorm) <- row.names(data.core.filter2)

  colnames(data.core.filter2.normlen.mednorm.log2) <- c("japo","brac")
  row.names(data.core.filter2.normlen.mednorm.log2) <- row.names(data.core.filter2)

  colnames(data.core.filter2.normlen.log2) <- c("japo","brac")
  row.names(data.core.filter2.normlen.log2) <- row.names(data.core.filter2)

   #output table
   write.table(x=data.core.filter2.normlen.mednorm.log2,row.names=T,col.names=T,sep='\t',file="data.core.filter2.normlen.mednorm.log2",quote=F)
   write.table(x=data.core.filter2,row.names=T,col.names=T,sep='\t',file="data.core.filter2",quote=F)

  #extract H1 knob genes
 
   knob.data <- data.core.filter2.normlen.log2[grep("B",row.names(data.core.filter2.normlen.log2),perl=T,value=F),] #21 pairs remain
   row.names(knob.data) <- c("JU5|BU5","JU6|BU6","JU7|BU7","JU8|BU8","JU9|BU9","JU10|BU10","J1|B1","J2|B2","J10|B10","J11|B11","J12|B12","J15|B15","J21|B21","J26-1|B26-1","J26-2|B26-2","J29|B29","JD1|BD1","JD2|BD2","JD3|BD3","JD4|BD4","JD5|BD5")

  #read and process chr04Eu 

  chr04Eu <- read.table("japo_brac.chr04Eu.orthGeneExpreLevel.mat.1.cut",header=T,row.names=1,comment.char="") 
  chr04Eu.filter2 = chr04Eu[ rowSums(chr04Eu[,1:2] > 2 ) == 2 ,] #each lib count MUST > 2 17 pairs remain 

  chr04Eu.filter2.normlen <- cbind(chr04Eu.filter2$cnt1/chr04Eu.filter2$len1,chr04Eu.filter2$cnt2/chr04Eu.filter2$len2)
  chr04Eu.filter2.normlen.log2 <- log2(chr04Eu.filter2.normlen)
 
  colnames(chr04Eu.filter2.normlen.log2) <- c("japo","brac")
  row.names(chr04Eu.filter2.normlen.log2) <- row.names(chr04Eu.filter2)
  row.names(chr04Eu.filter2.normlen.log2) <- c("J1|B1","J10|B10","J12|B12","J13|B13","J14|B14","J15|B15","J16-2|B16-2","J16-3|B16-3","J17-2|B17","J18|B18","J2-2|B2","J20|B20","J3|B3","J4|B4","J6|B6","J8|B8","J9|B9")
   chr04Eu.filter2.normlen.log2.sort <- chr04Eu.filter2.normlen.log2[c("J1|B1","J2-2|B2","J3|B3","J4|B4","J6|B6","J8|B8","J9|B9","J10|B10","J12|B12","J13|B13","J14|B14","J15|B15","J16-2|B16-2","J16-3|B16-3","J17-2|B17","J18|B18","J20|B20"),])

   write.table(x=knob.data,row.names=T,col.names=T,sep='\t',file="knob.data.filter2.normlen.log2",quote=F)
   write.table(x=chr04Eu.filter2.normlen.log2,row.names=T,col.names=T,sep='\t',file="chr04Eu.filter2.normlen.log2",quote=F)


#################### H1, chr04Eu plotting########################################

   par(mar=c(10,5,5,5),oma=c(1,0,0,0))
   plot(1:21,log10(knob.data[,1]/knob.data[,2]),type='l',ylim=c(-10,10), xlab="",ylab="",main="Orthologous gene expression variation \nin H1 region",xaxt='n',mgp=c(3,1,0))
   axis(side=1,at=1:21,labels=row.names(knob.data),las=2,mgp=c(0,2,1))
   title(xlab="Orthologous gene pairs",line=7,cex.lab=1.2,family="Calibri Light")
   title(ylab="Gene expression level Log Fold Change",line=3,cex.lab=1.2,family="Calibri Light")
   points(1:21,log10(knob.data[,1]/knob.data[,2]),col="orange",pch=17)
   
   result.bxstats <- boxplot.stats(log10(knob.data[,1]/knob.data[,2]))
   rect(xleft=0,xright=21,ytop=result.bxstats$stats[4],ybottom=result.bxstats$stats[2],col=rgb(1,0,0,0.2),border=F)
   text(x=21.3,y=result.bxstats$stats[4]+0.5,labels="Q3")
   text(x=21.3,y=result.bxstats$stats[1]-0.2,labels="Q1")
   
   res.ttest <- t.test(x=knob.data[,1],y=knob.data[,2])
   text(x=15,y=8,labels=bquote(t-test~p.value~"="~.(res.ttest$p.value)) )

   savePlot("orthologous_gene_expression_variation-extractH1.tiff",type='tiff')

 
   dev.new()
   par(mar=c(10,5,5,5),oma=c(1,0,0,0))
   plot(1:17,log10(chr04Eu.filter2.normlen.log2.sort[,1]/chr04Eu.filter2.normlen.log2.sort[,2]),type='l',ylim=c(-10,10), xlab="",ylab="",main="Orthologous gene expression variation \nin chr04Eu region",xaxt='n',mgp=c(3,1,0))
   axis(side=1,at=1:17,labels=row.names(chr04Eu.filter2.normlen.log2.sort),las=2,mgp=c(0,2,1))
   title(xlab="Orthologous gene pairs",line=7,cex.lab=1.2,family="Calibri Light")
   title(ylab="Gene expression level Log Fold Change",line=3,cex.lab=1.2,family="Calibri Light")
   points(1:17,log10(chr04Eu.filter2.normlen.log2.sort[,1]/chr04Eu.filter2.normlen.log2.sort[,2]),col="orange",pch=17)

   result.bxstats <- boxplot.stats(log10(chr04Eu.filter2.normlen.log2.sort[,1]/chr04Eu.filter2.normlen.log2.sort[,2]))
   rect(xleft=0,xright=17,ytop=result.bxstats$stats[4],ybottom=result.bxstats$stats[2],col=rgb(1,0,0,0.2),border=F)
   text(x=17.3,y=result.bxstats$stats[4]+0.5,labels="Q3")
   text(x=17.3,y=result.bxstats$stats[1]-0.2,labels="Q1")
   
   res.ttest <- t.test(x=chr04Eu.filter2.normlen.log2.sort[,1],y=chr04Eu.filter2.normlen.log2.sort[,2])
   text(x=12,y=8,labels=bquote(t-test~p.value~"="~.(res.ttest$p.value)) )


   savePlot("orthologous_gene_expression_variation-chr04Eu.tiff",type='tiff')
 




################genome wide plotting########################

  #description plots : histogram and scatterplot and correlation
  boxplot(data.core.filter2.normlen.mednorm,outline=F)

  plot(data.core.filter2.normlen.mednorm,cex=.5)
  plot(data.core.filter2.normlen.mednorm.log2,cex=.5)
  plot(data.core.filter2.normlen.mednorm.log2,cex=.5);points(data.core.filter2.normlen.mednorm.log2[grep("B",row.names(data.core.filter2.normlen.mednorm.log2),perl=T,value=F ),], col='red',cex=0.5);abline(0,1,col='grey')
  cor(data.core.filter2.normlen.mednorm.log2[,1],data.core.filter2.normlen.mednorm.log2[,2],method='pearson') #0.7567166
  cor(data.core.filter2.normlen.mednorm.log2[,1],data.core.filter2.normlen.mednorm.log2[,2],method='spearman') #0.7775936


  #lm regression and ANOVA
  data.core.filter2.normlen.mednorm.log2.df = as.data.frame(x=data.core.filter2.normlen.mednorm.log2)
  lm.sol <- lm(formula=japo~brac,data=data.core.filter2.normlen.mednorm.log2.df)
  summary(lm.sol)
  abline(coefficients(lm.sol))
  points(y=fitted(lm.sol),x=residuals(lm.sol),pch=20,col='navy') #residuals vs fitted values
  dev.new();par(mfrow=c(2,2));plot(lm.sol)

  anova(lm.sol)



   #density distribution histogram overlap plot
   hist(data.core.filter2.normlen.mednorm.log2[,1],breaks=60,col=rgb(0.8,0.8,0.8,0.3),xlim=c(-10,10),ylim=c(0,1500),ylab="count",xlab="normalized gene expression",main="Overall gene expression profile");hist(data.core.filter2.normlen.mednorm.log2[,2],col=rgb(0.1,0.1,0.1,0.3), breaks=60,xlim=c(-10,10),ylim=c(0,1500),add=T);legend(x=6,y=1500,legend=c("japo","brac"),fill=c(rgb(0.8,0.8,0.8,0.3),rgb(0.1,0.1,0.1,0.3)))


   #MAplot  
   plot(y=log2(data.core.filter2.normlen.mednorm[,1]/data.core.filter2.normlen.mednorm[,2]),x=rowMeans(data.core.filter2.normlen.mednorm.log2),cex=0.3  )
   abline(h=0,col='grey')
   knob.data <- data.core.filter2.normlen.mednorm.log2[grep("B",row.names(data.core.filter2.normlen.mednorm.log2),perl=T,value=F ),]
   points(y=knob.data[,1]-knob.data[,2],x=rowMeans(knob.data),cex=0.3,col='red'  )


  #pairwise test
  t.test(data.core.filter2.normlen.mednorm.log2[,1],data.core.filter2.normlen.mednorm.log2[,2],paired=T) #p-value = 4.968e-13
  wilcox.test(data.core.filter2.normlen.mednorm.log2[,1],data.core.filter2.normlen.mednorm.log2[,2],paired=T) #p-value = 0.07357

  knob.data.diff <- knob.data[,1] - knob.data[,2]
  chr04Eu.filter2.normlen.log2.diff <- chr04Eu.filter2.normlen.log2[,1] - chr04Eu.filter2.normlen.log2[,2]

  

 #QQ-plot: fit normal distribution: log transformation can yield approximate normal distribution
 x <- data.core.filter2.normlen.mednorm.log2[,1]

 library("fitdistrplus")
 fit <- fitdist(x,'norm')
 plot(fit) #plot method in this package

 para <- fit$estimate
 hist(x,prob=T)
 curve(expr=dnorm(x,para[1],para[2]),col=2,add=T)

 x.znorm <- (x-mean(x))/sd(x)
 qqnorm(x.znorm)
 abline(0,1,col='grey')









