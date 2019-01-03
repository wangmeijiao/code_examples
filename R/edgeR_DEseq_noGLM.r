

#################1 readin and filter#######################################

  data <- read.table("orthGene.fpkm.all",header=T,row.names=1,comment.char="")
  data.core <- data[,1:4]

#  data.core.filter2= data.core[ rowSums(data.core[,1:2] > 2 ) == 2 ,]  #each lib count MUST >2   lost genes like japo=0 brac=210 
#  n <- nrow(data.core.filter2) #13459 
#  lib1 <-  sum(data.core.filter2$cnt1) #13204060
#  lib2 <- sum(data.core.filter2$cnt2) #18166031

  all.HH <- read.table("all.HH.orthFpkm.out",header=F,row.names=1) #1087
#  all.EE <- read.table("all.EE.orthFpkm.out",header=F,row.names=1) #12819
#  all.EH <- read.table("all.EH.orthFpkm.out",header=F,row.names=1) #479
#  all.HE <- read.table("all.HE.orthFpkm.out",header=F,row.names=1) #1836
  all.AA <- read.table("all.AA.orthFpkm.out",header=F,row.names=1)  #16005

  #select count from data.core
  all.HH.count <- data.core[rownames(all.HH),] #1087
#  all.EE.count <- data.core[rownames(all.EE),]
#  all.EH.count <- data.core[rownames(all.EH),]
#  all.HE.count <- data.core[rownames(all.HE),]
  all.AA.count <- data.core[rownames(all.AA),]

  #filter out any of the raw count is less than 2
  all.HH.count.filter2 <- all.HH.count[rowSums(all.HH.count[,1:2] > 2) == 2 ,] #818
 # all.EE.count.filter2 <- all.EE.count[rowSums(all.EE.count[,1:2] > 2) == 2 ,] #10759
 # all.EH.count.filter2 <- all.EH.count[rowSums(all.EH.count[,1:2] > 2) == 2 ,] #373
 # all.HE.count.filter2 <- all.HE.count[rowSums(all.HE.count[,1:2] > 2) == 2 ,] #1413
  all.AA.count.filter2 <- all.AA.count[rowSums(all.AA.count[,1:2] > 2) == 2 ,] #13175


######2 edgeR############### 
  #edgeR use raw counts 
  library(edgeR)

  #test and plot for AA

  AA.counttable <- all.AA.count.filter2[,1:2]
  colnames(AA.counttable) <- c("japo","brac")  
  condition <- factor(c("japo","brac"))
  condition <- relevel(condition,ref='japo')

  e <- DGEList(count=AA.counttable,group=condition)
  e <- calcNormFactors(e)
  e.test <- exactTest(e,dispersion=0.1,pair=c("japo","brac"))
  e.test.table <- topTags(e.test,n=nrow(e), sort.by = 'p.value')$table
  e.test.table.sort <- e.test.table[order(e.test.table$FDR),]


  #e.HE <- e.test.table[rownames(all.HE.count.filter2),]
  #e.EH <- e.test.table[rownames(all.EH.count.filter2),]
  e.HH <- e.test.table[rownames(all.HH.count.filter2),]  #818
  #e.EE <- e.test.table[rownames(all.EE.count.filter2),]

  with(e.test.table,plot(logCPM,logFC,pch=20,main=expression("Genome-wide ortholog gene pair expression divergence \n     (O. sativa ssp. japonica vs O. brachyantha)" ),xlab='log Count per Million reads',ylab='log Fold Change',col=ifelse(FDR < 0.05, 'black','grey'), cex=0.5));
  points(e.HH$logCPM,e.HH$logFC, pch=20,col='red',cex=0.5); 
  #text(e.knob.sort$logCPM+0.2,e.knob.sort$logFC-0.1,col='black',label=1:nrow(e.knob.sort),cex=0.5); 
  legend("topright",legend=c("Ortholog gene pairs in core heterochromatin region","Diverged expression gene pairs "),col=c("red","black"),pch=20,bg='white');


 savePlot("edgeR-dispersion0.1-MAplot.HHinAA.tiff",type='tiff')


 e.HH.FDRbig0.05 <- e.HH[(e.HH$FDR > 0.05),] #650
# 1087 gene paired were identified located in the core heterochromatin with good synteny, of which 818 expressed gene pairs were expressed (raw reads > 2), 650 (79%) may be similarly expressed (edgeR, use difference expression threshold FDR < 0.05)







############################################################



###
#  with(e.test.table,plot(logCPM,logFC,pch=20,main="edgeR: Fold change vs abundance"))
#  with(subset(e.test.table,FDR<0.05),points(logCPM,logFC,pch=20,col='red')) #FDR < 0.05

#  e.knob <- e.test.table[grep("B",row.names(e.test.table),perl=T,value=F ),]
#  points(e.knob$logCPM,e.knob$logFC, pch=20,col='orange') #H1 points locate in middle part, CPM:counts-per-million
#  savePlot("edgeR-dispersion0.1-MAplot.tiff",type='tiff')

###

  e.knob <- e.test.table[grep("B",row.names(e.test.table),perl=T,value=F ),]

  e.knob.sort <- e.knob[c("Os04g03920|BU5","Os04g03980|BU6","Os04g03990|BU7","Os04g04000|BU8","Os04g04010|BU9","Os04g04020|BU10","Os04g04254|B1","Os04g04320|B2","Os04g05010|B10","Os04g05030|B11","Os04g05050|B12","Os04g05080|B15","Os04g06770|B21","Os04g08034|B26-1","Os04g08060|B26-2","Os04g08270|B29","Os04g08310|BD1","Os04g08320|BD2","Os04g08330|BD3","Os04g08340|BD4","Os04g08350|BD5"),]

#  with(e.test.table,plot(logCPM,logFC,pch=20,main=expression("Genome-wide ortholog gene pair expression divergence \n     (O. sativa ssp. japonica vs O. brachyantha)" ),xlab='log Count per Million reads',ylab='log Fold Change',col=ifelse(FDR < 0.05, 'black','grey'), cex=0.5));points(e.knob$logCPM,e.knob$logFC, pch=20,col='red',cex=0.5); text(e.knob$logCPM+0.2,e.knob$logFC-0.1,col='black',label=1:nrow(e.knob),cex=0.5); legend("topright",legend=c("Ortholog gene pairs in H1 region","Gene pairs with FDR <0.05"),col=c("red","black"),pch=20,bg='white');

  with(e.test.table,plot(logCPM,logFC,pch=20,main=expression("Genome-wide ortholog gene pair expression divergence \n     (O. sativa ssp. japonica vs O. brachyantha)" ),xlab='log Count per Million reads',ylab='log Fold Change',col=ifelse(FDR < 0.05, 'black','grey'), cex=0.5));points(e.knob.sort$logCPM,e.knob.sort$logFC, pch=20,col='red',cex=0.5); text(e.knob.sort$logCPM+0.2,e.knob.sort$logFC-0.1,col='black',label=1:nrow(e.knob.sort),cex=0.5); legend("topright",legend=c("Ortholog gene pairs in H1 region","Gene pairs with FDR <0.05"),col=c("red","black"),pch=20,bg='white');


  #1, Os04g03920|BU5; 2, Os04g03980|BU6; 3, Os04g03990|BU7; 4, Os04g04000|BU8; 5, Os04g04010|BU9; 6, Os04g04020|BU10; 7, Os04g04254|B1; 8, Os04g04320|B2; 9, Os04g05010|B10; 10, Os04g05030|B11; 11, Os04g05050|B12; 12, Os04g05080|B15; 13, Os04g06770|B21; 14, Os04g08034|B26-1; 15, Os04g08060|B26-2; 16, Os04g08270|B29; 17, Os04g08310|BD1; 18, Os04g08320|BD2; 19, Os04g08330|BD3; 20, Os04g08340|BD4; 21, Os04g08350|BD5;
  #1, JU5|BU5; 2, JU6|BU6; 3, JU7|BU7; 4, JU8|BU8; 5, JU9|BU9; 6, JU10|BU10; 7, J1|B1; 8, J2|B2; 9, J10|B10; 10, J11|B11; 11, J12|B12; 12, J15|B15; 13, J21|B21; 14, J26-1|B26-1; 15, J26-1|B26-2; 16, J29|B29; 17, JD1|BD1; 18, JD2|BD2; 19, JD3|BD3; 20, JD4|BD4; 21, JD5|BD5;

  #2,5,11 locate in FDR < 0.05   (not agree with Real time qPCR)
   #2, Os04g03980|BU6  5, Os04g04010|BU9  11, Os04g05050|B12
   #lost J11|B11 in genome-wide analysis
   #J3|B3 most variable in knob region   

savePlot("edgeR-dispersion0.1-MAplot.new.tiff",type='tiff')





#####DESeq#########
  data <- read.table("data.core.filter2",header=T,row.names=1,comment.char="")  #use raw filted count
  counttable <- data[,1:2]
  colnames(counttable) <- c("japo","brac")

  condition <- factor(c("japo","brac")) #no biological replicates
  condition <- relevel(condition,ref='japo')

  library('DESeq')
  countData <- newCountDataSet(counttable,condition)
  countData <- estimateSizeFactors( countData)
  sizeFactors(countData)

   ##estimate dispersion
   countData <- estimateDispersions(countData,method='blind',sharingMode='fit-only')
   plotDispEsts(countData, main="DESeq: Per-gene dispersion estimates")

   #do the nbinom test
   result <- nbinomTest(countData,"japo","brac")
   head(result)
   plotMA(result)
   hist(result$pval,breaks=100,col='skyblue',border='slateblue',main='')
   hist(result$padj,breaks=100,col='skyblue',border='slateblue',main='',ylim=c(0,40))
   plot(result$baseMean,result$log2FoldChange,pch=20,log='x',cex=0.3,col=ifelse(result$padj < 0.1, 'red','black'),xlab='log Base mean',ylab=expression(Log[2]~Fold~Change));abline(h=0,col='red',lwd=2)


   #scatter plot
   result.knob <- result[grep("B",result$id,perl=T,value=F ),]
   cor(log2(result$baseMeanA),log2(result$baseMeanB),method='pearson')
   #plot(log2(result$baseMeanA),log2(result$baseMeanB),pch=20,xlab="Normalized expression level(japo)",ylab="Normalized expression level(brac)",main="Scatter plot");points(log2(result.knob$baseMeanA),log2(result.knob$baseMeanB), pch=20,col='red');abline(0,1,col='grey');mtext(text="Pearson's r = 0.7789")
   plot(log2(result$baseMeanA),log2(result$baseMeanB),pch=20,xlab="Normalized expression level(japo)",ylab="Normalized expression level(brac)",main="Scatter plot",col='grey');points(log2(result.knob$baseMeanA),log2(result.knob$baseMeanB), pch=20,col='black');abline(0,1,col='grey',lty=2);mtext(text="Pearson's r = 0.7789")  


   savePlot("DESeq.scatterplot.logxy.knob.tiff",type='tiff')



