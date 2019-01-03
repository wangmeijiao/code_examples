catList <-
function(data=NULL,out="outfile"){
  if(!is.list(data)){stop("not a list\n")}
  for(name in names(data)){
    print (name)
    cat(name,data[[name]],"\n",file=out,append=T,sep="\t")
  }
  return("Your data has been written sucessfully")
}
drawDensiBar <-
function(data=NULL,title=NULL,color.list=NULL,ylim.list=NULL){
  par(mfrow=c(10,1),oma=c(0,0,3,0),mar=c(1,5,0,0))
  cnt <- 0
  for(cols in colnames(data)){
    cnt <- cnt+1
    if(cols=="k363"){next} 
    densi.mark<-data[,cols]
    range<-unlist(strsplit(ylim.list[cnt],","))
    min<- as.integer(range[1])
    max<-as.integer(range[2])
    if(max!=0){
      barplot(densi.mark,ylab=cols,border=color.list[cnt],space=c(0,0),names.arg="",ylim=c(min,max))
    }else{
          barplot(densi.mark,ylab=cols,border=color.list[cnt],space=c(0,0),names.arg="")
         }
  }
  title(main=title,outer=T)
}
drawDensiBarAll <-
function(data=NULL,title=NULL,color.list=NULL,ylim.list=NULL){
  par(mfrow=c(10,1),oma=c(0,0,3,0),mar=c(1,5,0,0))
  cnt <- 0
  for(cols in colnames(data)){
    cnt <- cnt+1
    #if(cols=="k363"){next} 
    densi.mark<-data[,cols]
    range<-unlist(strsplit(ylim.list[cnt],","))
    min<- as.integer(range[1])
    max<-as.integer(range[2])
    if(max!=0){
      barplot(densi.mark,ylab=cols,border=color.list[cnt],space=c(0,0),names.arg="",ylim=c(min,max))
    }else{
          barplot(densi.mark,ylab=cols,border=color.list[cnt],space=c(0,0),names.arg="")
         }
  }
  title(main=title,outer=T)
}
drawDensiHeatAll <-
function(data=NULL,title=NULL,ylim.list=NULL,ncols=100,type=NULL){
  par(mfrow=c(13,1),oma=c(0,0,6,0),mar=c(1,5,0,0))
  height<-rep(1,nrow(data))
  color.list<-colorRampPalette(c("white","black"))(ncols)
  result<-list(colors=color.list)
  cnt <- 0
  if(type=="tigr"){
    tags<-c("H3K4me2","H3K4me3","H3K36me3","H3K9ac","H4K12ac","H3K56ac","H3K92r1","H3K92r2","H3K27me1","H3K27me3","DNAmeth_CHG","mRNAWu_seedlingr1","mRNAWu_seedlingr2")
    orders<-c("k42","k43","k363","k9ac","h4k12ac","k56ac","k92r1","k92r2","k271","k273","DNAmeth_CHG","mRNAWu_seedlingr1","mRNAWu_seedlingr2")    
  }else if(type=="brachy"){
     tags<-c("H3K4me2","H3K27me2","H3K4me3","H3K9ac","H3K36me3","H3K9me2")
     orders<-c("k42","k272","k43","k9ac","k363","k92")
   }
  for(cols in orders){
    cnt <- cnt+1
    densi.mark<-data[,cols]
    min<-min(densi.mark)
    max<-max(densi.mark)
    range<-unlist(strsplit(ylim.list[cnt],","))
    low<- as.integer(range[1])
    high<-as.integer(range[2])
    print (c("min:",min,"max:",max,"low:",low,"high:",high))
    #filter
    if(high!=0){
      for( i in (1:length(densi.mark)) ){
       if(densi.mark[i]<low){densi.mark[i]<-0}else if(densi.mark[i]>low){densi.mark[i]<-densi.mark[i]-low}
      }
      min<-0
      max <- high-low      
      for( i in (1:length(densi.mark)) ){if(densi.mark[i]>max){densi.mark[i]<-max}}
    }
    #transform
    print(c(cols,"transformed min max:",min,max))
    index.color<-0
    for(j in (1:length(densi.mark))){
      n<-as.integer(ncols*(densi.mark[j]-min)/(max-min))
      if(n!=0){index.color[j]<-n}else{index.color[j]<-1}
    }
    result[[tags[cnt]]]<-index.color
    print(head(index.color,n=100))
    #draw Heatbox
    if(length(height)==length(index.color)){
      barplot(height,ylab="",col=color.list[index.color],border=NA,space=c(0,0),names.arg="",axes=F,cex.main=3)
      mtext(tags[cnt],side=2,cex=0.6,font=2,las=1)
    }else(stop("not equal"))
  }# mark for end
  title(main=title,outer=T)
  return(result)
}
drawDensiWinBar <-
function(data=NULL,win=500,norm=mean,title=NULL,color.list=NULL,scale=NULL){
  par(mfrow=c(10,1),oma=c(0,0,3,0),mar=c(1,5,0,0))
  cnt <- 0
  for(cols in colnames(data)){
    cnt <- cnt+1
    densi.mark<-data[,cols]
    len<-length(densi.mark)
    step<- as.integer(len/win)
    winStep <- gl(step,win)
    densiWinStep <- tapply(densi.mark[1:(win*step)],winStep,norm)
    barplot(densiWinStep,ylab=cols,col=color.list[cnt],border=color.list[cnt],space=c(0,0),names.arg="")
  }
  title(main=title,outer=T)
}
getDensiWin <-
function(data=NA,win=500,method=mean){
   long<-as.integer(length(rownames(data))/win)
   list<-data.frame(x=numeric(long)) 
   cnt<-0
   for(mark in colnames(data)){
      densi.mark<-data[,mark]
      len<-length(densi.mark)
      winN<-as.integer(len/win)
      winMode<-gl(winN,win)
      winDensi<-tapply(densi.mark[1:(win*winN)],winMode,method)
      list<-cbind(list,winDensi)
      cnt<-cnt+1
   }
   list<-list[,2:(cnt+1)]
   colnames(list)<-colnames(data)
   return (list)
}
getGrey <-
function(nCols=10){
  color.list<-0
  cnt<-0
  step=1/nCols
  for(i in seq(0+step,1,by=step)){cnt<-cnt+1;color.list[cnt]<-grey(i)}
  return(color.list)
}
medianSM <-
function(data=NA,win=3){ 
   if(win%%2==0||win<3){stop("input odd win number(>=3)\n")}
   mod=as.integer(win/2)
   list<-data.frame( x=numeric(nrow(data)) )
   for (mark in colnames(data)){
     len<-length(data[,mark])
     list[,mark]<-0
     for(i in (mod+1):(len-mod)){
       list[,mark][i]<-median(c(data[,mark][(i-mod):(i+mod)]))
     }
   }
   list<-list[,-1]
   return(list)
}
multiMax <-
function(data=NA){
  max.list<-0
  i<-0
  for(mark in colnames(data)){
   i<-i+1
   max.list[i]<-paste(mark,max(data[,mark]),sep=" ")
  }
  return (max.list)
}
testColor <-
function(color.list=NA,width=0.5,xmax=5){
  cnt<-0 
  barplot(height=rep(1,length(color.list)),width=rep(width,length(color.list)),col=color.list,border=NA,space=c(0,0),xlim=c(0,xmax),names.arg=color.list) 
  box()
  return(c(length(color.list),color.list))
}
triPolish <-
function(data=NULL){
 #Xn = 0.25X(n-1) + 0.5Xn + 0.25X(n+1)  
 list<-data.frame( x=numeric(nrow(data)) )
 for (mark in colnames(data)){
     densi<-data[,mark]
     len<-length(densi)
     list[,mark]<-0
     for(i in 2:(len-1)){
       list[,mark][i]<-0.25*densi[i-1]+0.5*densi[i]+0.25*densi[i+1]
     }
 }
 list<-list[,-1]
 return(list) 
}
callSignalMax <-
function(data=NA){   
   i<-0
   signal.max<-0
   for(mark in colnames(data)){
     i<-i+1
     points.del0<-data[,mark][data[,mark]!=0]
     result.background<-enhancedHist(points.del0)
     result.signal<-enhancedHist(result.background$out)
     signal.max[i]<-paste(result.signal$stats[5,1],mark,sep=" ")
   }
   return(signal.max)
}
enhancedHist <-
function(data=NA){
  layout(mat=matrix(c(1,2),2,1,byrow=TRUE),heights=c(1,3))
  par(mar=c(5.1,4.1,1.1,2.1))
  index.boxplot<-boxplot(data,horizontal=T,frame=F)
  index.hist<-hist(data,main="")
  return(c(index.boxplot,index.hist))
}
