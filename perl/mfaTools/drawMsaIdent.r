  args <- commandArgs(TRUE)
  file <- args[1]
  outfile <- paste(file,".png",sep="")
  data <- read.table(file)
  png(outfile,width=13, height=2.5,units = "in",res=100)
  max.data <- max(data$V1)
  perc75 <- max.data * 0.75 
  plot(data$V1,type="l",ylim=c(0,max.data),ylab="bp_identity",xlab="position") 
  abline(h=perc75)
  gabarge <- dev.off()



