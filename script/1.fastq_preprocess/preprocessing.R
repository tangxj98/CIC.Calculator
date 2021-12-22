# from the summary file, generate a count table for index, barcode , and index2 combination. Also simple QC for index1. 
args <- commandArgs(trailingOnly = T)
if(length(args)!=2){
  cat( "Usage: Rscript preprocessing.R <input> <output prefix>\n")
}else{
  
  name=args[2]
  print(name)
  # read in spike-in control. 
  spikein.contol <- read.table("../data/spikin.control.txt",header=T,sep="\t",stringsAsFactors = F)
  spikein.contol$Spike.in.BC.short <- sapply(spikein.contol$Spike.in.Barcode.Sequence, function(x) substr(x,1,33) ) 

  
  temp=read.table(args[1],header=F,sep=" ",stringsAsFactors = F)  
  indexes <- t(sapply(temp[,1], function(x) unname(unlist(strsplit(x,",")))))
  #colnames(spikein.contol)
  colnames(indexes) <- c("index.1","Barcode.Sequence","index.2")
  BC <- data.frame(indexes,count=temp[,2],stringsAsFactors=FALSE)

  # check the length of Barcode.Keep only those with 33 bp 
  BC.length <- sapply(BC$Barcode.Sequence, function(x) nchar(x) ) 
  table(BC.length)
  BC <- subset(BC, BC.length==33)
  BC$Is.Spikein <- ifelse(BC$Barcode.Sequence %in% spikein.contol$Spike.in.BC.short,"Yes","No")
  index1.Count<- aggregate(BC$count,by=list(BC$index.1),FUN=sum)
  colnames(index1.Count) <- c("index.1","Copies")
  index1.Count$pct <- signif(index1.Count$Copies*100/sum(index1.Count$Copies),digits = 3)

  pdf(paste0(name,".index1.QC.pdf") ) 

  barplot(index1.Count$Copies/10^6,names.arg = index1.Count$index.1,las=2,ylab="Copies (Million)", ylim=c(0,85),width = 1,space = 0)
  text(1:nrow(index1.Count)-0.5, index1.Count$Copies/10^6, labels = paste0( signif(index1.Count$pct,digits = 3),"%" ),pos = 3)
  box()
  dev.off()
  write.table(BC,paste0(name,".mergeinput.txt"),quote=F,sep="\t")
}






