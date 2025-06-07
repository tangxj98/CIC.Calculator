# from the summary file, generate a count table for index, barcode , and index2 combination. Also simple QC for index1. 
# Notice that as PE libraries, we need to filter only those not match with the known index2
args <- commandArgs(trailingOnly = T)
if(length(args)!=4){
  cat( "Usage: Rscript preprocessing.R <input> <sample name> <spikein control> <sample & index2 mapping> \n")
}else{
  
  library(dplyr)
  name=args[2]
  print(name)

  keep.singleton=FALSE # set this to TRUE to keep singleton (i.e., count ==1). By default, remove them. 
 
  # read in spike-in control. 
  spikein.contol <- read.table(args[3], header=T,sep="\t",stringsAsFactors = F)
  spikein.contol$Spike.in.BC.short <- sapply(spikein.contol$Spike.in.Barcode.Sequence, function(x) substr(x,1,33) ) 

  sample.index <- read.table(args[4],header=T,sep="\t",stringsAsFactors=F) 
  print(head(sample.index[,1]))
  if(! name %in% sample.index[,1]){
	print("Sample Name not matched. Please double check!")
	quit()	
  }else{
	print(paste0("Index matched: ",paste(sample.index[match(name,sample.index[,1]),],collapse=":")) )
	index2=sample.index[match(name,sample.index[,1]),3]
  } 	  
  test<- read.table(args[1],sep="?",header=F,nrows=1)
  if(grepl("index",test[1,],ignore.case=T)){header=T}else{header=F}
  if(grepl(" ",test[1,])){sep=" "}else{sep="\t"}

  temp=read.table(args[1],header=header,sep=sep,stringsAsFactors = F)  
 	
  print(head(temp))
  temp <- temp[ grep(paste0(index2,"$"),temp[,1]),]
  print(dim(temp))
  print(paste("We have",nrow(temp),"lines."))
  if( nrow(temp)<= 58063 ){

	print(" Can be processed all at once!")
	indexes <- t(sapply(temp[,1], function(x) unlist(strsplit(x,","))))

  }else{

	# the file is too huge to be handeled directly. Need to split into mulitple matrix
	nblock= floor(nrow(temp)/20000)
	print(paste("Need to have", nblock+1,"blocks") ) 
	indexes <- as.data.frame(t(sapply(temp[1:20000,1], function(x) unlist(strsplit(x,",")))))
        print(dim(indexes))
	for(i in 1:nblock){
		start=i*20000+1
		end=ifelse( (i+1)*20000>nrow(temp),nrow(temp),(i+1)*20000 )
		print(paste("start:end=",start,";",end))
		test <- array()
		test <-  as.data.frame( t(sapply(temp[start:end,1], function(x) unlist(strsplit(x,",")))))
		print(dim(test))
		indexes=bind_rows(indexes,test)
  		print(paste("Now we have", nrow(indexes),"lines"))
        }
  }
  print(head(indexes))
  colnames(indexes) <- c("index.1","Barcode.Sequence","index.2")
  BC <- data.frame(indexes,count=temp[,2],stringsAsFactors=FALSE)
  print(table(BC$index.1,BC$count==1))
  if(!keep.singleton){
	print("remove singletons...")
	BC <- BC[ BC$count>1, ] 
  }  
  print(dim(BC))
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






