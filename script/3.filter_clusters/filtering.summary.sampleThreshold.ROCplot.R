# This one is going to normalize each sample on it is own. 

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){
  cat("Usage:Rscript filtering.summary.R <clustered file> <output prefix>\n")
}else{

library(tidyr)
library(pheatmap)
library(investr,lib.loc="~/staff_analysis/R")
file=args[1]
filename=tail(unlist(strsplit(file,"/")),n=1)
name=args[2]
print(name)

spikein.contol <- read.table("spikein.control.txt",header=T,sep="\t",stringsAsFactors = F)
spikein.contol$Spike.in.BC.short <- sapply(spikein.contol$Spike.in.Barcode.Sequence, function(x) substr(x,1,33) ) 

 
# read clustered file  e.g. BC_Run1_04.cluster.txt

index1.coding <- read.table("index1.coding.txt",header=F,sep="\t",stringsAsFactors = F)
rownames(index1.coding) <- index1.coding$V2
index2.coding <- read.table("index2.coding.txt",header=F,sep="\t",stringsAsFactors = F)
rownames(index2.coding) <- index2.coding$V2

# debug
print(head(index1.coding))
print(head(index2.coding))

#allcluster <- read.table(file,header=F,sep="\t",stringsAsFactors = F)
test<- read.table(file,sep="?",header=F,nrows=1)
if(grepl("index",test[1,],ignore.case=T)){header=T}else{header=F}

allcluster <- read.table(file,header=header,sep="\t",stringsAsFactors = F)

colnames(allcluster) <- c("index.1","Barcode.Sequence","index.2","count","Is.Spikein","cluster","FirstBC")
allcluster$Is.Spikein2 <- ifelse(allcluster$FirstBC %in% spikein.contol$Spike.in.BC.short,"Yes","No")

#print(head(allcluster))

pdf(paste(name,"withinSample.filtering.pdf",sep = "."),width=9,height=9 )  
par(mfrow=c(3,3))
# before filtering the distribution and percentage: 
allcluster.counts.index.stat <- aggregate(allcluster$count,by=list(allcluster$index.1,allcluster$index.2), FUN=sum)
colnames(allcluster.counts.index.stat) <- c("Index.1","Index.2","count")
#print(head(allcluster.counts.index.stat))

allcluster.counts.index.stat.matrix <- spread(allcluster.counts.index.stat,key = "Index.2", value = count)

nas <- which(is.na(allcluster.counts.index.stat.matrix),arr.ind = T) 
print(paste("there are",length(nas),"NAs"))


if(nrow(nas)>0){
  cat("###################\n")
  cat("The following samples have no counts: \n")
  for(i  in 1:nrow(nas)){
    cat("\t",allcluster.counts.index.stat.matrix[nas[i,1],1],"\t",colnames(allcluster.counts.index.stat.matrix)[nas[i,2]],"\n")
  }
  cat("###################\n")
}
allcluster.counts.index.stat.matrix[ is.na(allcluster.counts.index.stat.matrix)] <- 0
rownames(allcluster.counts.index.stat.matrix) <- allcluster.counts.index.stat.matrix[,1]
#print(allcluster.counts.index.stat.matrix)
allcluster.counts.index.stat.matrix <- data.matrix(allcluster.counts.index.stat.matrix[,-1])
allcluster.FRV.index.stat.matrix <- round(apply(allcluster.counts.index.stat.matrix, 2, function(x) x/sum(x)),digits = 3)


write.table(allcluster.FRV.index.stat.matrix,paste0(name,".beforeFiletering.stat.matrix.txt"),quote=F,sep="\t")


index1 <- unique(allcluster$index.1)
index2 <- unique(allcluster$index.2)
print(index1)
print(index2)
#filtered FRV, filtered counts 
nonsp.id1=index1[index1!="GTCA"]
filtered.counts <- matrix(data=NA,nrow=length(nonsp.id1),ncol=length(index2))
rownames(filtered.counts) <- nonsp.id1 
colnames(filtered.counts) <- index2

filtered.FRV <- matrix(data=NA,nrow=length(nonsp.id1),ncol=length(index2))
rownames(filtered.FRV) <- index1[index1!="GTCA"]
colnames(filtered.FRV) <- index2

# final output of spikein data (will append every index2 in the same data frame)
spikein.final <- data.frame()

# final output of ns data ( will append every index2 as a list)
nonspikein.final <- list()
print(index2)
for ( id2 in index2) {
#id2="CGATC"
  cluster <- subset(allcluster, allcluster$index.2==id2)
  
  #head(cluster)
  cat("##############################################\n\n")
  cat("Start processing sample ",id2,"!\n")
  cat("index2:", id2,"are", round(100*nrow(cluster)/nrow(allcluster),digits = 2),"% of the whole cluster\n")
  cat("**********************************************\n\n")
  cat("Start processing the spike-in control index GTCA in sample",id2,"...\n")
  
  # 1. separate the spikein and non-spikein 
  spikein_cluster <- subset(cluster, index.1=="GTCA")
  cat("We have ",length(unique(spikein_cluster$FirstBC))," different barcodes in this spike-in index.\n")
  sp <- subset(spikein_cluster,Is.Spikein2=="Yes")
  ns <- subset(spikein_cluster,Is.Spikein2=="No")
  cat("We have ",length(unique(sp$FirstBC))," different barcodes in the spike-in index.\n")
  cat("We have ",length(unique(ns$FirstBC))," different barcodes in the non-spike-in index.\n")
  
  cat("percentge of spike-in in the GTCAs: ", round(100*sum(sp$count)/sum(spikein_cluster$count),digits = 2),"%\n")
  head(sp)
  
  # 2. get aggregated counts a) by the BC sequence b) by sample (index 2) for noralization (FRV)
  # a) BC seq
  sp.output <- aggregate(sp$count, by=list(sp$FirstBC), FUN=sum)
  colnames(sp.output) <- c("FirstBC","count")
  # head(sp.output)
  #rownames(sp.output)  <- paste0(sp.output$index.2,sp.output$output.bc.seq)
  
  sp.output$Input.Copies <- sapply(sp.output$FirstBC, function(x) spikein.contol[which(spikein.contol$Spike.in.BC.short==x), "Expected.Copies"])

  ns.output <- aggregate(ns$count, by=list(ns$FirstBC), FUN=sum)
  colnames(ns.output) <- c("FirstBC","count")

  cat ("max counts in non-spike-in clusters:", max(ns.output$count),"\n")
  cat ("length of unique counts for ROC performance",length(unique(ns.output$count)),"\n")
 
  
  # # b) by sample (index2)

  sp.output$FRV <- sp.output$count/sum(sp.output$count)
  sp.output.rm1 <- subset(sp.output,Input.Copies>1)
  
  ns.output$FRV <- ns.output$count/sum(ns.output$count)
  
  
  reg <- lm(data = sp.output, formula = log10(FRV)~log10(Input.Copies))
   plot(log10(sp.output$Input.Copies), log10(sp.output$FRV),xlab="log10(Input Copies)", ylab="log10(FRV)",pch=19,ylim=c(-6,2),main=id2)
   abline(reg, lty=2, lwd=2, col="red")
   mtext(paste("R squared=", round(summary(reg)$r.squared, digits=3), sep=" "), adj=1)
   mtext(paste("y=", round(summary(reg)$coefficients[2,1], digits=4),"x", round(summary(reg)$coefficients[1,1], digits=4),   sep=" "), adj=0)
  reg2 <- lm(data = sp.output.rm1, formula = log10(FRV)~log10(Input.Copies))
   abline(reg2, lty=2, lwd=2, col="blue")
   
   mtext(paste("R squared=", round(summary(reg2)$r.squared, digits=3), sep=" "), adj=1,line=1)
   mtext(paste("y=", round(summary(reg2)$coefficients[2,1], digits=4),"x", round(summary(reg2)$coefficients[1,1], digits=4),   sep=" "), adj=0,line=1)
  

  
  
  noise.range <- range(ns.output$FRV)
  ROC <- matrix(nrow=0,ncol=2)
  colnames(ROC) <- c("FPR", "TPR")
  thresholds<-sort(unique(ns.output$FRV))
  for(i in thresholds){
    
    TPR<-length(which(sp.output.rm1$FRV>i))/nrow(sp.output.rm1)
    FPR<-length(which(ns.output$FRV>i))/nrow(ns.output)
    ROC<-rbind(ROC, c(FPR,TPR))
  }
  
  # The object function is maximizing TPR and minimizing FPR. i.e. TPR-FPR. 
  ROC <- data.frame(ROC)
  
  if(nrow(ROC)<=3){
    cat ("WARNINGS: No spike-in cluster left. Poor quality.", id2, " Skipped\n")
    next
  } 
  
  
  ROC$threshold = thresholds
  ROC$obj=ROC$TPR-ROC$FPR
  best.cutoff <- ROC[ which.max(ROC$obj),"threshold"]
  best.cutoff
  #[1] 5.863534e-05
  sugg.threshold=ROC[which(ROC$threshold==best.cutoff),]
  plot(ROC[,1:2],type="b" , xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.4, ylab="True Positive Rate", xlab="False Positive Rate",main=id2)
  abline(coef=c(0,1), lty=2, lwd=0.5)
  points(ROC[which(ROC$threshold==best.cutoff),1:2], pch=16, cex=2, col="red")
  #mtext(paste("Suggested Threshold (FRV)=", signif(sugg.threshold[,"threshold"], digits=4)), line=0, adj=1)
  #mtext(paste("Suggested TPR=", signif(sugg.threshold[,"TPR"], digits=4)), line=1, adj=0)
  #mtext(paste("Suggested FPR=", signif(sugg.threshold[,"FPR"], digits=4)), line=0, adj=0)
  
  
  
  # evalutate filtering. 
  noise_remaining<- table(ns.output$FRV>best.cutoff)[2]
  sp.output.rm1.filtered <- subset(sp.output.rm1,FRV > best.cutoff)
  ns.output.filtered <- subset(ns.output, FRV> best.cutoff)
  
  
  if(nrow(sp.output.rm1.filtered)<=3){
    cat ("WARNINGS: No spike-in cluster left. Poor quality.", id2, " Skipped\n")
    next
  } 
  # redo the lm fitting. 
  reg3 <- lm(data = sp.output.rm1.filtered, formula = log10(FRV)~log10(Input.Copies))
  # predict the results. 
  sp.output.rm1.filtered$predicted.input <- round(sapply(log10(sp.output.rm1.filtered$FRV),function(x) invest(reg3,x,upper = 10000000,lower=0,interval = "none") ),digits = 1)

  if(nrow(ns.output.filtered)==0){
    ns.output.filtered <- data.frame(FirstBC="NA",count=0, Input.Copies=0,FRV=0, predicted.input=0)    
  }else{
    ns.output.filtered$predicted.input <- round(sapply(log10(ns.output.filtered$FRV),function(x) invest(reg3,x,upper = 10000000,lower=0,interval = "none") ),digits = 1)
  }
  # predict minumum input cells 
  minimum.input.cell=round(invest(reg3,log10(best.cutoff),upper = 10000000,lower=0,interval = "none") )

  
   hist(log10(ns.output$FRV), breaks = 100, xlab="log10(FRV)",main=paste0(id2, " unwanted barcode counts after filtering % =",round(sum(ns.output.filtered$count)*100/sum(spikein_cluster$count),digits = 2),"%") ) 
   abline(v=log10(best.cutoff),col="red",lty=4,lwd=2)
  

  
 
  # filter the rest 5 index.1 using the cutoff. 
  cat("\n\n******************************\n")
  cat("Start filtering the non-spike-in barcodes....\n")
  nonspikein_cluster <- subset(cluster, index.1!="GTCA")
 
  # calculate normalized FRV (devided by the aggregate counts of index2 no matter how many index1 since it is for each sample (index2) ) 
  # check the spike-in barcode portion. should be small. 
  cat("\nPercentage of spike-in barcode in other samples: ")
  print(sum(nonspikein_cluster[ nonspikein_cluster$FirstBC %in% spikein.contol$Spike.in.BC.short , "count"])/sum(nonspikein_cluster$count))
  cat("\n")
  
  # merge down to index1, index2, by cluster.  
  nonspikein_cluster <- subset(nonspikein_cluster, !(FirstBC %in% spikein.contol$Spike.in.BC.short))
  nonsp.output <- aggregate(nonspikein_cluster$count,by=list(nonspikein_cluster$index.1,nonspikein_cluster$index.2,nonspikein_cluster$FirstBC), FUN=sum)
  colnames(nonsp.output) <- c("index.1","index.2","output.bc.seq","count")
  head(nonsp.output)
  # calculate sample size and normalized FRV 
  
  nonsp.output$FRV <- nonsp.output$count/sum(nonsp.output$count)
  dim(nonsp.output)
  
  # filtering: remove small outputBC, keep only those FRV > threshold. 
  nonsp.output.filtered <- subset(nonsp.output, FRV>best.cutoff)
  
  head(nonsp.output.filtered)
  dim(nonsp.output.filtered)
  
  
  if(nrow(nonsp.output.filtered)>0){
    nonsp.output.filtered$predicted.input <- round(sapply(log10(nonsp.output.filtered$FRV),function(x) invest(reg3,x,upper = 10000000,lower=0,interval = "none") ),digits = 1)
  }else{
    nonsp.output.filtered <- data.frame(index.1="NA",index.2="NA",output.bc.seq="NA", count=0,FRV=0,predicted.input=0)
  }
  
  
  nonsp.output.filtered.clone.stat <- aggregate(nonsp.output.filtered$count,by=list(nonsp.output.filtered$index.1),FUN=sum)
  colnames(nonsp.output.filtered.clone.stat) <- c("Index.1","counts")
  nonsp.output.filtered.clone.stat$FRV <- nonsp.output.filtered.clone.stat$counts/sum(nonsp.output.filtered.clone.stat$counts)
  rownames(nonsp.output.filtered.clone.stat) <- nonsp.output.filtered.clone.stat[,1]

  
  threshold.cell <- round(best.cutoff*sum(sp.output$count))
  temp <- data.frame(index1="GTAC",index2=id2,sp.output.rm1.filtered,threshold=signif(best.cutoff,digits = 3),threshold.cell,minimum.input.cell )
  temp$index1 <- index1.coding[  as.character(temp$index1),"V1" ]
  temp$index2 <- index2.coding[  as.character(temp$index2),"V1" ]
  
  spikein.final <- rbind(spikein.final,temp[order(temp[,"predicted.input"],decreasing = T),])
  
  threshold.cell <- round(best.cutoff*sum(nonsp.output$count))
  temp <- data.frame(nonsp.output.filtered,threshold=signif(best.cutoff,digits = 3),threshold.cell,minimum.input.cell)
  temp$index.1 <- index1.coding[  as.character(temp$index.1),"V1" ]
  temp$index.2 <- index2.coding[  as.character(temp$index.2),"V1" ]
  
  nonspikein.final[[id2]] <- temp[order(temp[,"predicted.input"],decreasing = T),]
  
  write.table(nonspikein.final[[id2]],file=paste0(name,".",id2,".non-spikein.output.txt"),quote=F,sep="\t",row.names = F)
  
  filtered.FRV[,id2] <- nonsp.output.filtered.clone.stat[rownames(filtered.FRV),"FRV"]
  filtered.counts[,id2] <- nonsp.output.filtered.clone.stat[rownames(filtered.counts),"counts"]  
  
}



write.table(spikein.final,file=paste0(name,".spikein.output.txt"),quote=F,sep="\t",row.names = F)


output.final <- append(list(spikein.final),nonspikein.final)
names(output.final)[1] <- "spike-in"



library(openxlsx)
write.xlsx(output.final,file=paste0(name,".output.xlsx"))

filtered.counts[which(is.na(filtered.counts) )] <- 0
filtered.FRV[which(is.na(filtered.FRV) )] <- 0
filtered.FRV <- round(filtered.FRV,digits = 3)
#filtered.FRV[is.na(filtered.FRV)] <- 0

cat("pheatmap1...\n")
#pheatmap(allcluster.FRV.index.stat.matrix,display_numbers = T,cluster_rows = F,cluster_cols = F, main="before filtering", width = 6, height = 3)

cat("pheatmap2...\n")
#pheatmap(filtered.FRV, display_numbers = T,cluster_rows = F,cluster_cols = F,main = "after filtering", width = 6, height = 3)
dev.off()
write.table(filtered.FRV, paste(name,"filtered.FRV.txt",sep="."),quote=F,sep="\t",row.names = F) 
write.table(filtered.counts, paste(name,"filtered.counts.txt",sep="."),quote=F,sep="\t",row.names = F) 
}

