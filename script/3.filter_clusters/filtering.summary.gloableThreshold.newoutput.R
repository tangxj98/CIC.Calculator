# Global filtering 

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){
  cat("Usage:Rscript filtering.summary.R <clustered file> <output prefix>\n")
}else{

library(tidyr)
library(pheatmap)
library(investr)
library(scales)
library(openxlsx)
  
file=args[1]
name=args[2]
pdf(paste0(name,".acrossSample.filtering.pdf"))  
# read in spike-in control. 

spikein.contol <- read.table("../data/spikin.control.txt",header=T,sep="\t",stringsAsFactors = F)
spikein.contol$Spike.in.BC.short <- sapply(spikein.contol$Spike.in.Barcode.Sequence, function(x) substr(x,1,33) ) 

index1.coding <- read.table("../data/index1.coding.txt",header=F,sep="\t",stringsAsFactors = F)
rownames(index1.coding) <- index1.coding$V2
index2.coding <- read.table("../data/index2.coding.txt",header=F,sep="\t",stringsAsFactors = F)
rownames(index2.coding) <- index2.coding$V2



# read clustered file 
cluster <- read.table(file,header=F,sep="\t",stringsAsFactors = F)
 head(cluster)
colnames(cluster) <- c("index.1","Barcode.Sequence","index.2","count","Is.Spikein","cluster","FirstBC")
cluster$Is.Spikein2 <- ifelse(cluster$FirstBC %in% spikein.contol$Spike.in.BC.short,"Yes","No")

# Sample QC based on raw Sample Total Counts 
cluster.sample <- aggregate(cluster$count,by=list(cluster$index.2),FUN=sum)
cluster.sample <- cluster.sample[order(cluster.sample$x,decreasing = T),]
colnames(cluster.sample) <- c("index.2","counts")
rownames(cluster.sample) <- cluster.sample$index.2
cluster.sample$coding <- index2.coding[rownames(cluster.sample),1]

# plot the raw total counts 
temp <- barplot(cluster.sample[,2],names.arg = cluster.sample[,3],las=2,ylim=c(0,max(cluster.sample$counts)*1.1))
box()
text(x=temp,y= cluster.sample[,2]+max(cluster.sample$counts)*0.025 ,labels = cluster.sample[,2],cex = 0.8)

if(FALSE){
## do not remove clusters due to low counts. DEC2021
## threshold: median of the sample total counts * 0.1, i.e., at least 10% of the median counts. 
keep <- cluster.sample[cluster.sample$counts >= median(cluster.sample$counts)*0.1,"index.2"]
drop <- cluster.sample[cluster.sample$counts < median(cluster.sample$counts)*0.1,]
if(nrow(drop)>0){
    print("The following samples have been removed due to low total counts: ")
    print(cluster.sample[cluster.sample$counts < median(cluster.sample$counts)*0.1, ])
  
}

cluster <- cluster[ cluster$index.2 %in% keep, ]
index2.coding <- index2.coding[ index2.coding$V2 %in% keep, ]
}

cat("**********************************************\n\n")
cat("Global Threshoding Start processing the spike-in control index GTCA ...\n")

# 1. separate the spikein and non-spikein 
spikein_cluster <- subset(cluster, index.1=="GTCA")
sp <- subset(spikein_cluster,Is.Spikein2=="Yes")
ns <- subset(spikein_cluster,Is.Spikein2=="No")
cat("percentge of Barcode Sequence in the GTCAs: ")
print ( sum(sp$count)/sum(spikein_cluster$count))

# 2. get aggregated counts a) by the BC sequence b) by sample (index 2) for noralization (FRV)
# a) BC seq
sp.output <- aggregate(sp$count, by=list(sp$index.2,sp$FirstBC), FUN=sum)
colnames(sp.output) <- c("index.2","output.bc.seq","count")
rownames(sp.output)  <- paste0(sp.output$index.2,sp.output$output.bc.seq)
sp.output$Input.Copies <- sapply(sp.output$output.bc.seq, function(x) spikein.contol[which(spikein.contol$Spike.in.BC.short==x), "Expected.Copies"])
table( sapply(sp.output$output.bc.seq, function(x) nchar(x) ) )

ns.output <- aggregate(ns$count, by=list(ns$index.2,ns$FirstBC), FUN=sum)
colnames(ns.output) <- c("index.2", "output.bc.seq","count")
rownames(ns.output)  <- paste0(ns.output$index.2,ns.output$output.bc.seq)


# b) by sample (index2)
sp.sample <- aggregate(sp$count, by=list(sp$index.2),FUN=sum)
colnames(sp.sample) <- c("index2","samplesize")
rownames(sp.sample) <- sp.sample$index2
sp.output$samplesize <- sp.sample[sp.output$index.2,"samplesize"] 
sp.output$FRV <- sp.output$count/sp.output$samplesize
sp.output.rm1 <- subset(sp.output,Input.Copies>1)

ns.sample <- aggregate(ns$count, by=list(ns$index.2),FUN=sum)
colnames(ns.sample) <- c("index2","samplesize")
rownames(ns.sample) <- ns.sample$index2
ns.output$samplesize <- ns.sample[ns.output$index.2,"samplesize"] 
ns.output$FRV <- ns.output$count/ns.output$samplesize

# spikein.total = sum(sp.output$count)
# nonspikein.total=sum(ns.output$count)

# 96.9% counts are from spike-in without any merge. 
# spikein.total/(spikein.total + nonspikein.total)
# [1] 0.9689623

reg <- lm(data = sp.output, formula = log10(FRV)~log10(Input.Copies))
boxplot(log10(sp.output$FRV)~as.numeric(log10(sp.output$Input.Copies)),ylim=c(-6,1),boxwex=0.25)
points(log10(sp.output$Input.Copies), log10(sp.output$FRV),xlab="log10(Input Copies)", ylab="log10(FRV)")
range(log10(sp.output[which(sp.output$Input.Copies==1),"FRV" ] ))
plot(log10(sp.output$Input.Copies), log10(sp.output$FRV),xlab="log10(Input Copies)", ylab="log10(FRV)",ylim=c(-6,1),xlim=c(-1,5))
abline(reg, lty=2, lwd=2, col="red")
mtext(paste("R squared=", round(summary(reg)$r.squared, digits=3), sep=" "), adj=1)
mtext(paste("y=", round(summary(reg)$coefficients[2,1], digits=4),"x", round(summary(reg)$coefficients[1,1], digits=4),   sep=" "), adj=0)
reg2 <- lm(data = sp.output.rm1, formula = log10(FRV)~log10(Input.Copies))
abline(reg2, lty=2, lwd=2, col="blue")

mtext(paste("R squared=", round(summary(reg2)$r.squared, digits=3), sep=" "), adj=1,line=1)
mtext(paste("y=", round(summary(reg2)$coefficients[2,1], digits=4),"x", round(summary(reg2)$coefficients[1,1], digits=4),   sep=" "), adj=0,line=1)

#dev.off()

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
ROC$threshold = thresholds
ROC$obj=ROC$TPR-ROC$FPR
best.cutoff <- ROC[ which.max(ROC$obj),"threshold"]
best.cutoff
#[1] 5.863534e-05
sugg.threshold=ROC[which(ROC$threshold==best.cutoff),]
plot(ROC[,1:2],type="b" , xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.4, ylab="True Positive Rate", xlab="False Positive Rate")
abline(coef=c(0,1), lty=2, lwd=0.5)
points(ROC[which(ROC$threshold==best.cutoff),1:2], pch=16, cex=2, col="red")
mtext(paste("Suggested Threshold (FRV)=", signif(sugg.threshold[,"threshold"], digits=4)), line=0, adj=1)
mtext(paste("Suggested TPR=", signif(sugg.threshold[,"TPR"], digits=4)), line=1, adj=0)
mtext(paste("Suggested FPR=", signif(sugg.threshold[,"FPR"], digits=4)), line=0, adj=0)

# evalutate filtering. 
noise_remaining<- table(ns.output$FRV>best.cutoff)[2]
hist(log10(ns.output$FRV), breaks = 100, xlab="log10(FRV)",main=paste0("unwanted barcode after filtering % =",round(noise_remaining*100/nrow(ns.output),digits = 2),"%") ) 
abline(v=log10(best.cutoff),col="red",lty=4,lwd=2)

sp.output.rm1.filtered <- subset(sp.output.rm1,FRV > best.cutoff)
ns.output.filtered <- subset(ns.output, FRV> best.cutoff)



# how many effective sampes in this flow cells?  
# a) prediction to the input cells 

# redo the lm fitting. 
reg3 <- lm(data = sp.output.rm1.filtered, formula = log10(FRV)~log10(Input.Copies))
# predict the results. 
sp.output.rm1.filtered$predicted.input <- round(sapply(log10(sp.output.rm1.filtered$FRV),function(x) invest(reg3,x,upper = 10000000,lower=0,interval = "none") ),digits = 1)

plot(log10(sp.output.rm1.filtered$Input.Copies), log10(sp.output.rm1.filtered$FRV),xlab="log10(Input Copies)", ylab="log10(FRV)",ylim=c(-6,1),xlim=c(-1,5))
abline(reg, lty=2, lwd=2, col="red")
mtext(paste("R squared=", round(summary(reg)$r.squared, digits=3), sep=" "), adj=1)
mtext(paste("y=", round(summary(reg)$coefficients[2,1], digits=4),"x", round(summary(reg)$coefficients[1,1], digits=4),   sep=" "), adj=0)


#sample.input <- aggregate(sp.output.rm1.filtered$predicted.input,by=list(sp.output.rm1.filtered$index.2),FUN=sum)


# blacklist: identify barcodes presented in every non-spikein cluster and will survive after filtering. 
presence <- table(ns.output.filtered$output.bc.seq)
blacklist <- names(presence)[presence==max(presence)]
# blacklist <- names(presence)[presence==20]
if(max(presence)> nrow(index2.coding)*0.5){
  print(paste0("Max frequence of a barcode:",max(presence)))
  blacklist <- names(presence)[presence==max(presence)]
  print("blacklisted barcodes:")
  print(blacklist)
  write.table(ns.output.filtered[which(ns.output.filtered$output.bc.seq %in% blacklist),],file = paste0(name,".spikein_index.blacklist.txt"),quote=F,sep="\t", row.names =F )
}


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

# calculate sample size and normalized FRV 
nonsp.sample <- aggregate(nonsp.output$count, by=list(nonsp.output$index.2),FUN=sum)
colnames(nonsp.sample) <- c("index.2","samplesize")
rownames(nonsp.sample) <- nonsp.sample$index.2
nonsp.output$FRV <- nonsp.output$count/nonsp.sample[nonsp.output$index.2,"samplesize"]

# filtering: remove small outputBC, keep only those FRV > threshold. 
nonsp.output.filtered <- subset(nonsp.output, FRV>best.cutoff)
nonsp.output.filtered$sample <- paste(nonsp.output.filtered$index.1,nonsp.output.filtered$index.2,sep=".")
nonsp.output.filtered.sample.stat <- aggregate(nonsp.output.filtered$FRV,by=list(nonsp.output.filtered$index.1,nonsp.output.filtered$index.2),FUN=sum)
colnames(nonsp.output.filtered.sample.stat) <- c("Index.1","Index.2","FRV")


# predict input cells
nonsp.output.filtered$predicted.input <- round(sapply(log10(nonsp.output.filtered$FRV),function(x) invest(reg3,x,upper = 10000000,lower=0,interval = "none") ),digits = 1)

write.table(nonsp.output.filtered,paste0(name,".nonsp.filtered.output.txt"),quote=F,sep="\t")
# list if the blacklist sequence presented.
 temp <- nonsp.output.filtered[which(nonsp.output.filtered$output.bc.seq==blacklist),] 
if(nrow(temp)>0){
  print("The spike-in blacklist in each sample:")
  print(temp)
}

# output format: list of data.frame
temp <- sp.output.rm1.filtered 
temp <- data.frame(index.2.coding=index2.coding[sp.output.rm1.filtered$index.2,"V1"],temp,stringsAsFactors = F) 
temp <- temp[order(temp$index.2.coding,-temp$Input.Copies),]
output.final <- list(temp)


temp <- data.frame(index1.coding=index1.coding[nonsp.output.filtered$index.1,"V1"] ,index.2.coding=index2.coding[nonsp.output.filtered$index.2,"V1"],nonsp.output.filtered,stringsAsFactors = F)

for(i in sort(unique(temp$index.2.coding))){
  temp.id2 <- temp[temp$index.2.coding==i,]
  temp.id2 <- temp.id2[order(-temp.id2$count),]
  output.final <- append(x = output.final,values = list(temp.id2))
}
names(output.final) <- c("spike-in", as.character(sort(unique(temp$index.2.coding))) )

write.xlsx(output.final,file=paste0(name,".output.xlsx"),overwrite=T)



sample.index.stat <- spread(nonsp.output.filtered.sample.stat,key = "Index.2", value = FRV)
sample.index.stat[is.na(sample.index.stat)] <- 0
rownames(sample.index.stat) <- sample.index.stat[,1]
sample.index.stat <- sample.index.stat[,-1]
sample.index.stat <- round(sample.index.stat,digits = 4)


#
pheatmap(sample.index.stat,display_numbers = T,cluster_rows = F,cluster_cols = F)
write.table(sample.index.stat, paste0(name,".index.stat.matrix.txt"),quote=F,sep="\t") 
dev.off()

}
