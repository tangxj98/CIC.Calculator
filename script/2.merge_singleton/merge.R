# For each index1 & index2 combination, merge the barcode with minimal difference (which might be results of sequencing error). 
# Default mindist=2
# sp.index1 is the specific 1.6 that was not used in the experiments. 
mindist=2 
sp.index1="GTCA"
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
  cat("Usage:Rscript merge.R <input file> <index1> <index2> <output prefix>\n")
}else{

library(scales)
# input spike-in

print(args)
file=args[1]
index1=args[2]
index2=args[3]
name=args[4]

outputdir=paste0("./",name)
if(! dir.exists(outputdir) ){
	dir.create(outputdir)
}

# parameter, if the distance small than mindist, merge them

BC <- read.table(args[1],header=T,sep="\t",stringsAsFactors = F)
colnames(BC)[1:4] <- c("Index.1","BarcodeSeq","Index.2","Counts")
#index.1 BC.Barcode.Sequence        index.2 count   Is.Spikein
#iGTCA,AAGTAACAATCGTGATCGAAATGGGTCGAACTT,CTAGA    GTCA    AAGTAACAATCGTGATCGAAATGGGTCGAACTT       CTAGA   3888855 Yes

# remove clusters with only 1 count to reduce running time
dim(BC)
BC <- BC[ BC$Counts >1,]
cat("remove singletons...\n")
dim(BC)
#is.data.frame(BC)
head(BC)
BC.2 <- BC[ nchar(BC$BarcodeSeq)==33, ]

BC.subset <- BC.2[ BC.2$Index.1==index1 & BC.2$Index.2 == index2, ]

dist <- function(bc1, bc2){
	  t1 <- unlist(strsplit(bc1,split = ""))
  t2 <- unlist(strsplit(bc2,split = ""))
    return( length (which(t1!=t2)) )
}

# order decreasingly.
if(index1==sp.index1){
	# if this is spike-in, Barcode needs to be in the front. 
	BC.subset <- BC.subset[ order(BC.subset$Is.Spikein,BC.subset$Counts,decreasing = T),]
}else{
	BC.subset <- BC.subset[ order(BC.subset$Counts,decreasing = T),]
}
cat("BC.subset:")
dim(BC.subset)
head(BC.subset)


# merge
BC.subset$cluster <- rep(0,nrow(BC.subset))
BC.subset$FirstBC <- rep("",nrow(BC.subset))
total=sum(BC.subset$Counts)
n=1
merged.pct=0
while(sum(BC.subset$cluster==0)>0 ){
	  # start from the first row that's zero, i.e., the largest left BC.subset
	  pool=which(BC.subset$cluster==0)
  	  core=pool[1]
          BC.subset[core,"cluster"] <- n
          BC.subset[core,"FirstBC"] <- BC.subset[core,"BarcodeSeq"]
   
    pool=pool[-1]
      #calculate current dist
      cdist=unlist(sapply(BC.subset[pool,"BarcodeSeq"],function(x) dist(BC.subset[core,"BarcodeSeq"] ,x)) )
      # if any BC.subset can be merged in, get their pool id. 
      if(sum(cdist<=mindist)>0 ){
	          BC.subset[ pool[cdist <=mindist] ,"cluster" ] <- n
	          BC.subset[ pool[cdist <=mindist] ,"FirstBC"] <- BC.subset[core,"BarcodeSeq"]
        }
        
        merged.pct = sum(BC.subset[BC.subset$cluster!=0,"Counts"])/total
        cat(n,"\t",sum(BC.subset$cluster==0),"\t",merged.pct,"\n")
	  n=n+1
}
table(BC.subset$cluster)
write.table(BC.subset,file=paste0(outputdir,"/",paste(c(name,index1,index2,"cluster.txt"),collapse=".")),quote=F,sep="\t", row.names=F,col.names=F)
}
