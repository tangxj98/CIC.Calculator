# 
# This script is to merge the barcode counts from two flow cells that belong to the same sample.
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
 cat("Usage: Rscript merge.input.R <flowcell output1> <flowcell output2> <output prefix>\n")

}else{
	tmp <- list()
	for(i in 1:(length(args)-1) ){
		temp <- read.table(args[i],header=T,sep="\t",stringsAsFactors=F)
		tmp[[i]] <- temp
		print(args[i])
		print(head(temp))
	}
	
	output <- tmp[[1]]
	for(i in 2:(length(args)-1)){
		overlap <- intersect(rownames(tmp[[i]]),rownames(output))
		additional <- setdiff(rownames(tmp[[i]]),rownames(output))
		print(range( tmp[[i]][additional,"count"]))
		output[overlap,"count"] <- output[overlap,"count"] + tmp[[i]][overlap,"count"]
		output <- rbind(output, tmp[[i]][additional,]) 
	}
	write.table(output, file=paste0(args[3],".combined.mergeinput.txt"),quote=F,sep="\t")
}
