# calculate the CIC frequency
# merged filtered barcode files from all mice and generate CIC freq. 

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  cat("Usage:Rscript filtering.summary.R <folder that holds the filtered results> <filtred non-spikein count file pattern> <group info file> \n")
}else{


path=args[1]
pattern=args[2]
groupfile=args[3]
  
#groupinfo <- read.table(file = "/research/bsi/projects/PI/tertiary/Kannan_Nagarajan_m161624/s214639.CIC_calculator/processing/batch_Feb2023/BatchFeb2023.barcode.manifest.PE.txt",header = T,sep="\t",stringsAsFactors = F,colClasses = 'character')
groupinfo <- read.table(file = groupfile, header=T, sep="\t",stringsAsFactors = F,colClasses = 'character') 

# index1 and index2 coding. 
index2.coding <- read.table("index2.coding.txt",header=F,sep="\t",stringsAsFactors = F,colClasses = 'character')
colnames(index2.coding) <- c("Coding","Index.2")
rownames(index2.coding) <- index2.coding$Index.2

index1.coding <- read.table("index1.coding.txt",header=F,sep="\t",stringsAsFactors = F,colClasses = 'character')
colnames(index1.coding) <- c("Coding","Index.1")
rownames(index1.coding) <- index1.coding$Index.1


for(i in 1:nrow(groupinfo)){
  if(groupinfo[i,"Injected.Cell.Number"]==""){
    groupinfo[i,"Injected.Cell.Number"] <-groupinfo[i-1,"Injected.Cell.Number"]
  }
}
groupinfo$Organ <- gsub("PT","Primary Tumor",groupinfo$Organ)
if(length(grep("^Index",colnames(groupinfo)))!=1){
  print("there should only have one Index column available. please check the input in the group info!")
  quit()
}
rownames(groupinfo) <- groupinfo[,grep("^Index",colnames(groupinfo)) ]

#sampleIds <- paste0("BC_Seq", gsub("^0","", gsub("2.","",unique(groupinfo[,grep("^Index",colnames(groupinfo)) ]) )),".nonsp.filtered.output.txt")

#fnames <- list.files(path=path, pattern = "nonsp.filtered.output.txt") 
fnames <- list.files(path=path, pattern =pattern) 

# for global test, read in all filtered non-sp read count files. 
barcodes=data.frame()
for (i in fnames){
  temp <- read.table(paste(path,i,sep="/"),header=T,sep="\t",stringsAsFactors = F,colClasses = 'character')
  #temp$Sample <- rep(gsub(".nonsp.filtered.output.txt","",i),nrow(temp))
  temp$Sample <- rep(gsub(pattern,"",i),nrow(temp))
  barcodes=rbind(barcodes,temp)  
}


barcodes$index2.coding <- index2.coding[ barcodes$index.2, "Coding"]
barcodes$index1.coding <- index1.coding[ barcodes$index.1, "Coding"]

barcodes$Organ <- groupinfo[ barcodes$index2.coding, "Organ"]
if(length(grep("Annotation",colnames(groupinfo)))>0 ){
  barcodes$Annotation <- groupinfo[ barcodes$index2.coding, "Annotation"]
}else{
  barcodes$Annotation <- paste(barcodes$Sample,barcodes$Organ,sep="_")
  groupinfo$Annotation <- barcodes[match(groupinfo$Index,barcodes$index2.coding),"Annotation"]
}


# remove 2.01 which is not a mouse sample
barcodes <- barcodes[ !is.na(barcodes$Annotation), ] 

barcodes$Mouse <- groupinfo[ barcodes$index2.coding,"Mouse"] 

save(barcodes,file=paste(path,"Mouse.barcode.RData",sep = "/"))

#library(openxlsx)

# DNA.yield is in the group info this time...
# DNA.yield <- DNA.yield[ DNA.yield$`RC#` %in% groupinfo$RC, ]
DNA.yield <- groupinfo


#  DNA.yield related predicted input. Summary based on both index1.coding and index2.coding 
predicted.input <- aggregate(as.numeric(barcodes$predicted.input),by=list(barcodes$Annotation,barcodes$index1.coding,barcodes$index2.coding),FUN=sum)
colnames(predicted.input) <- c("Annotation","Index.1","Index.2","Predicted")
#DNA.yield 
# CIC frequency 
# CIC frequency = iutput cells / output cells 
# Input cells to DNA input
# first identify sample.ID

# summarize the sample specif uniq barcode counts. 
sample.unique.barcode <- data.matrix(unclass(table(barcodes[,c("Annotation","index1.coding")])))

# make sure that no missing samples.
# setdiff(groupinfo$Annotation, rownames(sample.unique.barcode))
# setdiff( rownames(sample.unique.barcode),groupinfo$Annotation)

setdiff(sort(unique(DNA.yield$Annotation)), sort(rownames((sample.unique.barcode))))
setdiff( sort(rownames((sample.unique.barcode))),sort(unique(DNA.yield$Annotation)))


DNA.yield$UniqueBarcodeClones.1.1 <- sample.unique.barcode[match(DNA.yield$Annotation,rownames(sample.unique.barcode)),"1.1"]
DNA.yield$UniqueBarcodeClones.1.2 <- sample.unique.barcode[match(DNA.yield$Annotation,rownames(sample.unique.barcode)),"1.2"]
DNA.yield$UniqueBarcodeClones.1.3 <- sample.unique.barcode[match(DNA.yield$Annotation,rownames(sample.unique.barcode)),"1.3"]
DNA.yield$UniqueBarcodeClones.1.4 <- sample.unique.barcode[match(DNA.yield$Annotation,rownames(sample.unique.barcode)),"1.4"]
DNA.yield$UniqueBarcodeClones.1.5 <- sample.unique.barcode[match(DNA.yield$Annotation,rownames(sample.unique.barcode)),"1.5"]

DNA.yield$UniqueBarcodeClones.total <- DNA.yield$UniqueBarcodeClones.1.1 + DNA.yield$UniqueBarcodeClones.1.2 + DNA.yield$UniqueBarcodeClones.1.3 + DNA.yield$UniqueBarcodeClones.1.4 + DNA.yield$UniqueBarcodeClones.1.5

#DNA.yield$UniqueBarcodeClones.tumortotal <- DNA.yield$UniqueBarcodeClones.total * DNA.yield$`DNA.yield.(ug)/organ`

DNA.yield$CIC.Frequency.Uniq.1.1 <- DNA.yield$UniqueBarcodeClones.1.1* DNA.yield$DNA.Yield..µg..Organ.mg. *5/ DNA.yield$Injected.Cell.Number
DNA.yield$CIC.Frequency.Uniq.1.2 <- DNA.yield$UniqueBarcodeClones.1.2* DNA.yield$DNA.Yield..µg..Organ.mg. *5/ DNA.yield$Injected.Cell.Number
DNA.yield$CIC.Frequency.Uniq.1.3 <- DNA.yield$UniqueBarcodeClones.1.3* DNA.yield$DNA.Yield..µg..Organ.mg. *5/ DNA.yield$Injected.Cell.Number
DNA.yield$CIC.Frequency.Uniq.1.4 <- DNA.yield$UniqueBarcodeClones.1.4* DNA.yield$DNA.Yield..µg..Organ.mg. *5/ DNA.yield$Injected.Cell.Number
DNA.yield$CIC.Frequency.Uniq.1.5 <- DNA.yield$UniqueBarcodeClones.1.5* DNA.yield$DNA.Yield..µg..Organ.mg. *5/ DNA.yield$Injected.Cell.Number

DNA.yield$CIC.Frequency.Uniq.total <- DNA.yield$UniqueBarcodeClones.total* DNA.yield$DNA.Yield..µg..Organ.mg. / DNA.yield$Injected.Cell.Number



write.xlsx(DNA.yield,file = "/research/bsi/projects/PI/tertiary/Kannan_Nagarajan_m161624/s214639.CIC_calculator/processing/batch_Feb2023/BatchFeb2023.CIC.Frequency.xlsx")



# #rownames(sample.total.barcode) <- sample.total.barcode$Annotation
# DNA.yield$Total.Barocde 


#namely a) CIC of indeterminate potential (CICIP, present only in primary site), b) CIC of mono-metastatic potential (CICMono, observed only in one site +/- primary), c) CIC of pluri-metastatic potential (CICPluri, observed in more than one but not in all sites), and d) CIC of toti-metastatic potential (CICToti, present in 'all' sites).



bar_sum <- function(x){
  #x is the array contain the barcodes organ.
  paste(sort(unique(x)),collapse = ";")
}


ibar_stat <- aggregate(barcodes$Organ,by=list(barcodes$index1.coding,barcodes$Mouse,barcodes$output.bc.seq),FUN=bar_sum)
colnames(ibar_stat) <- c("Index1.coding","Mouse","Barcode","Organs")
ibar_stat$Organs <- gsub("PT","Primary Tumor",ibar_stat$Organs)
library(stringr)
# number of organs presented
ibar_stat$FreqSites <- sapply(ibar_stat$Organs,str_count,";")+1

# group
ibar_stat$BarcodeType=""

ibar_stat$BarcodeType <- ifelse( ibar_stat$Organs=="Primary Tumor",  "CIC.IP", ifelse( ibar_stat$FreqSites==1 | (ibar_stat$FreqSites==2 & grepl("Primary",ibar_stat$Organs) ) ,"CIC.Mono" , ifelse(ibar_stat$FreqSites== length(unique(groupinfo$Organ)) ,"CIC.Toti" ,"CIC.Pluri" )   ) ) 


# mouse CIC categorize summary for each index1  
CIC.cat.summary <- table(ibar_stat$Mouse,ibar_stat$BarcodeType,ibar_stat$Index1.coding)

output <- data.frame(rbind(data.frame(Index1="1.1",as.data.frame.matrix(CIC.cat.summary[,,1])),data.frame(Index1="1.2",as.data.frame.matrix(CIC.cat.summary[,,2])),data.frame(Index1="1.3",as.data.frame.matrix(CIC.cat.summary[,,3])),data.frame(Index1="1.4",as.data.frame.matrix(CIC.cat.summary[,,4])),data.frame(Index1="1.5",as.data.frame.matrix(CIC.cat.summary[,,5])) ) )

output <- data.frame(ID=rep(rownames(output)[1:dim(CIC.cat.summary)[1]],5), output)
write.table(output,file = "BatchFeb2023.CIC.category.summary.updated.txt",quote=F,sep="\t",row.names = F)



}

