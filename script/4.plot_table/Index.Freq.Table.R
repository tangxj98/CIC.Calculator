# calculate the CIC frequency
# What did that guy said. 
#setwd("/research/bsi/projects/PI/tertiary/Kannan_Nagarajan_m161624/s214639.CIC_calculator/processing/R/new_testrun/")


groupinfo <- read.table(file = "Breast_newGroupInfo.txt",header = T,sep="\t",stringsAsFactors = F,colClasses = 'character')
rownames(groupinfo) <- groupinfo$Index.2



index2.coding <- read.table("index2.coding.txt",header=F,sep="\t",stringsAsFactors = F,colClasses = 'character')
colnames(index2.coding) <- c("Coding","Index.2")
rownames(index2.coding) <- index2.coding$Index.2


index1.coding <- read.table("index1.coding.txt",header=F,sep="\t",stringsAsFactors = F,colClasses = 'character')
colnames(index1.coding) <- c("Coding","Index.1")
rownames(index1.coding) <- index1.coding$Index.1


barcodes <- read.table("FMT_2.FCHNFVCBCX3.nonsp.filtered.output.txt",header=T,sep="\t",stringsAsFactors = F)
barcodes$index2.coding <- index2.coding[ barcodes$index.2, "Coding"]
barcodes$index1.coding <- index1.coding[ barcodes$index.1, "Coding"]

barcodes$Organ <- groupinfo[ barcodes$index2.coding, "Organ"]
barcodes$Annotation <- groupinfo[ barcodes$index2.coding, "Annotation"]
# remove 2.01 which is not a mouse sample
barcodes <- barcodes[ !is.na(barcodes$Annotation), ] 

barcodes$Mouse <- paste0("RC_",groupinfo[ barcodes$index2.coding,"RC"] )

save(barcodes,file="Breast_FMT2.barcodes.18mice.RData")
# sum of output cells 
# output <- aggregate(barcodes$count,by=list(barcodes$Annotation,barcodes$index1.coding),FUN=sum)
# output.1.2 <- subset(output,output$Group.2=="1.2")
# output.1.4 <- subset(output,output$Group.2=="1.4")

 
# predicted.input.1.2 <- subset(predicted.input,predicted.input$Group.2=="1.2")
# predicted.input.1.4 <- subset(predicted.input,predicted.input$Group.2=="1.4")
# 
# samplelist <- sort(unique(sapply(barcodes$Annotation, function(x) substr(x,nchar(x)-2, nchar(x)) ) ) )
# matrix.1.4 <- list()
# matrix.1.2 <- list()
# for(i in samplelist) {
#   temp <- barcodes[grepl(i,barcodes$Annotation) & barcodes$index1.coding=="1.4",c(6,10,12) ]
#   temp.2 <- aggregate(temp$predicted.input, by=list(temp$output.bc.seq,temp$organ),FUN=sum)
#   colnames(temp.2) <- c("Barcode","Organ","Predicted.Input")
#   temp.3 <- dcast(temp.2, Barcode~Organ)
#   temp.3[is.na(temp.3)] <- 0
#   temp.3 <- data.frame(index.1="1.4",temp.3)
#   matrix.1.4 <- append(x=matrix.1.4,values = list(temp.3) ) 
#                                                   
#   temp <- barcodes[grepl(i,barcodes$Annotation) & barcodes$index1.coding=="1.2",c(6,10,12) ]
#   if(nrow(temp)>0){
#     temp.2 <- aggregate(temp$predicted.input, by=list(temp$output.bc.seq,temp$organ),FUN=sum)
#     colnames(temp.2) <- c("Barcode","Organ","Predicted.Input")
#     temp.4 <- dcast(temp.2, Barcode~Organ)
#     temp.4[is.na(temp.4)] <- 0
#     temp.4 <- data.frame(index.1="1.2",temp.4)
#   } else{
#     temp.4 <- data.frame()
#     colnames(temp.4) <- colnames(matrix.1.2)
#   }
#   matrix.1.2 <- append(x=matrix.1.2,values = list(temp.4))
# }
# names(matrix.1.4) <- paste0("Mouse_RC",samplelist)
# names(matrix.1.2) <- paste0("Mouse_RC",samplelist)
# write.xlsx(matrix.1.4,file="Table1.1.4.xlsx",sep="\t")
# write.xlsx(matrix.1.2,file="Table1.1.2.xlsx",sep="\t")

library(openxlsx)
DNA.yield <- read.xlsx("/research/bsi/projects/PI/tertiary/Kannan_Nagarajan_m161624/s214639.CIC_calculator/processing/batch_202112/Groups and Annotations_Breast Study.update.xlsx",sheet=1)

# keep only information related to the FMT_2 group. 
DNA.yield <- DNA.yield[ DNA.yield$`RC#` %in% groupinfo$RC, ]

#  DNA.yield related predicted input. Summary based on both index1.coding and index2.coding 
predicted.input <- aggregate(barcodes$predicted.input,by=list(barcodes$Annotation,barcodes$index1.coding,barcodes$index2.coding),FUN=sum)
colnames(predicted.input) <- c("Annotation","Index.1","Index.2","Predicted")
#DNA.yield 
# CIC frequency 
# CIC frequency = iutput cells / output cells 
# Input cells to DNA input
# first identify sample.ID

# summarize the sample specif uniq barcode counts. 
sample.unique.barcode <- data.matrix(table(barcodes[,c("Annotation","index1.coding")]))

# make sure that no missing samples.
setdiff(groupinfo$Annotation, rownames(sample.unique.barcode))
setdiff( rownames(sample.unique.barcode),groupinfo$Annotation)

setdiff(sort(unique(DNA.yield$Annotation)), sort(rownames((sample.unique.barcode))))
setdiff( sort(rownames((sample.unique.barcode))),sort(unique(DNA.yield$Annotation)))


DNA.yield$UniqueBarcodeClones.1.1 <- sample.unique.barcode[DNA.yield$Annotation,"1.1"]
DNA.yield$UniqueBarcodeClones.1.2 <- sample.unique.barcode[DNA.yield$Annotation,"1.2"]
DNA.yield$UniqueBarcodeClones.1.3 <- sample.unique.barcode[DNA.yield$Annotation,"1.3"]
DNA.yield$UniqueBarcodeClones.1.4 <- sample.unique.barcode[DNA.yield$Annotation,"1.4"]
DNA.yield$UniqueBarcodeClones.1.5 <- sample.unique.barcode[DNA.yield$Annotation,"1.5"]

DNA.yield$UniqueBarcodeClones.total <- DNA.yield$UniqueBarcodeClones.1.1 + DNA.yield$UniqueBarcodeClones.1.2 + DNA.yield$UniqueBarcodeClones.1.3 + DNA.yield$UniqueBarcodeClones.1.4 + DNA.yield$UniqueBarcodeClones.1.5

#DNA.yield$UniqueBarcodeClones.tumortotal <- DNA.yield$UniqueBarcodeClones.total * DNA.yield$`DNA.yield.(ug)/organ`

DNA.yield$CIC.Frequency.Uniq.1.1 <- DNA.yield$UniqueBarcodeClones.1.1* DNA.yield$`DNA.yield.(ug)/organ` / DNA.yield$Index.1.1
DNA.yield$CIC.Frequency.Uniq.1.2 <- DNA.yield$UniqueBarcodeClones.1.2* DNA.yield$`DNA.yield.(ug)/organ` / DNA.yield$Index.1.2
DNA.yield$CIC.Frequency.Uniq.1.3 <- DNA.yield$UniqueBarcodeClones.1.3* DNA.yield$`DNA.yield.(ug)/organ` / DNA.yield$Index.1.3
DNA.yield$CIC.Frequency.Uniq.1.4 <- DNA.yield$UniqueBarcodeClones.1.4* DNA.yield$`DNA.yield.(ug)/organ` / DNA.yield$Index.1.4
DNA.yield$CIC.Frequency.Uniq.1.5 <- DNA.yield$UniqueBarcodeClones.1.5* DNA.yield$`DNA.yield.(ug)/organ` / DNA.yield$Index.1.5
DNA.yield$CIC.Frequency.Uniq.total <- DNA.yield$UniqueBarcodeClones.total* DNA.yield$`DNA.yield.(ug)/organ` / DNA.yield$Total.cells.injected



write.xlsx(DNA.yield,file = "BreastStudy.FMT_2.CIC.Frequency.xlsx")



# #rownames(sample.total.barcode) <- sample.total.barcode$Annotation
# DNA.yield$Total.Barocde 


#namely a) CIC of indeterminate potential (CICIP, present only in primary site), b) CIC of mono-metastatic potential (CICMono, observed only in one site +/- primary), c) CIC of pluri-metastatic potential (CICPluri, observed in more than one but not in all sites), and d) CIC of toti-metastatic potential (CICToti, present in 'all' sites).



bar_sum <- function(x){
  #x is the array contain the barcodes organ.
  paste(sort(unique(x)),collapse = ";")
}


ibar_stat <- aggregate(barcodes$Organ,by=list(barcodes$index1.coding,barcodes$Mouse,barcodes$output.bc.seq),FUN=bar_sum)
colnames(ibar_stat) <- c("Index1.coding","Mouse","Barcode","Organs")
library(stringr)
# number of organs presented
ibar_stat$FreqSites <- sapply(ibar_stat$Organs,str_count,";")+1

# group
ibar_stat$BarcodeType=""

ibar_stat$BarcodeType <- ifelse( ibar_stat$Organs=="Primary Tumor",  "CIC.IP", ifelse( ibar_stat$FreqSites==1 | (ibar_stat$FreqSites==2 & grepl("Primary",ibar_stat$Organs) ) ,"CIC.Mono" , ifelse(ibar_stat$FreqSites== length(unique(groupinfo$Organ)) ,"CIC.Toti" ,"CIC.Pluri" )   ) ) 

# length(unique(groupinfo$Organ))
#[1] 10
#  unique(groupinfo$Organ)
#

# rownames(ibar_stat) <- paste(ibar_stat$Barcode,ibar_stat$Mouse,sep=".")
# barcodes$BarcodeType <- ibar_stat[ paste(barcodes$output.bc.seq,barcodes$Mouse,sep="."),"BarcodeType" ]
# barcodeType_summary <- data.frame(table(ibar_stat$Mouse,ibar_stat$BarcodeType))
# colnames(barcodeType_summary) <- c("Mouse","BarcodeType","Frequncy")

# mouse CIC categorize summary for each index1  
CIC.cat.summary <- table(ibar_stat$Mouse,ibar_stat$BarcodeType,ibar_stat$Index1.coding)

output <- data.frame(rbind(as.data.frame.matrix(CIC.cat.summary[,,1]),as.data.frame.matrix(CIC.cat.summary[,,2]),as.data.frame.matrix(CIC.cat.summary[,,3]),as.data.frame.matrix(CIC.cat.summary[,,4]),as.data.frame.matrix(CIC.cat.summary[,,5]) ) )

output <- data.frame(ID=rep(rownames(output)[1:dim(CIC.cat.summary)[1]],5), output)
write.table(output,file = "FMT_2.CIC.category.summary.txt",quote=F,sep="\t",row.names = F)


