# generate MM plot from output. 
# Only works for this experimental design.
#

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=5){
  
  cat("Usage:Rscript Circos.Barcode.breast.R <group file>  <non-spikein output file> <mouse id> <group: one of the group from the group file> <options: 1.1 to 1.5 or all>\n")
  
}else{
  
  
  idp <- args[5]
  if(! idp %in% c("1.1","1.2","1.3","1.4","1.5","all")){
    cat("Please give either 1.1 ~ 1.5 or all!\n")
    quit()
  }
  groupinfo <- read.table(file = args[1],header = T,sep="\t",stringsAsFactors = F,colClasses = 'character')
  rownames(groupinfo) <- paste(groupinfo$Flowcell,groupinfo$Index.2,sep=".")
  gname=args[4]
  if(! gname %in% groupinfo$Groups){
    cat("Please give group in one of the group name in the group file!\n")
    quit()
  }
  
  library(circlize)
  library(RColorBrewer)
  library(gplots)
  library(scales)
  library(ComplexHeatmap)
  
  
  
  index2.coding <- read.table("index2.coding.txt",header=F,sep="\t",stringsAsFactors = F,colClasses = 'character')
  colnames(index2.coding) <- c("Coding","Index.2")
  rownames(index2.coding) <- index2.coding$Index.2
  
  
  index1.coding <- read.table("index1.coding.txt",header=F,sep="\t",stringsAsFactors = F,colClasses = 'character')
  colnames(index1.coding) <- c("Coding","Index.1")
  rownames(index1.coding) <- index1.coding$Index.1
  
  
  barcodes <- read.table(args[2],header=T,sep="\t",stringsAsFactors = F)
  barcodes$index2.coding <- index2.coding[ barcodes$index.2, "Coding"]
  barcodes$index1.coding <- index1.coding[ barcodes$index.1, "Coding"]
  barcodes$sample <- barcodes$index2.coding
  
  # color table
  colors <- read.table("BreastStudy.color.txt",header = F,sep="\t",stringsAsFactors = F,comment.char = "",row.names = 1)
  
  
  
  idx1.color <- read.table("BreastStudy.index1.color.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
  
 
    
    
    # plot for a specific mouse ### please note that this might need manual check. !!!!!!!
    mouse=args[3]
    plot.info <- groupinfo [ grep(mouse,groupinfo$Annotation),]
    rownames(plot.info) <- paste(plot.info$Flowcell,plot.info$Index.Used,sep=".")
    
    # generate plot data for this plot id. 
    plotid="RC_46"
    plot.data <- barcodes [ which(barcodes$sample %in% plot.info$Index.Used), ]
    plot.data$Organ <- groupinfo[ match(plot.data$sample,groupinfo$Index.Used),"Organ" ]
    plot.data$Organ <- factor(plot.data$Organ, levels = rownames(colors) )
    #
    df <- plot.data[,c("index1.coding","Organ","output.bc.seq","count")]
    # sectors: decided by plot.info$Organ
    df$sectors <- df$Organ
    
    df <- df[order(df$sectors,-df$count),]
    # what's x? using all barcode.
    # df$x <- unlist(sapply(table(df$sectors),function(x) seq(x)))
    df$x <- as.numeric(factor(df$output.bc.seq,levels = unique(barcodes$output.bc.seq)))
    df$y <- log10(df$count)
    freq <- table(df$output.bc.seq)
    maxfreq <- max(freq)
    df$barcodeFreq <- unlist(freq[df$output.bc.seq])
    df$barcodeFreq.norm <- unlist(freq[df$output.bc.seq])/maxfreq
    if(idp!="all"){
      df <- df[ df$index1.coding==idp,]
    }
    if(nrow(df)>0){
      print(paste0(gname,".",plotid,".",idp,".df.txt"))
      write.table(df,file=paste0(gname,".",plotid,".",idp,".df.txt"),quote=F,sep="\t",row.names = F)
    }else{
      print(paste0(gname," does not contain ",plotid,".",idp,"!"))	
      next;
    }	
    
    png(paste0(gname,".",plotid,".",idp,".png"),height = 400, width = 800)
    
    circos.clear()
    circos.par("start.degree"=90)
    
    allorgan= factor(rownames(colors),levels = rownames(colors))
    xlim.mat <- cbind(rep(0,nlevels(allorgan)), rep(max(df$x), nlevels(allorgan)))
    
    circos.initialize(factors = allorgan, xlim=xlim.mat )
    # track for just the name
    circos.track(df$sectors, ylim=c(0,1),track.margin=c(0.02,0.15),track.height = 0.1,bg.col=colors[levels(df$sectors),1],
                 panel.fun = function(x, y) {
                    circos.axis(labels.cex = 0.4)
                   
                 })
    # track for shared barcodes. 
    circos.track(df$sectors, ylim=c(0,1),track.height = 0.1, 
                 panel.fun = function(x, y) {
                   
                 }
    )
    
   
    
    if(nrow(df)>0){
      for(i in 1:nrow(df)){
        
        
          circos.segments(x0 = df[i,"x"],y0=0, x1=df[i,"x"], y1 = 1,straight=T,col=alpha(idx1.color[df[i,"index1.coding"], "color"],alpha = df[i,"barcodeFreq.norm"]/2),sector.index = df[i,"sectors"])
        
      } 
      
      
      # Logic: the barcode was identical for all sites with all possible barcodes (x is the barcode id sorted by characters). the color of the barcode was read and the opacity was decided by the normalized sharing frequency (barcodeFreq/MaxFreq). The links are from     
      
      # link
      # from primary to others. 
      k=0
      organs <- levels(df$Organ)
      for(l in 1:length(organs) ){
        for (i in which(df$Organ== organs[l] )  ) {
          for(j in which( df$x==df[i,"x"] & ! (df$Organ %in% organs[1:l] ) ) ){
            circos.link(df[i,"Organ"],df[i,"x"],df[j,"Organ"],df[j,"x"],col = alpha(idx1.color[df[i,"index1.coding"], "color"],0.25)) 
            k=k+1 
          }
        }
      }  
      lgd <- Legend(at = levels(df$sectors), type = "grid", 
                    legend_gp = gpar(fill = colors[levels(df$sectors),1]), title_position = "topleft", 
                    title = "Organs")
      draw(lgd, just = c("right", "bottom"),x = unit(1.6, "snpc"),y= unit(0.1,"snpc"))
    }
    dev.off()  
  }
  

