####
## Genome Loc Plot
###

## bedfile: bed file of locations of interes
# must be in the form of:
# Chr - Start - End - Gene
#  1    10000   100010 GeneX
#  1    100100  100200 GeneX
# All samples with the same gene name will be merged together


## chrAnnot: annotation file for genome of interest
# must have column called "sumdist" - if all chromosomes are put end to end, what is the base number
# similar deal for CEPlocation 

## function to plot intervals across the genome
#chrAnnot="../annotations/human_chromosome_summary.txt"
#bedfile=data$BedFile

GenomeLocationSeg=function(bedfile, chrAnnot, windsize=5000000){
  
  chrInfo=read.delim(chrAnnot, sep="\t")
  ## Set up the regions
  ExpectedNoWind=round(chrInfo$bp/windsize)
  cIdx2=sapply(1:23, function(x) seq(chrInfo$sumdist[x], chrInfo$cumdist2[x], length=ExpectedNoWind[x]+1))
  cIdx2N=lapply(1:length(cIdx2), function(x) cIdx2[[x]]-chrInfo$sumdist[x])
  cIdx2M=sapply(cIdx2N, function(x) cbind(x[1:(length(x)-1)], x[2:(length(x))]))
  cIdxL=sapply(cIdx2M, nrow)
  temp=do.call(rbind,cIdx2M)
  WindowReg=GRanges(seqnames = rep(chrInfo$Chromosome[1:23], times=cIdxL), ranges=IRanges(start= temp[,1], end=temp[,2]))
  ## Intersect the bed file
  BFile=read.delim(bedfile, header=F)
  t2=unique(BFile$V4)
  BFile2=sapply(t2, function(x) c(min(BFile$V2[which(BFile$V4==x)]), max(BFile$V3[which(BFile$V4==x)])))
  tx=GRanges(seqnames = BFile$V1[match(t2, BFile$V4)], ranges=IRanges(start=BFile2[1, ], end=BFile2[2, ]), Gene=t2)   
  ## Concatenate the gene names
  winannot=GenomicRanges::findOverlaps(tx,WindowReg)
  WindowReg$NoGenes=as.numeric(table(factor(subjectHits(winannot), c(1:length(WindowReg)))))
  Genes1=sapply(1:length(WindowReg), function(x) paste(tx$Gene[queryHits(winannot)[which(subjectHits(winannot)==x)]], sep=", ", collapse=", "))
  WindowReg$Genes=Genes1
  ## make the plot
  ## Output the data
  
  WindowReg$Site=c(1:length(WindowReg))
  WindowReg$Col=ifelse(as.numeric(seqnames(WindowReg))%%2==0, "Even Chr", "Odd Chr")
  
  ret=list(WindowReg=WindowReg, BedFile=tx)
  ret
  
}
 


WindowPlot=function(WindowM, Pvar){
  
  library(ggplot2)
  library(ggrepel)
  library(stringi)
  tx=which(WindowM$Genes=="")
  WindowM$Genes=sapply(WindowM$Genes, function(x) paste(stri_wrap(x, 20), collapse="\n "))
  WindowM$Genes[tx]=""
  WindowM$Yval=eval(WindowM[ ,Pvar])
  matidx=match(c(seq(1:22), "X"),WindowM$Chr)
  
  ggplot(WindowM, aes(x=Site, y=Yval, fill=Col, label=Genes))+geom_vline(xintercept=matidx-0.5, col="grey")+
    geom_bar(stat="identity")+  ylab(Pvar)+
    geom_label_repel(size=2.5, segment.size = 0.2, fill="white", fontface=3)+
    scale_fill_manual(values=c("coral", "dodgerblue", "grey"))+theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


MergeBed=function(singGR, sing2GR, sing3GR, lower=2, upper=100, namList=c("Temous", "CNVkit", "GATK")){

N4x=sing3GR[which(abs(sing3GR$mean)>lower & abs(sing3GR$mean)<upper) ]
N2x=sing2GR[which(abs(sing2GR$mean)>lower & abs(sing2GR$mean)<upper) ]
Nx=singGR[which(abs(singGR$mean)>lower & abs(singGR$mean)<upper) ]

colnames(elementMetadata(N4x))="A3"
colnames(elementMetadata(N2x))="A2"
colnames(elementMetadata(Nx))="A1"

AllMerge=c(Nx, N2x, N4x)
ax1=findOverlaps(AllMerge, BedList)
BedList$A1=NA
lx=which(!is.na(AllMerge$A1[queryHits(ax1)]))
BedList$A1[subjectHits(ax1)[lx]]=AllMerge$A1[queryHits(ax1)[lx]]
BedList$A2=NA
lx=which(!is.na(AllMerge$A2[queryHits(ax1)]))
BedList$A2[subjectHits(ax1)[lx]]=AllMerge$A2[queryHits(ax1)[lx]]
BedList$A3=NA
lx=which(!is.na(AllMerge$A3[queryHits(ax1)]))
BedList$A3[subjectHits(ax1)[lx]]=AllMerge$A3[queryHits(ax1)[lx]]

t2=unique(subjectHits(ax1))
#BedList[t2]
repThis=melt(BedList[t2])
colnames(repThis)=c("Chr", "Start", "End", "width", "strand", "Gene", namList)

 ret=list(report=repThis, GRanges=AllMerge)
 ret
}

#DT::datatable(repThis, rownames=F, class='cell-border stripe',
#              extensions="Buttons", options=list(dom="Bfrtip", buttons=c('csv', 'excel')))


