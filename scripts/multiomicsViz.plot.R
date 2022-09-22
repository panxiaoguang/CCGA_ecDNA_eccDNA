# Figure4A-----plot -------------------------------------------------------

EPMs<-readRDS("~/Project/Bladder/gene_drops.RDS")
mat_ecc<-EPMs%>%
  select(1:81)%>%
  setNames(stringr::str_remove(stringr::str_remove(names(EPMs)[1:81],"cBca_"),"T"))%>%
  tibble::column_to_rownames(var="gene")%>%
  as.data.frame()

CNVs<-read_tsv("~/Project/WGS/CNV/annotate/all_CNV.tsv")
mat_cnv<-CNVs%>%
  setNames(stringr::str_remove(stringr::str_remove(names(CNVs),"Bca_"),"T"))%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
#####process chrome location###

chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                        "12","13","14","15","16","17","18","19","20",
                        "21","22","X","Y")

chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                        "12","13","14","15","16","17","18","19","20",
                        "21","22","X","Y")

#######Extract sub list#########
genelocate<-read_tsv("~/Project/Bladder/dbs/hg38.coding.bed",col_names = F)
genelocate<-genelocate%>%
  select(4,1,2,3)%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("Symbol","chrom","start","end"))%>%
  as.data.frame()

genelocate_sourceOmics <- genelocate[genelocate[,2] %in%
                                       chrome_sourceOmics,]
genelocate_targetOmics <- genelocate[genelocate[,2] %in%
                                       chrome_targetOmics,]

intG <- intersect(rownames(targetOmics),genelocate_targetOmics[,1])

targetOmics <- targetOmics[intG,]

source_gene <- rownames(sourceOmics)
source_gene_locate <-intersect(unique(genelocate_sourceOmics[,1]),source_gene)
source_gene <- sourceOmics[source_gene_locate,]
genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[,1] %in% source_gene_locate,]
genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[,1] %in% intG,]
###Calculate the correlation between cna and other omics data######
corrArray <-calculateCorForTwoMatrices(source_gene,targetOmics,0.01)
## functions#############
calculateChromLength <- function(chromLength,selectedChrom,genelocate){
  chromLength <- chromLength[chromLength[,1] %in% selectedChrom,,drop=FALSE]
  
  if(length(selectedChrom)==1){
    x <- 0
  }else{
    x <- c(0,chromLength[1:(nrow(chromLength)-1),2])
  }
  chromLength[,3] <- cumsum(as.numeric(x))
  chromLength[,4] <- cumsum(as.numeric(chromLength[,2]))
  
  genelocate <- cbind(genelocate,0,0)
  
  colnames(genelocate)[5:6] <- c("finalstart","finalend")
  
  for(i in c(1:nrow(genelocate))){
    chr <- genelocate[i,2]
    s <- genelocate[i,3]
    e <- genelocate[i,4]
    cs <- chromLength[chromLength[,1]==chr,3]
    genelocate[i,5] <- s+cs
    genelocate[i,6] <- e+cs
  }
  re <- list(chromLength=chromLength,genelocate=genelocate)
  return(re)
}
plotHeatMap <- function(corrArray,genelocate_sourceOmics,
                        chromLength_sourceOmics,genelocate_targetOmics,chromLength_targetOmics,
                        sourceOmicsName,targetOmicsName,dim=1){
  
  allChromlen_sourceOmics <-
    chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]
  allChromlen_targetOmics <-
    chromLength_targetOmics[nrow(chromLength_targetOmics),4]
  
  if(dim==1){
    par(mar=c(4,4,4,0))
  }else{
    par(mar=c(0,4,4,0))
  }
  
  p <- which(corrArray!=0,arr.ind=TRUE)
  allcnagene <- rownames(corrArray)
  allovgene <- colnames(corrArray)
  
  la <- 1
  for(i in c(1:nrow(p))){
    
    cnag <- allcnagene[p[i,1]]
    ovg <- allovgene[p[i,2]]
    cnagp <- genelocate_sourceOmics[genelocate_sourceOmics[,1]==cnag,5]
    ovgp <- genelocate_targetOmics[genelocate_targetOmics[,1]==ovg,5]
    
    if(length(cnagp)==0 || length(ovgp)==0){
      next
    }
    
    cov <- corrArray[cnag,ovg]
    color <- ifelse(cov>0,"#E63126","#0932E3")
    
    if(la==1){
      if(dim==1){
        plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", targetOmicsName,
                                   " correlation",sep=""), xlim=c(0,allChromlen_sourceOmics),
             ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
             frame.plot=FALSE,xlab=paste(sourceOmicsName,
                                         " chromosomal location",sep=""),ylab=paste(targetOmicsName,
                                                                                    " chromosomal location",sep=""),pch=20,col=color,cex=0.2)
        axis(side=1,at=(chromLength_sourceOmics[,4]-
                          chromLength_sourceOmics[,2]/2),
             labels=chromLength_sourceOmics[,1])
      }else{
        plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", 
                                   targetOmicsName," correlation",sep=""),
             xlim=c(0,allChromlen_sourceOmics),
             ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
             frame.plot=FALSE,ylab=paste(targetOmicsName," chromosomal
location",sep=""),
             xlab="",pch=20,col=color,cex=0.2)
      }
      axis(side=2,at=(chromLength_targetOmics[,4]-
                        chromLength_targetOmics[,2]/2),
           labels=chromLength_targetOmics[,1])
      
      #abline(h=c(0,chromLength_targetOmics[,4]),v=c(0,chromLength_sourceOmics[,4]),
      #       col="gray",lty=3)
      la <- la+1
    }else{
      for(u in seq_len(length(cnagp))){
        for(v in seq_len(length(ovgp))){
          points(cnagp[u],ovgp[v],pch=20,col=color,cex=0.2)
        }
      }
    }
  }
}

#####################
##Calculate the location of genes in the heatmap
chromLength <- read_tsv("../../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
chromLength<-chromLength%>%
  filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("V1","V2"))
chromLength$V1<-factor(chromLength$V1,levels = c(seq(1,22),"X","Y"))
chromLength<-chromLength%>%
  arrange(V1)%>%
  as.data.frame()

re <-calculateChromLength(chromLength,chrome_sourceOmics,genelocate_sourceOmics)
genelocate_sourceOmics <- re$genelocate
chromLength_sourceOmics <- re$chromLength
re <-calculateChromLength(chromLength,chrome_targetOmics,genelocate_targetOmics)
genelocate_targetOmics <- re$genelocate
chromLength_targetOmics <- re$chromLength
###plot
plotHeatMap(corrArray,genelocate_sourceOmics,chromLength_sourceOmics,
            genelocate_targetOmics,chromLength_targetOmics,"eccDNA",
            "mRNA",dim=1)
