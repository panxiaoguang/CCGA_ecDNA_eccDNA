# Figure7F -------------plot---------------------------------------------------


##########################diffGene_with_eccDNA(higher and lower)#####################################################
counts<-read_tsv("RNAseq/allSample.count.txt")
counts<-counts%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
counts<-counts[,1:70]
groupInfo<-read_tsv("~/Project/Bladder/EPM and clinical variables.tsv")
groupInfo<-groupInfo%>%
  select(1,4)
colData<-tibble(Sample = colnames(counts))
colData<-colData%>%left_join(groupInfo,by="Sample")
colData<-colData%>%
  mutate(group=if_else(`EPM Group`=="Low","A","B"))%>%
  select(Sample,group)%>%
  as.data.frame()
colData$group<-factor(colData$group,levels = c("A","B"))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group) 
dds <- dds[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
res <- results(dds)
allresult<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="Geneid")%>%
  as_tibble()

allresult<-allresult%>%
  mutate(logp=-log10(padj))

allresult<-allresult%>%
  mutate(col=case_when((log2FoldChange<(-1))&(logp>2)~"downregular",
                       (log2FoldChange>1)&(logp>2)~"upregular",
                       TRUE ~ "none"))

geneNames = bitr(allresult$Geneid,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
names(geneNames)<-c("Geneid","id")
allresult<-allresult%>%inner_join(geneNames,by="Geneid")
geneList = allresult$log2FoldChange
names(geneList) = allresult$id
geneList = sort(geneList, decreasing = TRUE)

##download GSEA db
sub1<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:BIOCARTA")
sub2<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
sub3<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:PID")
sub4<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
C2_t2g<-rbind(sub1,sub2,sub3,sub4)%>%
  dplyr::select(gs_name, entrez_gene)
C2_t2n<-rbind(sub1,sub2,sub3,sub4)%>%
  dplyr::select(gs_name,gs_description)

em2 <- GSEA(geneList, TERM2GENE = C2_t2g,TERM2NAME=C2_t2n,minGSSize = 20,eps = 0)

rankvalues<-tibble(geneid=names(rev(geneList)),foldchange=rev(geneList),order=1:18680)
needID = c(6,9,11,20,22,26,30,33,34,38,43,58)
getPlotData<-function(i){
  tibble(geneid=em2@geneSets[[em2@result[needID[i],]$ID]])%>%
    left_join(rankvalues,by="geneid")%>%
    mutate(y1=i-1,y2=i-1+0.9,col=if_else(foldchange>0,"yes","no"))%>%
    tidyr::drop_na()
}


haha<-do.call("rbind",lapply(1:12, getPlotData)) 
ggplot(haha,aes(x=order,y=y1))+
  geom_segment(aes(xend=order,yend=y2,color=col))+
  scale_color_manual(values = c("no"="#5B799D","yes"="#CC726A"))+
  theme_void()

ggsave("GSEA_result.part3.png",width = 7.04,height = 2.37,dpi=300)

rankvalues<-rankvalues%>%
  mutate(col=if_else(foldchange>0,"yes","no"))

ggplot(rankvalues,aes(x=order,y=foldchange))+
  geom_area(aes(fill=col))+
  scale_fill_manual(values = c("no"="#5B799D","yes"="#CC726A"))+
  theme_pubr()
ggsave("GSEA_result.part2.pdf",width =7.04 ,height = 4.14)


ggplot(rankvalues,aes(x=foldchange))+
  geom_histogram(fill="#5B799D",color="black",bins = 100,size=0.2)+
  scale_x_continuous(limits = c(-2,2))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene Density")+
  xlab("Correlation")
ggsave("GSEA_result.part1.pdf",width =9.44,height = 4.14)