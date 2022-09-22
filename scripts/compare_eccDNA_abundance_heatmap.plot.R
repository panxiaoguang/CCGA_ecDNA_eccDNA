###first get abundance matrix
genes<-read_tsv("../qd-ECC4/S/ECC_report/FinallyData/bedFile/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)
get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("cgene_anno/{x}.startAnno.bed"),col_names = F)
  df<-df%>%
    select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    select(1,8)
  names(df)[2]<-x
  df
}
hebing<-function(x,y){full_join(x,y,by="gene")}
fin<-Reduce(hebing,lapply(c(cases,normals), get_genes))
fin<-fin%>%tidyr::gather(sample,value,-gene)%>%tidyr::replace_na(replace = list(value=0))%>%tidyr::pivot_wider(names_from = "sample",values_from = "value")
##calculate logfc
logfc<-fin%>%
  mutate(mean1=rowMeans(across(cBca_1T:cBca_88T)),
         mean2=rowMeans(across(cBca_1N:cBca_88N)))%>%
  select(gene,mean1,mean2)%>%
  mutate(logfc=log2(mean1)-log2(mean2))
### different analysis
dffs<-fin%>%
  tidyr::gather(sample,TPM,-gene)%>%
  tidyr::replace_na(replace = list(TPM=0))%>%
  mutate(gp=if_else(sample %in% cases,"Case","Normal"))%>%
  group_by(gene)%>%
  rstatix::pairwise_wilcox_test(TPM ~ gp,p.adjust.method = "BH")

dffs<-dffs%>%left_join(logfc,by="gene")
### prepare data for heatmap
plotData<-fin%>%tibble::column_to_rownames(var="gene")%>%as.matrix()
nima<-dffs%>%filter(p.adj<0.01,abs(logfc)>0.5)
nima<-nima%>%
  arrange(logfc)
plt<-plotData[nima$gene,]
plt[is.na(plt)]<-0
plt<-log2(plt+1)
col_fun<-colorRamp2(c(-2,0,2),c("#2166ac","#f7f7f7","#b2182b"))
ht<-Heatmap(plt2,
            show_row_names = F,
            show_column_names = F,
            cluster_rows = F,
            cluster_columns = F,
            #col = col_fun,
            top_annotation = HeatmapAnnotation(group=c(rep("Case",80),rep("Normal",80)),
                                               col = list(group=c("Case"="#e41a1c","Normal"="#377eb8"))),
            name = "Log normalized \neccDNA count",
            row_split = c(rep("down regular",326),rep("up regular",1019)))