# Figure3I------data ------------------------------------------------------
#proteinGenes<-read_tsv("~/Project/Bladder/dbs/hg38.coding.bed",col_names = F)%>%
#  pull(X4)
#newdb<-dbs%>%
#  tidyr::separate(X4,into = c("gene","type"),sep = ":")%>%
#  filter(gene%in%proteinGenes)%>%
#  tidyr::unite("X4",gene,type,sep = ":")
dbs<-read_tsv("~/Project/Bladder/dbs/hg38.genetic_elements_exceptCPG_fix.bed",col_names = F)
length_corr<-dbs%>%
  tidyr::separate(X4,into = c("gene","type"),sep = ":")%>%
  mutate(length=X3-X2)%>%
  group_by(type)%>%
  summarise(total_l=sum(length))
total_chrom_length<-read_tsv("~/Project/甲状腺癌/old/plasma/dbs/hg38.chromo.size",col_names = F)
total_chrom_length<-sum((total_chrom_length[1:24,])$X2)
length_corr<-length_corr%>%
  mutate(pct=total_l/total_chrom_length)

cal_element<-function(x){
  fs<-read_tsv(stringr::str_glue("~/Project/Bladder/elementAnno/{x}.startAnno.bed"),col_names = F)
  tongji<-fs%>%
    select(1,2,3,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(ecc)%>%
    distinct(X8,.keep_all = T)%>%
    tidyr::separate(X8,into=c("gene","type"),sep=":")%>%
    ungroup()%>%
    count(type)
  tongji$samples<-x
  tongji
}
fin<-do.call('rbind',lapply(samples,cal_element))
### pie plot data
pie_plot_data<-fin%>%
  group_by(type)%>%
  summarise(counts=sum(n))%>%
  mutate(ratio=counts/sum(counts))

## correct by element length for barplot

barplot_data<-pie_plot_data%>%
  left_join(length_corr,by="type")%>%
  mutate(enrichment=ratio/pct)