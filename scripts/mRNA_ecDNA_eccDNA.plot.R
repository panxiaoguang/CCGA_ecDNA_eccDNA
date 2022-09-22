# Figure4E----plot --------------------------------------------------------

TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
## should filter
goodGenes<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>5,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==126)%>%
  pull(Geneid)
TPMs<-TPMs%>%
  filter(Geneid%in%goodGenes)
TPM2<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
TPM3<-TPM2%>%filter(stringr::str_ends(sample,"T"))


## mean+sd
stst<-TPM3%>%
  group_by(Geneid)%>%
  summarise(mean=mean(TPM),sd=sd(TPM))
## zscore
zscores<-TPM3%>%
  left_join(stst,by="Geneid")%>%
  mutate(zscore=(TPM-mean)/sd)

zscores<-zscores%>%
  dplyr::select(Geneid,sample,zscore)%>%
  dplyr::rename(TPM=zscore)%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))

eccgenes<-read_tsv("~/Project/Bladder/allgene20.bed",col_names = F)
eccgenes<-eccgenes%>%
  mutate(sample=stringr::str_extract(X4,"cBca_\\d+[NT]"))%>%
  dplyr::select(sample,X8)%>%
  distinct(sample,X8,.keep_all = T)
names(eccgenes)[2]<-"Geneid"

part1<-eccgenes%>%
  filter(sample%in%c("cBca_19N","cBca_20N","cBca_39N"))
part2<-eccgenes%>%
  filter(stringr::str_ends(sample,"T"))%>%
  filter(!(sample%in%c("cBca_19T","cBca_20T","cBca_39T")))
eccgenes<-bind_rows(part1,part2)

eccgenes<-eccgenes%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))%>%
  left_join(zscores,by=c("sample","Geneid"))

ecgenes<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")%>%
  as_tibble()%>%
  filter(stringr::str_starts(feature,"ecDNA"))%>%
  dplyr::select(sample_name,gene)%>%
  distinct(sample_name,gene,.keep_all = T)%>%
  dplyr::rename(sample=sample_name,Geneid=gene)
ecgenes$sample<-as.character(ecgenes$sample)
ecgenes<-ecgenes%>%
  left_join(zscores,by=c("sample","Geneid"))%>%
  na.omit()
ecgenes<-ecgenes%>%
  mutate(ding=paste(sample,Geneid,sep="-"))
eccgenes<-eccgenes%>%
  mutate(ding=paste(sample,Geneid,sep="-"))
### use all focal-amp genes to filter
allGenes<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")%>%
  as_tibble()%>%
  dplyr::select(sample_name,gene)%>%
  distinct(sample_name,gene,.keep_all = T)%>%
  dplyr::rename(sample=sample_name,Geneid=gene)%>%
  mutate(ding=paste(sample,Geneid,sep="-"))
jiaoji<-intersect(allGenes$ding,eccgenes$ding)
eccgenes<-eccgenes%>%
  filter(!ding%in%jiaoji)
eccgenes$type<-"eccDNA"
ecgenes$type<-"ecDNA"

fin<-bind_rows(ecgenes,eccgenes)

fin<-fin%>%
  na.omit()

ggplot(fin,aes(x=TPM,color=type))+
  geom_density()+
  scale_color_manual(values = c("#9E3735","#48436D"))+
  theme_pubr()+
  xlab("RNA expression(Z-score)")

# Figure4F----plot --------------------------------------------------------

ecgenes<-ecgenes%>%filter(sample=="16")

plotData<-TPM3%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))%>%
  filter(sample=="16")%>%
  mutate(rank=dense_rank(desc(TPM)))%>%
  mutate(ecgene=if_else(Geneid%in%c("PABPC1","BRK1","YWHAZ","TADA3","CTNNB1","OXR1","UBR5","ANKRD46"),Geneid,""))%>%
  mutate(col=if_else(ecgene=="","no","yes"))


ggplot(plotData,aes(x=rank,y=TPM))+
  geom_line()+
  geom_point(aes(color=col))+
  scale_y_log10()+
  scale_color_manual(values = list(yes="red",no="transparent"))+
  geom_text_repel(aes(label=ecgene),max.overlaps = 100000000000)+
  theme(legend.position="none")+
  xlab("Genes ranked by expression level\n(TPM>5)")+
  ylab("TPM")+
  theme_pubr()