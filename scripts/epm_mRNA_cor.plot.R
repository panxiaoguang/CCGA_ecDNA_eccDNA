# Figure6A,6B,6C-----data -------------------------------------------------

##exps
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
TPM2<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
TPM3<-TPM2%>%filter(stringr::str_ends(sample,"T"))
TPM3$sample<-stringr::str_remove(TPM3$sample,"_T")

##group
groupInfo<-read_tsv("class_sample.tsv")
names(groupInfo)<-c("sample","type")
groupInfo<-groupInfo%>%
  mutate(type=if_else(type=="ecDNA","ecDNA","others"))
TPM3<-TPM3%>%
  left_join(groupInfo,by="sample")

EPMS<-openxlsx::read.xlsx("~/Project/Bladder/Bladder_circle_numbers.xlsx")%>%
  as_tibble()
EPMs<-EPMS%>%
  filter(group=="Tumour")%>%
  dplyr::select(sample,EPM)
EPMs$sample<-stringr::str_remove(EPMs$sample,"c")
EPMs$sample<-stringr::str_remove(EPMs$sample,"T")
EPMs$sample<-stringr::str_remove(EPMs$sample,"N")
TPM3<-TPM3%>%
  left_join(EPMs,by="sample")

TPM3<-TPM3%>%
  mutate(logT=log2(TPM+1),logE=log2(EPM+1))

goodGenes<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  filter(stringr::str_ends(sample,"T"))%>%
  mutate(type=if_else(TPM>1,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==70)%>%
  pull(Geneid)

diffGenes<-TPM3%>%
  dplyr::filter(Geneid%in%goodGenes)%>%
  dplyr::group_by(Geneid)%>%
  rstatix::wilcox_test(logT~type)

diffCors2<-TPM3%>%
  dplyr::group_by(Geneid)%>%
  rstatix::cor_test(vars = "logT",vars2 = "logE",method = "spearman")

diffGene1<-diffGenes%>%
  filter(p<0.05)
diffCor1<-diffCors%>%
  filter(p<0.05)

jiaoji2<-intersect(diffGene1$Geneid,diffCor1$Geneid)

zhuanhua<-bitr(jiaoji,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

ego<-enrichGO(zhuanhua2$ENTREZID,OrgDb =org.Hs.eg.db,pvalueCutoff = 0.05,ont = 'BP',readable = T )

eko<-enrichKEGG(zhuanhua2$ENTREZID,organism = "hsa",pvalueCutoff = 0.05)