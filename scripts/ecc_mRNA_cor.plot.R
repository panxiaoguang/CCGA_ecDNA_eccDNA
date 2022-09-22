# Figure4L ------data --------------------------------------------------


fin<-readRDS("~/Project/Bladder/gene_drops.RDS")
fin<-fin%>%
  dplyr::select(1:81)%>%
  tidyr::gather(sample,value,-gene)%>%
  dplyr::rename(Geneid=gene,EPM=value)%>%
  tidyr::replace_na(replace = list(EPM=0))

TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
goodGenes<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>5,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==126)%>%
  pull(Geneid)

TPMs<-TPMs%>%
  filter(Geneid%in%goodGenes)%>%
  dplyr::select(1:71)%>%
  tidyr::gather(sample,TPM,-Geneid)

fin$sample<-stringr::str_remove(stringr::str_remove(fin$sample,"cBca_"),"T")
TPMs$sample<-stringr::str_remove(stringr::str_remove(TPMs$sample,"Bca_"),"_T")
fin<-inner_join(fin,TPMs,by=c("Geneid","sample"))
rst<-fin%>%
  mutate(logT=log2(TPM+1),
         logE=log2(EPM+1))%>%
  group_by(Geneid)%>%
  rstatix::cor_test(vars = "logE",vars2 = "logT",method = "spearman")