# Figure4H,4J,4K----- data-----------------------------------------------------

AA_derived_CNV<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")
AA_derived_CNV$sample_name<-as.character(AA_derived_CNV$sample_name)

TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
allRNA_samples<-stringr::str_remove(stringr::str_remove(names(TPMs)[2:71],"Bca_"),"_T")

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
TPM_td<-TPM3%>%
  dplyr::rename(sample_name=sample,gene=Geneid)
TPM_td$sample_name<-stringr::str_remove(stringr::str_remove(TPM_td$sample_name,"Bca_"),"_T")

higher_AA<-AA_derived_CNV%>%
  filter(sample_name%in%allRNA_samples)

gene_in_sample<-higher_AA%>%
  group_by(gene)%>%
  summarise(nohave=paste0(setdiff(allRNA_samples,sample_name),collapse = ","))

higher_AA<-higher_AA%>%
  left_join(TPM_td,by=c("sample_name","gene"))
higher_AA<-na.omit(higher_AA)

ecc_amp<-higher_AA%>%
  filter(stringr::str_detect(feature,"ecDNA"))
non_ecc_amp<-higher_AA%>%
  filter(!stringr::str_detect(feature,"ecDNA"))

get_fold<-function(x){
  noneed<-higher_AA%>%
    filter(gene==x)%>%
    pull(sample_name)
  tmp<-TPM_td%>%
    filter(gene==x)%>%
    filter(!(sample_name%in%noneed))%>%
    summarise(noamp=mean(TPM))
  tmp$noamp
}

ecc_amp<-ecc_amp%>%
  mutate(noampTPM=purrr::map_dbl(gene,function(x){get_fold(x)}))
non_ecc_amp<-non_ecc_amp%>%
  mutate(noampTPM=purrr::map_dbl(gene,function(x){get_fold(x)}))
ecc_amp<-ecc_amp%>%
  mutate(FC=(TPM+1)/(noampTPM+1))
non_ecc_amp<-non_ecc_amp%>%
  mutate(FC=(TPM+1)/(noampTPM+1))

ecc_amp$type<-"ecDNA"
non_ecc_amp$type<-"others"
######for all protein coding
plotData<-bind_rows(ecc_amp,non_ecc_amp)
######for oncogenes
oncogenes<-openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx")
hao<-oncogenes%>%
  pull(Gene)%>%
  unique()
plotData2<-plotData%>%filter(gene%in%hao)
