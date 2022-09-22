# Figure4M-----plot -------------------------------------------------------

need_genes1<-openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx",sheet = 1)
need_genes1<-need_genes1%>%
  as_tibble()%>%
  filter(stringr::str_detect(Role,"oncogene"))

tmp=rst%>%mutate(rank=dense_rank(desc(cor)))%>%mutate(cosmic=if_else(Geneid%in%(need_genes1$Gene),"cosmic","non-cosmic"))%>%mutate(tp=if_else(cosmic=="cosmic","cosmic",tp))
tmp$tp<-factor(tmp$tp,levels = c("no","yes","cosmic"))

lizi<-tmp%>%
  filter(tp=="cosmic")%>%
  sample_frac(size=0.02)%>%
  arrange(rank)%>%
  pull(Geneid)

tmp<-tmp%>%
  mutate(label=if_else(Geneid%in%c("KDM5A",lizi,"PRIM2"),Geneid,""))

ggplot(tmp,aes(x=rank,y=cor))+
  geom_point(size=0.5)+
  geom_text_repel(aes(label=label),size=3.4,max.overlaps = 10000000000)+
  facet_wrap(~tp)+
  theme_pubr()
ggsave("cor.generank.pdf",width = 12.24,height = 3.73)