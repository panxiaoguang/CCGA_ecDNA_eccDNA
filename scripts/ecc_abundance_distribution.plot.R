# Figure3G-----plot -------------------------------------------------------

EPMs<-readRDS("../Bladder/gene_drops.RDS")
newEPM<-EPMs%>%tidyr::gather(samples,counts,-gene)%>%mutate(gp=if_else(stringr::str_ends(samples,"T"),"Tumor","Normal"))
newEPM<-newEPM%>%
  tidyr::replace_na(replace = list(counts=0))%>%
  mutate(logc=log2(counts+1))%>%
  group_by(samples)%>%
  mutate(rank=dense_rank(desc(logc)))%>%
  ungroup()

ggplot(newEPM,aes(x=rank,y=logc))+
  geom_line(aes(group=samples,color=gp))+
  scale_color_manual(values = c("Tumor"="#C8413B","Normal"="#364BBA"))+
  xlab("Genes ranks")+
  ylab("eccDNA density")+
  facet_wrap(.~gp)+
  theme_pubr()+
  theme(strip.background = element_blank())