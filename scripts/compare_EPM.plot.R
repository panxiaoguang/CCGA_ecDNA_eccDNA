# Figure3A-----plot -------------------------------------------------------

EPM<-openxlsx::read.xlsx("~/Project/Bladder/Bladder_circle_numbers.xlsx")%>%
  as_tibble()
EPM<-EPM%>%mutate(paired=stringr::str_extract(sample,"cBca_\\d+"))
EPM$paired<-factor(EPM$paired,levels = unique(stringr::str_sort(EPM$paired,numeric = T)))
NAT<-EPM%>%filter(group=="Normal")  
TUM<-EPM%>%filter(group=="Tumour")
NAT<-NAT%>%
  arrange(paired)%>%
  mutate(label=seq(1,80))
TUM<-TUM%>%
  arrange(paired)%>%
  mutate(label=seq(1,80))
plotData<-bind_rows(NAT,TUM)
## one line to plot figures
ggplot(plotData,aes(x=label,y=EPM))+
  geom_line(aes(group=paired),
            color="#bdbdbd",
            size=0.2)+
  geom_point(aes(color=group))+
  geom_smooth(aes(color=group,fill=group),
              linetype="longdash",
              size=0.7,
              method = "loess",
              alpha=0.5)+
  scale_color_manual(values = c("Normal"="#263272",
                                "Tumour"="#B83C3E"))+
  scale_fill_manual(values = c("Normal"="#263272",
                               "Tumour"="#B83C3E"
  ))+
  theme_pubr()+
  xlab("Samples")+
  ylab("Number of eccDNAs per \nmillion mapped reads")