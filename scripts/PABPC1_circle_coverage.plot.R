# Figure4C-----plot -------------------------------------------------------


df<-read_tsv("PABPC1_geneRegion_readcoverage3.txt")
plotData<-df%>%
  mutate(coord=c(seq(100680816,100727809,by=50)[2:940],100727809))%>%
  tidyr::gather(sample,value,-coord)%>%
  mutate(group=if_else(sample %in%c("5T","16T","17T","21T","50T","84T"),"PABPC1-amplified","non-PABPC1-amplified"))

plotData$value2<-log2(plotData$value+1)

plotData3<-plotData%>%
  group_by(coord,group)%>%
  summarise(value3=list(mean_ci(value2)))%>%
  tidyr::unnest()%>%
  ungroup()

ggplot(plotData3,aes(x=coord,y=y,group=group))+
  geom_ribbon(aes(ymin = ymin, ymax = ymax),fill="#D3D2D3")+
  geom_line(aes(color=group))+
  geom_vline(xintercept = 100685816,color="grey")+
  geom_vline(xintercept = 100722809,color="grey")+
  scale_color_manual(values = c("#34327F","#CF3430"))+
  xlab("Genomic range")+
  ylab("Mean read coverage(log2)")+
  theme_pubr()