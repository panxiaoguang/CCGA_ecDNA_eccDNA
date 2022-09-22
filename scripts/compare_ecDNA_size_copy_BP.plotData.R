# Figure2D,2E,2F -----data ------------------------------------------------

amp_class<-openxlsx::read.xlsx("all_sum.xlsx")%>%
  as_tibble(0)

newClass<-amp_class%>%
  as_tibble()%>%
  filter(AmpliconType!="No amp/Invalid")%>%
  select(1:3,7)
names(newClass)<-c("sample","amp","type","length")
breakPoints<-read_tsv("graphs/all.detective.breakpoints.tsv")

breakPoints<-breakPoints%>%
  left_join(newClass,by=c("sample","amp"))
wocao<-breakPoints%>%
  mutate(name=paste(sample,amp,sep = "_"),start=bp-1)%>%
  select(chrom,start,bp,type)
wocao<-na.omit(wocao)

plotData<-breakPoints%>%
  na.omit()%>%
  group_by(sample,amp,type,length)%>%
  summarise(count=n())%>%
  mutate(ratio=count/length*1000000)

plotData<-newClass%>%
  left_join(plotData,by=c("sample","amp","type","length"))%>%
  tidyr::replace_na(replace = list(count=0,ratio=0,logratio=0))

amp_class%>%
  filter(AmpliconType!="No amp/Invalid")%>%
  rstatix::wilcox_test(AmplifiedIntervalSize~AmpliconType)%>%
  filter(group1=="ecDNA"|group2=="ecDNA")

amp_class%>%
  filter(AmpliconType!="No amp/Invalid")%>%
  rstatix::wilcox_test(AverageAmplifiedCopyCount~AmpliconType)%>%
  filter(group1=="ecDNA"|group2=="ecDNA")

plotData%>%
  rstatix::wilcox_test(ratio~type)%>%
  filter(group1=="ecDNA"|group2=="ecDNA")
##ecDNA_linear = 1.64 e-16 , ecDNA_BFB = 3.11e-1,ecDNA_complex = 4.4e-2
##ecDNA_linear = 3.89 e-22 , ecDNA_BFB = 1.2e-1,ecDNA_complex = 7.15e-10
##ecDNA_linear = 4.87 e-13 , ecDNA_BFB = 2.7e-2,ecDNA_complex = 3.9e-1