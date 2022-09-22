library(Gviz)
# Figure4G-----plot -------------------------------------------------------

#######################WGS##################################

df<-read_tsv("PABPC1/PABPC1.WGS.allele.txt")
df<-df%>%
  filter(Good_depth>10)
wgs_alle<-df%>%
  tidyr::gather(alleleType,count,-`#CHR`,-POS,-Good_depth)%>%
  arrange(`#CHR`,POS,Good_depth,desc(count))%>%
  group_by(`#CHR`,POS,Good_depth)%>%
  slice_head(n=2)%>%
  ungroup()%>%
  mutate(tp=stringr::str_sub(alleleType,start=7,end=7))%>%
  mutate(AF=count/Good_depth)%>%
  mutate(tp2=rep(c("major","minor"),times=242))

data<-wgs_alle%>%
  select(`#CHR`,POS,AF,tp2)%>%
  tidyr::pivot_wider(names_from = "tp2",values_from = "AF")%>%
  arrange(POS)%>%
  as.data.frame()
gtrack <- GenomeAxisTrack(range=IRanges(start = 69253332,
                                        end = 70626244,
))

gr <- GRanges(seqnames = "chr8", strand = "*",
              ranges = IRanges(start = data$POS, width = 1),
              major=data[,3],minor=data[,4])

dtrack<-DataTrack(gr, name = "WGS AF",
                  groups=c("major","minor"),
                  col=c("#c5362c","#4475a7"),
                  type="p",
                  yTicksAt=c(0,0.5,1),
                  ylim=c(0,1),grid=T,lty.grid="dashed",lwd.grid=0.5,h=3)

wgsCOV<-read_tsv("PABPC1/WGS.cov.bdg",col_names = F)

cr <- GRanges(seqnames = "chr8", strand = "*",
              ranges = IRanges(start = wgsCOV$X2, end = wgsCOV$X3),
              count=wgsCOV$X4)
covtrack<-DataTrack(cr, name = "WGS cov",col="grey",type="histogram")

################################################################################
###############################RNAseq##########################################

df2<-read_tsv("PABPC1/PABPC1.RNA.allele.txt")
df2<-df2%>%
  filter(Good_depth>15)%>%
  mutate(Count_A=if_else(Count_A>2,Count_A,0),
         Count_C=if_else(Count_C>2,Count_C,0),
         Count_G=if_else(Count_G>2,Count_G,0),
         Count_T=if_else(Count_T>2,Count_T,0))%>%
  mutate(Good_depth=Count_A+Count_C+Count_G+Count_T)
RNA_alle<-df2%>%
  tidyr::gather(alleleType,count,-`#CHR`,-POS,-Good_depth)%>%
  mutate(tp=stringr::str_sub(alleleType,start=7,end=7),
         AF=count/Good_depth)%>%
  select(`#CHR`,POS,tp,AF)%>%
  dplyr::rename(RNAAF=AF)

data2<-wgs_alle%>%
  inner_join(RNA_alle,by=c("#CHR","POS","tp"))%>%
  select(`#CHR`,POS,RNAAF,tp2)%>%
  tidyr::pivot_wider(names_from = "tp2",values_from = "RNAAF")%>%
  arrange(POS)%>%
  as.data.frame()


gr2 <- GRanges(seqnames = "chr8", strand = "*",
               ranges = IRanges(start = data2$POS, width = 1),
               major=data2[,3],minor=data2[,4])

dtrack2<-DataTrack(gr2, 
                   name = "RNAcount",
                   groups=c("major","minor"),
                   col=c("#c5362c","#4475a7"),
                   yTicksAt=c(0,0.5,1),
                   ylim=c(0,1),grid=T,lty.grid="dashed",lwd.grid=0.5,h=3)

rnaCOV<-read_tsv("PABPC1/RNA.cov.bdg",col_names = F)

cr2 <- GRanges(seqnames = "chr8", strand = "*",
               ranges = IRanges(start = rnaCOV$X2, end = rnaCOV$X3),
               count=rnaCOV$X4)
covtrack2<-DataTrack(cr2, name = "RNA cov",col="grey",type="histogram",ylim=c(0,500))
##############################################################################
#####################################CIRCLE###################################
df3<-read_tsv("PABPC1/PABPC1.CIRCLE.allele.txt")
df3<-df3%>%
  filter(Good_depth>15)%>%
  mutate(Count_A=if_else(Count_A>2,Count_A,0),
         Count_C=if_else(Count_C>2,Count_C,0),
         Count_G=if_else(Count_G>2,Count_G,0),
         Count_T=if_else(Count_T>2,Count_T,0))%>%
  mutate(Good_depth=Count_A+Count_C+Count_G+Count_T)

circle_alle<-df3%>%
  tidyr::gather(alleleType,count,-`#CHR`,-POS,-Good_depth)%>%
  mutate(tp=stringr::str_sub(alleleType,start=7,end=7),
         AF=count/Good_depth)%>%
  select(`#CHR`,POS,tp,AF)%>%
  dplyr::rename(CIRCLEAF=AF)

data3<-wgs_alle%>%
  inner_join(circle_alle,by=c("#CHR","POS","tp"))%>%
  select(`#CHR`,POS,CIRCLEAF,tp2)%>%
  tidyr::pivot_wider(names_from = "tp2",values_from = "CIRCLEAF")%>%
  arrange(POS)%>%
  as.data.frame()

gr3 <- GRanges(seqnames = "chr8", strand = "*",
               ranges = IRanges(start = data3$POS, width = 1),
               major=data3[,3],minor=data3[,4])

dtrack3<-DataTrack(gr3, name = "CIRCLEcount",groups=c("major","minor"),col=c("#c5362c","#4475a7"),yTicksAt=c(0,0.5,1),
                   ylim=c(0,1),grid=T,lty.grid="dashed",lwd.grid=0.5,h=3)

circleCOV<-read_tsv("PABPC1/CIRCLE.cov.bdg",col_names = F)

cr3 <- GRanges(seqnames = "chr8", strand = "*",
               ranges = IRanges(start = circleCOV$X2, end = circleCOV$X3),
               count=circleCOV$X4)
covtrack3<-DataTrack(cr3, name = "CIRCLE cov",col="grey",type="histogram",ylim=c(0,1000))

png("test_region.png",width = 6.24,height =4.86,units = "in",res=300)
plotTracks(list(covtrack,dtrack,covtrack3,dtrack3,covtrack2,dtrack2,gtrack),
           from = 99497918,
           to=119476539,
           sizes=c(0.5,1,0.5,1,0.5,1,0.5),legend = F,
           background.title = "white",col.title="black",col.axis="black"
)
dev.off()
