library(ggplot2)
library(grid)
#library(Rmisc)
library(gridExtra)
pdf("Stat_Kataegis_WithinEvents_or_10k.pdf",w=6,h=5)
rawdata <- read.table("Stat_Kataegis_WithinEvents_or_10k.tsv",sep="\t",header=T)
data <- rawdata
#options(digits=3)

data$Group <- factor(data$Group,level=c("SV","SV (10kb)","ecDNA","ecDNA (10kb)","rearranged","rearranged (10kb)","BFB","Linear"))

#ggplot(data,aes(x=reorder(Group,-PCT),weight=PCT,fill=Group))+
ggplot(data,aes(x=Group,weight=PCT,fill=Group))+  
  geom_bar(position="dodge",color="black")+
#  scale_fill_manual(values=c("#660033", "#003366","#404040"))+
  ylim(c(0,25))+
	labs(title="",x="",y="Proportion of kataegic mutations (%)")+
  geom_text(aes(y=PCT,label=round(PCT,digits=2)),position=position_dodge(width=0.9),vjust=0,size=6)+
# facet_wrap(~Region,nrow=1,ncol=12)+
	theme(
	plot.margin=unit(c(1,1,-1,1),"lines"),
	panel.border=element_blank(),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	panel.spacing=unit(2,"lines"),
	strip.background=element_blank(),
	axis.line.x=element_line(colour="black",size=0.8),
	axis.line.y=element_line(colour="black",size=0.8),
	axis.ticks=element_line(size=0.8),
	axis.ticks.length=unit(0.3,"cm"),
#	legend.position=c(0.9,0.75),
  legend.position="none",
	legend.title=element_blank(),
	text=element_text(family="Helvetica",colour="black"),
	legend.text=element_text(size=18,colour="black"),
	axis.text=element_text(size=18,colour="black"),
#	axis.text.x=element_text(angle=0),
#	axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
  axis.text.x=element_text(angle=45,hjust=1),
	axis.title=element_text(size=18,colour="black"),
	strip.text=element_text(size=18)
	)
dev.off()
