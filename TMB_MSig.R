library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
pdf("TMB_MSig.pdf",w=4,h=7)
rawdata <- read.table("UBC_SomaticMutationsPerMb_Info.tsv",sep="\t",header=T)
data <- rawdata

dodge <- position_dodge(width = 0.4)
ggplot(data,aes(x=MSig,y=TMB_NonSyn,fill=MSig))+
  geom_boxplot(stat="boxplot",position="dodge")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF"))+
  ylim(c(0,52))+

  stat_compare_means(comparisons = list(c("MSig1","MSig2"),c("MSig2","MSig3"),c("MSig1","MSig3")),size=7)+
    labs(title="",x="",y="Mutations/Mb")+
	theme(
	plot.margin=unit(c(1,1,1,1),"lines"),
	plot.title = element_text(size = 28),
	panel.border=element_blank(),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	panel.spacing=unit(0.5,"lines"),
	strip.background=element_blank(),
	axis.line.x=element_line(colour="black",size=0.8),
	axis.line.y=element_line(colour="black",size=0.8),
	axis.ticks=element_line(size=0.8),
	axis.ticks.length=unit(0.3,"cm"),
	legend.title=element_blank(),
	legend.position="none",
	text=element_text(family="Helvetica",colour="black"),
	legend.text=element_text(size=20,colour="black"),
	axis.text=element_text(size=20,colour="black"),
	axis.title=element_text(size=20,colour="black"),
	strip.text=element_text(size=20)
	)
dev.off()
