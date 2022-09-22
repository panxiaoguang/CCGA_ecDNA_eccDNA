library(ggplot2)
library(grid)
library(gridExtra)
pdf("MutationFreq.Sel.pdf",w=7,h=14)
rawdata <- read.table("MutationFreq.Sel.tsv",sep="\t",header=T)
data <- rawdata
data$Study <- factor(data$Study, levels = c("This study","TCGA 2017","CornellTrento 2016","DFCI MSKCC 2014","BGI 2013"))
data$Symbol <- factor(data$Symbol,levels = c("STAG2","RB1","EP300","KDM6A","KMT2C","KMT2D","MUC4","MUC16","TP53","TTN"))
ggplot(data,aes(x=Symbol,weight=PCT,fill=Study))+
  geom_bar(position="dodge")+
  scale_fill_manual(values=c("#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF","#79AF97FF"))+
  	ylim(c(0,100))+
	labs(title="",x="",y="Mutation Frequency(%)")+
  coord_flip()+
	theme(
	plot.margin=unit(c(0,2,1,0),"lines"),
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
	legend.position=c(0.7,0.2),
	legend.title=element_blank(),
	text=element_text(family="Helvetica",colour="black"),
	legend.text=element_text(size=22,colour="black"),
	axis.text=element_text(size=25,colour="black"),
	axis.title=element_text(size=25,colour="black"),
	strip.text=element_text(size=18+2)
	)
dev.off()
