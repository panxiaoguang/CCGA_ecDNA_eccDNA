library(ggplot2)
pdf("TumorVAF_Kataegis_Kyklonas.pdf",w=6,h=3)
rawdata <- read.table("UBC_Combine_GW_with_Kataegis_kyklonas_EC.maf.tsv",sep="\t",header=T)
data <- rawdata
ggplot(data,aes(x=TumorVAF,fill=Kataegis,color=Kataegis))+
  geom_density(alpha=0.75, colour = NA) +
	labs(title="",x="Tumor VAF",y="Density") +
  xlim(c(0,1))+
  scale_fill_manual(values=c("#374E55FF","#B24745FF"))+
 	theme(
	plot.margin=unit(c(0,1,1,1),"lines"),
	panel.border=element_blank(),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	panel.spacing=unit(1,"lines"),
	strip.background=element_blank(),
	axis.line.x=element_line(colour="black",size=0.8),
	axis.line.y=element_line(colour="black",size=0.8),
	axis.ticks=element_line(size=0.8),
	axis.ticks.length=unit(0.3,"cm"),
	legend.title=element_blank(),
	legend.position=c(0.85,0.85),
	text=element_text(family="Helvetica",colour="black"),
	legend.text=element_text(size=18,colour="black"),
	axis.text=element_text(size=18,colour="black"),
	axis.title=element_text(size=18,colour="black"),
	strip.text=element_text(size=18)
	)
dev.off()
