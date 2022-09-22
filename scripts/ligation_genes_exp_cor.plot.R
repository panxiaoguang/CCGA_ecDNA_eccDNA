# Figure6D----plot --------------------------------------------------------

genes<-c("LIG3","LIG4","POLM","POLQ","PRKDC","BRCA1","BRCA2","MSH3")
exp_plot<-TPM3%>%
  filter(Geneid%in%genes)

p1<-ggplot(exp_plot%>%filter(Geneid==genes[1]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[1]} expression value"))+
  theme(legend.position = "none")

p2<-ggplot(exp_plot%>%filter(Geneid==genes[2]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[2]} expression value"))+
  theme(legend.position = "none")

p3<-ggplot(exp_plot%>%filter(Geneid==genes[3]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[3]} expression value"))+
  theme(legend.position = "none")

p4<-ggplot(exp_plot%>%filter(Geneid==genes[4]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[4]} expression value"))+
  theme(legend.position = "none")

p5<-ggplot(exp_plot%>%filter(Geneid==genes[5]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[5]} expression value"))+
  theme(legend.position = "none")

p6<-ggplot(exp_plot%>%filter(Geneid==genes[6]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[6]} expression value"))+
  theme(legend.position = "none")

p7<-ggplot(exp_plot%>%filter(Geneid==genes[7]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[7]} expression value"))+
  theme(legend.position = "none")

p8<-ggplot(exp_plot%>%filter(Geneid==genes[8]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[8]} expression value"))+
  theme(legend.position = "none")

ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow=2,align = "v")

ggsave("somegenes.plot1.pdf",width = 6.85,height = 5.3)

# Figure6E-----plot -------------------------------------------------------


p1<-ggplot(exp_plot%>%filter(Geneid==genes[1]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[1]} eccDNA value"))

p2<-ggplot(exp_plot%>%filter(Geneid==genes[2]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[2]} eccDNA value"))

p3<-ggplot(exp_plot%>%filter(Geneid==genes[3]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[3]} eccDNA value"))

p4<-ggplot(exp_plot%>%filter(Geneid==genes[4]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[4]} eccDNA value"))

p5<-ggplot(exp_plot%>%filter(Geneid==genes[5]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[5]} eccDNA value"))

p6<-ggplot(exp_plot%>%filter(Geneid==genes[6]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[6]} eccDNA value"))

p7<-ggplot(exp_plot%>%filter(Geneid==genes[7]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[7]} eccDNA value"))

p8<-ggplot(exp_plot%>%filter(Geneid==genes[8]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[8]} eccDNA value"))

ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow=2,align = "v")