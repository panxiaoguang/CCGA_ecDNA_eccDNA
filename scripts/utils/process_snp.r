library(readr)
library(dplyr)
args<-commandArgs()
args<-args[6:length(args)]
df<-read_tsv(args[1],col_names=F)
df<-df%>%select(1,2,14)
write_tsv(df,args[2])
df<-df%>%select(1,2)%>%distinct(X1,X2,.keep_all=T)
write_tsv(df,args[3],col_names=F)

