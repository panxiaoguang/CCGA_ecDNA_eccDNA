# Figure2A--data ----------------------------------------------------------

df<-openxlsx::read.xlsx("~/Project/WGS/allsamples.ecDNA_det.xlsx")%>%
  as_tibble()
### define amplicon type
amp_type<-df%>%
  mutate(type=amplicon_decomposition_class,
         type=if_else(`ecDNA+`=="Positive","ecDNA",type),
         type=if_else(`BFB+`=="Positive","BFB",type))

get_sample_numbers<-function(x){
  amp_type%>%
  filter(type==x)%>%
  pull(sample_name)%>%
  unique()%>%
  length()
}
five_class<-c("ecDNA","BFB","Complex non-cyclic",
              "Linear amplification","No amp/Invalid")
amp_type$type<-factor(amp_type$type,levels = five_class)
## define sample type

  
sample_tables<-tibble(class=five_class,numbers=sapply(five_class,get_sample_numbers))
sample_tables