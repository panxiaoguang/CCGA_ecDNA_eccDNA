### Rscripts for getBAF
library(ASCAT)
get_BAF<-function(x){
    ascat.prepareHTS(stringr::str_glue("alignData/sorted_cBca_{x}T_circle.bam"),
                 stringr::str_glue("alignData/sorted_cBca_{x}N_circle.bam"),
                 stringr::str_glue("{x}T"),
                 stringr::str_glue("{x}N"),
                 "/dellfsqd2/ST_LBI/USER/panxiaoguang/app/miniconda3/envs/biosoft/bin/alleleCounter",
                 "/dellfsqd2/ST_LBI/USER/panxiaoguang/thyroid_cancer/WGS/CNV/G1000_alleles_hg38/G1000_alleles_hg38_chr",
                 "/dellfsqd2/ST_LBI/USER/panxiaoguang/thyroid_cancer/WGS/CNV/G1000_loci_hg38/G1000_loci_hg38_chr",
                 "XX",
                 "hg38",
                 nthreads = 10)
}

args <- commandArgs(trailingOnly = TRUE)

get_BAF(args[1])
###########################################################################

# Figure3L----plot --------------------------------------------------------


get_BAF<-function(x){
  df1<-read_tsv(stringr::str_glue("BAF/circleseq/AlleleCNV/{x}T_tumourBAF.txt"))%>%
    select(1,4)%>%
    setNames(c("pos","Circle-seq"))
  df2<-read_tsv(stringr::str_glue("BAF/wgs/WGS_BAF/{x}T_tumourBAF.txt"))%>%
    select(1,4)%>%
    setNames(c("pos","WGS"))
  its<-df1%>%
    inner_join(df2,by="pos")%>%
    filter(WGS!=0)%>%
    filter(WGS!=1)%>%
    tidyr::gather(type,BAF,-pos)
  its
}

samples<-c("1","2","3","4","5","6","7",
           "8","9","10","11","12","13",
           "14","15","16","17","18","19",
           "20","21","23","24","25","26",
           "27","28","29","30","31","32",
           "33","34","35","36","37","38",
           "39","40","41","42","43","44",
           "45","46","47","48","50",
           "51","52","54","55","56","57",
           "58","59","60","62","63","64",
           "65","66","67","68","69","70",
           "71","72","74","75","76","77",
           "78","79","80","84","85","86",
           "87","88")

plotData<-do.call("bind_rows",lapply(samples, get_BAF))

plotData<-do.call("bind_rows",lapply(c("11","12"), get_BAF))

ggplot(plotData,aes(x=BAF,
                    color=type,
                    fill=type,
                    y=..count../sum(..count..)/2))+
  geom_histogram(alpha=0.7,position = "identity")+
  scale_color_manual(values = c("Circle-seq"="#393A80","WGS"="#D1352B"))+
  scale_fill_manual(values = c("Circle-seq"="#393A80","WGS"="#D1352B"))+
  xlab("B-allele frequency")+
  ylab("Fraction of heterozygous\nSNPs")+
  theme_pubr()