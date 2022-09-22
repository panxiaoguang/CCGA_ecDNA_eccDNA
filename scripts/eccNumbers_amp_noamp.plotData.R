# Figure4B-----plot -------------------------------------------------------


chromSize<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
total_chrom_size<-chromSize%>%
  filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
  summarise(total=sum(X2))%>%
  pull(total)


get_stat<-function(x,total_chrom_size=3088269832){
  seed_region<-read_tsv(stringr::str_glue("~/Project/Bladder/seeds/{x}_AA_CNV_SEEDS.bed"),col_names = F)
  seed_length<-seed_region%>%
    mutate(length=X3-X2)%>%
    summarise(total=sum(length))%>%
    pull(total)
  inters<-read_tsv(stringr::str_glue("~/Project/Bladder/seeds/{x}.inters.bed"),col_names = F)
  eccs<-read_tsv(stringr::str_glue("~/Project/Bladder/bedFile/cBca_{x}T.ecc.bed"),col_names = F)
  if(nrow(inters)!=nrow(eccs)){
    cat(x,"may be have some problems!")
  }
  inters<-inters%>%
    filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
    group_by(X5)%>%
    summarise(count=n())%>%
    filter(X5<2)%>%
    arrange(X5)%>%
    mutate(length=c(total_chrom_size-seed_length,seed_length))%>%
    mutate(count2=count/length,
           count3=count2/sum(count2))%>%
    select(X5,count3)%>%
    setNames(c("type","counts"))
  inters$sample<-x
  inters
}

samples<-c("1","2","5","8","9","13","14","15","16","17","21","23","24","25","28","29","30","31","32","33","34","36","37","39","40","41","42","44","46","47","48","50","51","54","56","57","60","62","63","65","67","68","69","70","71","72","74","75","76","77","79","80","84","85","86","87","88")

fin<-do.call("bind_rows",lapply(samples, function(x) get_stat(x)))

plotData<-fin%>%
  mutate(type=if_else(type==0,"non-amp","amp"))