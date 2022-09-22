# Figure3M ------------data---------------------------------------------------


samples<-c("5N", "5T", "13N", "13T", "15N", "15T", "16N", "16T", "17N", "17T", "21N", "21T", "29N", "29T", "37N", "37T", "50N", "50T")

get_nsegment<-function(x){
  test<-read_tsv(stringr::str_glue("~/Project/Bladder/longreads/eccDNAs/{x}.info.txt"))
  test<-test%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)
  test$label<-seq(1,nrow(test))
  haha<-test%>%
    select(label,fragments)%>%
    tidyr::separate_rows(fragments,sep="\\|")%>%
    mutate(orign=sapply(stringr::str_split(fragments,":"),function(x) x[[1]]))%>%
    group_by(label)%>%
    summarise(count=n())%>%
    group_by(count)%>%
    summarise(num=n())
  haha$sample<-x
  haha
}

plotData<-do.call("bind_rows",lapply(samples, get_nsegment))
plotData<-plotData%>%
  tidyr::pivot_wider(names_from = count,
                     values_from = num,values_fill = 0)