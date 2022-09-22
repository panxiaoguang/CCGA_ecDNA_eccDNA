# Figure2B---data ---------------------------------------------------------

sample_tables<-amp_type%>%
  arrange(sample_name,type)%>%
  mutate(type=as.character(type))%>%
  group_by(sample_name)%>%
  tidyr::nest()%>%
  mutate(sampletype=purrr::map(data,function(x){(x$type)[[1]]}))%>%
  select(sample_name,sampletype)%>%
  ungroup()%>%
  group_by(sampletype)%>%
  summarise(count=n())

Pvalue<-fisher.test(matrix(c(28,45,29,3,10,5,4,4,29,23),nrow = 2))

Pvalue$p.value