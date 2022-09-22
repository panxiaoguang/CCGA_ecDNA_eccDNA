# Figure3B,C-----data -------------------------------------------------------

normals<-c("cBca_1N", "cBca_2N", "cBca_3N", "cBca_4N", "cBca_5N", "cBca_6N", "cBca_7N", "cBca_8N", "cBca_9N", "cBca_10N", "cBca_11N", "cBca_12N", "cBca_13N", "cBca_14N", "cBca_15N", "cBca_16N", "cBca_17N", "cBca_18N", "cBca_19T", "cBca_20T", "cBca_21N", "cBca_23N", "cBca_24N", "cBca_25N", "cBca_26N", "cBca_27N", "cBca_28N", "cBca_29N", "cBca_30N", "cBca_31N", "cBca_32N", "cBca_33N", "cBca_34N", "cBca_35N", "cBca_36N", "cBca_37N", "cBca_38N", "cBca_39T", "cBca_40N", "cBca_41N", "cBca_42N", "cBca_43N", "cBca_44N", "cBca_45N", "cBca_46N", "cBca_47N", "cBca_48N", "cBca_50N", "cBca_51N", "cBca_52N", "cBca_54N", "cBca_55N", "cBca_56N", "cBca_57N", "cBca_58N", "cBca_59N", "cBca_60N", "cBca_62N", "cBca_63N", "cBca_64N", "cBca_65N", "cBca_66N", "cBca_67N", "cBca_68N", "cBca_69N", "cBca_70N", "cBca_71N", "cBca_72N", "cBca_74N", "cBca_75N", "cBca_76N", "cBca_77N", "cBca_78N", "cBca_79N", "cBca_80N", "cBca_84N", "cBca_85N", "cBca_86N", "cBca_87N", "cBca_88N")

cases<-c("cBca_1T", "cBca_2T", "cBca_3T", "cBca_4T", "cBca_5T", "cBca_6T", "cBca_7T", "cBca_8T", "cBca_9T", "cBca_10T", "cBca_11T", "cBca_12T", "cBca_13T", "cBca_14T", "cBca_15T", "cBca_16T", "cBca_17T", "cBca_18T", "cBca_19N", "cBca_20N", "cBca_21T", "cBca_23T", "cBca_24T", "cBca_25T", "cBca_26T", "cBca_27T", "cBca_28T", "cBca_29T", "cBca_30T", "cBca_31T", "cBca_32T", "cBca_33T", "cBca_34T", "cBca_35T", "cBca_36T", "cBca_37T", "cBca_38T", "cBca_39N", "cBca_40T", "cBca_41T", "cBca_42T", "cBca_43T", "cBca_44T", "cBca_45T", "cBca_46T", "cBca_47T", "cBca_48T", "cBca_50T", "cBca_51T", "cBca_52T", "cBca_54T", "cBca_55T", "cBca_56T", "cBca_57T", "cBca_58T", "cBca_59T", "cBca_60T", "cBca_62T", "cBca_63T", "cBca_64T", "cBca_65T", "cBca_66T", "cBca_67T", "cBca_68T", "cBca_69T", "cBca_70T", "cBca_71T", "cBca_72T", "cBca_74T", "cBca_75T", "cBca_76T", "cBca_77T", "cBca_78T", "cBca_79T", "cBca_80T", "cBca_84T", "cBca_85T", "cBca_86T", "cBca_87T", "cBca_88T")
cal_total_eccs<-function(x){
  read_tsv(stringr::str_glue("~/Project/Bladder/filter_eccs/{x}_circle_site.filter.tsv"))%>%
    nrow()
}

eccnumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),cal_total_eccs))

cal_ingenes<-function(x){
  read_tsv(stringr::str_glue("~/Project/Bladder/gene_anno/{x}.startAnno.bed"),col_names = F)%>%
    select(X4)%>%
    distinct(X4)%>%
    nrow()
}
genenumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),cal_ingenes))
finally<-full_join(eccnumbers,genenumbers,by="sample")
finally<-finally%>%
  mutate(geneP=eccs.y/eccs.x*100)

calrepeat<-function(x){
  df<-read_tsv(stringr::str_glue("~/Project/Bladder/repeatStat/repStats.{x}.tsv"))
  sum(df$ratio)
}
repnumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),calrepeat))

all<-finally%>%full_join(repnumbers,by="sample")



