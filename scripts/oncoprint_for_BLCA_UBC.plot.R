########CCGA PART#####################################################################################
df<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")
need_genes1<-openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx",sheet = 1)
need_genes1<-need_genes1%>%
  as_tibble()%>%
  filter(stringr::str_detect(Role,"oncogene"))
df1<-df%>%
  filter(gene%in%need_genes1$Gene)
mat<-df1%>%as_tibble()%>%
  select(sample_name,feature,gene)%>%
  mutate(feature=stringr::str_remove(feature,"_\\d"))%>%
  tidyr::pivot_wider(names_from = "sample_name",
                     values_from = "feature",
                     values_fn =function(x){paste(x,collapse = ";")},values_fill = "")%>%
  tibble::column_to_rownames(var="gene")%>%as.matrix()

needs<-df1%>%group_by(gene)%>%summarise(ct=n())%>%arrange(desc(ct))%>%filter(ct>2)%>%pull(gene)
needs<-needs[1:20]

mat2<-mat[needs,]
#############################################################################################
########################################TCGA PART################################################
TCGA<-openxlsx::read.xlsx("TCGA_UBC.xlsx")
TCGAbed<-TCGA%>%
  as_tibble()%>%
  mutate(info=paste(sample_barcode,
                    amplicon_index,
                    amplicon_classification,sep = "&"))%>%
  select(amplicon_intervals,info)%>%
  tidyr::separate_rows(amplicon_intervals,sep=",")%>%
  tidyr::separate(amplicon_intervals,into = c("chrom","start","end"),sep = "[:-]")%>%
  mutate(chrom=paste0("chr",chrom))

tcga_anno<-read_tsv("TCGA_data.anno.bed",col_names = F)
mat3<-tcga_anno%>%
  select(X4,X8)%>%
  tidyr::separate(X4,into = c("sample","amp","feature"),sep="&")%>%
  mutate(feature=case_when(feature == "BFB" ~ "BFB",
                           feature == "Circular" ~ "ecDNA",
                           feature == "Heavily-rearranged" ~ "Complex non-cyclic",
                           feature == "Linear" ~ "Linear amplification"))%>%
  filter(X8%in%(need_genes1$Gene))%>%
  select(sample,feature,X8)%>%
  dplyr::rename(gene=X8)%>%
  tidyr::pivot_wider(names_from = "sample",
                     values_from = "feature",
                     values_fn =function(x){paste(x,collapse = ";")},values_fill = "")%>%
  tibble::column_to_rownames(var="gene")%>%as.matrix()

mat4<-mat3[needs,]
#############################################################################################
final_data<-cbind(mat2,mat4)
nima=tibble(genes=rownames(final_data),num=as.character(apply(final_data,1,function(x){sum(x=="ecDNA")})))

col = c(BFB = "#436693", `Complex non-cyclic` = "#C58D65",
        ecDNA="#BC4137",`Linear amplification`="#A1B2C0",
        unknown="#ACABAC")
##############################################anno part#############################################
anno_table<-openxlsx::read.xlsx("TCGA-CCGA.xlsx",sheet = 2)
anno_table<-anno_table%>%
  mutate(label=stringr::str_extract(Sample,"\\d+"))
anno_table<-anno_table%>%
  filter(label%in%(colnames(mat2)))%>%
  tidyr::replace_na(replace = list("N"="No avalible","M"="No avalible"))%>%
  select(label,Survival,Gender,Age,N,M,`NMIBC/MIBC`,Grade)

anno_table2<-openxlsx::read.xlsx("TCGA-CCGA.xlsx",sheet = 1)
colnames(anno_table2)[6]<-"NMIBC/MIBC"
anno_table2<-anno_table2%>%
  filter(sample_barcode%in%(colnames(mat4)))%>%
  tidyr::replace_na(replace = list("NMIBC/MIBC"="No avalible","N"="No avalible","M"="No avalible"))%>%
  select(`sample_barcode`,Survival,Gender,Age,N,M,`NMIBC/MIBC`,Grade)%>%
  dplyr::rename(label="sample_barcode")%>%
  mutate(Grade=stringr::str_remove(Grade," "))

fin_anno<-bind_rows(anno_table,anno_table2)
fin_anno<-fin_anno%>%tibble::column_to_rownames(var="label")%>%
  as.data.frame()

fin_anno<-fin_anno[colnames(final_data),]
###############################################################################################
############################################plot oncoprint#######################################
ht<-oncoPrint(final_data,
              alter_fun = function(x, y, w, h, v) {
                n = sum(v)
                h = h*0.9
                if(n){grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
                                gp = gpar(fill = col[names(which(v))],lwd=0.5,col="white"), just = "top")
                }else{grid.rect(x,y,w,h,gp = gpar(fill = "#E4E4E4",lwd=0.5,col="white"))}}, col = col,
              row_names_gp = gpar(fontface="italic",fontsize=8.5),
              column_split = c(rep("CCGA",46),rep("TCGA",68)),
              row_order = order(nima$num,decreasing = T),
              top_annotation = 
                HeatmapAnnotation(#cbar = anno_oncoprint_barplot(),
                  Age = fin_anno$Age,
                  N = fin_anno$N,
                  M = fin_anno$M,
                  Survival = fin_anno$Survival,
                  Gender = fin_anno$Gender,
                  `NMIBC/MIBC`=fin_anno$`NMIBC/MIBC`,
                  Grade=fin_anno$Grade,
                  simple_anno_size = unit(0.3, "cm"),
                  annotation_name_gp = gpar(fontsize=8.5),
                  col = list(Age=c("<=65"="#B9DFFB",">65"="#68C84D"),
                             N = c("N0"="#F3F3F4","N1"="#ABDAE4","N2"="#4B95E9","N3"="#123294","No avalible"="#747070"),
                             M = c("M0"="#F3F3F4","M1"="#F09F37","No avalible" = "#747070"),
                             Survival =c("Alive"="#F3F3F4","Death"="#010101"),
                             Gender=c("FEMALE"="#E93420","MALE"="#316DBB"),
                             `NMIBC/MIBC`=c("MIBC"="#AE2417","NMIBC"="#F3F3F4","No avalible"="#747070"),
                             Grade=c("High"="#4FADEB","Low"="#F3F3F4")
                  )
                ))
pdf("oncoprint5.pdf",width = 12.55,height =4.37)
draw(ht)
dev.off()