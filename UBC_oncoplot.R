set.seed(80)
library(maftools)
pdf("UBC_oncoplot.pdf",w=10,h=6)

ubc.maf <- "UBC_Combine.maf.tsv"
ubc.clin <- "UBC-clinical.tsv"

ubc = read.maf(maf = ubc.maf,
                clinicalData = ubc.clin,
                verbose = FALSE)
#By default the function plots top20 mutated genes

col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Multi_Hit','Nonsense_Mutation','Splice_Site','Frame_Shift_Ins','In_Frame_Ins',  'In_Frame_Del','Translation_Start_Site','Nonstop_Mutation')

colors_MSig = c("#E64B35FF","#4DBBD5FF","#00A087FF")
colors_Type = c("#9a1159","#E3F9A6")
colors_N = c("#CFECEC","#8EEBEC","#0AFFFF","#C0C0C0")
colors_M = c("#659EC7","#368BC1","#C0C0C0")
colors_Grade = c("#98cacd","#ff9700")
colors_Survival = c("#4E9258","#C0C0C0")
colors_Age = c("#34A56F","#CD853F")
colors_Gender = c("#0041C2","#E4287C")

names(colors_MSig) = c('MSig1','MSig2','MSig3')
names(colors_Type) = c('MIBC','NMIBC')
names(colors_N) = c('N0','N1','N2','Unknown')
names(colors_M) = c('M0','M1','Unknown')
names(colors_Grade) = c('Low','High')
names(colors_Survival) = c('Alive','Dead')
names(colors_Age) = c('<=65','>65')
names(colors_Gender) = c('Male','Female')

color2 = list(MSig = colors_MSig, Type = colors_Type, N=colors_N, M=colors_M, Grade = colors_Grade, Survival=colors_Survival, Age=colors_Age, Gender = colors_Gender)

SMG=c("KDM6A","KIAA0040","STC1","TP53","PLEKHO1","CBX3","WIZ","CALY","RBMX","RETSAT","CENPB","STAG2","CDKN1A","AK2","FOXQ1","ZFP36L1","RB1","EP300","CHD1L","FANK1","GPATCH2L","KRT76","KMT2D","HLA-A")

oncoplot(maf = ubc,
         colors = col,
         annotationColor = color2,
         genes = SMG,
         clinicalFeatures = c("MSig","Type","N","M","Grade","Survival","Age","Gender"),
         sortByAnnotation = TRUE,
         fontSize = 0.75,
         gene_mar = 6,
         showTumorSampleBarcodes = FALSE,
         removeNonMutated= FALSE,
         showTitle = FALSE,
         draw_titv = TRUE,
          legend_height = 5,
         anno_height =4)
dev.off()