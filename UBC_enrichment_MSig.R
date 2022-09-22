set.seed(80)
library(maftools)

pdf("UBC_enrichment_MSig.pdf",w=7,h=3)

ubc.maf <- "UBC_Combine.maf.tsv"
ubc.clin <- "UBC-clinical.tsv"

ubc = read.maf(maf = ubc.maf, clinicalData = ubc.clin, verbose = FALSE)

#SMG=c("KDM6A","KIAA0040","STC1","TP53","PLEKHO1","CBX3","WIZ","CALY","RBMX","RETSAT","CENPB","STAG2","CDKN1A","AK2","FOXQ1","ZFP36L1","RB1","EP300","CHD1L","FANK1","GPATCH2L","KRT76","KMT2D","HLA-A")

MSig.ce = clinicalEnrichment(maf = ubc, clinicalFeature = 'MSig', minMut = 5, pathways = FALSE)
write.table(MSig.ce$groupwise_comparision, file="UBC_enrichment_MSig.tsv", sep ="\t", col.names = NA)
plotEnrichmentResults(enrich_res = MSig.ce, pVal = 0.01, geneFontSize = 0.9, annoFontSize = 0.9, legendFontSize = 0.9, cols = c("#E64B35FF","#4DBBD5FF","#00A087FF"))

#c("#E64B35FF","#4DBBD5FF","#00A087FF")

dev.off()