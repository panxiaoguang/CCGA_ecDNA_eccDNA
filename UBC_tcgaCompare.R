set.seed(80)
library(maftools)
pdf("UBC_tcgaCompare.pdf",w=8,h=4)
ubc.maf <- "UBC_Combine.maf.tsv"
ubc = read.maf(maf = ubc.maf,verbose = FALSE)
tcgaCompare(maf = ubc, cohortName = 'UBC')
dev.off()